#!/usr/bin/env python

"""
Renders the cladebreaker analyze results as a PDF: the tree, then each test.

Page 1 carries the rooted tree the statistics were actually computed on, with
tips coloured or marked by group, above the report's conclusion. The pages after
it give one section per test.

DEPENDENCIES: matplotlib only. The Newick reader here is deliberately small
because it never sees arbitrary input -- it reads the normalised tree written by
monophyly.py, which dendropy produced. It is checked against dendropy in the
test suite.

COLOUR: group colours come from a validated categorical palette, used for at
most three groups. A fourth slot puts yellow and orange on the page together and
that pair falls below the normal-vision separation floor, so beyond three groups
identity moves to marker shape and colour is dropped rather than invented. Tip
labels are always drawn, so identity is never carried by colour alone.

Example usage:
    analyze_pdf.py --tree rooted.nwk --out-prefix cladebreaker *.json
"""

import argparse
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")                      # no display in a pipeline task
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

# ---------------------------------------------------------------------------
# Palette (light surface; a PDF is printed on white)
# ---------------------------------------------------------------------------

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_FAINT = "#8a8a86"
RULE = "#d8d8d4"

# Validated categorical slots 1-3: all-pairs CVD dE 9.2, normal-vision dE 24.0.
GROUP_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a"]
MAX_COLOURED_GROUPS = len(GROUP_COLOURS)
# Secondary encoding, always applied so colour is never the only identity channel.
GROUP_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

WITHIN_COLOUR = "#2a78d6"    # slot 1
BETWEEN_COLOUR = "#eb6834"   # slot 2

SUFFIXES = {
    ".monophyly.json": "monophyly",
    ".rosenberg.json": "rosenberg",
    ".slatkin_maddison.json": "slatkin_maddison",
    ".gsi.json": "gsi",
    ".snp_separation.json": "snp_separation",
    ".analyze_report.json": "report",
}


# ---------------------------------------------------------------------------
# Newick
# ---------------------------------------------------------------------------


class Node(object):
    __slots__ = ("label", "length", "children", "x", "y")

    def __init__(self):
        self.label = None
        self.length = 0.0
        self.children = []
        self.x = 0.0
        self.y = 0.0

    @property
    def is_leaf(self):
        return not self.children


def parse_newick(text):
    """Read one Newick tree. Handles quoted labels, branch lengths, comments."""
    text = text.strip()
    if not text:
        raise ValueError("empty Newick string")
    text = text.split(";")[0] + ";"

    pos = [0]

    def error(message):
        raise ValueError("{} at character {}".format(message, pos[0]))

    def skip():
        while pos[0] < len(text):
            char = text[pos[0]]
            if char.isspace():
                pos[0] += 1
            elif char == "[":                       # a comment, e.g. support values
                depth = 0
                while pos[0] < len(text):
                    if text[pos[0]] == "[":
                        depth += 1
                    elif text[pos[0]] == "]":
                        depth -= 1
                        if depth == 0:
                            pos[0] += 1
                            break
                    pos[0] += 1
            else:
                return

    def read_label():
        skip()
        if pos[0] >= len(text):
            return None
        if text[pos[0]] in "'\"":
            quote = text[pos[0]]
            pos[0] += 1
            chars = []
            while pos[0] < len(text):
                if text[pos[0]] == quote:
                    # a doubled quote inside a quoted label is an escaped quote
                    if pos[0] + 1 < len(text) and text[pos[0] + 1] == quote:
                        chars.append(quote)
                        pos[0] += 2
                        continue
                    pos[0] += 1
                    break
                chars.append(text[pos[0]])
                pos[0] += 1
            return "".join(chars)
        start = pos[0]
        while pos[0] < len(text) and text[pos[0]] not in "(),:;[]":
            pos[0] += 1
        return text[start:pos[0]].strip()

    def read_length():
        skip()
        if pos[0] < len(text) and text[pos[0]] == ":":
            pos[0] += 1
            skip()
            start = pos[0]
            while pos[0] < len(text) and (text[pos[0]].isdigit() or text[pos[0]] in ".eE+-"):
                pos[0] += 1
            try:
                return float(text[start:pos[0]])
            except ValueError:
                error("unreadable branch length '{}'".format(text[start:pos[0]]))
        return 0.0

    def read_node():
        skip()
        node = Node()
        if pos[0] < len(text) and text[pos[0]] == "(":
            pos[0] += 1
            while True:
                node.children.append(read_node())
                skip()
                if pos[0] >= len(text):
                    error("unbalanced parentheses")
                if text[pos[0]] == ",":
                    pos[0] += 1
                    continue
                if text[pos[0]] == ")":
                    pos[0] += 1
                    break
                error("expected ',' or ')' but found '{}'".format(text[pos[0]]))
        label = read_label()
        node.label = label or None
        node.length = read_length()
        return node

    root = read_node()
    skip()
    if pos[0] >= len(text) or text[pos[0]] != ";":
        error("tree does not end with ';'")
    return root


def leaves(node):
    found = []
    stack = [node]
    while stack:
        current = stack.pop()
        if current.is_leaf:
            found.append(current)
        else:
            stack.extend(reversed(current.children))
    return found


def postorder(node):
    order, stack = [], [node]
    while stack:
        current = stack.pop()
        order.append(current)
        stack.extend(current.children)
    return list(reversed(order))


def layout(root, use_lengths=True):
    """Rectangular phylogram coordinates.

    Leaves are spread evenly on y in the order they appear; each internal node
    sits at the midpoint of its children. x is cumulative branch length, or
    depth when the tree has no lengths to speak of.
    """
    order = postorder(root)
    if use_lengths:
        total = sum(node.length or 0.0 for node in order)
        use_lengths = total > 0

    root.x = 0.0
    for node in reversed(order):                 # parents before children
        for child in node.children:
            child.x = node.x + ((child.length or 0.0) if use_lengths else 1.0)

    tips = leaves(root)
    for index, tip in enumerate(tips):
        tip.y = float(index)
    for node in order:
        if node.children:
            node.y = sum(child.y for child in node.children) / float(len(node.children))
    return tips, use_lengths


# ---------------------------------------------------------------------------
# Loading the tool outputs
# ---------------------------------------------------------------------------


def load(paths):
    loaded = {}
    for path in paths:
        name = os.path.basename(path)
        for suffix, tool in SUFFIXES.items():
            if name.endswith(suffix):
                with open(path) as handle:
                    loaded[tool] = json.load(handle)
                break
    return loaded


def fmt(value, places=4):
    return "NA" if value is None else "{:.{}f}".format(value, places)


def group_order(loaded):
    order = []
    for tool in ("report", "monophyly", "rosenberg", "gsi", "snp_separation"):
        for row in (loaded.get(tool) or {}).get("results", []):
            if row["group"] not in order:
                order.append(row["group"])
    return order


def read_groups(path):
    """tip label -> group. The groups file is the authoritative mapping; the
    JSON outputs only ever name the clade-breaking tips, never every tip."""
    mapping = {}
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split(",")
            if len(fields) < 2:
                continue
            tip, group = fields[0].strip(), fields[1].strip()
            if number == 1 and tip.lower() in ("tip", "tip_label", "taxon", "sample", "label"):
                continue
            mapping[tip] = group
    return mapping


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


def draw_tree(ax, root, tip_group, groups):
    """Rectangular phylogram, tips carrying their group's colour and marker."""
    tips, used_lengths = layout(root)
    coloured = len(groups) <= MAX_COLOURED_GROUPS
    colour_of = {g: (GROUP_COLOURS[i] if coloured else INK_SOFT) for i, g in enumerate(groups)}
    marker_of = {g: GROUP_MARKERS[i % len(GROUP_MARKERS)] for i, g in enumerate(groups)}

    for node in postorder(root):
        if node.children:
            ys = [child.y for child in node.children]
            ax.plot([node.x, node.x], [min(ys), max(ys)], color=INK_SOFT, linewidth=0.8, zorder=1)
            for child in node.children:
                ax.plot([node.x, child.x], [child.y, child.y], color=INK_SOFT, linewidth=0.8, zorder=1)

    span = max((tip.x for tip in tips), default=1.0) or 1.0
    offset = span * 0.012
    n = len(tips)
    label_size = 7.0 if n <= 40 else (5.5 if n <= 80 else 4.0)
    show_labels = n <= 120

    for tip in tips:
        group = tip_group.get(tip.label)
        colour = colour_of.get(group, INK_FAINT)
        ax.plot([tip.x], [tip.y], marker=marker_of.get(group, "o"), markersize=3.2,
                color=colour, markeredgecolor=SURFACE, markeredgewidth=0.4, zorder=3)
        if show_labels:
            ax.text(tip.x + offset, tip.y, tip.label, fontsize=label_size,
                    va="center", ha="left", color=INK, zorder=3)

    # Scale bar below the tips, so branch lengths mean something.
    bottom = n
    if used_lengths:
        bar = span / 5.0
        magnitude = 10 ** int(round(math.log10(bar))) if bar > 0 else 1
        bar = max(magnitude, bar - bar % magnitude) if magnitude else bar
        ax.plot([0, bar], [n + 1.2, n + 1.2], color=INK_SOFT, linewidth=1.2)
        ax.text(bar / 2.0, n + 2.2, "{:g} substitutions/site".format(bar),
                fontsize=6, ha="center", va="top", color=INK_SOFT)
        bottom = n + 3

    ax.set_xlim(-span * 0.02, span * (1.30 if show_labels else 1.02))
    ax.set_ylim(bottom, -1)          # inverted once: tip 0 at the top
    ax.axis("off")

    handles = []
    for group in groups:
        handles.append(plt.Line2D([], [], linestyle="none", marker=marker_of[group],
                                  markersize=5, color=colour_of[group], label=group))
    if handles:
        ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=7,
                  handletextpad=0.4, borderpad=0.2)
    return coloured


def draw_table(ax, headers, rows, title, note=None, widths=None):
    """A plain text table. Beyond a handful of columns a table beats a chart.

    Everything is placed in axes coordinates (transAxes) with the limits pinned:
    ax.text defaults to DATA coordinates, and a single plotted line is enough to
    autoscale them, which throws the rows across the page.
    """
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=6)
    if not rows:
        ax.text(0, 0.5, "no result", fontsize=8, color=INK_FAINT, va="center",
                transform=ax.transAxes)
        return

    columns = len(headers)
    widths = widths or [1.0 / columns] * columns
    xs, running = [], 0.0
    for width in widths:
        xs.append(running)
        running += width

    show_head = any(str(head).strip() for head in headers)
    top = 0.88
    step = min(0.17, 0.74 / max(len(rows), 1))
    if show_head:
        for x, head in zip(xs, headers):
            ax.text(x, top, head, fontsize=7.5, color=INK_SOFT, va="top", ha="left",
                    weight="bold", transform=ax.transAxes)
        ax.axhline(top - step * 0.55, color=RULE, linewidth=0.6)
    else:
        top += step * 0.8
    for index, row in enumerate(rows):
        y = top - step * (index + 1)
        for x, cell in zip(xs, row):
            ax.text(x, y, str(cell), fontsize=7.5, color=INK, va="top", ha="left",
                    transform=ax.transAxes)
    if note:
        y = top - step * (len(rows) + 1) - 0.02
        for line in wrap(note, 118):
            ax.text(0, y, line, fontsize=6.5, color=INK_FAINT, va="top", ha="left",
                    transform=ax.transAxes)
            y -= 0.075


def draw_snp_ranges(ax, rows):
    """Within- vs between-group SNP distance ranges, one row per group.

    Two series, so the legend plus a direct min-max label on every bar carry
    identity; the visible gap between the two bars is the statistic the reader
    wants.

    Within- and between-group distances are the same measure but routinely
    differ by orders of magnitude -- a handful of SNPs inside an outbreak
    against thousands between lineages. On a linear axis the within bar
    collapses to nothing, so the axis goes logarithmic once the range is wide
    enough to warrant it. Zeros (identical genomes) cannot go on a log axis, so
    those fall back to linear.
    """
    usable = [row for row in rows if row["within"]["n"] and row["between"]["n"]]
    if not usable:
        ax.axis("off")
        ax.text(0, 0.5, "no SNP distances available", fontsize=8, color=INK_FAINT, va="center")
        return

    values = [v for row in usable for key in ("within", "between")
              for v in (row[key]["min"], row[key]["max"])]
    smallest = min(values)
    largest = max(values)
    log_scale = smallest > 0 and largest / smallest > 30
    if log_scale:
        ax.set_xscale("log")

    for index, row in enumerate(usable):
        y = len(usable) - index - 1
        for offset, key, colour in ((0.15, "within", WITHIN_COLOUR),
                                    (-0.15, "between", BETWEEN_COLOUR)):
            stats = row[key]
            low, high = stats["min"], stats["max"]
            # A zero-width range still needs to be visible as a mark.
            if low == high:
                ax.plot([low], [y + offset], marker="o", markersize=5, color=colour,
                        markeredgecolor=SURFACE, markeredgewidth=0.8, zorder=3)
            else:
                ax.plot([low, high], [y + offset] * 2, color=colour, linewidth=4,
                        solid_capstyle="round", zorder=2)
                ax.plot([stats["mean"]], [y + offset], marker="|", markersize=8,
                        color=SURFACE, markeredgewidth=1.2, zorder=3)
            label = "{:g}".format(low) if low == high else "{:g}-{:g}".format(low, high)
            gap = (high * 1.12) if log_scale else (high + (largest - smallest) * 0.015)
            ax.text(gap, y + offset, label, fontsize=6, color=INK_SOFT,
                    ha="left", va="center")

    ax.set_yticks(range(len(usable)))
    ax.set_yticklabels([row["group"] for row in reversed(usable)], fontsize=7.5, color=INK)
    ax.set_xlabel("SNP distance" + ("  (log scale)" if log_scale else ""),
                  fontsize=7.5, color=INK_SOFT)
    ax.tick_params(axis="x", labelsize=7, colors=INK_SOFT, length=2)
    ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.grid(axis="x", color=RULE, linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_ylim(-0.6, len(usable) - 0.4)
    if log_scale:
        ax.set_xlim(smallest * 0.55, largest * 2.6)
    else:
        span = (largest - smallest) or 1.0
        ax.set_xlim(max(0, smallest - span * 0.06), largest + span * 0.22)
    # Legend above the plot, where it cannot land on a bar.
    ax.legend(handles=[
        plt.Line2D([], [], color=WITHIN_COLOUR, linewidth=4, label="within group"),
        plt.Line2D([], [], color=BETWEEN_COLOUR, linewidth=4, label="between groups"),
    ], loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, frameon=False, fontsize=7)


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


def page_header(fig, title, subtitle=None):
    fig.text(0.06, 0.965, title, fontsize=13, color=INK, weight="bold", va="top")
    if subtitle:
        fig.text(0.06, 0.935, subtitle, fontsize=8, color=INK_SOFT, va="top")
    fig.patch.set_facecolor(SURFACE)


def wrap(text, width):
    words, line, lines = text.split(), "", []
    for word in words:
        if line and len(line) + 1 + len(word) > width:
            lines.append(line)
            line = word
        else:
            line = word if not line else line + " " + word
    if line:
        lines.append(line)
    return lines


def build(loaded, root, tip_group, out_path):
    groups = group_order(loaded)
    report = loaded.get("report") or {}
    monophyly_rows = (loaded.get("monophyly") or {}).get("results", [])

    # Any group only the groups file knows about still deserves a legend entry.
    for group in sorted(set(tip_group.values())):
        if group not in groups:
            groups.append(group)

    with PdfPages(out_path) as pdf:
        # ---- page 1: the tree and the conclusion --------------------------
        fig = plt.figure(figsize=(8.27, 11.69))          # A4 portrait
        page_header(fig, "Cladebreaker analyze report",
                    "Tools: {}   ·   rooting: {}   ·   alpha = {}".format(
                        ", ".join(report.get("tools", sorted(loaded))) or "none",
                        (loaded.get("monophyly") or loaded.get("gsi") or {}).get("rooting", "?"),
                        report.get("alpha", "?")))

        conclusion = report.get("conclusion") or []
        text_height = 0.035 + 0.018 * len(conclusion)
        if root is not None:
            fig.text(0.06, 0.905, "Rooted tree the statistics were computed on",
                     fontsize=9.5, color=INK, va="top")
            ax = fig.add_axes([0.06, 0.06 + text_height, 0.88, 0.825 - text_height])
            draw_tree(ax, root, tip_group, groups)
        y = 0.055 + 0.018 * (len(conclusion) - 1)
        fig.text(0.06, y + 0.028, "CONCLUSION", fontsize=9, color=INK_SOFT, weight="bold")
        for line in conclusion:
            fig.text(0.06, y, line.strip(), fontsize=8.5, color=INK, va="top")
            y -= 0.018
        pdf.savefig(fig)
        plt.close(fig)

        # ---- page 2: per-test results -------------------------------------
        fig = plt.figure(figsize=(8.27, 11.69))
        page_header(fig, "Results by test")
        slots = [0.74, 0.50, 0.26, 0.03]

        rows = [[r["group"], r["n_tips"], "yes" if r["monophyletic"] else "no",
                 r.get("n_clusters", "-"), r.get("n_clade_breaking", 0),
                 fmt(r.get("monophyly_fraction"), 2)] for r in monophyly_rows]
        ax = fig.add_axes([0.06, slots[0], 0.88, 0.17])
        draw_table(ax, ["group", "tips", "clade", "clusters", "breaking", "support"], rows,
                   "1 · Monophyly",
                   "Clusters is how many separate blocks a broken group falls into; 1 means monophyletic.",
                   widths=[0.30, 0.10, 0.11, 0.14, 0.15, 0.20])

        ros = (loaded.get("rosenberg") or {})
        rows = [[r["group"], r["n_tips"], fmt(r.get("p_a"), 6), fmt(r.get("log10_p_a"), 2),
                 fmt(r.get("p_a_holm"), 6),
                 "-" if r.get("min_b_for_alpha") is None else r["min_b_for_alpha"]]
                for r in ros.get("results", [])]
        summary = ros.get("summary") or {}
        note = "P_A is withheld where a group is not a clade: the test is conditional on observing monophyly."
        joint = (summary.get("joint") or {}).get("p")
        if joint is not None:
            note = "Joint monophyly of all groups: p = {} (10^{}).  ".format(
                fmt(joint, 8), fmt(summary["joint"]["log10_p"], 1)) + note
        ax = fig.add_axes([0.06, slots[1], 0.88, 0.17])
        draw_table(ax, ["group", "tips", "P_A", "log10", "P_A (Holm)", "min b"], rows,
                   "2 · Rosenberg — could the monophyly be chance?", note,
                   widths=[0.30, 0.10, 0.16, 0.12, 0.18, 0.14])

        gsi_rows = (loaded.get("gsi") or {}).get("results", [])
        rows = [[r["group"], r["n_tips"], fmt(r.get("gsi"), 4), fmt(r.get("min_gs"), 4),
                 fmt(r.get("p_value"), 6)] for r in gsi_rows]
        ax = fig.add_axes([0.06, slots[2], 0.88, 0.17])
        draw_table(ax, ["group", "tips", "gsi", "min(gs)", "p-value"], rows,
                   "3 · Genealogical sorting index — how far towards exclusivity?",
                   "0 is ancestry fully mixed with the other groups; 1 is monophyly.",
                   widths=[0.30, 0.12, 0.18, 0.18, 0.22])

        sm = loaded.get("slatkin_maddison") or {}
        if sm:
            rows = [["migration events (s)", fmt(sm.get("s_observed"), 1)],
                    ["expected by chance", "{}{}".format(
                        fmt(sm.get("null_mean"), 1),
                        "" if sm.get("null_sd") is None else " +/- " + fmt(sm["null_sd"], 1))],
                    ["minimum possible", sm.get("s_minimum_possible", "-")],
                    ["p-value", fmt(sm.get("p_value"), 6)]]
        else:
            rows = []
        ax = fig.add_axes([0.06, slots[3], 0.88, 0.17])
        draw_table(ax, ["", ""], rows,
                   "4 · Slatkin-Maddison — is the residual mixing less than random?",
                   "Fewer changes than chance demands means structure even without monophyly. "
                   "Independent of rooting.", widths=[0.40, 0.60])
        pdf.savefig(fig)
        plt.close(fig)

        # ---- page 3: SNP separation and caveats ---------------------------
        snp = (loaded.get("snp_separation") or {})
        fig = plt.figure(figsize=(8.27, 11.69))
        page_header(fig, "SNP separation and caveats" if snp.get("results") else "Caveats")

        snp_rows = snp.get("results", [])
        if snp_rows:
            ax = fig.add_axes([0.14, 0.62, 0.80, 0.26])
            draw_snp_ranges(ax, snp_rows)
            rows = [[r["group"], r["n_tips"], fmt(r["within"]["max"], 0), fmt(r["between"]["min"], 0),
                     fmt(r.get("gap"), 2), fmt(r.get("ratio_of_means"), 2), fmt(r.get("p_value"), 4)]
                    for r in snp_rows]
            ax = fig.add_axes([0.06, 0.40, 0.88, 0.17])
            draw_table(ax, ["group", "tips", "max within", "min between", "gap", "ratio", "p-value"], rows,
                       "5 · SNP separation",
                       "Gap above 1 means every pair inside the group is closer than any pair leaving it.",
                       widths=[0.24, 0.09, 0.16, 0.17, 0.11, 0.11, 0.12])
            caveat_top = 0.34
        else:
            fig.text(0.06, 0.90, "5 · SNP separation", fontsize=9.5, color=INK, va="top")
            fig.text(0.06, 0.875, "Not run: no alignment or distance matrix was supplied.",
                     fontsize=8, color=INK_FAINT, va="top")
            caveat_top = 0.83

        caveats = [payload["caveat"] for key, payload in sorted(loaded.items())
                   if isinstance(payload, dict) and payload.get("caveat")]
        y = caveat_top
        if caveats:
            fig.text(0.06, y, "CAVEATS", fontsize=9, color=INK_SOFT, weight="bold")
            y -= 0.022
            for caveat in caveats:
                for line in wrap(caveat, 108):
                    fig.text(0.06, y, line, fontsize=7, color=INK_SOFT, va="top")
                    y -= 0.014
                y -= 0.010
        pdf.savefig(fig)
        plt.close(fig)


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Render the cladebreaker analyze results as a PDF.",
        epilog="Example usage: analyze_pdf.py --tree rooted.nwk --out-prefix cladebreaker *.json",
    )
    parser.add_argument("inputs", nargs="+", help="JSON outputs from the analyze tools.")
    parser.add_argument("--tree", help="The rooted Newick tree to draw (monophyly.py writes it).")
    parser.add_argument("--groups", help="Groups file; supplies the tip-to-group map for the figure.")
    parser.add_argument("--out-prefix", default="cladebreaker", help="Prefix for the PDF.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    loaded = load(args.inputs)
    if not loaded:
        sys.exit("ERROR: none of the given files are recognised analyze tool outputs.")

    root = None
    if args.tree:
        with open(args.tree) as handle:
            root = parse_newick(handle.read())

    tip_group = read_groups(args.groups) if args.groups else {}

    out_path = "{}.analyze_report.pdf".format(args.out_prefix)
    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)
    build(loaded, root, tip_group, out_path)
    print("Wrote {}".format(out_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
