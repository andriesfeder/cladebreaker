#!/usr/bin/env python

"""
Combines the analyze tools into one verdict, following the decision path.

Each tool answers a different question and none of them answers it alone:

    monophyly          is the group a clade on the rooted tree?
    rosenberg          if it is, could that be chance? (sample size adequacy)
    gsi                if it is not, how far along the way to exclusivity is it?
    slatkin_maddison   if it is not, is the residual mixing still non-random?
    snp_separation     how far apart are the groups in SNP distance?

The path taken here is:

    monophyletic  -> Rosenberg is the formal test, reported with SNP separation.
                     A non-significant result means the sample is too small, not
                     that the group is undistinguished, so the report says how
                     many more outside genomes would settle it.

    not monophyletic
                  -> Slatkin-Maddison asks whether the residual mixing is still
                     significantly less than random, with the gsi giving the
                     magnitude per group and monophyly naming the clade-breaking
                     tips. This is a weaker but real positive result, NOT a
                     failure -- which is the whole reason for the fallback.

Every tool is run regardless and all outputs are kept; the decision only chooses
which result leads the report. Inputs are the JSON files the tools write, matched
by their filename suffix, and any subset may be supplied.

Example usage:
    analyze_report.py --out-prefix cladebreaker *.json
"""

import argparse
import json
import os
import sys

SUFFIXES = {
    ".monophyly.json": "monophyly",
    ".rosenberg.json": "rosenberg",
    ".slatkin_maddison.json": "slatkin_maddison",
    ".gsi.json": "gsi",
    ".snp_separation.json": "snp_separation",
}


# This script's own output, which a glob of *.json on a re-run will pick up.
OWN_OUTPUT = ".analyze_report.json"


def load(paths):
    """Match each JSON to the tool that wrote it, by filename suffix."""
    loaded = {}
    for path in paths:
        name = os.path.basename(path)
        if name.endswith(OWN_OUTPUT):
            continue          # re-running over a glob is normal, not worth a warning
        for suffix, tool in SUFFIXES.items():
            if name.endswith(suffix):
                with open(path) as handle:
                    loaded[tool] = json.load(handle)
                break
        else:
            sys.stderr.write("WARNING: ignoring {}, which is not a recognised tool output\n".format(name))
    return loaded


def fmt(value, places=4):
    return "NA" if value is None else "{:.{}f}".format(value, places)


def by_group(payload, key="results"):
    if not payload:
        return {}
    return {row["group"]: row for row in payload.get(key, [])}


def verdict(group, mono, ros, gsi, alpha):
    """The one-line conclusion for a single group."""
    if mono is None:
        return "no monophyly result", None

    if mono["monophyletic"]:
        if ros is None or ros.get("p_a") is None:
            return "monophyletic (chance not tested)", None
        if ros["p_a"] <= alpha:
            return "distinct: monophyletic, and too large to be chance", True
        needed = ros.get("min_b_for_alpha")
        detail = "" if not needed else "; {} outside genomes would settle it".format(needed)
        return ("monophyletic, but the sample is too small to exclude chance" + detail), False

    broken = "{} clade-breaking tip(s), {} cluster(s)".format(
        mono.get("n_clade_breaking", 0), mono.get("n_clusters", "?"))
    if gsi and gsi.get("gsi") is not None and gsi.get("p_value") is not None and gsi["p_value"] <= alpha:
        return ("partially sorted: gsi={} (P={}), {}".format(
            fmt(gsi["gsi"], 3), fmt(gsi["p_value"], 4), broken)), True
    if gsi and gsi.get("gsi") is not None:
        return ("not distinguishable from random sorting: gsi={} (P={}), {}".format(
            fmt(gsi["gsi"], 3), fmt(gsi.get("p_value"), 4), broken)), False
    return "not monophyletic: " + broken, None


def headline(loaded, groups, alpha):
    """The leading statement: Rosenberg if monophyly held, otherwise Slatkin-Maddison."""
    mono = by_group(loaded.get("monophyly"))
    ros_summary = (loaded.get("rosenberg") or {}).get("summary", {})
    sm = loaded.get("slatkin_maddison")

    monophyletic = [g for g in groups if mono.get(g, {}).get("monophyletic")]
    all_mono = bool(groups) and bool(mono) and len(monophyletic) == len(groups)

    lines = []
    if not mono:
        # Without a monophyly result there is no branch to take; saying "no group
        # is monophyletic" here would assert something never tested.
        lines.append("Monophyly was not tested, so no decision path applies.")
        if sm and sm.get("p_value") is not None:
            lines.append("  Slatkin-Maddison: {} migration event(s) against {} expected by chance, P={}.".format(
                fmt(sm["s_observed"], 1), fmt(sm.get("null_mean"), 1), fmt(sm["p_value"], 6)))
        return lines, False

    if all_mono:
        lines.append("All {} group(s) are monophyletic.".format(len(groups)))
        joint = ros_summary.get("joint") or {}
        recip = ros_summary.get("reciprocal") or {}
        if joint.get("p") is not None:
            lines.append(
                "  Joint monophyly of all groups under random branching: p={} (10^{}).".format(
                    fmt(joint["p"], 8), fmt(joint["log10_p"], 1)))
        elif recip.get("p_ab") is not None:
            lines.append(
                "  Reciprocal monophyly of {} and {}: P_AB={} (10^{}).".format(
                    recip["groups"][0], recip["groups"][1],
                    fmt(recip["p_ab"], 8), fmt(recip["log10_p_ab"], 1)))
        elif joint.get("skipped"):
            lines.append("  Joint probability not computed: {}.".format(joint["skipped"]))
        if not ros_summary.get("exhaustive", True):
            lines.append("  Groups do not cover every tip, so only per-group P_A applies.")
    else:
        if monophyletic:
            lines.append("{} of {} group(s) are monophyletic: {}.".format(
                len(monophyletic), len(groups), ", ".join(monophyletic)))
        else:
            lines.append("No group is monophyletic, so Rosenberg's test does not apply.")
        if sm and sm.get("p_value") is not None:
            significant = sm["p_value"] <= alpha
            lines.append(
                "  Slatkin-Maddison: {} migration event(s) against {} expected by chance, P={}.".format(
                    fmt(sm["s_observed"], 1), fmt(sm.get("null_mean"), 1), fmt(sm["p_value"], 6)))
            lines.append(
                "  The groups are structured despite not being clades."
                if significant else
                "  The mixing is not distinguishable from random.")
        elif sm:
            lines.append("  Slatkin-Maddison: {} migration event(s) (no permutation test run).".format(
                fmt(sm["s_observed"], 1)))
    return lines, all_mono


COLUMNS = ["group", "n_tips", "monophyletic", "n_clusters", "gsi", "gsi_p",
           "p_a", "p_a_holm", "snp_gap", "snp_ratio", "snp_p", "verdict"]


def build_rows(loaded, groups, alpha):
    mono = by_group(loaded.get("monophyly"))
    ros = by_group(loaded.get("rosenberg"))
    gsi = by_group(loaded.get("gsi"))
    snp = by_group(loaded.get("snp_separation"))

    rows = []
    for group in groups:
        m, r, g, s = mono.get(group), ros.get(group), gsi.get(group), snp.get(group)
        text, significant = verdict(group, m, r, g, alpha)
        rows.append({
            "group": group,
            "n_tips": (m or r or g or s or {}).get("n_tips"),
            "monophyletic": None if m is None else m["monophyletic"],
            "n_clusters": (m or {}).get("n_clusters"),
            "gsi": (g or {}).get("gsi"),
            "gsi_p": (g or {}).get("p_value"),
            "p_a": (r or {}).get("p_a"),
            "p_a_holm": (r or {}).get("p_a_holm"),
            "snp_gap": (s or {}).get("gap"),
            "snp_ratio": (s or {}).get("ratio_of_means"),
            "snp_p": (s or {}).get("p_value"),
            "verdict": text,
            "significant": significant,
        })
    return rows


def write_table(path, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join([
                row["group"],
                "NA" if row["n_tips"] is None else str(row["n_tips"]),
                "NA" if row["monophyletic"] is None else ("yes" if row["monophyletic"] else "no"),
                "NA" if row["n_clusters"] is None else str(row["n_clusters"]),
                fmt(row["gsi"], 4), fmt(row["gsi_p"], 6),
                fmt(row["p_a"], 8), fmt(row["p_a_holm"], 8),
                fmt(row["snp_gap"], 3), fmt(row["snp_ratio"], 3), fmt(row["snp_p"], 6),
                row["verdict"],
            ]) + "\n")


def render(loaded, groups, rows, lines, alpha):
    out = []
    out.append("=" * 72)
    out.append(" CLADEBREAKER ANALYZE REPORT")
    out.append("=" * 72)

    tools = ", ".join(sorted(loaded)) or "none"
    out.append("Tools run    : {}".format(tools))
    source = loaded.get("monophyly") or loaded.get("gsi") or loaded.get("rosenberg") or {}
    if source.get("tree_file"):
        out.append("Tree         : {}".format(source["tree_file"]))
    if source.get("rooting"):
        note = "" if source["rooting"] != "as-is" else "   <- trusted as written; check the tree really is rooted"
        out.append("Rooting      : {}{}".format(source["rooting"], note))
    if source.get("n_trees", 1) > 1:
        out.append("Trees        : {} (ensemble)".format(source["n_trees"]))
    out.append("Significance : alpha = {}".format(alpha))
    out.append("")

    out.append("CONCLUSION")
    out.extend("  " + line for line in lines)
    out.append("")

    out.append("BY GROUP")
    for row in rows:
        out.append("  {:<20} n={:<5} {}".format(row["group"], row["n_tips"] or "?", row["verdict"]))
        if row["snp_gap"] is not None or row["snp_ratio"] is not None:
            out.append("  {:<20}       SNP separation: gap={} ratio={} (P={})".format(
                "", fmt(row["snp_gap"], 2), fmt(row["snp_ratio"], 2), fmt(row["snp_p"], 4)))
    out.append("")

    caveats = [payload["caveat"] for key, payload in sorted(loaded.items()) if payload.get("caveat")]
    if caveats:
        out.append("CAVEATS")
        for caveat in caveats:
            for chunk in wrap(caveat, 68):
                out.append("  " + chunk)
            out.append("")
    out.append("=" * 72)
    return "\n".join(out) + "\n"


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


def write_multiqc(path, rows, lines):
    with open(path, "w") as handle:
        handle.write("# id: 'analyze_report'\n")
        handle.write("# section_name: 'Cladebreaker analyze summary'\n")
        handle.write("# description: '{}'\n".format(" ".join(lines).replace("'", "")))
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    namespace: 'analyze_report'\n")
        handle.write("Group\tTips\tMonophyletic\tgsi\tP_A\tSNP gap\tVerdict\n")
        for row in rows:
            handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["group"], row["n_tips"] if row["n_tips"] is not None else "NA",
                "NA" if row["monophyletic"] is None else ("yes" if row["monophyletic"] else "no"),
                fmt(row["gsi"], 3), fmt(row["p_a"], 6), fmt(row["snp_gap"], 2), row["verdict"]))


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Combine the cladebreaker analyze tool outputs into a single verdict.",
        epilog="Example usage: analyze_report.py --out-prefix cladebreaker *.json",
    )
    parser.add_argument("inputs", nargs="+", help="JSON outputs from the analyze tools.")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level (default: 0.05).")
    parser.add_argument("--out-prefix", default="cladebreaker", help="Prefix for output files.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if not 0.0 < args.alpha < 1.0:
        sys.exit("ERROR: --alpha must be between 0 and 1.")

    loaded = load(args.inputs)
    if not loaded:
        sys.exit("ERROR: none of the given files are recognised analyze tool outputs.")

    # Group order follows whichever tool ran, preferring the tree-based ones.
    groups = []
    for tool in ("monophyly", "rosenberg", "gsi", "snp_separation"):
        for row in (loaded.get(tool) or {}).get("results", []):
            if row["group"] not in groups:
                groups.append(row["group"])

    rows = build_rows(loaded, groups, args.alpha)
    lines, all_mono = headline(loaded, groups, args.alpha)

    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)

    text = render(loaded, groups, rows, lines, args.alpha)
    with open("{}.analyze_report.txt".format(args.out_prefix), "w") as handle:
        handle.write(text)
    write_table("{}.analyze_report.tsv".format(args.out_prefix), rows)
    write_multiqc("{}.analyze_report_mqc.tsv".format(args.out_prefix), rows, lines)
    with open("{}.analyze_report.json".format(args.out_prefix), "w") as handle:
        json.dump({
            "alpha": args.alpha, "tools": sorted(loaded), "all_monophyletic": all_mono,
            "conclusion": lines, "results": rows,
        }, handle, indent=2)
        handle.write("\n")

    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
