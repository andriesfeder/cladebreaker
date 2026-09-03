#!/usr/bin/env python

"""
Genealogical sorting index (gsi) for labeled groups on a rooted tree.

Implements the statistic of Cummings, Neel & Shaw (2008) Evolution 62:2411-2422,
doi:10.1111/j.1558-5646.2008.00442.x

    gs      = n / sum_{u in U} (d_u - 2)        (eq. 1)  n = group size - 1
    max(gs) = 1                                 (eq. 2)
    min(gs) = n / sum_{i in I} (d_i - 2)        (eq. 3)
    gsi     = (gs - min(gs)) / (1 - min(gs))    (eq. 4)
    gsi_T   = sum_t gsi_t * P_t                 (eq. 5)

U is the set of internal nodes of the smallest subtree containing every member of
the group (its MRCA included); I is every internal node in the tree; d is the node
degree, i.e. the number of branch ends at the node.

NOTE ON DEGREE: the root of a rooted genealogy is treated as carrying a stem
branch, so every internal node has degree = children + 1, the root included.
Without the root stem the worked example in Fig. 1 of the paper does not
reproduce (min(gs) there is 3/7, which requires the root to count as degree 3).

Significance is assessed by permuting group labels across the tips while holding
the topology fixed, per Maddison & Slatkin (1991):

    P = (1 + #{gsi_permuted >= gsi_observed}) / (permutations + 1)

Example usage:
    gsi.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix gsi
"""

import argparse
import json
import os
import random
import sys

import dendropy


# ---------------------------------------------------------------------------
# Tree indexing
#
# The dendropy tree is flattened into plain arrays once, then all of the
# permutation work runs over those arrays. Internal nodes are stored in
# post-order, so a node always appears after its children.
# ---------------------------------------------------------------------------


class TreeIndex(object):
    """Flattened rooted topology: internal nodes in post-order plus degrees."""

    def __init__(self, tree, canonical_labels):
        label_to_canonical = {label: i for i, label in enumerate(canonical_labels)}

        self.children = []       # per internal node: list of (is_leaf, index)
        self.degree = []         # per internal node: branch ends at the node
        self.leaf_canonical = []  # per local leaf: index into canonical_labels

        leaf_index = {}
        internal_index = {}

        for node in tree.postorder_node_iter():
            if node.is_leaf():
                label = node.taxon.label if node.taxon is not None else node.label
                if label not in label_to_canonical:
                    raise KeyError(label)
                leaf_index[id(node)] = len(self.leaf_canonical)
                self.leaf_canonical.append(label_to_canonical[label])
                continue

            kids = []
            for child in node.child_node_iter():
                if child.is_leaf():
                    kids.append((True, leaf_index[id(child)]))
                else:
                    kids.append((False, internal_index[id(child)]))

            internal_index[id(node)] = len(self.children)
            self.children.append(kids)
            # + 1 for the parent branch; the root carries a stem branch too.
            self.degree.append(len(kids) + 1)

        self.n_internal = len(self.children)
        self.n_leaves = len(self.leaf_canonical)
        # Denominator of eq. 3: every internal node in the tree.
        self.total_degree_sum = sum(d - 2 for d in self.degree)


def index_trees(trees, canonical_labels):
    return [TreeIndex(tree, canonical_labels) for tree in trees]


# ---------------------------------------------------------------------------
# The statistic
# ---------------------------------------------------------------------------


def gs_values(tree_index, assignment, group_sizes):
    """gs for every group on one tree.

    assignment maps each canonical tip index to a group index, or -1 for tips
    that carry no group label. Returns a list of gs values (None where the
    statistic is undefined).
    """
    n_groups = len(group_sizes)
    counts = [[0] * n_groups for _ in range(tree_index.n_internal)]

    for i, kids in enumerate(tree_index.children):
        row = counts[i]
        for is_leaf, j in kids:
            if is_leaf:
                group = assignment[tree_index.leaf_canonical[j]]
                if group >= 0:
                    row[group] += 1
            else:
                child_row = counts[j]
                for g in range(n_groups):
                    row[g] += child_row[g]

    results = []
    for g, size in enumerate(group_sizes):
        if size < 2:
            # n = 0; a single tip has no ancestry to share, gsi is undefined.
            results.append(None)
            continue

        # Post-order puts descendants first, so the first node holding the whole
        # group is the MRCA; every later full-count node is one of its ancestors.
        mrca = None
        for i in range(tree_index.n_internal):
            if counts[i][g] == size:
                mrca = i
                break
        if mrca is None:
            # Group members are not all in this tree.
            results.append(None)
            continue

        # U: nodes at or below the MRCA that subtend at least one group member.
        denominator = 0
        for i in range(mrca + 1):
            if counts[i][g] > 0:
                denominator += tree_index.degree[i] - 2

        results.append((size - 1) / denominator if denominator > 0 else None)

    return results


def min_gs(tree_index, size):
    """Equation 3: the smallest gs attainable by a group of this size."""
    if size < 2 or tree_index.total_degree_sum <= 0:
        return None
    return (size - 1) / tree_index.total_degree_sum


def normalize(gs, minimum):
    """Equation 4."""
    if gs is None or minimum is None or minimum >= 1.0:
        return None
    return (gs - minimum) / (1.0 - minimum)


def ensemble_gsi(per_tree_gsi):
    """Equation 5, with each sampled topology weighted by its frequency.

    Trees are taken as drawn from the ensemble, so P_t is 1/T for each.
    """
    usable = [value for value in per_tree_gsi if value is not None]
    if not usable:
        return None
    return sum(usable) / len(usable)


def analyse(trees, mapping, groups, permutations=0, seed=42):
    """Run the whole statistic over an ensemble of trees.

    Returns one result row per group, in the order of `groups`, plus the
    per-tree gsi values.
    """
    canonical_labels = sorted(set(tip_labels(trees[0])))
    group_index = {group: g for g, group in enumerate(groups)}
    assignment = [group_index.get(mapping.get(label, None), -1) for label in canonical_labels]
    group_sizes = [sum(1 for g in assignment if g == index) for index in range(len(groups))]

    tree_indices = index_trees(trees, canonical_labels)
    single = len(tree_indices) == 1

    gsi_values, per_tree = observed_gsi(tree_indices, assignment, group_sizes)
    p_values, n_permutations = permutation_test(
        tree_indices, assignment, group_sizes, gsi_values, permutations, seed
    )
    single_gs = gs_values(tree_indices[0], assignment, group_sizes) if single else None

    rows = []
    for g, group in enumerate(groups):
        rows.append(
            {
                "group": group,
                "n_tips": group_sizes[g],
                "n_trees": len(tree_indices),
                "gs": single_gs[g] if single else None,
                "min_gs": min_gs(tree_indices[0], group_sizes[g]) if single else None,
                "gsi": gsi_values[g],
                "p_value": p_values[g],
                "n_permutations": n_permutations,
            }
        )
    return rows, per_tree


def observed_gsi(tree_indices, assignment, group_sizes):
    """Per-group gsi across the ensemble, plus the per-tree values."""
    per_tree = [[] for _ in group_sizes]
    for tree_index in tree_indices:
        gs_list = gs_values(tree_index, assignment, group_sizes)
        for g, gs in enumerate(gs_list):
            per_tree[g].append(normalize(gs, min_gs(tree_index, group_sizes[g])))
    return [ensemble_gsi(values) for values in per_tree], per_tree


def permutation_test(tree_indices, assignment, group_sizes, observed, permutations, seed):
    """Shuffle group labels across the labeled tips, holding topologies fixed.

    One shuffle per replicate is applied to every tree in the ensemble, which
    preserves any correlated structure across loci (the paper notes both
    permutation schemes converge for the general case).
    """
    if permutations <= 0:
        return [None] * len(group_sizes), 0

    labeled = [i for i, g in enumerate(assignment) if g >= 0]
    labels = [assignment[i] for i in labeled]

    rng = random.Random(seed)
    shuffled = list(assignment)
    at_least = [0] * len(group_sizes)

    for _ in range(permutations):
        rng.shuffle(labels)
        for position, label in zip(labeled, labels):
            shuffled[position] = label

        permuted, _ = observed_gsi(tree_indices, shuffled, group_sizes)
        for g, value in enumerate(permuted):
            if observed[g] is None or value is None:
                continue
            if value >= observed[g]:
                at_least[g] += 1

    # The observed labeling counts as one of the replicates.
    p_values = [
        None if observed[g] is None else (1 + at_least[g]) / (permutations + 1)
        for g in range(len(group_sizes))
    ]
    return p_values, permutations


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------


def read_groups(path):
    """Read a two-column tip-to-group table (TSV or CSV, optional header)."""
    mapping = {}
    order = []
    with open(path, "r") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split(",")
            fields = [field.strip() for field in fields if field.strip()]
            if len(fields) < 2:
                sys.exit(
                    "ERROR: Please check the groups file -> line {} does not have two "
                    "columns:\n{}".format(number, line)
                )
            tip, group = fields[0], fields[1]
            if number == 1 and tip.lower() in ("tip", "tip_label", "taxon", "sample", "label"):
                continue
            if tip in mapping:
                sys.exit("ERROR: Please check the groups file -> tip '{}' listed twice.".format(tip))
            mapping[tip] = group
            if group not in order:
                order.append(group)

    if not mapping:
        sys.exit("ERROR: Please check the groups file -> no tip/group assignments found.")
    if len(order) < 2:
        sys.exit(
            "ERROR: Please check the groups file -> gsi compares labeled groups, but only "
            "one group ('{}') was found.".format(order[0])
        )
    return mapping, order


def reroot_at_midpoint(tree):
    """Root a tree at the midpoint of its longest leaf-to-leaf path.

    dendropy has reroot_at_midpoint(), but as of 5.0.13 it walks up from only
    one end of the diameter and asserts if it reaches the MRCA with distance
    left over, which floating-point drift makes happen on larger trees (its
    _tree.py, "ugly, ugly, ugly code to find two nodes that span the
    midpoint"). Rerooting is not the interesting part of this script, so we do
    it here rather than crash on somebody's tree.
    """
    # Deepest leaf below each node, and the child that path goes through.
    depth = {}
    descend = {}
    diameter = 0.0
    meeting_node = None

    for node in tree.postorder_node_iter():
        best = []
        for child in node.child_node_iter():
            best.append((depth[id(child)] + (child.edge.length or 0.0), child))
        best.sort(key=lambda pair: pair[0], reverse=True)

        depth[id(node)] = best[0][0] if best else 0.0
        descend[id(node)] = best[0][1] if best else None

        # The longest path meeting at this node runs down its two deepest branches.
        span = (best[0][0] if best else 0.0) + (best[1][0] if len(best) > 1 else 0.0)
        if span > diameter:
            diameter = span
            meeting_node = node

    if meeting_node is None or diameter <= 0.0:
        sys.exit(
            "ERROR: The tree has no branch lengths to midpoint-root on. "
            "Use --root outgroup or --root as-is."
        )

    # The midpoint sits on the deeper of the two branches, this far below the
    # meeting node. depth >= diameter/2 there, so this never goes negative.
    remaining = depth[id(meeting_node)] - diameter / 2.0

    node = meeting_node
    while True:
        child = descend[id(node)]
        if child is None:
            break
        length = child.edge.length or 0.0
        if remaining <= length:
            break
        remaining -= length
        node = child

    child = descend[id(node)]
    if child is None or remaining <= 0.0:
        tree.reseed_at(node, update_bipartitions=False, suppress_unifurcations=True)
        return

    length = child.edge.length or 0.0
    tree.reroot_at_edge(
        child.edge,
        length1=remaining,
        length2=max(length - remaining, 0.0),
        update_bipartitions=False,
        suppress_unifurcations=True,
    )


def read_trees(path, root, outgroup):
    """Load one or more Newick trees and root them as requested."""
    trees = dendropy.TreeList.get(
        path=path,
        schema="newick",
        rooting="force-rooted",
        preserve_underscores=True,
    )
    if not trees:
        sys.exit("ERROR: No trees could be read from {}".format(path))

    for number, tree in enumerate(trees, start=1):
        if root == "midpoint":
            if any(edge.length is None for edge in tree.postorder_edge_iter() if edge.tail_node):
                sys.exit(
                    "ERROR: Tree {} has no branch lengths, so --root midpoint cannot be "
                    "used. Use --root outgroup or --root as-is.".format(number)
                )
            reroot_at_midpoint(tree)
        elif root == "outgroup":
            present = set(tip_labels(tree))
            absent = [tip for tip in outgroup if tip not in present]
            if absent:
                sys.exit(
                    "ERROR: Please check --outgroup -> {} tip(s) are not in tree {}: {}".format(
                        len(absent), number, ", ".join(absent)
                    )
                )
            node = tree.mrca(taxon_labels=outgroup)
            if node is tree.seed_node:
                sys.exit(
                    "ERROR: The outgroup spans the whole of tree {}, so it cannot root it. "
                    "Check that the outgroup is monophyletic.".format(number)
                )
            tree.reroot_at_edge(node.edge, update_bipartitions=False, suppress_unifurcations=True)

    return trees


def tip_labels(tree):
    return [
        node.taxon.label if node.taxon is not None else node.label
        for node in tree.leaf_node_iter()
    ]


def check_tips(trees, mapping, ignore_unlabeled):
    """Reconcile the groups file against the tips present in the trees."""
    first = set(tip_labels(trees[0]))
    for number, tree in enumerate(trees[1:], start=2):
        if set(tip_labels(tree)) != first:
            sys.exit(
                "ERROR: Tree {} does not have the same tips as tree 1. All trees in an "
                "ensemble must share one set of tips.".format(number)
            )

    missing = sorted(tip for tip in mapping if tip not in first)
    if missing:
        sys.exit(
            "ERROR: Please check the groups file -> {} tip(s) are not in the tree: {}".format(
                len(missing), ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
            )
        )

    unlabeled = sorted(tip for tip in first if tip not in mapping)
    if unlabeled and not ignore_unlabeled:
        sys.exit(
            "ERROR: {} tip(s) in the tree have no group in the groups file: {}\n"
            "Assign them a group, or pass --ignore-unlabeled to leave them out of the "
            "statistic.".format(
                len(unlabeled), ", ".join(unlabeled[:10]) + (" ..." if len(unlabeled) > 10 else "")
            )
        )
    return unlabeled


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def format_value(value, places=4):
    return "NA" if value is None else "{:.{}f}".format(value, places)


COLUMNS = [
    "group",
    "n_tips",
    "n_trees",
    "gs",
    "min_gs",
    "gsi",
    "p_value",
    "n_permutations",
]


def write_results(path, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        row["group"],
                        str(row["n_tips"]),
                        str(row["n_trees"]),
                        format_value(row["gs"]),
                        format_value(row["min_gs"]),
                        format_value(row["gsi"]),
                        format_value(row["p_value"]),
                        str(row["n_permutations"]),
                    ]
                )
                + "\n"
            )


def write_per_tree(path, groups, per_tree):
    with open(path, "w") as handle:
        handle.write("tree\tgroup\tgsi\n")
        for g, group in enumerate(groups):
            for number, value in enumerate(per_tree[g], start=1):
                handle.write("{}\t{}\t{}\n".format(number, group, format_value(value)))


def write_multiqc(path, rows, n_trees):
    """Custom-content table for the existing MultiQC report."""
    statistic = "gsi_T" if n_trees > 1 else "gsi"
    with open(path, "w") as handle:
        handle.write("# id: 'gsi'\n")
        handle.write("# section_name: 'Genealogical sorting index'\n")
        handle.write(
            "# description: 'Degree of exclusive ancestry of each labeled group on the "
            "phylogeny, after Cummings, Neel & Shaw (2008). 1 is monophyly; the p-value "
            "comes from permuting group labels across the tips.'\n"
        )
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    namespace: 'gsi'\n")
        handle.write("Group\tTips\t{}\tp-value\n".format(statistic))
        for row in rows:
            handle.write(
                "{}\t{}\t{}\t{}\n".format(
                    row["group"], row["n_tips"], format_value(row["gsi"]), format_value(row["p_value"])
                )
            )


# ---------------------------------------------------------------------------


def parse_args(args=None):
    description = "Calculate the genealogical sorting index (gsi) for labeled groups on a rooted tree."
    epilog = "Example usage: gsi.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix gsi"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("--tree", required=True, help="Newick file holding one tree, or an ensemble of trees.")
    parser.add_argument("--groups", required=True, help="Two-column table mapping tip label to group.")
    parser.add_argument(
        "--root",
        required=True,
        choices=["as-is", "midpoint", "outgroup"],
        help="How to root the tree. gsi is undefined on an unrooted tree, so this is "
        "explicit: 'as-is' trusts the Newick as written, 'midpoint' needs branch "
        "lengths, 'outgroup' needs --outgroup.",
    )
    parser.add_argument("--outgroup", help="Comma-separated tip labels to root on, for --root outgroup.")
    parser.add_argument(
        "--permutations", type=int, default=9999,
        help="Label permutations for the significance test (default: 9999, i.e. 10^4 replicates "
        "including the observed labeling). 0 skips the test.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for the permutation test (default: 42).")
    parser.add_argument("--groups-to-test", help="Comma-separated subset of groups to report (default: all).")
    parser.add_argument(
        "--ignore-unlabeled", action="store_true",
        help="Allow tips that are absent from the groups file; they stay on the tree but take no group label.",
    )
    parser.add_argument("--out-prefix", default="gsi", help="Prefix for output files (default: gsi).")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if args.root == "outgroup" and not args.outgroup:
        sys.exit("ERROR: --root outgroup requires --outgroup <tip[,tip...]>")
    if args.permutations < 0:
        sys.exit("ERROR: --permutations cannot be negative.")

    outgroup = [tip.strip() for tip in args.outgroup.split(",") if tip.strip()] if args.outgroup else []

    mapping, groups = read_groups(args.groups)
    trees = read_trees(args.tree, args.root, outgroup)
    unlabeled = check_tips(trees, mapping, args.ignore_unlabeled)

    if args.groups_to_test:
        wanted = [group.strip() for group in args.groups_to_test.split(",") if group.strip()]
        unknown = [group for group in wanted if group not in groups]
        if unknown:
            sys.exit(
                "ERROR: --groups-to-test names group(s) that are not in the groups file: {}".format(
                    ", ".join(unknown)
                )
            )
    else:
        wanted = groups

    all_rows, per_tree = analyse(trees, mapping, groups, args.permutations, args.seed)
    rows = [row for row in all_rows if row["group"] in wanted]
    n_trees = all_rows[0]["n_trees"]

    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)

    write_results("{}.gsi.tsv".format(args.out_prefix), rows)
    write_multiqc("{}.gsi_mqc.tsv".format(args.out_prefix), rows, n_trees)
    if n_trees > 1:
        write_per_tree("{}.gsi_per_tree.tsv".format(args.out_prefix), groups, per_tree)

    with open("{}.gsi.json".format(args.out_prefix), "w") as handle:
        json.dump(
            {
                "tree_file": args.tree,
                "n_trees": n_trees,
                "rooting": args.root,
                "outgroup": outgroup,
                "permutations": rows[0]["n_permutations"],
                "seed": args.seed,
                "unlabeled_tips": unlabeled,
                "results": rows,
            },
            handle,
            indent=2,
        )
        handle.write("\n")

    statistic = "gsi_T" if n_trees > 1 else "gsi"
    print(
        "Genealogical sorting index ({} tree(s), {} permutations)".format(
            n_trees, rows[0]["n_permutations"]
        )
    )
    for row in rows:
        print(
            "  {:<20} n={:<5} {}={:<8} P={}".format(
                row["group"], row["n_tips"], statistic, format_value(row["gsi"]), format_value(row["p_value"])
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
