#!/usr/bin/env python

"""
Is each labeled group a clade, and if not, exactly what breaks it?

A group is monophyletic when the smallest subtree containing all of its tips
contains nothing else. Reported per group:

  * monophyletic          the group is a clade on the tree
  * monophyly_fraction    across an ensemble, the fraction of trees in which it
                          is a clade -- clade support for the grouping itself
  * clade_breaking        the named tips sitting inside the group's MRCA that do
                          not belong to it. This is the actionable half of a
                          negative answer: it names the genomes to look at.
  * n_clusters            how many separate blocks the group is broken into,
                          counting maximal clades made only of its members. One
                          means monophyletic; three means the group sits in three
                          places on the tree.

ROOTING: monophyly is a property of a ROOTED tree, and the same caveat applies as
for the genealogical sorting index -- a mis-rooted tree still gives a confident
looking answer. RAxML-NG output is unrooted, so pass --root midpoint or
--root outgroup rather than relying on --root as-is. (Slatkin-Maddison, by
contrast, is rooting independent, which is why it is the safer fallback when
monophyly fails.)

Example usage:
    monophyly.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix monophyly
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cladebreaker_phylo import (  # noqa: E402
    add_tree_arguments,
    format_value,
    group_counts,
    load_inputs,
    monophyly,
    write_tree,
)


def count_clusters(tree_index, assignment, group, counts):
    """Maximal clades made up only of this group's tips.

    A node is 'pure' when every tip below it belongs to the group. The blocks we
    want are the pure nodes whose parent is not pure, plus any group tip whose
    parent is not pure. Counting them says how badly a broken group is scattered:
    1 is monophyletic, 2 means it sits in two places, and so on.
    """
    pure = [counts[i][group] > 0 and counts[i][group] == tree_index.leaf_total[i]
            for i in range(tree_index.n_internal)]

    # Mark every node that has a pure ancestor, so only maximal blocks survive.
    covered = [False] * tree_index.n_internal
    leaf_covered = [False] * tree_index.n_leaves
    blocks = 0

    # Walk parents-before-children by reversing the post-order.
    for i in range(tree_index.n_internal - 1, -1, -1):
        if pure[i] and not covered[i]:
            blocks += 1
            covered[i] = True
        if covered[i]:
            for is_leaf, j in tree_index.children[i]:
                if is_leaf:
                    leaf_covered[j] = True
                else:
                    covered[j] = True

    # Group tips not inside any pure block are blocks of their own.
    for j in range(tree_index.n_leaves):
        if not leaf_covered[j] and assignment[tree_index.leaf_canonical[j]] == group:
            blocks += 1
    return blocks


def analyse(tree_indices, assignment, groups, group_sizes):
    n_groups = len(group_sizes)
    rows = [
        {
            "group": group,
            "n_tips": group_sizes[g],
            "monophyletic": None,
            "monophyly_fraction": 0.0,
            "mrca_clade_size": None,
            "n_clade_breaking": 0,
            "clade_breaking": [],
            "n_clusters": None,
        }
        for g, group in enumerate(groups)
    ]

    mono_counts = [0] * n_groups
    for t, tree_index in enumerate(tree_indices):
        counts = group_counts(tree_index, assignment, n_groups)
        for g in range(n_groups):
            is_mono, mrca, breaking = monophyly(tree_index, assignment, g, group_sizes[g], counts)
            mono_counts[g] += bool(is_mono)
            if t == 0:
                rows[g]["clade_breaking"] = breaking
                rows[g]["n_clade_breaking"] = len(breaking)
                rows[g]["n_clusters"] = count_clusters(tree_index, assignment, g, counts)
                rows[g]["mrca_clade_size"] = (
                    tree_index.leaf_total[mrca] if mrca is not None else group_sizes[g]
                )

    for g, row in enumerate(rows):
        row["monophyletic"] = mono_counts[g] == len(tree_indices)
        row["monophyly_fraction"] = mono_counts[g] / len(tree_indices)
    return rows


COLUMNS = ["group", "n_tips", "monophyletic", "monophyly_fraction",
           "mrca_clade_size", "n_clade_breaking", "n_clusters"]


def write_results(path, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join([
                row["group"], str(row["n_tips"]),
                "yes" if row["monophyletic"] else "no",
                format_value(row["monophyly_fraction"], 3),
                "NA" if row["mrca_clade_size"] is None else str(row["mrca_clade_size"]),
                str(row["n_clade_breaking"]),
                "NA" if row["n_clusters"] is None else str(row["n_clusters"]),
            ]) + "\n")


def write_breaking(path, rows):
    with open(path, "w") as handle:
        handle.write("group\tclade_breaking_tip\ttip_group\n")
        for row in rows:
            for tip, owner in row["clade_breaking"]:
                handle.write("{}\t{}\t{}\n".format(row["group"], tip, owner))


def write_multiqc(path, rows):
    with open(path, "w") as handle:
        handle.write("# id: 'monophyly'\n")
        handle.write("# section_name: 'Monophyly'\n")
        handle.write(
            "# description: 'Whether each labeled group forms a clade on the rooted tree. "
            "Clusters counts how many separate blocks a broken group falls into; 1 is "
            "monophyletic. Monophyly is undefined on an unrooted tree, so check the rooting "
            "used before reading these.'\n"
        )
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    namespace: 'monophyly'\n")
        handle.write("Group\tTips\tMonophyletic\tSupport\tClade-breaking\tClusters\n")
        for row in rows:
            handle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["group"], row["n_tips"], "yes" if row["monophyletic"] else "no",
                format_value(row["monophyly_fraction"], 3), row["n_clade_breaking"],
                "NA" if row["n_clusters"] is None else row["n_clusters"]))


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Report whether each labeled group is a clade, and name the tips that break it.",
        epilog="Example usage: monophyly.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix monophyly",
    )
    add_tree_arguments(parser)
    parser.add_argument("--out-prefix", default="monophyly", help="Prefix for output files (default: monophyly).")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    data = load_inputs(args)
    rows = analyse(data["tree_indices"], data["assignment"], data["groups"], data["group_sizes"])

    labels = data["canonical_labels"]
    group_of = {}
    for g, group in enumerate(data["groups"]):
        for i, value in enumerate(data["assignment"]):
            if value == g:
                group_of[i] = group
    for row in rows:
        row["clade_breaking"] = [(labels[i], group_of.get(i, "unlabeled")) for i in row["clade_breaking"]]

    reported = [row for row in rows if row["group"] in data["wanted"]]

    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)

    # The rooted tree the tests actually ran on, for the report figure.
    write_tree(data["trees"], "{}.rooted.nwk".format(args.out_prefix))

    write_results("{}.monophyly.tsv".format(args.out_prefix), reported)
    write_breaking("{}.monophyly_clade_breaking.tsv".format(args.out_prefix), reported)
    write_multiqc("{}.monophyly_mqc.tsv".format(args.out_prefix), reported)
    with open("{}.monophyly.json".format(args.out_prefix), "w") as handle:
        json.dump({
            "tree_file": args.tree, "n_trees": len(data["tree_indices"]),
            "rooting": args.root, "outgroup": data["outgroup"],
            "unlabeled_tips": data["unlabeled"],
            "all_monophyletic": all(row["monophyletic"] for row in rows),
            "results": [dict(row, clade_breaking=[list(t) for t in row["clade_breaking"]]) for row in reported],
        }, handle, indent=2)
        handle.write("\n")

    print("Monophyly ({} tree(s), rooting: {})".format(len(data["tree_indices"]), args.root))
    for row in reported:
        if row["monophyletic"]:
            support = "" if len(data["tree_indices"]) == 1 else "  (clade in {} of trees)".format(
                format_value(row["monophyly_fraction"], 2))
            print("  {:<20} n={:<5} monophyletic{}".format(row["group"], row["n_tips"], support))
        else:
            names = ", ".join(tip for tip, _ in row["clade_breaking"][:6])
            more = " ..." if row["n_clade_breaking"] > 6 else ""
            print("  {:<20} n={:<5} NOT monophyletic - {} cluster(s), {} clade-breaking tip(s): {}{}".format(
                row["group"], row["n_tips"], row["n_clusters"], row["n_clade_breaking"], names, more))
    return 0


if __name__ == "__main__":
    sys.exit(main())
