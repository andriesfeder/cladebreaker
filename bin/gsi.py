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
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cladebreaker_phylo import (  # noqa: E402
    build_assignment,
    check_tips,
    format_value,
    group_counts,
    index_trees,
    permutation_test,
    read_groups,
    read_trees,
    reroot_at_midpoint,  # noqa: F401  (re-exported for the test suite)
    tip_labels,
    TreeIndex,           # noqa: F401  (re-exported for the test suite)
)


# ---------------------------------------------------------------------------
# The statistic
# ---------------------------------------------------------------------------


def gs_values(tree_index, assignment, group_sizes):
    """gs for every group on one tree.

    assignment maps each canonical tip index to a group index, or -1 for tips
    that carry no group label. Returns a list of gs values (None where the
    statistic is undefined).
    """
    counts = group_counts(tree_index, assignment, len(group_sizes))

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
    canonical_labels, assignment, group_sizes = build_assignment(trees, mapping, groups)

    tree_indices = index_trees(trees, canonical_labels)
    single = len(tree_indices) == 1

    gsi_values, per_tree = observed_gsi(tree_indices, assignment, group_sizes)
    p_values, n_permutations = permutation_test(
        lambda shuffled: observed_gsi(tree_indices, shuffled, group_sizes)[0],
        assignment,
        len(group_sizes),
        gsi_values,
        permutations,
        seed,
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


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


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
