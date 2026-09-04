#!/usr/bin/env python

"""
Slatkin-Maddison test: is the mixing of labeled groups on a tree less than random?

Slatkin & Maddison (1989) Genetics 123:603-613.

The group label is treated as an unordered multistate character on the tips, and
Fitch parsimony gives s, the fewest state changes ("migration events") that can
explain the observed labeling on the tree. Freely mixing groups are interleaved
across the genealogy and need many changes; structured groups need few. The null
is built by shuffling the labels across the tips with the topology held fixed:

    P = (1 + #{s_permuted <= s_observed}) / (permutations + 1)

so a small P means the observed labeling requires significantly fewer migrations
than chance would demand -- evidence of structure short of monophyly. This is the
test to reach for when a group is NOT strictly monophyletic: it asks whether the
residual mixing is still less than random, rather than treating anything short of
a clean clade as a failure.

ROOTING: unlike monophyly, the genealogical sorting index and Rosenberg's tests,
the Fitch change count is a property of the UNROOTED topology -- rerooting cannot
change the minimum number of changes. So --root is accepted for consistency with
the other tools but does not affect the result, which makes this the safer test
when the correct rooting is exactly what you are unsure about.

Tips carrying no group label are treated as missing data (a wildcard state), so
they constrain nothing and are never permuted.

Example usage:
    slatkin_maddison.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix sm
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cladebreaker_phylo import (  # noqa: E402
    add_tree_arguments,
    format_value,
    load_inputs,
    permutation_test,
)

WILDCARD = -1  # an unlabeled tip constrains nothing


def fitch_changes(tree_index, assignment, n_groups):
    """Minimum number of state changes on the tree, by Fitch parsimony.

    Post-order over the flattened tree. At each internal node, take the states
    appearing in the largest number of children; the cost is the number of
    children that did not carry one of them. On a binary node that reduces to the
    familiar rule -- intersect if you can, else union and add one.

    States are held as bitmasks. A tip with no group label gets the full mask, so
    it is compatible with anything and forces no change.
    """
    full = (1 << n_groups) - 1
    states = [0] * tree_index.n_internal
    changes = 0

    for i, kids in enumerate(tree_index.children):
        child_states = []
        for is_leaf, j in kids:
            if is_leaf:
                group = assignment[tree_index.leaf_canonical[j]]
                child_states.append(full if group == WILDCARD else (1 << group))
            else:
                child_states.append(states[j])

        # How many children admit each state?
        tally = [0] * n_groups
        for mask in child_states:
            for g in range(n_groups):
                if mask >> g & 1:
                    tally[g] += 1
        best = max(tally) if tally else 0

        if best == 0:                      # no shared state at all
            states[i] = full
            changes += len(child_states) - 1
            continue

        combined = 0
        for g in range(n_groups):
            if tally[g] == best:
                combined |= 1 << g
        states[i] = combined
        changes += len(child_states) - best

    return changes


def observed_changes(tree_indices, assignment, n_groups):
    """Mean s across an ensemble; the plain count for a single tree."""
    values = [fitch_changes(tree_index, assignment, n_groups) for tree_index in tree_indices]
    return sum(values) / len(values), values


def analyse(tree_indices, assignment, groups, group_sizes, permutations=9999, seed=42):
    n_groups = len(group_sizes)
    observed, per_tree = observed_changes(tree_indices, assignment, n_groups)

    null_values = []

    def statistic(shuffled):
        return [observed_changes(tree_indices, shuffled, n_groups)[0]]

    p_values, n_permutations = permutation_test(
        statistic, assignment, 1, [observed], permutations, seed,
        larger_is_extreme=False,                     # FEWER changes means structure
        observer=lambda values: null_values.append(values[0]),
    )

    null_mean = sum(null_values) / len(null_values) if null_values else None
    null_sd = None
    if len(null_values) > 1:
        variance = sum((v - null_mean) ** 2 for v in null_values) / (len(null_values) - 1)
        null_sd = variance ** 0.5

    # A labeling can never need fewer than one change per group beyond the first.
    present = sum(1 for size in group_sizes if size > 0)
    return {
        "n_tips": len(assignment),
        "n_labeled_tips": sum(group_sizes),
        "n_groups": present,
        "n_trees": len(tree_indices),
        "s_observed": observed,
        "s_per_tree": per_tree,
        "s_minimum_possible": max(present - 1, 0),
        "null_mean": null_mean,
        "null_sd": null_sd,
        "p_value": p_values[0],
        "n_permutations": n_permutations,
        "groups": list(groups),
    }


COLUMNS = ["n_tips", "n_labeled_tips", "n_groups", "n_trees", "s_observed",
           "s_minimum_possible", "null_mean", "null_sd", "p_value", "n_permutations"]


def write_results(path, result):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        handle.write("\t".join([
            str(result["n_tips"]), str(result["n_labeled_tips"]), str(result["n_groups"]),
            str(result["n_trees"]), format_value(result["s_observed"], 2),
            str(result["s_minimum_possible"]), format_value(result["null_mean"], 2),
            format_value(result["null_sd"], 2), format_value(result["p_value"], 6),
            str(result["n_permutations"]),
        ]) + "\n")


def write_multiqc(path, result):
    with open(path, "w") as handle:
        handle.write("# id: 'slatkin_maddison'\n")
        handle.write("# section_name: 'Slatkin-Maddison migration test'\n")
        handle.write(
            "# description: 'Fewest label changes (migration events) needed to explain the "
            "grouping on the tree, against a null from shuffling the labels. Fewer changes "
            "than chance demands means the groups are structured even where they are not "
            "strictly monophyletic. The statistic is independent of rooting.'\n"
        )
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    namespace: 'slatkin_maddison'\n")
        handle.write("Metric\tValue\n")
        handle.write("Migration events (s)\t{}\n".format(format_value(result["s_observed"], 2)))
        handle.write("Expected under null\t{}\n".format(format_value(result["null_mean"], 2)))
        handle.write("Minimum possible\t{}\n".format(result["s_minimum_possible"]))
        handle.write("p-value\t{}\n".format(format_value(result["p_value"], 6)))


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Slatkin-Maddison test for population structure on a tree.",
        epilog="Example usage: slatkin_maddison.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix sm",
    )
    add_tree_arguments(parser, permutations_default=9999)
    parser.add_argument("--out-prefix", default="slatkin_maddison", help="Prefix for output files.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if args.permutations < 0:
        sys.exit("ERROR: --permutations cannot be negative.")

    data = load_inputs(args)
    result = analyse(
        data["tree_indices"], data["assignment"], data["groups"], data["group_sizes"],
        permutations=args.permutations, seed=args.seed,
    )
    result.update({
        "tree_file": args.tree, "rooting": args.root, "outgroup": data["outgroup"],
        "unlabeled_tips": data["unlabeled"],
    })

    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)

    write_results("{}.slatkin_maddison.tsv".format(args.out_prefix), result)
    write_multiqc("{}.slatkin_maddison_mqc.tsv".format(args.out_prefix), result)
    with open("{}.slatkin_maddison.json".format(args.out_prefix), "w") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")

    print("Slatkin-Maddison migration test ({} tree(s), {} permutations)".format(
        result["n_trees"], result["n_permutations"]))
    print("  groups            : {}".format(", ".join(result["groups"])))
    print("  migration events  : {}   (minimum possible {})".format(
        format_value(result["s_observed"], 2), result["s_minimum_possible"]))
    print("  expected by chance: {}{}".format(
        format_value(result["null_mean"], 2),
        "" if result["null_sd"] is None else " +/- {}".format(format_value(result["null_sd"], 2))))
    print("  p-value           : {}".format(format_value(result["p_value"], 6)))
    if result["p_value"] is not None and result["p_value"] <= 0.05:
        print("\n  The groups need significantly fewer migrations than chance: structured,")
        print("  even though this says nothing about whether any group is a clade.")
    elif result["p_value"] is not None:
        print("\n  Not distinguishable from random mixing of the group labels.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
