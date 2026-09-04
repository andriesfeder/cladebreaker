#!/usr/bin/env python

"""
Tests whether observed monophyly could be a chance outcome of random branching.

Monophyly on its own is weak evidence when samples are small: a specified set of
lineages can form a clade purely by chance. These tests ask how improbable the
observed monophyly would be under a null model in which all lineages come from a
single taxon branching at random (the Yule / Yule-Harding-Kingman process, whose
topology distribution is identical to Kingman's coalescent).

Rosenberg (2007) Evolution 61:317-323, doi:10.1111/j.1558-5646.2007.00023.x

    P_A(a,b)  = [2 / C(a+b,a)] * [(a+b) / (a(a+1))]     (eq. 1)
    P_AB(a,b) = 2 / [C(a+b,a) * (a+b-1)]                (eq. 5)

P_A is monophyly of one group of size a among b other lineages of arbitrary
origin; P_AB is reciprocal monophyly of two groups that between them account for
every tip.

Zhu, Degnan & Steel (2011) Theoretical Population Biology 80:118-124,
arXiv:1101.1311, generalise this to k groups that partition the tips:

    p(a_1..a_k; T_k) = 2^(k-1) (prod a_i!) / n! * prod_{v in I(T_k)} 1/(A_v - 1)
    p(a_1..a_k)      = sum over rooted binary trees T_k on k leaves    (Thm 5.1)

where A_v is the number of tips below interior vertex v. That sum is evaluated
here by recursion over subsets rather than by enumerating the (2k-3)!! trees,
which is what makes k beyond about 8 tractable.

WHAT THESE TESTS DO NOT DO: the probability depends only on the group sizes, not
on the tree. It is not a measure of how well supported the clade is, nor of how
divergent the groups are. It answers one question only -- is the sample large
enough that chance can be excluded as the cause of the monophyly you observed.

Example usage:
    rosenberg.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix rosenberg
"""

import argparse
import json
import os
import sys
from math import comb, exp, lgamma, log, log10

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cladebreaker_phylo import (  # noqa: E402
    add_tree_arguments,
    format_value,
    group_counts,
    load_inputs,
    monophyly,
)

LN10 = log(10.0)

CAVEAT = (
    "The null model assumes the sampled lineages are a random draw from one taxon "
    "branching at random. In a cladebreaker run they are not: WhatsGNU selects the "
    "reference genomes precisely because they are the closest relatives of your "
    "isolates, and groups defined after inspecting the tree bias these p-values "
    "downwards. Treat a significant result as evidence that the sample is large "
    "enough for monophyly to be meaningful, not as proof of taxonomic distinctness."
)


# ---------------------------------------------------------------------------
# The probabilities
#
# Everything is computed in natural logs via lgamma so that the astronomically
# small values reached at realistic sample sizes neither overflow the binomial
# nor underflow the float.
# ---------------------------------------------------------------------------


def log_p_a(a, b):
    """log P_A: monophyly of a group of size a among b other lineages (eq. 1)."""
    if a < 1 or b < 1:
        return None
    n = a + b
    # 1 / C(n,a) == a! b! / n!
    return log(2.0) + lgamma(a + 1) + lgamma(b + 1) - lgamma(n + 1) + log(n) - log(a) - log(a + 1.0)


def log_p_ab(a, b):
    """log P_AB: reciprocal monophyly of two groups of sizes a and b (eq. 5)."""
    if a < 1 or b < 1 or a + b < 2:
        return None
    n = a + b
    return log(2.0) + lgamma(a + 1) + lgamma(b + 1) - lgamma(n + 1) - log(n - 1.0)


def log_p_joint(sizes):
    """log p(a_1..a_k): all k groups simultaneously clades (Zhu et al. Thm 5.1).

    Sums over rooted binary trees on k leaves by recursion over subsets. Each
    unordered split is counted once by pinning the lowest-indexed group to the
    left block. O(3^k), against O((2k-3)!!) for enumerating the trees.
    """
    k = len(sizes)
    if k < 2 or any(size < 1 for size in sizes):
        return None
    n = sum(sizes)

    below = [0] * (1 << k)
    for mask in range(1, 1 << k):
        low = mask & -mask
        below[mask] = below[mask ^ low] + sizes[low.bit_length() - 1]

    # f[S] = sum over rooted binary trees on S of prod 1/(A_v - 1) over interior v
    f = [0.0] * (1 << k)
    for i in range(k):
        f[1 << i] = 1.0

    for mask in range(1, 1 << k):
        if bin(mask).count("1") < 2:
            continue
        low = mask & -mask          # pin the lowest member to the left block
        rest = mask ^ low
        total = 0.0
        sub = rest
        while True:
            left = low | sub
            right = mask ^ left
            if right:
                total += f[left] * f[right]
            if sub == 0:
                break
            sub = (sub - 1) & rest
        f[mask] = total / (below[mask] - 1.0)

    full = (1 << k) - 1
    if f[full] <= 0.0:
        return None
    lead = (k - 1) * log(2.0) + sum(lgamma(size + 1) for size in sizes) - lgamma(n + 1)
    return lead + log(f[full])


def probability(log_value):
    """exp() that degrades to 0.0 rather than raising, for extreme logs."""
    if log_value is None:
        return None
    try:
        return exp(log_value)
    except OverflowError:
        return 0.0


def minimum_b(a, alpha, reciprocal=False, limit=100000):
    """Smallest number of outside lineages that would make the test significant.

    This is the actionable output when a result is non-significant: it says how
    many more reference genomes the run needs.
    """
    if a < 2:
        return None
    logarithm = log_p_ab if reciprocal else log_p_a
    target = log(alpha)
    b = 2 if reciprocal else 1
    while b <= limit:
        value = logarithm(a, b)
        if value is not None and value <= target:
            return b
        b += 1
    return None


def holm(p_values):
    """Holm-Bonferroni adjusted p-values, preserving input order.

    Clade events for disjoint groups are positively correlated under this null
    (Zhu et al. Corollary 4.6), so Holm is conservative here rather than invalid.
    """
    indexed = [(p, i) for i, p in enumerate(p_values) if p is not None]
    if not indexed:
        return [None] * len(p_values)
    indexed.sort()
    m = len(indexed)
    adjusted = [None] * len(p_values)
    running = 0.0
    for rank, (p, i) in enumerate(indexed):
        running = max(running, min(1.0, (m - rank) * p))
        adjusted[i] = running
    return adjusted


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def analyse(tree_indices, assignment, groups, group_sizes, alpha=0.05, max_joint_groups=14):
    """Monophyly status per group, then the Rosenberg probabilities that apply."""
    n_groups = len(group_sizes)
    n_tips = len(assignment)
    labeled_tips = sum(group_sizes)

    # Monophyly is a property of each tree; sizes are not, so the probabilities
    # below are identical across an ensemble while the status may vary.
    mono_counts = [0] * n_groups
    intruders = [[] for _ in range(n_groups)]
    for t, tree_index in enumerate(tree_indices):
        counts = group_counts(tree_index, assignment, n_groups)
        for g in range(n_groups):
            is_mono, _, breaking = monophyly(tree_index, assignment, g, group_sizes[g], counts)
            mono_counts[g] += bool(is_mono)
            if t == 0:
                intruders[g] = breaking

    rows = []
    for g, group in enumerate(groups):
        a = group_sizes[g]
        b = labeled_tips - a if labeled_tips == n_tips else n_tips - a
        fraction = mono_counts[g] / len(tree_indices)
        # The test is conditional on having observed monophyly. Reporting what
        # P_A "would have been" for a group that is not a clade invites the
        # number to be quoted as if it meant something, so withhold it.
        is_clade = mono_counts[g] == len(tree_indices)
        log_value = log_p_a(a, b) if (a >= 2 and is_clade) else None
        rows.append(
            {
                "group": group,
                "n_tips": a,
                "n_other": b,
                "monophyletic": mono_counts[g] == len(tree_indices),
                "monophyly_fraction": fraction,
                "n_clade_breaking": len(intruders[g]),
                "clade_breaking": [i for i in intruders[g]],
                "p_a": probability(log_value),
                "log10_p_a": None if log_value is None else log_value / LN10,
                "min_b_for_alpha": minimum_b(a, alpha) if (a >= 2 and is_clade) else None,
            }
        )

    monophyletic = [row for row in rows if row["monophyletic"] and row["n_tips"] >= 2]

    # Holm across only the groups the test actually applies to.
    adjusted = holm([row["p_a"] for row in monophyletic])
    for row, value in zip(monophyletic, adjusted):
        row["p_a_holm"] = value
    for row in rows:
        row.setdefault("p_a_holm", None)

    summary = {
        "n_groups": n_groups,
        "n_tips": n_tips,
        "exhaustive": labeled_tips == n_tips,
        "all_monophyletic": len(monophyletic) == n_groups,
        "alpha": alpha,
        "reciprocal": None,
        "joint": None,
    }

    # P_AB applies only when there are exactly two groups, both monophyletic,
    # and between them they account for every tip.
    if n_groups == 2 and summary["exhaustive"] and len(monophyletic) == 2:
        log_value = log_p_ab(group_sizes[0], group_sizes[1])
        summary["reciprocal"] = {
            "groups": list(groups),
            "p_ab": probability(log_value),
            "log10_p_ab": None if log_value is None else log_value / LN10,
            "min_b_for_alpha": minimum_b(group_sizes[0], alpha, reciprocal=True),
        }

    # The joint probability needs the groups to partition the tips.
    if n_groups >= 2 and summary["exhaustive"] and len(monophyletic) == n_groups:
        if n_groups <= max_joint_groups:
            log_value = log_p_joint(group_sizes)
            summary["joint"] = {
                "p": probability(log_value),
                "log10_p": None if log_value is None else log_value / LN10,
            }
        else:
            summary["joint"] = {
                "p": None,
                "log10_p": None,
                "skipped": "k={} exceeds --max-joint-groups={}".format(n_groups, max_joint_groups),
            }

    return rows, summary


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


COLUMNS = [
    "group", "n_tips", "n_other", "monophyletic", "monophyly_fraction",
    "n_clade_breaking", "p_a", "log10_p_a", "p_a_holm", "min_b_for_alpha",
]


def write_results(path, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join([
                    row["group"],
                    str(row["n_tips"]),
                    str(row["n_other"]),
                    "yes" if row["monophyletic"] else "no",
                    format_value(row["monophyly_fraction"], 3),
                    str(row["n_clade_breaking"]),
                    format_value(row["p_a"], 6),
                    format_value(row["log10_p_a"], 2),
                    format_value(row["p_a_holm"], 6),
                    "NA" if row["min_b_for_alpha"] is None else str(row["min_b_for_alpha"]),
                ]) + "\n"
            )


def write_multiqc(path, rows, summary):
    with open(path, "w") as handle:
        handle.write("# id: 'rosenberg'\n")
        handle.write("# section_name: 'Chance monophyly (Rosenberg)'\n")
        handle.write(
            "# description: 'Probability that the observed monophyly of each group is a chance "
            "outcome of random branching, after Rosenberg (2007). Depends on group sizes only, "
            "so it tests whether the sample is large enough for monophyly to mean anything. "
            + CAVEAT.replace("'", "") + "'\n"
        )
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    namespace: 'rosenberg'\n")
        handle.write("Group\tTips\tMonophyletic\tP_A\tlog10 P_A\tP_A (Holm)\n")
        for row in rows:
            handle.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    row["group"], row["n_tips"], "yes" if row["monophyletic"] else "no",
                    format_value(row["p_a"], 6), format_value(row["log10_p_a"], 2),
                    format_value(row["p_a_holm"], 6),
                )
            )


def parse_args(args=None):
    description = "Test whether observed monophyly could arise by chance under random branching."
    epilog = "Example usage: rosenberg.py --tree tree.nwk --groups groups.tsv --root midpoint --out-prefix rosenberg"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    add_tree_arguments(parser)
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level (default: 0.05).")
    parser.add_argument(
        "--max-joint-groups", type=int, default=14,
        help="Skip the joint probability beyond this many groups; the sum is O(3^k) (default: 14).",
    )
    parser.add_argument("--out-prefix", default="rosenberg", help="Prefix for output files (default: rosenberg).")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if not 0.0 < args.alpha < 1.0:
        sys.exit("ERROR: --alpha must be between 0 and 1.")

    data = load_inputs(args)
    rows, summary = analyse(
        data["tree_indices"], data["assignment"], data["groups"], data["group_sizes"],
        alpha=args.alpha, max_joint_groups=args.max_joint_groups,
    )
    labels = data["canonical_labels"]
    for row in rows:
        row["clade_breaking"] = [labels[i] for i in row["clade_breaking"]]

    reported = [row for row in rows if row["group"] in data["wanted"]]

    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)

    write_results("{}.rosenberg.tsv".format(args.out_prefix), reported)
    write_multiqc("{}.rosenberg_mqc.tsv".format(args.out_prefix), reported, summary)
    with open("{}.rosenberg.json".format(args.out_prefix), "w") as handle:
        json.dump(
            {
                "tree_file": args.tree, "n_trees": len(data["tree_indices"]),
                "rooting": args.root, "outgroup": data["outgroup"],
                "unlabeled_tips": data["unlabeled"], "summary": summary,
                "results": reported, "caveat": CAVEAT,
            },
            handle, indent=2,
        )
        handle.write("\n")

    print("Chance-monophyly test ({} tree(s), {} groups)".format(len(data["tree_indices"]), summary["n_groups"]))
    for row in reported:
        if not row["monophyletic"]:
            print("  {:<20} n={:<5} not monophyletic ({} clade-breaking tip(s)) - P_A does not apply".format(
                row["group"], row["n_tips"], row["n_clade_breaking"]))
        elif row["n_tips"] < 2:
            print("  {:<20} n=1     monophyly is trivial for a single tip".format(row["group"]))
        else:
            verdict = "" if row["p_a"] is None or row["p_a"] <= args.alpha else \
                "  <- not significant; needs b>={}".format(row["min_b_for_alpha"])
            print("  {:<20} n={:<5} monophyletic  P_A={:<12} (10^{}){}".format(
                row["group"], row["n_tips"], format_value(row["p_a"], 6),
                format_value(row["log10_p_a"], 1), verdict))

    if summary["reciprocal"]:
        r = summary["reciprocal"]
        print("\n  Reciprocal monophyly of {} and {}: P_AB={} (10^{})".format(
            r["groups"][0], r["groups"][1], format_value(r["p_ab"], 8), format_value(r["log10_p_ab"], 1)))
    if summary["joint"] and summary["joint"].get("p") is not None:
        print("  Joint monophyly of all {} groups: p={} (10^{})".format(
            summary["n_groups"], format_value(summary["joint"]["p"], 8),
            format_value(summary["joint"]["log10_p"], 1)))
    elif summary["joint"] and summary["joint"].get("skipped"):
        print("  Joint monophyly not computed: {}".format(summary["joint"]["skipped"]))

    print("\nNOTE: " + CAVEAT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
