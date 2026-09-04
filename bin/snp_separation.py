#!/usr/bin/env python

"""
How well separated are the labeled groups in SNP distance, not just in topology?

Monophyly, the gsi and Slatkin-Maddison all read the tree. This reads the
alignment underneath it, which catches the case where a group is a clean clade
sitting only a handful of SNPs away from its neighbours -- topologically tidy,
epidemiologically indistinguishable.

For every group two summaries are computed from the pairwise distance matrix:

    within    distances between two members of the group
    between   distances from a member of the group to any tip outside it

and from those, two ratios:

    ratio_of_means = mean(between) / mean(within)
    gap            = min(between) / max(within)

The gap is the one outbreak work usually turns on: above 1 means every pair
inside the group is closer than any pair crossing out of it, so the group is
cleanly separated in SNP space with no overlap at all.

NOT A PUBLISHED STATISTIC. Unlike the gsi (Cummings et al. 2008), Rosenberg's
P_A/P_AB (2007) or Slatkin-Maddison (1989), these ratios have no literature
behind them and no distribution of their own; the thresholds that matter are
organism- and outbreak-specific. The permutation p-value reported here is
therefore a test of a locally defined statistic: it asks whether this labeling
separates the samples better than a random labeling of the same sizes would.
Report the underlying distances alongside it, never the ratio on its own.

Distances come from a snp-dists matrix. Pass --distances for one already
computed, or --alignment to have snp-dists run over an alignment.

Example usage:
    snp_separation.py --distances core.dists.tsv --groups groups.tsv --out-prefix snpsep
"""

import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cladebreaker_phylo import (  # noqa: E402
    format_value,
    permutation_test,
    read_groups,
)

CAVEAT = (
    "The separation ratios are defined by this pipeline, not drawn from the "
    "literature, and the SNP thresholds that matter are organism- and "
    "outbreak-specific. Read them together with the raw within- and "
    "between-group distances rather than on their own."
)


# ---------------------------------------------------------------------------
# Distances
# ---------------------------------------------------------------------------


def run_snp_dists(alignment):
    """Shell out to snp-dists over an alignment."""
    try:
        completed = subprocess.run(
            ["snp-dists", "-q", alignment], check=True, stdout=subprocess.PIPE, universal_newlines=True
        )
    except FileNotFoundError:
        sys.exit(
            "ERROR: snp-dists was not found on PATH. Install it (bioconda::snp-dists) or "
            "pass an already-computed matrix with --distances."
        )
    except subprocess.CalledProcessError as error:
        sys.exit("ERROR: snp-dists failed on {} (exit {}).".format(alignment, error.returncode))
    return completed.stdout


def parse_matrix(text, source="matrix"):
    """Parse a snp-dists matrix into (sample names, list of rows).

    snp-dists puts its own name and version in the first header cell, so the
    sample names start at the second column.
    """
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        sys.exit("ERROR: {} does not look like a snp-dists matrix (fewer than two lines).".format(source))

    names = lines[0].split("\t")[1:]
    if not names:
        sys.exit("ERROR: {} has no sample columns.".format(source))

    rows = []
    seen = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(names) + 1:
            sys.exit(
                "ERROR: {} is ragged -- row '{}' has {} value(s) for {} sample(s).".format(
                    source, fields[0], len(fields) - 1, len(names)
                )
            )
        seen.append(fields[0])
        try:
            rows.append([float(value) for value in fields[1:]])
        except ValueError:
            sys.exit("ERROR: {} has a non-numeric distance in row '{}'.".format(source, fields[0]))

    if seen != names:
        sys.exit(
            "ERROR: {} is not square -- row labels do not match the column order.".format(source)
        )
    return names, rows


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------


def summarise(values):
    """min / max / mean / median over a list of distances."""
    if not values:
        return {"n": 0, "min": None, "max": None, "mean": None, "median": None}
    ordered = sorted(values)
    n = len(ordered)
    middle = n // 2
    median = ordered[middle] if n % 2 else (ordered[middle - 1] + ordered[middle]) / 2.0
    return {
        "n": n,
        "min": ordered[0],
        "max": ordered[-1],
        "mean": sum(ordered) / n,
        "median": median,
    }


def split_distances(rows, assignment, group):
    """Within-group and out-of-group distances for one group."""
    members = [i for i, g in enumerate(assignment) if g == group]
    member_set = set(members)
    within, between = [], []
    for a_index, i in enumerate(members):
        for j in members[a_index + 1:]:
            within.append(rows[i][j])
        for j in range(len(assignment)):
            # Only labeled tips outside the group count as "between".
            if j not in member_set and assignment[j] >= 0:
                between.append(rows[i][j])
    return within, between


def ratios(within, between):
    """ratio_of_means and gap, or None where they are undefined."""
    if not within or not between:
        return None, None
    mean_within = sum(within) / len(within)
    mean_between = sum(between) / len(between)
    max_within = max(within)
    ratio = mean_between / mean_within if mean_within > 0 else None
    gap = min(between) / max_within if max_within > 0 else None
    return ratio, gap


def statistic_for(rows, n_groups, key="ratio"):
    """A statistic function over an assignment, for the permutation driver."""
    def statistic(assignment):
        values = []
        for g in range(n_groups):
            within, between = split_distances(rows, assignment, g)
            ratio, gap = ratios(within, between)
            values.append(ratio if key == "ratio" else gap)
        return values
    return statistic


def analyse(names, rows, assignment, groups, permutations=9999, seed=42):
    n_groups = len(groups)
    results = []
    for g, group in enumerate(groups):
        within, between = split_distances(rows, assignment, g)
        ratio, gap = ratios(within, between)
        results.append({
            "group": group,
            "n_tips": sum(1 for value in assignment if value == g),
            "within": summarise(within),
            "between": summarise(between),
            "ratio_of_means": ratio,
            "gap": gap,
            "zero_within_diversity": bool(within) and max(within) == 0,
            "cleanly_separated": gap is not None and gap > 1.0,
        })

    observed = [row["ratio_of_means"] for row in results]
    p_values, n_permutations = permutation_test(
        statistic_for(rows, n_groups), assignment, n_groups, observed,
        permutations, seed, larger_is_extreme=True,
    )
    for row, p_value in zip(results, p_values):
        row["p_value"] = p_value
        row["n_permutations"] = n_permutations

    pairs = []
    for i in range(n_groups):
        for j in range(i + 1, n_groups):
            values = [
                rows[a][b]
                for a in range(len(assignment)) if assignment[a] == i
                for b in range(len(assignment)) if assignment[b] == j
            ]
            pairs.append({"group_a": groups[i], "group_b": groups[j], "distances": summarise(values)})

    return results, pairs


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


COLUMNS = ["group", "n_tips", "within_n", "within_min", "within_max", "within_mean", "within_median",
           "between_n", "between_min", "between_max", "between_mean", "between_median",
           "ratio_of_means", "gap", "cleanly_separated", "p_value", "n_permutations"]


def write_results(path, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            w, b = row["within"], row["between"]
            handle.write("\t".join([
                row["group"], str(row["n_tips"]),
                str(w["n"]), format_value(w["min"], 1), format_value(w["max"], 1),
                format_value(w["mean"], 2), format_value(w["median"], 1),
                str(b["n"]), format_value(b["min"], 1), format_value(b["max"], 1),
                format_value(b["mean"], 2), format_value(b["median"], 1),
                format_value(row["ratio_of_means"], 3), format_value(row["gap"], 3),
                "yes" if row["cleanly_separated"] else "no",
                format_value(row["p_value"], 6), str(row["n_permutations"]),
            ]) + "\n")


def write_pairs(path, pairs):
    with open(path, "w") as handle:
        handle.write("group_a\tgroup_b\tn_pairs\tmin\tmax\tmean\tmedian\n")
        for pair in pairs:
            d = pair["distances"]
            handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                pair["group_a"], pair["group_b"], d["n"], format_value(d["min"], 1),
                format_value(d["max"], 1), format_value(d["mean"], 2), format_value(d["median"], 1)))


def write_multiqc(path, rows):
    with open(path, "w") as handle:
        handle.write("# id: 'snp_separation'\n")
        handle.write("# section_name: 'SNP separation'\n")
        handle.write(
            "# description: 'Within- and between-group SNP distances. Gap is the smallest "
            "between-group distance over the largest within-group distance; above 1 means the "
            "group is cleanly separated with no overlap. " + CAVEAT.replace("'", "") + "'\n"
        )
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    namespace: 'snp_separation'\n")
        handle.write("Group\tTips\tMax within\tMin between\tGap\tMean ratio\tp-value\n")
        for row in rows:
            handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["group"], row["n_tips"], format_value(row["within"]["max"], 1),
                format_value(row["between"]["min"], 1), format_value(row["gap"], 3),
                format_value(row["ratio_of_means"], 3), format_value(row["p_value"], 6)))


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Summarise within- and between-group SNP distances and their separation.",
        epilog="Example usage: snp_separation.py --distances core.dists.tsv --groups groups.tsv --out-prefix snpsep",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--distances", help="A snp-dists matrix (TSV).")
    source.add_argument("--alignment", help="An alignment; snp-dists is run over it to build the matrix.")
    parser.add_argument("--groups", required=True, help="Two-column table mapping sample name to group.")
    parser.add_argument("--groups-to-test", help="Comma-separated subset of groups to report (default: all).")
    parser.add_argument(
        "--ignore-unlabeled", action="store_true",
        help="Allow samples absent from the groups file; they take no group label.",
    )
    parser.add_argument("--permutations", type=int, default=9999, help="Label permutations (default: 9999). 0 skips.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    parser.add_argument("--out-prefix", default="snp_separation", help="Prefix for output files.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if args.permutations < 0:
        sys.exit("ERROR: --permutations cannot be negative.")

    if args.distances:
        with open(args.distances) as handle:
            names, rows = parse_matrix(handle.read(), args.distances)
    else:
        names, rows = parse_matrix(run_snp_dists(args.alignment), "the snp-dists output")

    mapping, groups = read_groups(args.groups)

    missing = sorted(name for name in mapping if name not in set(names))
    if missing:
        sys.exit(
            "ERROR: Please check the groups file -> {} sample(s) are not in the distance "
            "matrix: {}".format(len(missing), ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))
        )
    unlabeled = sorted(name for name in names if name not in mapping)
    if unlabeled and not args.ignore_unlabeled:
        sys.exit(
            "ERROR: {} sample(s) in the distance matrix have no group: {}\n"
            "Assign them a group, or pass --ignore-unlabeled to leave them out.".format(
                len(unlabeled), ", ".join(unlabeled[:10]) + (" ..." if len(unlabeled) > 10 else "")
            )
        )

    group_index = {group: g for g, group in enumerate(groups)}
    assignment = [group_index.get(mapping.get(name, None), -1) for name in names]

    wanted = groups
    if args.groups_to_test:
        wanted = [group.strip() for group in args.groups_to_test.split(",") if group.strip()]
        unknown = [group for group in wanted if group not in groups]
        if unknown:
            sys.exit("ERROR: --groups-to-test names group(s) not in the groups file: {}".format(", ".join(unknown)))

    results, pairs = analyse(names, rows, assignment, groups, args.permutations, args.seed)
    reported = [row for row in results if row["group"] in wanted]

    directory = os.path.dirname(args.out_prefix)
    if directory:
        os.makedirs(directory, exist_ok=True)

    write_results("{}.snp_separation.tsv".format(args.out_prefix), reported)
    write_pairs("{}.snp_separation_pairs.tsv".format(args.out_prefix), pairs)
    write_multiqc("{}.snp_separation_mqc.tsv".format(args.out_prefix), reported)
    with open("{}.snp_separation.json".format(args.out_prefix), "w") as handle:
        json.dump({
            "distances": args.distances, "alignment": args.alignment,
            "n_samples": len(names), "unlabeled_samples": unlabeled,
            "results": reported, "pairs": pairs, "caveat": CAVEAT,
        }, handle, indent=2)
        handle.write("\n")

    print("SNP separation ({} samples, {} permutations)".format(len(names), args.permutations))
    for row in reported:
        w, b = row["within"], row["between"]
        if w["n"] == 0:
            print("  {:<20} n={:<5} only one member - no within-group distances".format(
                row["group"], row["n_tips"]))
            continue
        print("  {:<20} n={:<5} within {}-{} (mean {})  between {}-{} (mean {})".format(
            row["group"], row["n_tips"], format_value(w["min"], 0), format_value(w["max"], 0),
            format_value(w["mean"], 1), format_value(b["min"], 0), format_value(b["max"], 0),
            format_value(b["mean"], 1)))
        verdict = "cleanly separated" if row["cleanly_separated"] else "overlapping"
        print("  {:<20}       gap={}  ratio={}  P={}  ({})".format(
            "", format_value(row["gap"], 2), format_value(row["ratio_of_means"], 2),
            format_value(row["p_value"], 4), verdict))

    print("\nNOTE: " + CAVEAT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
