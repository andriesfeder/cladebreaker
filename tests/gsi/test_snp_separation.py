#!/usr/bin/env python

"""
Tests for bin/snp_separation.py.

The ratios here are defined by this pipeline rather than taken from a paper, so
there are no published values to check against. Instead the tests pin the
arithmetic against hand-computable matrices, cover the degenerate cases that a
real outbreak throws up (identical genomes, singleton groups, overlapping
distance ranges), and check the snp-dists parser against real output format.

Run with pytest, or directly:  python tests/gsi/test_snp_separation.py
"""

import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "bin"))

spec = importlib.util.spec_from_file_location("snp_separation", os.path.join(ROOT, "bin", "snp_separation.py"))
snpsep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(snpsep)


def close(observed, expected, tolerance=1e-9):
    assert observed is not None, "expected {}, got None".format(expected)
    assert abs(observed - expected) <= tolerance, "expected {}, got {}".format(expected, observed)


def matrix(names, rows):
    """Build a snp-dists style matrix string."""
    out = ["snp-dists 0.8.2\t" + "\t".join(names)]
    for name, row in zip(names, rows):
        out.append(name + "\t" + "\t".join(str(v) for v in row))
    return "\n".join(out) + "\n"


# Two tight clusters, far apart: a = {s1,s2,s3}, b = {s4,s5,s6}
NAMES = ["s1", "s2", "s3", "s4", "s5", "s6"]
ROWS = [
    [0, 2, 4, 100, 102, 104],
    [2, 0, 3, 101, 103, 105],
    [4, 3, 0, 100, 100, 102],
    [100, 101, 100, 0, 5, 6],
    [102, 103, 100, 5, 0, 4],
    [104, 105, 102, 6, 4, 0],
]
ASSIGN = [0, 0, 0, 1, 1, 1]
GROUPS = ["a", "b"]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def test_parses_a_snp_dists_matrix():
    names, rows = snpsep.parse_matrix(matrix(NAMES, ROWS))
    assert names == NAMES
    assert rows[0][3] == 100.0
    assert rows[5][0] == 104.0


def test_rejects_a_ragged_matrix():
    broken = "snp-dists 0.8.2\ts1\ts2\ns1\t0\t3\ns2\t3\n"
    try:
        snpsep.parse_matrix(broken)
    except SystemExit as error:
        assert "ragged" in str(error)
    else:
        raise AssertionError("expected a ragged matrix to be refused")


def test_rejects_a_matrix_whose_rows_are_reordered():
    scrambled = "snp-dists 0.8.2\ts1\ts2\ns2\t0\t3\ns1\t3\t0\n"
    try:
        snpsep.parse_matrix(scrambled)
    except SystemExit as error:
        assert "not square" in str(error)
    else:
        raise AssertionError("expected mismatched row labels to be refused")


def test_rejects_non_numeric_distances():
    try:
        snpsep.parse_matrix("snp-dists 0.8.2\ts1\ts2\ns1\t0\tNA\ns2\tNA\t0\n")
    except SystemExit as error:
        assert "non-numeric" in str(error)
    else:
        raise AssertionError("expected a non-numeric distance to be refused")


# ---------------------------------------------------------------------------
# The summaries
# ---------------------------------------------------------------------------


def test_summarise():
    s = snpsep.summarise([4, 1, 3, 2])
    assert s["n"] == 4
    close(s["min"], 1); close(s["max"], 4)
    close(s["mean"], 2.5); close(s["median"], 2.5)      # even count -> midpoint
    close(snpsep.summarise([5, 1, 3])["median"], 3)     # odd count -> middle
    assert snpsep.summarise([])["n"] == 0
    assert snpsep.summarise([])["mean"] is None


def test_within_and_between_are_split_correctly():
    within, between = snpsep.split_distances(ROWS, ASSIGN, 0)
    assert sorted(within) == [2, 3, 4]                   # the three a-a pairs
    assert len(between) == 9                             # 3 a tips x 3 b tips
    close(min(between), 100)


def test_ratios_on_a_hand_computed_matrix():
    within, between = snpsep.split_distances(ROWS, ASSIGN, 0)
    ratio, gap = snpsep.ratios(within, between)
    close(ratio, (sum(between) / 9.0) / 3.0)             # mean within = (2+3+4)/3 = 3
    close(gap, 100.0 / 4.0)                              # min between 100, max within 4
    assert gap > 1


def test_gap_below_one_when_the_ranges_overlap():
    rows = [
        [0, 30, 5, 6],
        [30, 0, 7, 8],
        [5, 7, 0, 25],
        [6, 8, 25, 0],
    ]
    within, between = snpsep.split_distances(rows, [0, 0, 1, 1], 0)
    ratio, gap = snpsep.ratios(within, between)
    close(gap, 5.0 / 30.0)                               # closest outsider beats the furthest member
    assert gap < 1
    assert ratio < 1


def test_identical_genomes_leave_the_gap_undefined():
    """Zero within-group diversity would divide by zero, so the gap is withheld."""
    rows = [
        [0, 0, 50],
        [0, 0, 50],
        [50, 50, 0],
    ]
    within, between = snpsep.split_distances(rows, [0, 0, 1], 0)
    ratio, gap = snpsep.ratios(within, between)
    assert gap is None
    assert ratio is None
    results, _ = snpsep.analyse(["s1", "s2", "s3"], rows, [0, 0, 1], ["a", "b"], permutations=0)
    assert results[0]["zero_within_diversity"] is True
    assert results[0]["cleanly_separated"] is False


def test_a_singleton_group_has_no_within_distances():
    results, _ = snpsep.analyse(NAMES, ROWS, [0, 1, 1, 1, 1, 1], ["a", "b"], permutations=0)
    assert results[0]["within"]["n"] == 0
    assert results[0]["ratio_of_means"] is None
    assert results[0]["gap"] is None


def test_unlabeled_samples_are_left_out_of_between():
    """A sample with no group must not count as an outsider."""
    labeled = snpsep.split_distances(ROWS, [0, 0, 0, 1, 1, 1], 0)[1]
    dropped = snpsep.split_distances(ROWS, [0, 0, 0, 1, 1, -1], 0)[1]
    assert len(labeled) == 9
    assert len(dropped) == 6


# ---------------------------------------------------------------------------
# The test itself
# ---------------------------------------------------------------------------


def clustered_matrix(sizes, within=5, between=100):
    """Tight clusters far apart, with slight jitter so ties do not dominate."""
    labels = [g for g, size in enumerate(sizes) for _ in range(size)]
    n = len(labels)
    rows = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            base = within if labels[i] == labels[j] else between
            rows[i][j] = rows[j][i] = base + ((i * 7 + j * 3) % 5)
    return ["s%d" % i for i in range(n)], rows, labels


def test_separated_groups_are_significant():
    names, rows, assignment = clustered_matrix([6, 6])
    results, pairs = snpsep.analyse(names, rows, assignment, GROUPS, permutations=999, seed=3)
    for row in results:
        assert row["cleanly_separated"] is True
        assert row["p_value"] < 0.01, row["p_value"]
    assert len(pairs) == 1
    assert pairs[0]["distances"]["n"] == 36


def test_the_p_value_floor_is_set_by_sample_size():
    """Six samples split 3/3 give only C(6,3)=20 labelings, so no amount of
    separation can push p below about 1/20. Both groups here are cleanly
    separated -- gap well above 1 -- yet neither reaches a convincing p-value.
    The same sample-size limit applies to Slatkin-Maddison and to Rosenberg's
    tables, and it is why a null result on a small set is not evidence of
    mixing."""
    results, pairs = snpsep.analyse(NAMES, ROWS, ASSIGN, GROUPS, permutations=999, seed=3)
    for row in results:
        assert row["cleanly_separated"] is True     # descriptively unambiguous
        assert row["p_value"] >= 0.04, row["p_value"]   # yet cannot beat the floor
    assert len(pairs) == 1
    close(pairs[0]["distances"]["min"], 100)
    assert pairs[0]["distances"]["n"] == 9


def test_a_meaningless_labeling_is_not_significant():
    """Splitting each real cluster across both labels destroys the signal."""
    results, _ = snpsep.analyse(NAMES, ROWS, [0, 1, 0, 1, 0, 1], GROUPS, permutations=999, seed=3)
    for row in results:
        assert row["p_value"] > 0.05, row["p_value"]
        assert row["cleanly_separated"] is False


def test_permutation_is_reproducible_and_skippable():
    first = snpsep.analyse(NAMES, ROWS, ASSIGN, GROUPS, permutations=499, seed=11)[0]
    second = snpsep.analyse(NAMES, ROWS, ASSIGN, GROUPS, permutations=499, seed=11)[0]
    assert [r["p_value"] for r in first] == [r["p_value"] for r in second]
    none = snpsep.analyse(NAMES, ROWS, ASSIGN, GROUPS, permutations=0)[0]
    assert none[0]["p_value"] is None
    assert none[0]["gap"] is not None      # the descriptive part still works


def test_three_groups_give_every_pairwise_summary():
    results, pairs = snpsep.analyse(NAMES, ROWS, [0, 0, 1, 1, 2, 2], ["a", "b", "c"], permutations=0)
    assert len(results) == 3
    assert {(p["group_a"], p["group_b"]) for p in pairs} == {("a", "b"), ("a", "c"), ("b", "c")}


# ---------------------------------------------------------------------------
# The fixture the test_analyze profile runs on
# ---------------------------------------------------------------------------


def read_fasta(path):
    names, seqs, current = [], [], []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                names.append(line[1:])
                current = []
            elif line:
                current.append(line)
    seqs.append("".join(current))
    return names, seqs


def snp_matrix(seqs):
    """Pairwise differing-site counts, the same thing snp-dists computes."""
    n = len(seqs)
    rows = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = sum(1 for x, y in zip(seqs[i], seqs[j]) if x != y)
            rows[i][j] = rows[j][i] = d
    return rows


def test_the_test_analyze_fixture_has_the_distances_it_claims():
    """conf/test_analyze.config documents this alignment as 2 within a, 4 within
    b and 9 between. Pin it, so editing the fixture cannot quietly invalidate
    the profile's documented behaviour."""
    names, seqs = read_fasta(os.path.join(HERE, "cummings_fig1.aln"))
    assert names == ["a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4"]
    rows = snp_matrix(seqs)
    assignment = [0, 0, 0, 0, 1, 1, 1, 1]

    results, pairs = snpsep.analyse(names, rows, assignment, ["a", "b"], permutations=0)
    a, b = results

    assert (a["within"]["min"], a["within"]["max"]) == (2, 2)
    assert (b["within"]["min"], b["within"]["max"]) == (4, 4)
    assert (a["between"]["min"], a["between"]["max"]) == (9, 9)
    close(a["gap"], 9 / 2.0)          # 4.5
    close(b["gap"], 9 / 4.0)          # 2.25
    close(a["ratio_of_means"], 9 / 2.0)
    close(b["ratio_of_means"], 9 / 4.0)
    assert a["cleanly_separated"] and b["cleanly_separated"]
    assert pairs[0]["distances"]["n"] == 16


def main():
    tests = [v for n, v in sorted(globals().items()) if n.startswith("test_")]
    failures = 0
    for test in tests:
        try:
            test()
        except AssertionError as error:
            failures += 1
            print("FAIL  {}: {}".format(test.__name__, error))
        else:
            print("ok    {}".format(test.__name__))
    print("\n{} passed, {} failed".format(len(tests) - failures, failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
