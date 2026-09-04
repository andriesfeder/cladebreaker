#!/usr/bin/env python

"""
Tests for bin/rosenberg.py.

The reference values are published, which makes this unusually checkable:

  * Rosenberg (2007) Tables 1 and 6, plus the worked examples in its text.
  * Zhu, Degnan & Steel (2011) Theorem 5.1, whose stated p(2,2,2) = 2/225 and
    whose k=3 closed form both have to fall out of the general subset recursion,
    and which must reduce to Rosenberg's P_AB at k=2.

Run with pytest, or directly:  python tests/gsi/test_rosenberg.py
"""

import importlib.util
import os
import sys
from fractions import Fraction
from math import comb, exp, factorial

import dendropy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "bin"))

spec = importlib.util.spec_from_file_location("rosenberg", os.path.join(ROOT, "bin", "rosenberg.py"))
rosenberg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rosenberg)

import cladebreaker_phylo as phylo  # noqa: E402


def P_A(a, b):
    return exp(rosenberg.log_p_a(a, b))


def P_AB(a, b):
    return exp(rosenberg.log_p_ab(a, b))


def p_joint(sizes):
    return exp(rosenberg.log_p_joint(sizes))


def close(observed, expected, tolerance=1e-9):
    assert observed is not None, "expected {}, got None".format(expected)
    assert abs(observed - expected) <= tolerance, "expected {}, got {}".format(expected, observed)


# ---------------------------------------------------------------------------
# Rosenberg (2007)
# ---------------------------------------------------------------------------

TABLE_1 = {  # P_A, significance levels for monophyly of group A
    (2, 1): 0.333, (2, 2): 0.222, (2, 3): 0.167, (2, 4): 0.133, (2, 5): 0.111,
    (2, 10): 0.061, (2, 20): 0.032, (2, 50): 0.013, (3, 1): 0.167, (3, 2): 0.083,
    (3, 3): 0.050, (3, 4): 0.033, (3, 5): 0.024, (4, 1): 0.100, (4, 2): 0.040,
    (4, 3): 0.020, (4, 4): 0.011, (5, 1): 0.067, (5, 2): 0.022, (10, 1): 0.018,
}

TABLE_6 = {  # P_AB, reciprocal monophyly
    (2, 2): 0.111, (2, 3): 0.050, (2, 4): 0.027, (2, 5): 0.016,
    (3, 2): 0.050, (3, 3): 0.020, (4, 2): 0.027, (5, 2): 0.016,
}


def test_table_1_probability_of_monophyly():
    for (a, b), expected in TABLE_1.items():
        close(P_A(a, b), expected, tolerance=0.0006)


def test_table_6_probability_of_reciprocal_monophyly():
    for (a, b), expected in TABLE_6.items():
        close(P_AB(a, b), expected, tolerance=0.0006)


def test_worked_examples_from_the_text():
    close(P_A(6, 6), 6.18e-4, tolerance=1e-6)
    close(P_AB(6, 6), 1.97e-4, tolerance=1e-6)
    close(P_AB(6, 6) / P_A(6, 6), 7.0 / 22.0, tolerance=1e-9)
    close(P_A(986, 1), 2.06e-6, tolerance=1e-8)
    close(P_A(4, 1), 0.1, tolerance=1e-9)


def test_closed_forms_match_the_exact_rationals():
    """Guards the lgamma implementation against drift from the exact formulas."""
    for a in range(2, 12):
        for b in range(1, 12):
            close(P_A(a, b), float(Fraction(2, comb(a + b, a)) * Fraction(a + b, a * (a + 1))), 1e-12)
            close(P_AB(a, b), float(Fraction(2, comb(a + b, a) * (a + b - 1))), 1e-12)


def test_monophyly_of_a_single_lineage_is_certain():
    """The paper notes a=1 is degenerate: monophyly is guaranteed."""
    for b in (1, 5, 50):
        close(P_A(1, b), 1.0)


def test_logs_survive_sample_sizes_that_would_overflow_a_binomial():
    # C(2000,1000) has ~600 digits; the value itself is far below float minimum.
    assert rosenberg.log_p_ab(1000, 1000) / rosenberg.LN10 < -500
    assert rosenberg.probability(rosenberg.log_p_ab(1000, 1000)) == 0.0


def test_minimum_b_matches_rosenbergs_table_3():
    """Table 3: smallest b making P_A <= alpha, for a fixed a."""
    assert rosenberg.minimum_b(2, 1e-1) == 6
    assert rosenberg.minimum_b(3, 1e-1) == 2
    assert rosenberg.minimum_b(4, 1e-1) == 1
    assert rosenberg.minimum_b(2, 1e-2) == 66
    assert rosenberg.minimum_b(3, 1e-2) == 9


# ---------------------------------------------------------------------------
# Zhu, Degnan & Steel (2011) Theorem 5.1
# ---------------------------------------------------------------------------


def test_joint_probability_matches_the_papers_example():
    close(p_joint([2, 2, 2]), 2.0 / 225.0, tolerance=1e-12)


def test_joint_probability_reduces_to_rosenberg_at_k_equals_two():
    for a, b in [(2, 2), (3, 3), (4, 2), (6, 6), (5, 9), (10, 7), (20, 13)]:
        close(p_joint([a, b]), P_AB(a, b), tolerance=1e-12)


def test_joint_probability_matches_the_k_equals_three_closed_form():
    def closed(a1, a2, a3):
        n = a1 + a2 + a3
        s = sum(Fraction(1, n - ai - 1) for ai in (a1, a2, a3))
        return float(Fraction(4 * factorial(a1) * factorial(a2) * factorial(a3),
                              factorial(n) * (n - 1)) * s)

    for triple in [(2, 2, 2), (3, 4, 5), (5, 5, 5), (2, 7, 4), (10, 3, 6)]:
        close(p_joint(list(triple)), closed(*triple), tolerance=1e-15)


def test_joint_probability_is_order_independent():
    """Exchangeability: only the multiset of group sizes matters."""
    close(p_joint([2, 5, 3]), p_joint([5, 3, 2]), tolerance=1e-15)
    close(p_joint([1, 4, 2, 6]), p_joint([6, 2, 4, 1]), tolerance=1e-15)


def test_joint_probability_falls_as_groups_are_added():
    previous = 1.0
    for k in range(2, 9):
        value = p_joint([3] * k)
        assert value < previous
        previous = value


# ---------------------------------------------------------------------------
# Holm
# ---------------------------------------------------------------------------


def test_holm_adjustment():
    close(rosenberg.holm([0.01])[0], 0.01)
    # m=3: 0.01*3=0.03, 0.02*2=0.04, 0.03*1=0.03 -> monotone: 0.03, 0.04, 0.04
    adjusted = rosenberg.holm([0.01, 0.02, 0.03])
    close(adjusted[0], 0.03)
    close(adjusted[1], 0.04)
    close(adjusted[2], 0.04)
    # Order is preserved and Nones pass through. Only the two real tests count
    # towards m, so this is Holm at m=2: 2*0.01=0.02, then max(0.02, 1*0.03)=0.03.
    adjusted = rosenberg.holm([0.03, None, 0.01])
    close(adjusted[0], 0.03)
    assert adjusted[1] is None
    close(adjusted[2], 0.02)
    assert all(v <= 1.0 for v in rosenberg.holm([0.9, 0.8, 0.7]) if v is not None)


# ---------------------------------------------------------------------------
# End to end on the paper's tree
# ---------------------------------------------------------------------------

PAPER_GROUPS = {"a1": "a", "a2": "a", "a3": "a", "a4": "a",
                "b1": "b", "b2": "b", "b3": "b", "b4": "b"}


def load_fixture(name):
    with open(os.path.join(HERE, name)) as handle:
        return dendropy.TreeList.get(data=handle.read(), schema="newick",
                                     rooting="force-rooted", preserve_underscores=True)


def analyse_fixture(name, mapping=PAPER_GROUPS, groups=("a", "b")):
    trees = load_fixture(name)
    labels, assignment, sizes = phylo.build_assignment(trees, mapping, list(groups))
    indices = phylo.index_trees(trees, labels)
    rows, summary = rosenberg.analyse(indices, assignment, list(groups), sizes)
    return {row["group"]: row for row in rows}, summary, labels


def test_cummings_figure_1_only_group_a_is_monophyletic():
    """Group a is a clade; group b is paraphyletic, so P_A applies only to a."""
    rows, summary, labels = analyse_fixture("cummings_fig1.nwk")

    assert rows["a"]["monophyletic"] is True
    close(rows["a"]["p_a"], P_A(4, 4), tolerance=1e-12)
    assert rows["b"]["monophyletic"] is False
    assert rows["b"]["n_clade_breaking"] == 4      # all four a tips sit inside b's MRCA
    # P_A is conditional on observing monophyly, so it is withheld for group b
    # rather than reported as a number someone could quote.
    assert rows["b"]["p_a"] is None
    assert rows["b"]["log10_p_a"] is None
    assert rows["b"]["min_b_for_alpha"] is None

    # Not all groups are clades, so neither reciprocal nor joint applies.
    assert summary["reciprocal"] is None
    assert summary["joint"] is None
    assert summary["all_monophyletic"] is False


def test_reciprocal_and_joint_fire_when_both_groups_are_clades():
    trees = dendropy.TreeList.get(data="(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));",
                                  schema="newick", rooting="force-rooted")
    labels, assignment, sizes = phylo.build_assignment(trees, PAPER_GROUPS, ["a", "b"])
    rows, summary = rosenberg.analyse(phylo.index_trees(trees, labels), assignment, ["a", "b"], sizes)

    assert summary["all_monophyletic"] is True
    close(summary["reciprocal"]["p_ab"], P_AB(4, 4), tolerance=1e-12)
    close(summary["joint"]["p"], P_AB(4, 4), tolerance=1e-12)   # k=2 joint == P_AB


def test_three_monophyletic_groups_give_the_joint_probability():
    newick = "(((x1,x2),(x3,x4)),((y1,y2),(z1,z2)));"
    mapping = {"x1": "x", "x2": "x", "x3": "x", "x4": "x",
               "y1": "y", "y2": "y", "z1": "z", "z2": "z"}
    trees = dendropy.TreeList.get(data=newick, schema="newick", rooting="force-rooted")
    labels, assignment, sizes = phylo.build_assignment(trees, mapping, ["x", "y", "z"])
    rows, summary = rosenberg.analyse(phylo.index_trees(trees, labels), assignment, ["x", "y", "z"], sizes)

    assert summary["all_monophyletic"] is True
    close(summary["joint"]["p"], p_joint([4, 2, 2]), tolerance=1e-15)
    # Every group gets its own P_A against all other tips, Holm adjusted.
    close(rows[0]["p_a"], P_A(4, 4), tolerance=1e-12)
    close(rows[1]["p_a"], P_A(2, 6), tolerance=1e-12)
    assert all(row["p_a_holm"] >= row["p_a"] for row in rows)


def test_joint_requires_an_exhaustive_grouping():
    """With an unlabeled tip the groups no longer partition the tree."""
    mapping = {k: v for k, v in PAPER_GROUPS.items() if k != "b4"}
    trees = dendropy.TreeList.get(data="(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));",
                                  schema="newick", rooting="force-rooted")
    labels, assignment, sizes = phylo.build_assignment(trees, mapping, ["a", "b"])
    rows, summary = rosenberg.analyse(phylo.index_trees(trees, labels), assignment, ["a", "b"], sizes)

    assert summary["exhaustive"] is False
    assert summary["joint"] is None
    assert summary["reciprocal"] is None
    # P_A still applies: its b lineages may be of arbitrary origin.
    close(rows[0]["p_a"], P_A(4, 4), tolerance=1e-12)


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
