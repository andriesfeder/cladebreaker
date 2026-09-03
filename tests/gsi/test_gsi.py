#!/usr/bin/env python

"""
Tests for bin/gsi.py.

The two worked examples in Cummings, Neel & Shaw (2008) are the reference
values: Fig. 1 (fully resolved) and Fig. 2 (with a polytomy) both give
gsi = 0.563 for the paraphyletic group b and gsi = 1 for the monophyletic
group a, with min(gs) = 3/7 on either tree.

Run with pytest, or directly:  python tests/gsi/test_gsi.py
"""

import importlib.util
import itertools
import os
import sys

import dendropy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

spec = importlib.util.spec_from_file_location("gsi", os.path.join(ROOT, "bin", "gsi.py"))
gsi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gsi)


PAPER_GROUPS = {
    "a1": "a", "a2": "a", "a3": "a", "a4": "a",
    "b1": "b", "b2": "b", "b3": "b", "b4": "b",
}


def load(newick):
    return dendropy.TreeList.get(
        data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True
    )


def load_fixture(name):
    with open(os.path.join(HERE, name)) as handle:
        return load(handle.read())


def rows_by_group(rows):
    return {row["group"]: row for row in rows}


def close(observed, expected, tolerance=1e-9):
    assert observed is not None, "expected {}, got NA".format(expected)
    assert abs(observed - expected) < tolerance, "expected {}, got {}".format(expected, observed)


# ---------------------------------------------------------------------------
# The paper's worked examples
# ---------------------------------------------------------------------------


def test_figure_1_fully_resolved():
    """Fig. 1: monophyletic a, paraphyletic b, on a fully resolved tree."""
    rows, _ = gsi.analyse(load_fixture("cummings_fig1.nwk"), PAPER_GROUPS, ["a", "b"])
    result = rows_by_group(rows)

    close(result["a"]["gs"], 1.0)
    close(result["a"]["gsi"], 1.0)

    close(result["b"]["gs"], 3.0 / 4.0)
    close(result["b"]["min_gs"], 3.0 / 7.0)
    close(result["b"]["gsi"], 0.5625)


def test_figure_2_polytomy():
    """Fig. 2: same gsi values once a polytomy replaces resolved structure."""
    rows, _ = gsi.analyse(load_fixture("cummings_fig2.nwk"), PAPER_GROUPS, ["a", "b"])
    result = rows_by_group(rows)

    close(result["a"]["gs"], 1.0)
    close(result["a"]["gsi"], 1.0)

    close(result["b"]["gs"], 3.0 / 4.0)
    close(result["b"]["min_gs"], 3.0 / 7.0)
    close(result["b"]["gsi"], 0.5625)


def test_root_carries_a_stem_branch():
    """min(gs) = 3/7 only holds if the root counts as degree children + 1.

    Treating the root as degree 2 (what a plain Newick reading gives) would
    make the denominator of eq. 3 six, not seven, and quietly inflate every
    gsi value on the tree.
    """
    index = gsi.TreeIndex(load_fixture("cummings_fig1.nwk")[0], sorted(PAPER_GROUPS))
    assert index.n_internal == 7
    assert index.total_degree_sum == 7
    close(gsi.min_gs(index, 4), 3.0 / 7.0)

    polytomous = gsi.TreeIndex(load_fixture("cummings_fig2.nwk")[0], sorted(PAPER_GROUPS))
    assert polytomous.n_internal == 4
    assert sorted(polytomous.degree) == [3, 3, 4, 5]
    assert polytomous.total_degree_sum == 7


def exact_null_p_value(tree_index, size, observed, n_tips):
    """Exact P by enumerating every labeling, for checking the sampler.

    Feasible only for a toy tree: there are C(n_tips, size) labelings.
    """
    at_least = 0
    total = 0
    for combination in itertools.combinations(range(n_tips), size):
        assignment = [0 if i in combination else 1 for i in range(n_tips)]
        gs = gsi.gs_values(tree_index, assignment, [size, n_tips - size])[0]
        value = gsi.normalize(gs, gsi.min_gs(tree_index, size))
        total += 1
        if value is not None and value >= observed - 1e-12:
            at_least += 1
    return at_least / total


def test_permutation_converges_to_the_exact_null():
    """The sampled P-value should land on the exhaustively enumerated one.

    On the Fig. 1 tree there are only C(8,4) = 70 distinct labelings, so the
    null can be enumerated outright: exactly one reaches gsi = 1 (group a's own
    labeling) and eight reach gsi >= 0.5625.

    The paper reports P = 0.0121 for group a, which agrees with 1/70 = 0.0143
    within Monte Carlo error. Its P = 0.0130 for group b does not reproduce:
    both groups share one null distribution here, so P(b) >= P(a) must hold,
    and the monophyletic labeling alone already clears group b's lower
    threshold. 8/70 is what the statistic as defined in the paper gives.
    """
    trees = load_fixture("cummings_fig1.nwk")
    index = gsi.TreeIndex(trees[0], sorted(PAPER_GROUPS))

    exact_a = exact_null_p_value(index, 4, 1.0, 8)
    exact_b = exact_null_p_value(index, 4, 0.5625, 8)
    close(exact_a, 1.0 / 70.0)
    close(exact_b, 8.0 / 70.0)

    rows, _ = gsi.analyse(trees, PAPER_GROUPS, ["a", "b"], permutations=9999, seed=1)
    result = rows_by_group(rows)

    close(result["a"]["p_value"], exact_a, tolerance=0.005)
    close(result["b"]["p_value"], exact_b, tolerance=0.015)


# ---------------------------------------------------------------------------
# Behaviour
# ---------------------------------------------------------------------------


def balanced_newick(depth, labels):
    """Balanced rooted tree over 2**depth tips, taking labels in order."""
    def build(level, offset):
        if level == 0:
            return labels[offset], offset + 1
        left, offset = build(level - 1, offset)
        right, offset = build(level - 1, offset)
        return "({},{})".format(left, right), offset

    newick, _ = build(depth, 0)
    return newick + ";"


def test_monophyly_scores_one():
    labels = ["t{}".format(i) for i in range(16)]
    groups = {label: ("a" if i < 8 else "b") for i, label in enumerate(labels)}
    rows, _ = gsi.analyse(load(balanced_newick(4, labels)), groups, ["a", "b"], permutations=999, seed=7)
    result = rows_by_group(rows)

    for group in ("a", "b"):
        close(result[group]["gsi"], 1.0)
        assert result[group]["p_value"] < 0.01


def test_fully_interleaved_labels_score_near_zero():
    """Maximally mixed ancestry is the floor of the statistic, so P is high."""
    labels = ["t{}".format(i) for i in range(16)]
    groups = {label: ("a" if i % 2 == 0 else "b") for i, label in enumerate(labels)}
    rows, _ = gsi.analyse(load(balanced_newick(4, labels)), groups, ["a", "b"], permutations=999, seed=7)
    result = rows_by_group(rows)

    for group in ("a", "b"):
        assert result[group]["gsi"] < 0.05, result[group]["gsi"]
        assert result[group]["p_value"] > 0.5, result[group]["p_value"]


def test_single_tip_group_is_not_a_number():
    groups = dict(PAPER_GROUPS)
    groups["b1"] = "c"
    groups["b2"] = "c"
    groups["b3"] = "c"
    groups["b4"] = "c"
    groups["a4"] = "solo"
    rows, _ = gsi.analyse(load_fixture("cummings_fig1.nwk"), groups, ["a", "c", "solo"])
    result = rows_by_group(rows)

    assert result["solo"]["n_tips"] == 1
    assert result["solo"]["gsi"] is None
    assert result["solo"]["p_value"] is None


def test_unlabeled_tips_are_excluded():
    """A tip left out of the groups file stays on the tree but takes no label."""
    groups = {tip: group for tip, group in PAPER_GROUPS.items() if tip != "b4"}
    rows, _ = gsi.analyse(load_fixture("cummings_fig1.nwk"), groups, ["a", "b"])
    result = rows_by_group(rows)

    assert result["b"]["n_tips"] == 3
    # b1-b3 unite through three nodes of degree 3, so gs = 2/3.
    close(result["b"]["gs"], 2.0 / 3.0)
    close(result["b"]["min_gs"], 2.0 / 7.0)


def test_ensemble_of_identical_trees_matches_the_single_tree():
    single = gsi.analyse(load_fixture("cummings_fig1.nwk"), PAPER_GROUPS, ["a", "b"])[0]
    with open(os.path.join(HERE, "cummings_fig1.nwk")) as handle:
        newick = handle.read().strip()
    ensemble = gsi.analyse(load(newick + "\n" + newick + "\n"), PAPER_GROUPS, ["a", "b"])[0]

    assert ensemble[0]["n_trees"] == 2
    for one, many in zip(single, ensemble):
        close(many["gsi"], one["gsi"])


def test_ensemble_averages_over_topologies():
    """gsi_T is the mean of the per-tree values when topologies are equally weighted."""
    with open(os.path.join(HERE, "cummings_fig1.nwk")) as handle:
        resolved = handle.read().strip()
    monophyletic = "(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));"

    rows, per_tree = gsi.analyse(load(resolved + "\n" + monophyletic + "\n"), PAPER_GROUPS, ["a", "b"])
    result = rows_by_group(rows)

    close(result["b"]["gsi"], (0.5625 + 1.0) / 2.0)
    assert result["b"]["gs"] is None, "gs is per-tree only; the ensemble reports gsi_T"
    assert len(per_tree[1]) == 2


def deepest_leaf_from_root(tree):
    """Distance to the furthest leaf, excluding the root's own stem branch."""
    stem = tree.seed_node.edge.length or 0.0
    return max(leaf.distance_from_root() for leaf in tree.leaf_node_iter()) - stem


def test_midpoint_rooting_halves_the_longest_path():
    """The furthest leaf must sit at exactly half the tree's diameter."""
    tree = load("((A:1,B:1):1,(C:1,D:9):1);")[0]
    gsi.reroot_at_midpoint(tree)
    # Longest path is D..A = 9 + 1 + 1 + 1 = 12, so the root lands 6 along it,
    # which is 3 up D's own branch.
    close(deepest_leaf_from_root(tree), 6.0)


def test_midpoint_rooting_on_an_ultrametric_tree():
    """The midpoint is an existing node here, which dendropy's own version trips on."""
    tree = load("(((A:1,B:1):1,C:2):1,((D:1,E:1):1,F:2):1);")[0]
    gsi.reroot_at_midpoint(tree)
    close(deepest_leaf_from_root(tree), 3.0)


def test_midpoint_rooting_without_branch_lengths_is_rejected():
    tree = load_fixture("cummings_fig1.nwk")[0]
    try:
        gsi.reroot_at_midpoint(tree)
    except SystemExit as error:
        assert "branch lengths" in str(error)
    else:
        raise AssertionError("expected midpoint rooting to be refused without branch lengths")


def test_permutation_is_reproducible_for_a_seed():
    first = gsi.analyse(load_fixture("cummings_fig1.nwk"), PAPER_GROUPS, ["a", "b"], 499, seed=11)[0]
    second = gsi.analyse(load_fixture("cummings_fig1.nwk"), PAPER_GROUPS, ["a", "b"], 499, seed=11)[0]
    assert [row["p_value"] for row in first] == [row["p_value"] for row in second]


# ---------------------------------------------------------------------------


def main():
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
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
