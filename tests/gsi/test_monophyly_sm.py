#!/usr/bin/env python

"""
Tests for bin/monophyly.py and bin/slatkin_maddison.py.

The Slatkin-Maddison tests include the claim both scripts' docstrings make: that
the Fitch change count is a property of the unrooted topology, so rerooting
cannot alter it. That is what makes it the safe fallback when rooting is the
uncertain part, so it is asserted rather than assumed -- by rerooting a tree at
every single edge and checking s never moves.

Run with pytest, or directly:  python tests/gsi/test_monophyly_sm.py
"""

import importlib.util
import os
import sys

import dendropy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "bin"))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "bin", name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mono = load("monophyly")
sm = load("slatkin_maddison")
import cladebreaker_phylo as phylo  # noqa: E402

PAPER_GROUPS = {"a1": "a", "a2": "a", "a3": "a", "a4": "a",
                "b1": "b", "b2": "b", "b3": "b", "b4": "b"}


def trees_from(data):
    return dendropy.TreeList.get(data=data, schema="newick", rooting="force-rooted",
                                 preserve_underscores=True)


def fixture(name):
    with open(os.path.join(HERE, name)) as handle:
        return trees_from(handle.read())


def setup(trees, mapping, groups):
    labels, assignment, sizes = phylo.build_assignment(trees, mapping, list(groups))
    return phylo.index_trees(trees, labels), assignment, list(groups), sizes, labels


# ---------------------------------------------------------------------------
# monophyly.py
# ---------------------------------------------------------------------------


def test_figure_1_group_a_is_a_clade_and_b_is_not():
    indices, assignment, groups, sizes, labels = setup(fixture("cummings_fig1.nwk"), PAPER_GROUPS, ("a", "b"))
    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}

    assert rows["a"]["monophyletic"] is True
    assert rows["a"]["n_clusters"] == 1
    assert rows["a"]["n_clade_breaking"] == 0
    assert rows["a"]["mrca_clade_size"] == 4

    assert rows["b"]["monophyletic"] is False
    # b's MRCA is the whole tree, so all four a tips break it.
    assert rows["b"]["mrca_clade_size"] == 8
    assert sorted(labels[i] for i in rows["b"]["clade_breaking"]) == ["a1", "a2", "a3", "a4"]
    # b1..b4 attach one after another along the spine: four separate blocks.
    assert rows["b"]["n_clusters"] == 4


def test_both_groups_clades():
    indices, assignment, groups, sizes, _ = setup(
        trees_from("(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));"), PAPER_GROUPS, ("a", "b"))
    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}
    for group in ("a", "b"):
        assert rows[group]["monophyletic"] is True
        assert rows[group]["n_clusters"] == 1


def test_cluster_count_finds_scattered_blocks():
    """Group a sits in two places, so two clusters even though four tips."""
    indices, assignment, groups, sizes, _ = setup(
        trees_from("(((a1,a2),(b1,b2)),((a3,a4),(b3,b4)));"), PAPER_GROUPS, ("a", "b"))
    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}
    assert rows["a"]["monophyletic"] is False
    assert rows["a"]["n_clusters"] == 2
    assert rows["b"]["n_clusters"] == 2


def test_fully_interleaved_group_is_all_singletons():
    labels = ["t%d" % i for i in range(8)]
    mapping = {t: ("a" if i % 2 == 0 else "b") for i, t in enumerate(labels)}
    newick = "(((t0,t1),(t2,t3)),((t4,t5),(t6,t7)));"
    indices, assignment, groups, sizes, _ = setup(trees_from(newick), mapping, ("a", "b"))
    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}
    assert rows["a"]["n_clusters"] == 4      # no two a tips are ever sisters
    assert rows["b"]["n_clusters"] == 4


def test_monophyly_fraction_across_an_ensemble():
    """One tree of two has group a as a clade, so support is 0.5."""
    both = "(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));"
    broken = "((((a1,a2),(a3,b1)),a4),((b2,b3),b4));"
    indices, assignment, groups, sizes, _ = setup(trees_from(both + "\n" + broken), PAPER_GROUPS, ("a", "b"))
    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}
    assert rows["a"]["monophyletic"] is False       # not a clade in every tree
    assert abs(rows["a"]["monophyly_fraction"] - 0.5) < 1e-12


def test_single_tip_group_is_trivially_a_clade():
    mapping = dict(PAPER_GROUPS)
    mapping["a4"] = "solo"
    indices, assignment, groups, sizes, _ = setup(
        fixture("cummings_fig1.nwk"), mapping, ("a", "b", "solo"))
    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}
    assert rows["solo"]["monophyletic"] is True
    assert rows["solo"]["n_clusters"] == 1


# ---------------------------------------------------------------------------
# slatkin_maddison.py -- Fitch counts
# ---------------------------------------------------------------------------


def changes(newick, mapping, groups):
    indices, assignment, _, sizes, _ = setup(trees_from(newick), mapping, groups)
    return sm.fitch_changes(indices[0], assignment, len(sizes))


def test_two_clean_clades_need_one_change():
    assert changes("(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));", PAPER_GROUPS, ("a", "b")) == 1


def test_the_caterpillar_from_figure_1():
    """a is a clade, b1..b4 hang off the spine: still a single change."""
    with open(os.path.join(HERE, "cummings_fig1.nwk")) as handle:
        assert changes(handle.read(), PAPER_GROUPS, ("a", "b")) == 1


def test_two_blocks_per_group_need_two_changes():
    """Four alternating blocks, but only two changes: paint the tree a, then pay
    once for each of the two b blocks. Block count minus one would be wrong."""
    assert changes("(((a1,a2),(b1,b2)),((a3,a4),(b3,b4)));", PAPER_GROUPS, ("a", "b")) == 2


def test_fully_interleaved_labels_need_the_most_changes():
    labels = ["t%d" % i for i in range(8)]
    mapping = {t: ("a" if i % 2 == 0 else "b") for i, t in enumerate(labels)}
    assert changes("(((t0,t1),(t2,t3)),((t4,t5),(t6,t7)));", mapping, ("a", "b")) == 4


def test_three_groups_need_two_changes_when_each_is_a_clade():
    mapping = {"x1": "x", "x2": "x", "y1": "y", "y2": "y", "z1": "z", "z2": "z"}
    assert changes("(((x1,x2),(y1,y2)),(z1,z2));", mapping, ("x", "y", "z")) == 2


def test_polytomies_are_handled():
    """A star of one group plus one outsider is a single change."""
    mapping = {"a1": "a", "a2": "a", "a3": "a", "b1": "b"}
    assert changes("(a1,a2,a3,b1);", mapping, ("a", "b")) == 1
    # A star split 2/2 costs two: the single internal node holds one state, so
    # whichever pair it does not match each needs a change.
    mapping2 = {"a1": "a", "a2": "a", "b1": "b", "b2": "b"}
    assert changes("(a1,a2,b1,b2);", mapping2, ("a", "b")) == 2


def test_unlabeled_tips_are_free():
    """A wildcard tip must not force a change."""
    mapping = {"a1": "a", "a2": "a", "b1": "b", "b2": "b"}      # x1 unlabeled
    newick = "(((a1,a2),x1),(b1,b2));"
    trees = trees_from(newick)
    labels, assignment, sizes = phylo.build_assignment(trees, mapping, ["a", "b"])
    indices = phylo.index_trees(trees, labels)
    assert sm.fitch_changes(indices[0], assignment, 2) == 1


def test_fitch_matches_an_exhaustive_minimum():
    """Cross-check Fitch against brute force over every internal-node labeling.

    Independent of any hand arithmetic, and it covers the polytomy and
    three-state generalisations where the rule is easiest to get wrong.
    """
    from itertools import product

    def brute_force(tree_index, assignment, n_groups):
        best = None
        for states in product(range(n_groups), repeat=tree_index.n_internal):
            cost = 0
            for i, kids in enumerate(tree_index.children):
                for is_leaf, j in kids:
                    if is_leaf:
                        g = assignment[tree_index.leaf_canonical[j]]
                        if g >= 0 and g != states[i]:
                            cost += 1
                    elif states[j] != states[i]:
                        cost += 1
            best = cost if best is None else min(best, cost)
        return best

    three = {"x1": "x", "x2": "x", "y1": "y", "y2": "y", "z1": "z", "z2": "z"}
    interleaved = {"t%d" % i: ("a" if i % 2 == 0 else "b") for i in range(8)}
    cases = [
        ("(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));", PAPER_GROUPS, ["a", "b"]),
        ("(((a1,a2),(b1,b2)),((a3,a4),(b3,b4)));", PAPER_GROUPS, ["a", "b"]),
        ("(((a1,a2),(a3,(a4,b1))),((b2,b3),b4));", PAPER_GROUPS, ["a", "b"]),
        ("(a1,a2,b1,b2);", {"a1": "a", "a2": "a", "b1": "b", "b2": "b"}, ["a", "b"]),
        ("(a1,a2,a3,b1);", {"a1": "a", "a2": "a", "a3": "a", "b1": "b"}, ["a", "b"]),
        ("(((x1,x2),(y1,y2)),(z1,z2));", three, ["x", "y", "z"]),
        ("((x1,y1,z1),(x2,y2,z2));", three, ["x", "y", "z"]),
        ("(((t0,t1),(t2,t3)),((t4,t5),(t6,t7)));", interleaved, ["a", "b"]),
    ]
    for newick, mapping, groups in cases:
        indices, assignment, _, sizes, _ = setup(trees_from(newick), mapping, groups)
        assert sm.fitch_changes(indices[0], assignment, len(sizes)) == \
            brute_force(indices[0], assignment, len(sizes)), newick


def test_fitch_count_is_independent_of_rooting():
    """The claim both docstrings make, checked at every edge of a real tree."""
    with open(os.path.join(HERE, "cummings_fig1.nwk")) as handle:
        newick = handle.read()
    baseline = changes(newick, PAPER_GROUPS, ("a", "b"))

    seen = set()
    source = trees_from(newick)[0]
    for index in range(len(list(source.postorder_edge_iter()))):
        tree = trees_from(newick)[0]
        edges = [e for e in tree.postorder_edge_iter() if e.tail_node is not None]
        if index >= len(edges):
            break
        tree.reroot_at_edge(edges[index], update_bipartitions=False, suppress_unifurcations=True)
        rerooted = dendropy.TreeList([tree])
        labels, assignment, sizes = phylo.build_assignment(rerooted, PAPER_GROUPS, ["a", "b"])
        indices = phylo.index_trees(rerooted, labels)
        seen.add(sm.fitch_changes(indices[0], assignment, 2))

    assert seen == {baseline}, "rerooting changed s: saw {}".format(sorted(seen))


# ---------------------------------------------------------------------------
# slatkin_maddison.py -- the test itself
# ---------------------------------------------------------------------------


def test_structured_groups_are_significant():
    indices, assignment, groups, sizes, _ = setup(
        trees_from("(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));"), PAPER_GROUPS, ("a", "b"))
    result = sm.analyse(indices, assignment, groups, sizes, permutations=999, seed=3)
    assert result["s_observed"] == 1
    assert result["s_minimum_possible"] == 1
    assert result["p_value"] < 0.05
    assert result["null_mean"] > result["s_observed"]


def test_interleaved_groups_are_not_significant():
    labels = ["t%d" % i for i in range(8)]
    mapping = {t: ("a" if i % 2 == 0 else "b") for i, t in enumerate(labels)}
    indices, assignment, groups, sizes, _ = setup(
        trees_from("(((t0,t1),(t2,t3)),((t4,t5),(t6,t7)));"), mapping, ("a", "b"))
    result = sm.analyse(indices, assignment, groups, sizes, permutations=999, seed=3)
    assert result["p_value"] > 0.5, result["p_value"]


def balanced_newick(depth, labels):
    def build(level, offset):
        if level == 0:
            return labels[offset], offset + 1
        left, offset = build(level - 1, offset)
        right, offset = build(level - 1, offset)
        return "({},{})".format(left, right), offset
    newick, _ = build(depth, 0)
    return newick + ";"


def test_the_broken_case_rosenberg_cannot_handle():
    """Neither group is a clade, yet the mixing is far from random.

    This is the branch the decision path exists for: one misplaced tip destroys
    monophyly for both groups, so Rosenberg declines, but Slatkin-Maddison still
    detects the structure.
    """
    labels = ["t%d" % i for i in range(32)]
    mapping = {t: ("a" if i < 16 else "b") for i, t in enumerate(labels)}
    mapping["t0"] = "b"                       # one intruder inside the a half
    indices, assignment, groups, sizes, _ = setup(
        trees_from(balanced_newick(5, labels)), mapping, ("a", "b"))

    rows = {r["group"]: r for r in mono.analyse(indices, assignment, groups, sizes)}
    assert rows["a"]["monophyletic"] is False
    assert rows["b"]["monophyletic"] is False

    result = sm.analyse(indices, assignment, groups, sizes, permutations=999, seed=5)
    assert result["s_observed"] == 2
    assert result["p_value"] < 0.01, result["p_value"]
    assert result["null_mean"] > 8


def test_small_trees_lack_power():
    """An 8-tip tree cannot reject even when the labeling looks structured.

    Worth pinning: it is the same sample-size limitation Rosenberg's tables make
    explicit, and it stops a non-significant result on a tiny tree being read as
    evidence of mixing.
    """
    indices, assignment, groups, sizes, _ = setup(
        trees_from("(((a1,a2),(a3,(a4,b1))),((b2,b3),b4));"), PAPER_GROUPS, ("a", "b"))
    result = sm.analyse(indices, assignment, groups, sizes, permutations=9999, seed=5)
    assert result["s_observed"] == 2
    assert result["p_value"] > 0.05, result["p_value"]


def test_permutation_is_reproducible_and_skippable():
    indices, assignment, groups, sizes, _ = setup(
        trees_from("(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));"), PAPER_GROUPS, ("a", "b"))
    first = sm.analyse(indices, assignment, groups, sizes, permutations=499, seed=11)
    second = sm.analyse(indices, assignment, groups, sizes, permutations=499, seed=11)
    assert first["p_value"] == second["p_value"]

    none = sm.analyse(indices, assignment, groups, sizes, permutations=0)
    assert none["p_value"] is None
    assert none["s_observed"] == 1


def test_ensemble_averages_the_change_count():
    clean = "(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));"
    broken = "(((a1,a2),(b1,b2)),((a3,a4),(b3,b4)));"
    indices, assignment, groups, sizes, _ = setup(trees_from(clean + "\n" + broken), PAPER_GROUPS, ("a", "b"))
    result = sm.analyse(indices, assignment, groups, sizes, permutations=0)
    assert result["s_per_tree"] == [1, 2]
    assert abs(result["s_observed"] - 1.5) < 1e-12


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
