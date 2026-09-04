#!/usr/bin/env python

"""
Tests for bin/analyze_pdf.py.

The Newick reader is the part worth pinning: it is small on purpose, because it
only ever reads the normalised tree monophyly.py writes. These tests hold it to
dendropy's reading of the same strings -- tip set, tip order and total branch
length -- across the fixtures, a real 77-tip tree, quoted labels, comments and
scientific notation.

Run with pytest, or directly:  python tests/gsi/test_analyze_pdf.py
"""

import importlib.util
import os
import sys

import dendropy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "bin"))

spec = importlib.util.spec_from_file_location("analyze_pdf", os.path.join(ROOT, "bin", "analyze_pdf.py"))
apdf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(apdf)


def fixture(name):
    with open(os.path.join(HERE, name)) as handle:
        return handle.read()


def dendropy_tips(newick):
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted",
                             preserve_underscores=True)
    tips = [n.taxon.label if n.taxon is not None else n.label for n in tree.leaf_node_iter()]
    length = sum(e.length or 0.0 for e in tree.postorder_edge_iter())
    return tips, length


def check_against_dendropy(newick):
    root = apdf.parse_newick(newick)
    mine = [tip.label for tip in apdf.leaves(root)]
    mine_length = sum(node.length or 0.0 for node in apdf.postorder(root))
    theirs, their_length = dendropy_tips(newick)
    assert mine == theirs, "tip order differs:\n  mine   {}\n  theirs {}".format(mine[:6], theirs[:6])
    assert abs(mine_length - their_length) < 1e-9, "branch length total differs"
    return root


def test_matches_dendropy_on_the_fixtures():
    check_against_dendropy(fixture("cummings_fig1.nwk"))
    check_against_dendropy(fixture("cummings_fig2.nwk"))       # a polytomy


def test_matches_dendropy_on_quoted_labels_and_comments():
    check_against_dendropy("(('a b':0.1,'c,d':0.2):0.3,e:0.4);")
    check_against_dendropy("((a:0.1,b:0.2)95:0.3[100],c:0.4);")  # support as label and as comment


def test_matches_dendropy_without_branch_lengths():
    check_against_dendropy("((a,b),(c,d));")


def test_matches_dendropy_on_scientific_notation():
    check_against_dendropy("((a:1.5e-3,b:2E-4):1e-2,c:0.5);")


def test_underscores_in_labels_are_preserved():
    """GCA accessions and Path_Reg style IDs must survive verbatim: plain Newick
    would read an underscore as a space."""
    root = apdf.parse_newick("((Path_Reg_22:0.1,GCA_000012345.1:0.2):0.3,SAMPLE_1_T1:0.4);")
    assert sorted(tip.label for tip in apdf.leaves(root)) == \
        ["GCA_000012345.1", "Path_Reg_22", "SAMPLE_1_T1"]


def test_malformed_newick_is_rejected():
    for bad in ["((a,b)", "", "   ", "(a,b))"]:
        try:
            apdf.parse_newick(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("expected {!r} to be rejected".format(bad))


def test_a_missing_terminating_semicolon_is_tolerated():
    """Plenty of tools emit a tree without the closing semicolon, and refusing
    to draw a figure over punctuation would be unhelpful."""
    root = apdf.parse_newick("((a:1,b:1):1,c:2)")
    assert sorted(tip.label for tip in apdf.leaves(root)) == ["a", "b", "c"]


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def test_layout_places_leaves_evenly_and_parents_at_their_midpoint():
    root = apdf.parse_newick("(((a:1,b:1):1,c:2):1,d:3);")
    tips, used_lengths = apdf.layout(root)
    assert used_lengths is True
    assert [tip.y for tip in tips] == [0.0, 1.0, 2.0, 3.0]
    # x is cumulative branch length from the root
    depth = {tip.label: tip.x for tip in tips}
    assert abs(depth["a"] - 3.0) < 1e-9      # 1 + 1 + 1
    assert abs(depth["d"] - 3.0) < 1e-9
    ab = [n for n in apdf.postorder(root) if not n.is_leaf and len(n.children) == 2
          and {c.label for c in n.children} == {"a", "b"}][0]
    assert abs(ab.y - 0.5) < 1e-9            # midway between its two tips


def test_layout_falls_back_to_depth_without_branch_lengths():
    root = apdf.parse_newick("((a,b),(c,d));")
    tips, used_lengths = apdf.layout(root)
    assert used_lengths is False
    assert all(abs(tip.x - 2.0) < 1e-9 for tip in tips)   # every tip two steps deep


def test_layout_handles_a_polytomy():
    root = apdf.parse_newick(fixture("cummings_fig2.nwk"))
    tips, _ = apdf.layout(root)
    assert len(tips) == 8
    star = [n for n in apdf.postorder(root) if len(n.children) == 4][0]
    assert abs(star.y - sum(c.y for c in star.children) / 4.0) < 1e-9


# ---------------------------------------------------------------------------
# Groups and colour policy
# ---------------------------------------------------------------------------


def test_reads_the_groups_file():
    mapping = apdf.read_groups(os.path.join(HERE, "cummings_groups.tsv"))
    assert mapping["a1"] == "a" and mapping["b4"] == "b"
    assert len(mapping) == 8          # the header row is not an entry


def test_colour_is_capped_at_three_groups():
    """A fourth categorical slot puts yellow and orange together, which falls
    below the normal-vision separation floor. Past three groups the renderer
    drops colour rather than inventing a hue; marker shape and the tip labels
    carry identity instead."""
    assert apdf.MAX_COLOURED_GROUPS == 3
    assert len(apdf.GROUP_COLOURS) == 3
    assert len(apdf.GROUP_MARKERS) >= 8      # enough shapes to outlast the colours


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
