#!/usr/bin/env python

"""
Tests for bin/gsi_groups.py, which resolves tree tip labels against the sample
IDs the pipeline knows about.

Run with pytest, or directly:  python tests/gsi/test_gsi_groups.py
"""

import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

spec = importlib.util.spec_from_file_location("gsi_groups", os.path.join(ROOT, "bin", "gsi_groups.py"))
gsi_groups = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gsi_groups)


def assign(tips, query, reference, allow_unmatched=False):
    return gsi_groups.assign(
        tips, [("query", query), ("reference", reference)], allow_unmatched
    )


def test_exact_labels():
    assignments, unmatched = assign(
        ["SAMPLE_1", "GCA_000012345.1"], ["SAMPLE_1"], ["GCA_000012345.1"]
    )
    assert assignments == {"SAMPLE_1": "query", "GCA_000012345.1": "reference"}
    assert unmatched == []


def test_pipeline_suffixes_are_tolerated():
    """Annotation and pangenome steps append to the sample name."""
    assignments, _ = assign(
        ["SAMPLE_1", "GCA_000012345.1_genomic"],
        ["SAMPLE_1_T1"],
        ["GCA_000012345.1"],
    )
    assert assignments == {"SAMPLE_1": "query", "GCA_000012345.1_genomic": "reference"}


def test_numeric_ids_do_not_bleed_into_each_other():
    """Path_Reg_5 must not claim Path_Reg_50's tip: 0 is not a separator."""
    assignments, unmatched = assign(
        ["Path_Reg_5", "Path_Reg_50"], ["Path_Reg_5", "Path_Reg_50"], [], allow_unmatched=True
    )
    assert assignments == {"Path_Reg_5": "query", "Path_Reg_50": "query"}

    # The same boundary rule, across groups, where a wrong match would matter.
    assignments, _ = assign(["Path_Reg_50"], ["Path_Reg_5"], ["Path_Reg_50"])
    assert assignments == {"Path_Reg_50": "reference"}


def test_a_tip_matching_both_groups_is_refused():
    try:
        assign(["SAMPLE_1_T1"], ["SAMPLE_1"], ["SAMPLE_1_T1_extra"])
    except SystemExit as error:
        assert "more than one group" in str(error)
    else:
        raise AssertionError("expected an ambiguous tip to be refused")


def test_exact_match_wins_over_a_suffix_match():
    """An exact hit settles it, even when another ID also extends the tip."""
    assignments, _ = assign(["SAMPLE_1"], ["SAMPLE_1"], ["SAMPLE_1_backup"])
    assert assignments == {"SAMPLE_1": "query"}


def test_unmatched_tips_are_refused_by_default():
    try:
        assign(["SAMPLE_1", "MYSTERY"], ["SAMPLE_1"], [])
    except SystemExit as error:
        assert "match neither" in str(error)
    else:
        raise AssertionError("expected an unmatched tip to be refused")


def test_unmatched_tips_can_be_allowed():
    assignments, unmatched = assign(["SAMPLE_1", "MYSTERY"], ["SAMPLE_1"], [], allow_unmatched=True)
    assert assignments == {"SAMPLE_1": "query"}
    assert unmatched == ["MYSTERY"]


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
