#!/usr/bin/env python

"""
Build a gsi groups file by matching tree tips against lists of sample IDs.

Inside the pipeline the query isolates and the reference genomes pulled from
NCBI are already known, but the tip labels on the tree come from whatever the
pangenome or SNP caller wrote out, which may carry a suffix (a .gff basename, a
_genomic tag, Prokka's locus prefix). This resolves the two against each other
and refuses to guess when a tip is ambiguous.

Example usage:
    gsi_groups.py --tree tree.nwk --query-ids queries.txt \\
        --reference-ids references.txt --out groups.tsv
"""

import argparse
import os
import sys

import dendropy


def read_ids(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def tip_labels(path):
    tree = dendropy.Tree.get(
        path=path, schema="newick", rooting="force-rooted", preserve_underscores=True
    )
    return [
        node.taxon.label if node.taxon is not None else node.label
        for node in tree.leaf_node_iter()
    ]


SEPARATORS = ("_", ".", "-", "|")


def extends(longer, shorter):
    """True if `longer` is `shorter` plus a suffix starting at a separator.

    The separator matters: sample IDs in this pipeline are numbered, so a bare
    prefix test would let Path_Reg_5 claim Path_Reg_50's tip.
    """
    return longer.startswith(shorter) and longer[len(shorter):].startswith(SEPARATORS)


def candidates(tip, identifiers):
    """Identifiers that could name this tip, best match first.

    Exact first; then either side extending the other, which covers the
    suffixes the annotation and pangenome steps add (_T1, _genomic, a locus tag).
    """
    exact = [name for name in identifiers if name == tip]
    if exact:
        return exact
    return [name for name in identifiers if extends(tip, name) or extends(name, tip)]


def assign(tips, groups, allow_unmatched):
    """Map each tip to a group label, or report why it could not be mapped."""
    assignments = {}
    unmatched = []
    ambiguous = []

    for tip in tips:
        hits = []
        for label, identifiers in groups:
            for name in candidates(tip, identifiers):
                hits.append((label, name))

        if not hits:
            unmatched.append(tip)
            continue

        exact = [hit for hit in hits if hit[1] == tip]
        if exact:
            hits = exact

        labels = {label for label, _ in hits}
        if len(labels) > 1:
            ambiguous.append((tip, sorted(labels)))
            continue

        assignments[tip] = hits[0][0]

    if ambiguous:
        sys.exit(
            "ERROR: {} tip(s) match more than one group, so the grouping is not "
            "well defined:\n{}".format(
                len(ambiguous),
                "\n".join("  {} -> {}".format(tip, ", ".join(labels)) for tip, labels in ambiguous[:10]),
            )
        )

    if unmatched and not allow_unmatched:
        sys.exit(
            "ERROR: {} tree tip(s) match neither the query nor the reference IDs: {}\n"
            "Pass --allow-unmatched to leave them out of the statistic, or check that "
            "the tip labels correspond to the sample IDs.".format(
                len(unmatched), ", ".join(unmatched[:10]) + (" ..." if len(unmatched) > 10 else "")
            )
        )

    return assignments, unmatched


def parse_args(args=None):
    description = "Write a gsi groups file by matching tree tips to query and reference sample IDs."
    epilog = "Example usage: gsi_groups.py --tree tree.nwk --query-ids q.txt --reference-ids r.txt --out groups.tsv"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("--tree", required=True, help="Newick tree whose tips are to be labeled.")
    parser.add_argument("--query-ids", required=True, help="File of query sample IDs, one per line.")
    parser.add_argument("--reference-ids", help="File of reference genome IDs, one per line.")
    parser.add_argument("--query-label", default="query", help="Group name for the query isolates (default: query).")
    parser.add_argument(
        "--reference-label", default="reference", help="Group name for the reference genomes (default: reference)."
    )
    parser.add_argument(
        "--allow-unmatched", action="store_true", help="Leave tips that match no ID out of the groups file."
    )
    parser.add_argument("--out", default="groups.tsv", help="Output groups file (default: groups.tsv).")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    query = read_ids(args.query_ids)
    reference = read_ids(args.reference_ids)
    if not query:
        sys.exit("ERROR: No query sample IDs were found in {}".format(args.query_ids))

    tips = tip_labels(args.tree)
    groups = [(args.query_label, query), (args.reference_label, reference)]
    assignments, unmatched = assign(tips, groups, args.allow_unmatched)

    counts = {}
    with open(args.out, "w") as handle:
        handle.write("tip_label\tgroup\n")
        for tip in tips:
            if tip in assignments:
                handle.write("{}\t{}\n".format(tip, assignments[tip]))
                counts[assignments[tip]] = counts.get(assignments[tip], 0) + 1

    if len(counts) < 2:
        sys.exit(
            "ERROR: gsi compares labeled groups, but the tips resolved to only "
            "{}. Check that the reference genomes are present on the tree.".format(
                "one group ({})".format(list(counts)[0]) if counts else "no groups"
            )
        )

    print("Wrote {} with:".format(args.out))
    for label, count in sorted(counts.items()):
        print("  {:<20} {} tip(s)".format(label, count))
    if unmatched:
        print("  {} tip(s) left unlabeled: {}".format(
            len(unmatched), ", ".join(unmatched[:10]) + (" ..." if len(unmatched) > 10 else "")
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
