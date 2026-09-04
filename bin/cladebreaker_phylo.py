#!/usr/bin/env python

"""
Shared phylogenetic plumbing for the cladebreaker analyze tools.

Tree loading, the three rooting modes, the groups file, the flattened tree
index, monophyly, and the label-permutation driver. The statistics themselves
live in the tool scripts that import this: gsi.py, rosenberg.py, and so on.

Nextflow puts bin/ on PATH, not PYTHONPATH, so each tool needs

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

before importing this module. The whole bin/ directory is staged into the task
working directory, so the import resolves there.
"""

import random
import sys

import dendropy


# ---------------------------------------------------------------------------
# Tree indexing
#
# The dendropy tree is flattened into plain arrays once, then all of the
# permutation work runs over those arrays. Internal nodes are stored in
# post-order, so a node always appears after its children.
# ---------------------------------------------------------------------------


class TreeIndex(object):
    """Flattened rooted topology: internal nodes in post-order plus degrees."""

    def __init__(self, tree, canonical_labels):
        label_to_canonical = {label: i for i, label in enumerate(canonical_labels)}

        self.children = []       # per internal node: list of (is_leaf, index)
        self.degree = []         # per internal node: branch ends at the node
        self.leaf_canonical = []  # per local leaf: index into canonical_labels

        leaf_index = {}
        internal_index = {}

        for node in tree.postorder_node_iter():
            if node.is_leaf():
                label = node.taxon.label if node.taxon is not None else node.label
                if label not in label_to_canonical:
                    raise KeyError(label)
                leaf_index[id(node)] = len(self.leaf_canonical)
                self.leaf_canonical.append(label_to_canonical[label])
                continue

            kids = []
            for child in node.child_node_iter():
                if child.is_leaf():
                    kids.append((True, leaf_index[id(child)]))
                else:
                    kids.append((False, internal_index[id(child)]))

            internal_index[id(node)] = len(self.children)
            self.children.append(kids)
            # + 1 for the parent branch; the root carries a stem branch too.
            self.degree.append(len(kids) + 1)

        self.n_internal = len(self.children)
        self.n_leaves = len(self.leaf_canonical)
        # Denominator of eq. 3: every internal node in the tree.
        self.total_degree_sum = sum(d - 2 for d in self.degree)

        # Tips below each internal node, in the same post-order. Lets a caller
        # tell a clade made only of group members from one that is mixed.
        self.leaf_total = [0] * self.n_internal
        for i, kids in enumerate(self.children):
            self.leaf_total[i] = sum(1 if is_leaf else self.leaf_total[j] for is_leaf, j in kids)


def index_trees(trees, canonical_labels):
    return [TreeIndex(tree, canonical_labels) for tree in trees]

# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------


def read_groups(path):
    """Read a two-column tip-to-group table (TSV or CSV, optional header)."""
    mapping = {}
    order = []
    with open(path, "r") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split(",")
            fields = [field.strip() for field in fields if field.strip()]
            if len(fields) < 2:
                sys.exit(
                    "ERROR: Please check the groups file -> line {} does not have two "
                    "columns:\n{}".format(number, line)
                )
            tip, group = fields[0], fields[1]
            if number == 1 and tip.lower() in ("tip", "tip_label", "taxon", "sample", "label"):
                continue
            if tip in mapping:
                sys.exit("ERROR: Please check the groups file -> tip '{}' listed twice.".format(tip))
            mapping[tip] = group
            if group not in order:
                order.append(group)

    if not mapping:
        sys.exit("ERROR: Please check the groups file -> no tip/group assignments found.")
    if len(order) < 2:
        sys.exit(
            "ERROR: Please check the groups file -> gsi compares labeled groups, but only "
            "one group ('{}') was found.".format(order[0])
        )
    return mapping, order


def reroot_at_midpoint(tree):
    """Root a tree at the midpoint of its longest leaf-to-leaf path.

    dendropy has reroot_at_midpoint(), but as of 5.0.13 it walks up from only
    one end of the diameter and asserts if it reaches the MRCA with distance
    left over, which floating-point drift makes happen on larger trees (its
    _tree.py, "ugly, ugly, ugly code to find two nodes that span the
    midpoint"). Rerooting is not the interesting part of this script, so we do
    it here rather than crash on somebody's tree.
    """
    # Deepest leaf below each node, and the child that path goes through.
    depth = {}
    descend = {}
    diameter = 0.0
    meeting_node = None

    for node in tree.postorder_node_iter():
        best = []
        for child in node.child_node_iter():
            best.append((depth[id(child)] + (child.edge.length or 0.0), child))
        best.sort(key=lambda pair: pair[0], reverse=True)

        depth[id(node)] = best[0][0] if best else 0.0
        descend[id(node)] = best[0][1] if best else None

        # The longest path meeting at this node runs down its two deepest branches.
        span = (best[0][0] if best else 0.0) + (best[1][0] if len(best) > 1 else 0.0)
        if span > diameter:
            diameter = span
            meeting_node = node

    if meeting_node is None or diameter <= 0.0:
        sys.exit(
            "ERROR: The tree has no branch lengths to midpoint-root on. "
            "Use --root outgroup or --root as-is."
        )

    # The midpoint sits on the deeper of the two branches, this far below the
    # meeting node. depth >= diameter/2 there, so this never goes negative.
    remaining = depth[id(meeting_node)] - diameter / 2.0

    node = meeting_node
    while True:
        child = descend[id(node)]
        if child is None:
            break
        length = child.edge.length or 0.0
        if remaining <= length:
            break
        remaining -= length
        node = child

    child = descend[id(node)]
    if child is None or remaining <= 0.0:
        tree.reseed_at(node, update_bipartitions=False, suppress_unifurcations=True)
        return

    length = child.edge.length or 0.0
    tree.reroot_at_edge(
        child.edge,
        length1=remaining,
        length2=max(length - remaining, 0.0),
        update_bipartitions=False,
        suppress_unifurcations=True,
    )


def read_trees(path, root, outgroup):
    """Load one or more Newick trees and root them as requested."""
    trees = dendropy.TreeList.get(
        path=path,
        schema="newick",
        rooting="force-rooted",
        preserve_underscores=True,
    )
    if not trees:
        sys.exit("ERROR: No trees could be read from {}".format(path))

    for number, tree in enumerate(trees, start=1):
        if root == "midpoint":
            if any(edge.length is None for edge in tree.postorder_edge_iter() if edge.tail_node):
                sys.exit(
                    "ERROR: Tree {} has no branch lengths, so --root midpoint cannot be "
                    "used. Use --root outgroup or --root as-is.".format(number)
                )
            reroot_at_midpoint(tree)
        elif root == "outgroup":
            present = set(tip_labels(tree))
            absent = [tip for tip in outgroup if tip not in present]
            if absent:
                sys.exit(
                    "ERROR: Please check --outgroup -> {} tip(s) are not in tree {}: {}".format(
                        len(absent), number, ", ".join(absent)
                    )
                )
            node = tree.mrca(taxon_labels=outgroup)
            if node is tree.seed_node:
                sys.exit(
                    "ERROR: The outgroup spans the whole of tree {}, so it cannot root it. "
                    "Check that the outgroup is monophyletic.".format(number)
                )
            tree.reroot_at_edge(node.edge, update_bipartitions=False, suppress_unifurcations=True)

    return trees


def tip_labels(tree):
    return [
        node.taxon.label if node.taxon is not None else node.label
        for node in tree.leaf_node_iter()
    ]


def check_tips(trees, mapping, ignore_unlabeled):
    """Reconcile the groups file against the tips present in the trees."""
    first = set(tip_labels(trees[0]))
    for number, tree in enumerate(trees[1:], start=2):
        if set(tip_labels(tree)) != first:
            sys.exit(
                "ERROR: Tree {} does not have the same tips as tree 1. All trees in an "
                "ensemble must share one set of tips.".format(number)
            )

    missing = sorted(tip for tip in mapping if tip not in first)
    if missing:
        sys.exit(
            "ERROR: Please check the groups file -> {} tip(s) are not in the tree: {}".format(
                len(missing), ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
            )
        )

    unlabeled = sorted(tip for tip in first if tip not in mapping)
    if unlabeled and not ignore_unlabeled:
        sys.exit(
            "ERROR: {} tip(s) in the tree have no group in the groups file: {}\n"
            "Assign them a group, or pass --ignore-unlabeled to leave them out of the "
            "statistic.".format(
                len(unlabeled), ", ".join(unlabeled[:10]) + (" ..." if len(unlabeled) > 10 else "")
            )
        )
    return unlabeled


def write_tree(trees, path):
    """Write the rooted tree(s) back out, exactly as the tests saw them.

    Downstream rendering must draw the tree the statistics were computed on, not
    the file the user passed in -- those differ whenever --root did anything.
    """
    trees.write(path=path, schema="newick", suppress_rooting=True,
                unquoted_underscores=True, suppress_internal_node_labels=True)


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def format_value(value, places=4):
    return "NA" if value is None else "{:.{}f}".format(value, places)


# ---------------------------------------------------------------------------
# Group assignment
# ---------------------------------------------------------------------------


def build_assignment(trees, mapping, groups):
    """Canonical tip order plus the group index of each tip.

    Tips carrying no group label are assigned -1: they still shape the topology
    but belong to no group and take no part in permutations.
    """
    canonical_labels = sorted(set(tip_labels(trees[0])))
    group_index = {group: g for g, group in enumerate(groups)}
    assignment = [group_index.get(mapping.get(label, None), -1) for label in canonical_labels]
    group_sizes = [sum(1 for g in assignment if g == index) for index in range(len(groups))]
    return canonical_labels, assignment, group_sizes


def group_counts(tree_index, assignment, n_groups):
    """Tips of each group below every internal node, in post-order."""
    counts = [[0] * n_groups for _ in range(tree_index.n_internal)]
    for i, kids in enumerate(tree_index.children):
        row = counts[i]
        for is_leaf, j in kids:
            if is_leaf:
                group = assignment[tree_index.leaf_canonical[j]]
                if group >= 0:
                    row[group] += 1
            else:
                child_row = counts[j]
                for g in range(n_groups):
                    row[g] += child_row[g]
    return counts


def leaves_below(tree_index, node):
    """Canonical indices of every tip descended from an internal node."""
    found = []
    stack = [node]
    while stack:
        for is_leaf, j in tree_index.children[stack.pop()]:
            if is_leaf:
                found.append(tree_index.leaf_canonical[j])
            else:
                stack.append(j)
    return found


# ---------------------------------------------------------------------------
# Monophyly
# ---------------------------------------------------------------------------


def monophyly(tree_index, assignment, group, size, counts=None):
    """Is this group a clade, and if not, which tips break it?

    A group is monophyletic when the smallest subtree containing all of its tips
    contains nothing else. The tips that are in that subtree but not in the group
    are the clade-breaking taxa, which is the actionable part of a negative
    answer.

    Returns (is_monophyletic, mrca_node_or_None, intruder_canonical_indices).
    A group of fewer than two tips is trivially a clade; that is reported as
    monophyletic, but callers should treat it as uninformative.
    """
    if size < 1:
        return False, None, []
    if counts is None:
        counts = group_counts(tree_index, assignment, len(set(g for g in assignment if g >= 0)) or 1)

    if size == 1:
        return True, None, []

    # Post-order puts descendants first, so the first node holding the whole
    # group is its MRCA.
    mrca = None
    for i in range(tree_index.n_internal):
        if counts[i][group] == size:
            mrca = i
            break
    if mrca is None:
        return False, None, []

    intruders = [leaf for leaf in leaves_below(tree_index, mrca) if assignment[leaf] != group]
    return (not intruders), mrca, sorted(intruders)


# ---------------------------------------------------------------------------
# Permutation testing
# ---------------------------------------------------------------------------


def permutation_test(statistic, assignment, n_groups, observed, permutations, seed,
                     larger_is_extreme=True, observer=None):
    """Shuffle group labels across the labeled tips, holding the topology fixed.

    `statistic` takes an assignment and returns one value per group (None where
    undefined). One shuffle per replicate is applied to every tree in an
    ensemble, which preserves any correlated structure across loci.

    `larger_is_extreme` picks the tail: True for statistics where a high value
    means structure (gsi), False for statistics where a low value does (the
    Slatkin-Maddison migration count).

    `observer`, if given, is called with each replicate's values, which is how a
    caller recovers the null distribution rather than only the p-value.

    P = (1 + #{permuted at least as extreme}) / (permutations + 1), so the
    observed labeling counts as one of the replicates.
    """
    if permutations <= 0:
        return [None] * n_groups, 0

    labeled = [i for i, g in enumerate(assignment) if g >= 0]
    labels = [assignment[i] for i in labeled]

    rng = random.Random(seed)
    shuffled = list(assignment)
    at_least = [0] * n_groups

    for _ in range(permutations):
        rng.shuffle(labels)
        for position, label in zip(labeled, labels):
            shuffled[position] = label

        values = statistic(shuffled)
        if observer is not None:
            observer(values)
        for g, value in enumerate(values):
            if observed[g] is None or value is None:
                continue
            if (value >= observed[g]) if larger_is_extreme else (value <= observed[g]):
                at_least[g] += 1

    p_values = [
        None if observed[g] is None else (1 + at_least[g]) / (permutations + 1)
        for g in range(n_groups)
    ]
    return p_values, permutations


# ---------------------------------------------------------------------------
# Shared command-line surface
# ---------------------------------------------------------------------------


def add_tree_arguments(parser, permutations_default=None):
    """The --tree/--groups/--root options every analyze tool shares."""
    parser.add_argument("--tree", required=True, help="Newick file holding one tree, or an ensemble of trees.")
    parser.add_argument("--groups", required=True, help="Two-column table mapping tip label to group.")
    parser.add_argument(
        "--root",
        required=True,
        choices=["as-is", "midpoint", "outgroup"],
        help="How to root the tree. Monophyly is undefined on an unrooted tree, so this is "
        "explicit: 'as-is' trusts the Newick as written, 'midpoint' needs branch "
        "lengths, 'outgroup' needs --outgroup.",
    )
    parser.add_argument("--outgroup", help="Comma-separated tip labels to root on, for --root outgroup.")
    parser.add_argument(
        "--ignore-unlabeled", action="store_true",
        help="Allow tips that are absent from the groups file; they stay on the tree but take no group label.",
    )
    parser.add_argument("--groups-to-test", help="Comma-separated subset of groups to report (default: all).")
    if permutations_default is not None:
        parser.add_argument(
            "--permutations", type=int, default=permutations_default,
            help="Label permutations for the significance test (default: {}). 0 skips the test.".format(
                permutations_default
            ),
        )
        parser.add_argument("--seed", type=int, default=42, help="Random seed for the permutation test (default: 42).")
    return parser


def load_inputs(args):
    """Validate the shared arguments and load everything the tools need."""
    if args.root == "outgroup" and not args.outgroup:
        sys.exit("ERROR: --root outgroup requires --outgroup <tip[,tip...]>")

    outgroup = [tip.strip() for tip in args.outgroup.split(",") if tip.strip()] if args.outgroup else []

    mapping, groups = read_groups(args.groups)
    trees = read_trees(args.tree, args.root, outgroup)
    unlabeled = check_tips(trees, mapping, args.ignore_unlabeled)

    wanted = groups
    if getattr(args, "groups_to_test", None):
        wanted = [group.strip() for group in args.groups_to_test.split(",") if group.strip()]
        unknown = [group for group in wanted if group not in groups]
        if unknown:
            sys.exit(
                "ERROR: --groups-to-test names group(s) that are not in the groups file: {}".format(
                    ", ".join(unknown)
                )
            )

    canonical_labels, assignment, group_sizes = build_assignment(trees, mapping, groups)
    return {
        "trees": trees,
        "tree_indices": index_trees(trees, canonical_labels),
        "mapping": mapping,
        "groups": groups,
        "wanted": wanted,
        "canonical_labels": canonical_labels,
        "assignment": assignment,
        "group_sizes": group_sizes,
        "unlabeled": unlabeled,
        "outgroup": outgroup,
    }
