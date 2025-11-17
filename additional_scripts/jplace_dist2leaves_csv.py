#!/usr/bin/env python

import json
import csv
import re
from io import StringIO

from Bio import Phylo


def parse_tree_with_edge_nums(tree_str):
    """
    Parse the jplace tree string with Biopython, preserving tip names and
    mapping each jplace edge_num to the distal clade on that branch.
    """

    # Convert ":len{edge}" into ":len[&edge_num=edge]" so it becomes a comment
    annotated = re.sub(
        r":([0-9.eE+-]+)\{(\d+)\}",
        r":\1[&edge_num=\2]",
        tree_str,
    )

    tree = Phylo.read(StringIO(annotated), "newick")

    # Annotate parent pointers (Biopython Clade doesn’t track these by default)
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            child.parent = parent

    # Build edge_num -> clade (distal node of that edge)
    edge_to_clade = {}
    for clade in tree.find_clades(order="level"):
        comment = getattr(clade, "comment", None)
        if comment:
            m = re.search(r"edge_num=(\d+)", comment)
            if m:
                edge_to_clade[int(m.group(1))] = clade

    return tree, edge_to_clade


def dist_from_pos_on_edge(tree, branch_node, distal_length, leaf):
    """
    Distance from a point P located `distal_length` away from the **distal end**
    of branch (parent -> branch_node) to `leaf`.

    jplace: distal_length is from the distal (away from root) side of the
    reference edge to the placement attachment point.
    """

    L = branch_node.branch_length or 0.0
    parent = getattr(branch_node, "parent", None)

    # Leaf is on the distal side (the subtree rooted at branch_node),
    # including the branch_node itself.
    if (leaf is branch_node) or branch_node.is_parent_of(leaf):
        return distal_length + tree.distance(branch_node, leaf)

    # Leaf is on the proximal side (the rest of the tree)
    if parent is None:
        # Should not happen: root should not have an edge_num
        raise ValueError("Branch node has no parent; cannot determine proximal side.")

    return (L - distal_length) + tree.distance(parent, leaf)


def placement_leaf_distances(jplace_path, out_path):
    """
    For each placement in a jplace file, compute its distance to every leaf on
    the reference tree and write a long-format TSV.

    Columns:
      pquery_index, placement_index, query_names, edge_num,
      like_weight_ratio, distal_length, pendant_length,
      leaf_name, distance
    """
    with open(jplace_path) as f:
        jplace = json.load(f)

    tree, edge_to_clade = parse_tree_with_edge_nums(jplace["tree"])
    leaves = tree.get_terminals()

    fields = jplace["fields"]

    try:
        edge_idx = fields.index("edge_num")
        distal_idx = fields.index("distal_length")
        pendant_idx = fields.index("pendant_length")
    except ValueError as e:
        raise SystemExit(f"Required jplace field missing: {e}")

    lwr_idx = fields.index("like_weight_ratio") if "like_weight_ratio" in fields else None

    with open(out_path, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(
            [
                "pquery_index",
                "placement_index",
                "query_names",
                "edge_num",
                "like_weight_ratio",
                "distal_length",
                "pendant_length",
                "leaf_name",
                "distance",
            ]
        )

        for pq_idx, pquery in enumerate(jplace["placements"]):
            # Names: "n" (just names) or "nm" (name, multiplicity)
            names = pquery.get("n")
            if names is None and "nm" in pquery:
                names = [nm[0] for nm in pquery["nm"]]
            if names is None:
                names = []
            query_names = ",".join(names)

            placements = pquery["p"]
            # EPA sometimes uses a single list instead of list-of-lists
            if placements and not isinstance(placements[0], (list, tuple)):
                placements = [placements]

            for pl_idx, pl in enumerate(placements):
                edge_num = int(pl[edge_idx])
                distal = float(pl[distal_idx])
                pendant = float(pl[pendant_idx])
                lwr = float(pl[lwr_idx]) if lwr_idx is not None else None

                branch_node = edge_to_clade.get(edge_num)
                if branch_node is None:
                    raise ValueError(f"Edge number {edge_num} not found in tree")

                for leaf in leaves:
                    d = dist_from_pos_on_edge(tree, branch_node, distal, leaf)
                    writer.writerow(
                        [
                            pq_idx,
                            pl_idx,
                            query_names,
                            edge_num,
                            lwr,
                            distal,
                            pendant,
                            leaf.name,
                            d,
                        ]
                    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compute distances from each placement to every leaf in a jplace file."
    )
    parser.add_argument("jplace", help="Input .jplace file")
    parser.add_argument(
        "-o",
        "--out",
        default="placement_leaf_distances.tsv",
        help="Output TSV file (default: placement_leaf_distances.tsv)",
    )
    args = parser.parse_args()

    placement_leaf_distances(args.jplace, args.out)
