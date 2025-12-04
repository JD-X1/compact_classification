# test_placement_distances.py
import csv
import json
from jplace_dist2leaves_csv import placement_leaf_distances

def write_test_jplace(path):
    data = {
        "tree": "((A:1{1},B:1{2})n1:1{3},C:2{4}):0{0};",
        "fields": [
            "edge_num",
            "likelihood",
            "like_weight_ratio",
            "distal_length",
            "pendant_length"
        ],
        "placements": [
            {
                "p": [
                    [3, -1.0, 1.0, 0.2, 0.5]
                ],
                "n": ["query1"]
            }
        ]
    }
    with open(path, "w") as f:
        json.dump(data, f)

def test_simple_tree(tmp_path):
    jplace_path = tmp_path / "test_simple.jplace"
    out_path = tmp_path / "test_simple_distances.tsv"

    write_test_jplace(jplace_path)
    placement_leaf_distances(str(jplace_path), str(out_path))

    rows = []
    with open(out_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    dist = {row["leaf_name"]: float(row["distance"]) for row in rows}
    assert dist["A"] == 1.2
    assert dist["B"] == 1.2
    assert dist["C"] == 2.8
