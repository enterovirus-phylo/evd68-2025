import argparse
import json
import pandas as pd
from Bio import Phylo


def main():
    parser = argparse.ArgumentParser(description="Build branch-length node data from TreeTime output")
    parser.add_argument("--tree", required=True, help="Timetree in nexus format")
    parser.add_argument("--divergence-tree", required=True, help="Divergence tree in nexus format")
    parser.add_argument("--dates", required=True, help="Dates TSV from TreeTime")
    parser.add_argument("--rate", default=1.0, type=float, help="Rate used in TreeTime (default: 1.0)")
    parser.add_argument("--outliers", required=True, help="Outliers TSV from TreeTime")
    parser.add_argument("--output-tree", required=True, help="Output tree in newick format")
    parser.add_argument("--output-node-data", required=True, help="Output node data JSON")
    args = parser.parse_args()

    T = list(Phylo.parse(args.tree, "newick"))[-2]
    T_div = list(Phylo.parse(args.divergence_tree, "newick"))[-2]
    outliers = pd.read_csv(args.outliers, sep="\t", index_col=0)
    for outlier in outliers.index:
        T.prune(outlier)
        T_div.prune(outlier)

    for node in T.find_clades():
        node.confidence = None  # remove confidence values if they exist, as augur doesn't expect them
        node.branch_length = node.branch_length*args.rate if node.branch_length is not None else 0.0

    # write back as newick to ensure branch lengths are in the format augur expects
    Phylo.write(T, args.output_tree, "newick")



    nodes_data = {}
    for node in T.find_clades():
        if node.branch_length is not None:
            nodes_data[node.name] = {"branch_length": node.branch_length,
                                     "clock_length": node.branch_length}
        else:
            nodes_data[node.name] = {"branch_length": None,
                                     "clock_length": None}

    for node in T_div.find_clades():
        if node.branch_length is not None:
            if node.name not in nodes_data:
                nodes_data[node.name] = {}
            nodes_data[node.name]["mutation_length"] = node.branch_length

    d = pd.read_csv(args.dates, sep="\t")

    for _, row in d.iterrows():
        strain = row["#node"]
        if strain in outliers.index:
            continue
        if strain not in nodes_data:
            continue
        if row["date"] != '--':
            nodes_data[strain]["numdate"] = float(row["numeric date"])
            nodes_data[strain]["date"] = row["date"]
        else:
            nodes_data[strain]["numdate"] = None
            nodes_data[strain]["date"] = None

    with open(args.output_node_data, "w") as f:
        json.dump({"nodes": nodes_data}, f, indent=2)


if __name__ == "__main__":
    main()
