#!/usr/bin/env python

from Bio import Phylo
from matplotlib import pyplot as plt
import csv
import argparse
import os
import io

colors = dict(Increase='green', Decrease='red', Rapid='blue')

labels = dict()

def display(tree_type):
    if tree_type == 'Rapid':
        return 'Rapidly evolving families'

    else:
        return tree_type


def label(n):
    family_count = n.info
    if n.name:
        return "%s (%s)" % (n.name, family_count)

    else:
        return family_count

def prettify_tree(fig):
    plt.ylabel('')
    plt.tick_params(
        axis='y',  # changes apply to the y-axis
        which='both',  # both major and minor ticks are affected
        left='off',  # ticks along the left edge are off
        right='off',  # ticks along the right edge are off
        labelleft='off')  # labels along the left edge are off

    [fig.gca().spines[border].set_visible(False) for border in ['left', 'right', 'top']]


def get_clades(node):
    clades = []
    clades.append(node)
    if node.clades:
        for clade in node.clades:
            clades += get_clades(clade)
    return clades


def draw_tree(datafile, tree_type, node_ids, output_file):
    with open(datafile) as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            labels[row["#Taxon_ID"]] = row

    with open(node_ids) as file:
        for line in file:
            if line.startswith("# IDs of nodes:"):
                newick_string = line.replace("# IDs of nodes:", "").strip()
                ids = Phylo.read(io.StringIO(newick_string), "newick")
            if line.startswith("Tree:"):
                newick_string = line.replace("# IDs of nodes:", "").strip()
                tree = Phylo.read(io.StringIO(newick_string), "newick")

#    for clade, c_id in zip(tree.find_clades(), id_tree.find_clades()):
    for clade, id in zip(get_clades(tree.root), get_clades(ids.root)):
       clade.info = labels[id.name][display(tree_type)] if id.name in labels else ""

    tree.ladderize()   # Flip branches so deeper clades are displayed at top

    plt.ion()
    fig = plt.figure(frameon=False)
    Phylo.draw(tree, axes=fig.gca(), do_show=False, label_func=label, label_colors = lambda n: colors[tree_type])
    plt.title(display(tree_type) + " (count)")
    prettify_tree(fig)
    plt.ioff()

    if output_file:
        fig.savefig(output_file, format='png', bbox_inches='tight', dpi=300)

    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, prog="draw_expansion_tree.py")
    parser.add_argument("-i", "--input-file", action="store", dest="input_file", required=True, type=str, help="full path to CAFE5's clade_results.txt ")
    parser.add_argument("-d", "--ids", action="store", dest="id_tree", required=True, help="Full path to CAFE5 report.cafe file")
    parser.add_argument("-y", "--tree-type", action="store", dest="tree_type", default="Increase", required=False, type=str, help="Expansions, Contractions, or Rapid")
    parser.add_argument("-o", "--output-file", action="store", dest="output_file", required=False, type=str, help="output PNG file name")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        exit("Could not find input file. Exiting...\n")

    draw_tree(args.input_file, args.tree_type, args.id_tree, args.output_file)
