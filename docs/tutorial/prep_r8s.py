"""
This script defines and runs the r8s_prep function, which reads in input from the user, and prepares a r8s control file (so r8s can be run).
"""

__author__ = "Fabio H. K. Mendes"

import os
import argparse

def r8s_prep(path_tree_file, sites_n, list_of_spp_tuples, list_of_spp_cal_points, path_output_file):
    """
    Prepare r8s input file

    [str] path_tree_file: path to .txt file containing tree in NEWICK format
    [str] n_sites: number of sites in alignment that was used to infer species tree
    [list] list_of_spp_tuples: list of tuples (each tuple being two species IDs whose mrca's age we are constraining; e.g., [('ENSG00','ENSPTR'),('ENSFCA','ENSECA')]
    [list] list_of_spp_cal_points: list of flats, one for each tuple in list_of_spp_tuples (e.g., [6.4,80])
    """

    tree_str = str()
    with open(path_tree_file, "r") as tree_file:
        tree_str = tree_file.readline().rstrip()

    fixage_lines = list()
    with open(path_output_file, "w") as output_file:
        output_file.write(
            "#NEXUS\nbegin trees;\ntree nj_tree = [&R] " + tree_str + "\nEnd;\nbegin rates;\nblformat nsites=" + \
            sites_n + " lengths=persite ultrametric=no;\ncollapse;\n")

        for idx, pair in enumerate(list_of_spp_tuples):
            node_name = pair[0][-3:] + pair[1][-3:]
            output_file.write(
                "mrca " + node_name + " " + pair[0] + " " + pair[1] + ";\n")
            fixage_lines.append("fixage taxon=" + node_name + " age=" + list_of_cal_points[idx] + ";\n")

        for line in fixage_lines:
            output_file.write(line)

        output_file.write(
            "divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;\ndescribe plot=chronogram;\ndescribe plot=tree_description;\nend;"
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, prog="mcl2rawcafe.py")
    parser.add_argument("-i", "--input-file", action="store", dest="input_file", required=True, type=str, help="full path to .txt file containing tree in NEWICK format")
    parser.add_argument("-o", "--output-file", action="store", dest="output_file", required=True, type=str, help="full path to file to be written (r8s input file)")
    parser.add_argument("-s", "--sites-n", action="store", dest="sites_n", required=True, type=str, help="number of sites in alignment used to estimate species tree")
    parser.add_argument("-p", "--pairs-species", action="store", dest="spp_pairs", required=True, type=str, help="")
    parser.add_argument("-c", "--calibration-points", action="store", dest="cal_points", required=True, type=str, help="")

    args = parser.parse_args()

    list_of_spp_tuples = list()
    for pair in args.spp_pairs.split(" "):
        list_of_spp_tuples.append(tuple(pair.split(",")))

    list_of_cal_points = args.cal_points.split(",")

    print "\nRunning cafetutorial_clade_and_size_filter.py as a standalone...\n"
    
    r8s_prep(args.input_file, args.sites_n, list_of_spp_tuples, list_of_cal_points, args.output_file)


