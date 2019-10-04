import os
import sys
import itertools
import argparse

def list_fa_files(path_to_dir):
    """ Return list of .fa files """
    for f in os.listdir(path_to_dir):
        print f
    return []
#    return [f for f in os.listdir(path_to_dir) if f.endswith(".fa") and not f.startswith("longest")]
    
def fasta_iter(fasta_filehandle):
    """ Yield tuples of header, sequence """

    faiter = (x[1] for x in \
              itertools.groupby(fasta_filehandle, lambda line: line[0] == ">")) # tuple, items from ">" to ">"

    for header in faiter:
        header = next(header)[1:].strip() # drop the ">"
        seq = "".join(s.strip() for s in next(faiter))

        yield header, seq

def seq_seventify(seq):
    """ Pretty print sequence """

    return "\n".join(seq[i:i+70] for i in range(0, len(seq), 70))

def find_longest_iso(path_fa_file):
    """ Return dict with only longest isoforms from .fa files """

    longest_iso_dict = dict()
    with open(path_fa_file, "r") as fa_file:
        for header, seq in fasta_iter(fa_file):

            # ignoring unavailable sequences
            if header.find("|") == -1 or seq.find("unavailable") != -1:
                continue

            seq_id, seq_len_str = header.split("|")
            seq_len = int(seq_len_str)

            try:
                if seq_len > longest_iso_dict[seq_id]:
                    longest_iso_dict[seq_id] = seq_len
            except:
                longest_iso_dict[seq_id] = seq_len

    return longest_iso_dict

def longest_iso_printer(path_fa_dir):
    """ Print .fa files with only longest isoforms """

    list_fa_filenames = list_fa_files(path_fa_dir)
    for fa_file_name in list_fa_filenames:
        longest_iso_dict = find_longest_iso(path_fa_dir + fa_file_name)

        with open(path_fa_dir + "longest_" + fa_file_name, "w") as longest_fa_file:
            with open(path_fa_dir + fa_file_name, "r") as fa_file:
                for header, seq in fasta_iter(fa_file):

                    try:
                        seq_id, seq_len_str = header.split("|")
                        seq_len = int(seq_len_str)

                        if seq_len == longest_iso_dict[seq_id]:
                            longest_fa_file.write(">"+seq_id + "\n" + "".join(seq_seventify(seq)) + "\n")

                    except:
                        continue

        print fa_file_name + "...done!"  
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, prog="cafetutorial_longest_iso.py")
    parser.add_argument("-d", "--data-directory", action="store", dest="dir_path", required=True, type=str, help="full path to directory containing .fa files")

    args = parser.parse_args()

    if not args.dir_path.endswith("/"):
        exit("\nPlease make sure the path to the directory containing the .fa files ends with \"/\"\n")
    
    longest_iso_printer(args.dir_path)
