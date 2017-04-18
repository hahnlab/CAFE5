#include "io.h"
#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <sstream>

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "prefix", optional_argument, NULL, 'p' },
  { "simulate", optional_argument, NULL, 's' },
  { "nsims", optional_argument, NULL, 'n' },
  { 0, 0, 0, 0 }
};


/* Populates famdist_map with root family distribution read from famdist_file_path */
map<int, int>* read_famdist(string famdist_file_path) {

  map<int, int> *p_famdist_map = new map<int, int>();
  ifstream famdist_file(famdist_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
  string line;
  while (getline(famdist_file, line)) {
    istringstream ist(line);
    int fam_size, fam_count;
    ist >> fam_size >> fam_count;
    (*p_famdist_map)[fam_size] = fam_count;
  }
  return p_famdist_map;
}
