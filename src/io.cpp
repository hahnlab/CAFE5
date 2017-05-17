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
  { "fixed_lambda", optional_argument, NULL, 'k' },
  { 0, 0, 0, 0 }
};


//! Populate famdist_map with root family distribution read from famdist_file_path
/*!
  This function is called by CAFExp's main function when "simulate" is specified 
*/
map<int, int>* read_rootdist(string rootdist_file_path) {

    map<int, int> *p_rootdist_map = new map<int, int>();
    ifstream rootdist_file(rootdist_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    string line;
    while (getline(rootdist_file, line)) {
        istringstream ist(line);
        int fam_size, fam_count;
        ist >> fam_size >> fam_count;
        (*p_rootdist_map)[fam_size] = fam_count;
    }
  
  return p_rootdist_map;
}

/* START: Printing functions for simulation engine */

//! Print simulations from provided root family distribution
void print_simulation(std::vector<vector<trial *> > &sim, std::ostream& ost) {

    trial::iterator names_it = sim[0][0]->begin();
    
    // printing header with '#' followed by species names, in the order they will appear below
    for (; names_it != sim[0][0]->end(); ++names_it) {
	    ost << "#" << names_it->first->get_taxon_name() << endl;
    }
    
    // printing gene family sizes
    for (int i = 0; i < sim.size(); ++i) { // iterating over gene family sizes     
        for (int j = 0; j < sim[i].size(); ++j) { // iterating over trial of ith gene family size
            for (trial::iterator jth_trial_it = sim[i][j]->begin(); jth_trial_it != sim[i][j]->end(); ++jth_trial_it) { // accessing jth trial inside ith gene family
                ost << jth_trial_it->second << "\t";
            }
            
            ost << endl;
        }
    }
}

/* END: Printing functions for simulation engine */