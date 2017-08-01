#include "io.h"
#include "utils.h"
#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <sstream>

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "tree", required_argument, NULL, 't' },
  { "lambda_tree", optional_argument, NULL, 'y' },
  { "prefix", optional_argument, NULL, 'p' },
  { "simulate", optional_argument, NULL, 's' },
  { "nsims", optional_argument, NULL, 'n' },
  { "fixed_lambda", optional_argument, NULL, 'k' },
  { 0, 0, 0, 0 }
};

//! Read tree from user-provided tree file
/*!
  This function is called by CAFExp's main function when "--tree"/"-t" is specified
*/
clade* read_tree(string tree_file_path, bool lambda_tree) {
    newick_parser parser(false);
    
    // if this function is used to read lambda tree instead of phylogenetic tree
    if (lambda_tree) {
       parser.parse_to_lambdas = true;
    }
    
    ifstream tree_file(tree_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    string line;
    
    if (tree_file.good()) {
        getline(tree_file, line);
    }
    tree_file.close();
    
    parser.newick_string = line;
    clade *p_tree = parser.parse_newick();
    return p_tree;
}

//! Read gene family data from user-provided tab-delimited file
/*!
  This function is called by CAFExp's main function when "--infile"/"-i" is specified  
*/
void read_gene_families(string input_file_path) {
    vector<gene_family> gene_families;
    ifstream input_file(input_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    string line;
    bool is_first_line = true;
    while (getline(input_file, line)) {
        std::vector<std::string> tokens;
        
        if (is_first_line) {
            tokens = tokenize_str(line, '\t');
            is_first_line = false;
        }
        
        
    }
}

//! Populate famdist_map with root family distribution read from famdist_file_path
/*!
  This function is called by CAFExp's main function when "--simulate"/"-s" is specified 
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

/* START: Reading in gene family data */

//! Return highest gene family count. 
/*!
  CAFE had a not_root_max (which we use; see below) and a root_max = MAX(30, rint(max*1.25));
*/
int gene_family::max_family_size() const {
    int max_size = max_element(species_size_map.begin(), species_size_map.end(), max_value<string, int>)->second;
    return max_size + std::max(50, max_size / 5);
}

//! Return first element of pair
template <typename type1, typename type2>
type1 do_get_species(const std::pair<type1, type2> & p1) {
    return p1.first;
}

//! Return vector of species names
vector<string> gene_family::get_species() const {
    vector<string> result(species_size_map.size());
    transform(species_size_map.begin(), species_size_map.end(), result.begin(), do_get_species<string,int>); // transform performs an operation on all elements of a container (here, the operation is the template)
    
    return result;
}

/* END: Reading in gene family data */