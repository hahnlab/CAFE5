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
  { "output_prefix", required_argument, NULL, 'o'}, 
  { "tree", required_argument, NULL, 't' },
  { "fixed_lambda", required_argument, NULL, 'l' },
  { "fixed_multiple_lambdas", required_argument, NULL, 'm' },
  { "lambda_tree", required_argument, NULL, 'y' },
  { "n_gamma_cats", required_argument, NULL, 'k' },
  { "fixed_alpha", required_argument, NULL, 'a' },
  { "rootdist", required_argument, NULL, 'f'},
  { "poisson", optional_argument, NULL, 'p' },
  { "estimate", optional_argument, NULL, 'e' },
  { "simulate", optional_argument, NULL, 's' },
  { "nsims", optional_argument, NULL, 'n' },
  { "log", optional_argument, NULL, 'g'},
  { 0, 0, 0, 0 }
};

void input_parameters::check_input() {
    //! The user cannot specify both -e and -l
    if (estimate && fixed_lambda > 0.0) {
        throw runtime_error("You cannot both estimate (-e) and fix the lambda(s) value(s) (-l). Exiting...");
    }
        
    //! The user cannot specify both -l and -y
    if (fixed_lambda > 0.0 && !fixed_multiple_lambdas.empty()) {
        throw runtime_error("You cannot fix one lambda value (-l) and many lambda values (-m). Exiting...");
    }

    //! The user cannot specify both -i and -s
    if (!input_file_path.empty() && simulate) {
        throw runtime_error("You cannot simulate and try to infer simultaneously. Exiting...");
    }
}

/* START: Reading in gene family data */

gene_family::gene_family(trial *t)
{
	for (auto it = t->begin(); it != t->end(); ++it)
	{
		if (it->first->is_leaf())
			this->_species_size_map[it->first->get_taxon_name()] = it->second;
	}
	find_max_size();
}

//! Find highest gene family count. 
/*!
  CAFE had a not_root_max (which we use = _parsed_max_family_size; see below) and a root_max = MAX(30, rint(max*1.25));
*/
void gene_family::find_max_size() {
    // Max family size can only be found if there is data inside the object in the first place
    if (!_species_size_map.empty()) {
        _max_family_size = (*max_element(_species_size_map.begin(), _species_size_map.end(), max_value<string, int>)).second;
        _parsed_max_family_size = _max_family_size + std::max(50, _max_family_size/5);
    }
}

//! Mainly for debugging: In case one want to grab the gene count for a given species
int gene_family::get_species_size(std::string species) const {
    // First checks if species data has been entered (i.e., is key in map?)
    if (_species_size_map.find(species) == _species_size_map.end()) {
        throw std::runtime_error(species + " was not found in gene family " + _id);
    }
      
    return _species_size_map.at(species);
}
/* END: Reading in gene family data */

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
void read_gene_families(std::istream& input_file, clade *p_tree, std::vector<gene_family> *p_gene_families) {
    map<int, std::string> sp_col_map; // {col_idx: sp_name} 
    std::string line;
    bool is_header = true;
    map<int, string> leaf_indices;
    int index = 0;
    
    while (getline(input_file, line)) {
        std::vector<std::string> tokens = tokenize_str(line, '\t');
        if (!leaf_indices.empty() && line[0] != '#')
            is_header = false;

        // header
        if (is_header) {
            if (line[0] == '#')
            {
                if (p_tree == NULL)
                    throw std::invalid_argument("No matching tree provided for simulation");
                string taxon_name = line.substr(1);
                clade *p_descendant = p_tree->find_descendant(taxon_name);
                if (p_descendant == NULL)
                    throw std::runtime_error(taxon_name + " not located in tree");
                if (p_descendant->is_leaf())
                {
                    leaf_indices[index] = taxon_name;
                }
                index++;
            }
            else
            {
                is_header = false;
                if (leaf_indices.empty())
                {
                    for (int i = 0; i < tokens.size(); ++i) {
                        if (i == 0 || i == 1) {} // ignores description and ID cols
                        sp_col_map[i] = tokens[i];
                    }
                }
            }
        }
        // not header
        else {
            gene_family genfam; 
            
            for (int i = 0; i < tokens.size(); ++i) {
                if (leaf_indices.empty())
                {
                    if (i == 0)
                        genfam.set_desc(tokens[i]);
                    if (i == 1)
                        genfam.set_id(tokens[i]);
                    else {
                        std::string sp_name = sp_col_map[i];
                        cout << sp_name << " " << tokens[i] << endl;
                        genfam.set_species_size(sp_name, atoi(tokens[i].c_str()));
                    }
                }
                else
                {
                    if (leaf_indices.find(i) != leaf_indices.end())
                    {
                        std::string sp_name = leaf_indices[i];
                        cout << sp_name << " " << tokens[i] << endl;
                        genfam.set_species_size(sp_name, atoi(tokens[i].c_str()));
                    }
                }
            }
            
            genfam.find_max_size(); // getting max gene family size for this gene family
            p_gene_families->push_back(genfam);
        }
    }
    
    // return p_gene_families;
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

//! Return first element of pair
template <typename type1, typename type2>
type1 do_get_species(const std::pair<type1, type2> & p1) {
    return p1.first;
}

//! Return vector of species names
vector<string> gene_family::get_species() const {
    vector<string> result(_species_size_map.size());
    transform(_species_size_map.begin(), _species_size_map.end(), result.begin(), do_get_species<string,int>); // transform performs an operation on all elements of a container (here, the operation is the template)
    
    return result;
}

/* END: Reading in gene family data */