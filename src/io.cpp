#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <sstream>
#include <cstring>
#include <set>
#include <numeric>
#include <cassert>
#include <cmath>

#include "io.h"
#include "utils.h"
#include "gene_family.h"

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "error_model", required_argument, NULL, 'e' },
  { "output_prefix", required_argument, NULL, 'o'}, 
  { "tree", required_argument, NULL, 't' },
  { "fixed_lambda", required_argument, NULL, 'l' },
  { "fixed_multiple_lambdas", required_argument, NULL, 'm' },
  { "lambda_tree", required_argument, NULL, 'y' },
  { "n_gamma_cats", required_argument, NULL, 'k' },
  { "fixed_alpha", required_argument, NULL, 'a' },
  { "rootdist", required_argument, NULL, 'f'},
  { "poisson", optional_argument, NULL, 'p' },
  { "simulate", optional_argument, NULL, 's' },
  { "chisquare_compare", required_argument, NULL, 'r' },
  { "log", optional_argument, NULL, 'g'},
  { 0, 0, 0, 0 }
};

void input_parameters::check_input() {
    //! Options -l and -m cannot both specified.
    if (fixed_lambda > 0.0 && !fixed_multiple_lambdas.empty()) {
        throw runtime_error("Options -l and -m are mutually exclusive.");
    }

    //! Option -m requires a lambda tree (-y)
    if (!fixed_multiple_lambdas.empty() && lambda_tree_file_path.empty()) {
        throw runtime_error("Multiple lambda values (-m) specified with no lambda tree (-y)");
    }
    
    //! Options -l and -i have to be both specified (if estimating and not simulating).
    if (fixed_lambda > 0.0 && input_file_path.empty() && !is_simulating) {
        throw runtime_error("Options -l and -i must both be provided an argument.");
    }
    
    //! The number of simulated families is specified either through -s, or through -f. Cannot be both. 
    if (is_simulating && !rootdist.empty()) {
        throw runtime_error("Option -s cannot be provided an argument if -f is specified.");
    }
    
    //! Options -i and -f cannot be both specified. Either one or the other is used to specify the root eq freq distr'n.
    if (!input_file_path.empty() && !rootdist.empty()) {
        throw runtime_error("Options -i and -f are mutually exclusive.");
    }

    if (fixed_alpha >= 0.0 && n_gamma_cats == 1) {
        throw runtime_error("Alpha specified with 1 gamma category.");
    }
}

/* START: Reading in tree data */
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
    if (!tree_file.is_open())
    {
        throw std::runtime_error("Failed to open " + tree_file_path);
    }

    string line;
    
    if (tree_file.good()) {
        getline(tree_file, line);
    }
    tree_file.close();
    
    parser.newick_string = line;
    clade *p_tree = parser.parse_newick();
    
    return p_tree;
}
/* END: Reading in tree data */

//! Read gene family data from user-provided tab-delimited file
/*!
  This function is called by execute::read_gene_family_data, which is itself called by CAFExp's main function when "--infile"/"-i" is specified  
*/
void read_gene_families(std::istream& input_file, clade *p_tree, std::vector<gene_family> *p_gene_families) {
    map<int, std::string> sp_col_map; // For dealing with CAFE input format, {col_idx: sp_name} 
    map<int, string> leaf_indices; // For dealing with CAFExp input format, {idx: sp_name}, idx goes from 0 to number of species
    std::string line;
    bool is_header = true;
    int index = 0;
    
    while (getline(input_file, line)) {
        std::vector<std::string> tokens = tokenize_str(line, '\t');
        // Check if we are done reading the headers
        if (!leaf_indices.empty() && line[0] != '#') { is_header = false; }

        // If still reading header
        if (is_header) {
            
            // Reading header lines, CAFExp input format
            if (line[0] == '#') {
                if (p_tree == NULL) { throw std::runtime_error("No tree was provided."); }
                string taxon_name = line.substr(1); // Grabs from character 1 and on
                if (taxon_name.back() == '\r')
                    taxon_name.pop_back();

                auto p_descendant = p_tree->find_descendant(taxon_name); // Searches (from root down) and grabs clade "taxon_name" root
                
                if (p_descendant == NULL) { throw std::runtime_error(taxon_name + " not located in tree"); }
                if (p_descendant->is_leaf()) { leaf_indices[index] = taxon_name; } // Only leaves matter for computation or estimation
                index++;
            }
            
            // Reading single header line, CAFE input format
            else {
                is_header = false;
                
                // If reading CAFE input format, leaf_indices (which is for CAFExp input format) should be empty
                if (leaf_indices.empty()) {
                    for (int i = 0; i < tokens.size(); ++i) {
                        if (i == 0 || i == 1) {} // Ignoring description and ID columns
                        sp_col_map[i] = tokens[i];
                    }
                }
            }
        }
        
        // Header has ended, reading gene family counts
        else {
            gene_family genfam; 
            
            for (int i = 0; i < tokens.size(); ++i) {                
                // If reading CAFE input format, leaf_indices (which is for CAFExp input format) should be empty
                if (leaf_indices.empty()) {
                    if (i == 0) { genfam.set_desc(tokens[i]); }
                    else if (i == 1) { genfam.set_id(tokens[i]); }
                    else {
                        std::string sp_name = sp_col_map[i];
                        genfam.set_species_size(sp_name, atoi(tokens[i].c_str()));
                        // cout << "Species " << sp_name << " has " << tokens[i] << "gene members." << endl;
                    }
                }
                
                // If reading CAFExp input format
                else {
                    // If index i is in leaf_indices
                    if (leaf_indices.find(i) != leaf_indices.end()) { // This should always be true...
                        std::string sp_name = leaf_indices[i];
                        genfam.set_species_size(sp_name, atoi(tokens[i].c_str()));
                        // cout << "Species " << sp_name << " has " << tokens[i] << "gene members." << endl;
                    }
                    else
                    {
                        if (i == tokens.size() - 1)
                            genfam.set_id(tokens[i]);
                    }
                }
            }
            
            p_gene_families->push_back(genfam);
        }
    }
}
/* END: Reading in gene family data */

/* START: Reading in error model data */
void error_model::set_max_cnt(int max_cnt) {
    _max_cnt = max_cnt;
}

void error_model::set_deviations(std::vector<std::string> deviations) {
    for (std::vector<std::string>::iterator it = deviations.begin(); it != deviations.end(); ++it) {
        _deviations.push_back(std::stoi(*it));
    }
}

inline bool is_nearly_equal(double x, double y)
{
    const double epsilon = 0.0000000000001;
    return std::abs(x - y) <= epsilon * std::abs(x);
}

void error_model::set_probs(int fam_size, std::vector<double> probs_deviation) {
    if (fam_size == 0 && !is_nearly_equal(probs_deviation[0], 0.0))
    {
        throw std::runtime_error("Cannot have a non-zero probability for family size 0 for negative deviation");
    }

    if (!is_nearly_equal(accumulate(probs_deviation.begin(), probs_deviation.end(), 0.0), 1.0))
    {
        throw std::runtime_error("Sum of probabilities must be equal to one");
    }

    if (_error_dists.empty())
        _error_dists.push_back(vector<double>(probs_deviation.size()));

    if (_error_dists.size() <= fam_size)
    {
        _error_dists.resize(fam_size + 1, _error_dists.back());
    }
    _error_dists[fam_size] = probs_deviation; // fam_size starts at 0 at tips, so fam_size = index of vector
}

std::vector<double> error_model::get_probs(int fam_size) const {
    return _error_dists[fam_size];
}

std::vector<double> error_model::get_epsilons() const {
    set<double> unique_values;
    for (auto& vec : _error_dists)
        unique_values.insert(vec.back());

    vector<double> result(unique_values.size());
    copy(unique_values.begin(), unique_values.end(), result.begin());
    return result;
}

// simple case where we have a single epsilon value in the tree
void error_model::update_single_epsilon(double new_epsilon)
{
    auto epsilons = get_epsilons();
    assert(epsilons.size() == 1);
    map<double, double> replacements;
    replacements[epsilons[0]] = new_epsilon;
    replace_epsilons(&replacements);
}

error_model* error_model::replace_epsilons(std::map<double, double> *new_epsilons)
{
    vector<double> vec = _error_dists[0];
    assert(vec.size() == 3);
    for (auto kv : *new_epsilons)
    {
        if (is_nearly_equal(kv.first, vec.back()))
        {
            vec.back() = kv.second;
            vec[1] = 1 - kv.second;
            set_probs(0, vec);
        }
    }

    for (size_t i = 1; i < _error_dists.size(); ++i)
    {
        vector<double> vec = _error_dists[i];
        assert(vec.size() == 3);

        for (auto kv : *new_epsilons)
        {
            if (is_nearly_equal(kv.first, vec.back()))
            {
                vec.back() = kv.second;
                vec.front() = kv.second;
                vec[1] = 1 - (kv.second * 2);
                set_probs(i, vec);
            }
        }
    }
}

double to_double(string s)
{
    return std::stod(s);
}
//! Read user-provided error model
/*!
  This function is called by execute::read_error_model, which is itself called by CAFExp's main function when "--error_model"/"-e" is specified  
*/
void read_error_model_file(std::istream& error_model_file, error_model *p_error_model) {
    std::string line;
    std::string max_header = "max";
    std::string cnt_diff_header = "cnt";
    
    while (getline(error_model_file, line)) {
        std::vector<std::string> tokens;
        
        // maxcnt line
        if (strncmp(line.c_str(), max_header.c_str(), max_header.size()) == 0) { 
            tokens = tokenize_str(line, ':');
            tokens[1].erase(remove_if(tokens[1].begin(), tokens[1].end(), ::isspace), tokens[1].end()); // removing whitespace
            int max_cnt = std::stoi(tokens[1]);
            p_error_model->set_max_cnt(max_cnt);
        }
        
        // cntdiff line
        else if (strncmp(line.c_str(), cnt_diff_header.c_str(), cnt_diff_header.size()) == 0) { 
            tokens = tokenize_str(line, ' ');
            
            if (tokens.size() % 2 != 0) { 
                throw std::runtime_error("Number of different count differences in the error model (including 0) is not an odd number. Exiting...");
            }
            
            std::vector<std::string> cnt_diffs(tokens.begin()+1, tokens.end());
            
//            cout << "Count diffs are" << endl;
//            for (auto it = cnt_diffs.begin(); it != cnt_diffs.end(); ++it) {
//                cout << *it << " ";
//            }
//            cout << endl;
            p_error_model->set_deviations(cnt_diffs);
        }
        else
        {
            tokens = tokenize_str(line, ' ');
            if (tokens.size() > 0)
            {
                int sz = std::stoi(tokens[0]);
                vector<double> values(tokens.size() - 1);
                transform(tokens.begin() + 1, tokens.end(), values.begin(), to_double);
                p_error_model->set_probs(sz, values);
            }
        }
    }
        // std::vector<std::string> tokens = tokenize_str(line, '\t'); 
}
/* END: Reading in error model data */

