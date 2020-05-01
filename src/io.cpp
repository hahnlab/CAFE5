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
#include <algorithm>
#include <iomanip>

#include <sys/stat.h>

#include "io.h"
#include "gene_family.h"
#include "error_model.h"
#include "clade.h"

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "error_model", optional_argument, NULL, 'e' },
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
  { "pvalue", required_argument, NULL, 'P' },
  { "log", optional_argument, NULL, 'g'},
  { "zero_root", no_argument, NULL, 'z' },
  { "lambda_per_family", no_argument, NULL, 'b' },
  { "log_config", required_argument, NULL, 'L' },
  { "optimizer_expansion", optional_argument, NULL, 'E' },
  { "optimizer_reflection", optional_argument, NULL, 'R' },
  { "optimizer_iterations", optional_argument, NULL, 'I' },
  { "help", no_argument, NULL, 'h'},
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
    
    if (is_simulating)
    {
        // Must specify a lambda
        if (fixed_lambda <= 0.0 && fixed_multiple_lambdas.empty()) {
            throw runtime_error("Cannot simulate without initial lambda values");
        }

        if (fixed_alpha <= 0.0 && this->n_gamma_cats > 1) {
            throw runtime_error("Cannot simulate gamma clusters without an alpha value");
        }
    }
    else
    {
        if (fixed_alpha >= 0.0 && n_gamma_cats == 1) {
            throw runtime_error("Alpha specified with 1 gamma category.");
    }


    if (lambda_per_family)
    {
        if (input_file_path.empty())
            throw runtime_error("No family file provided");
        if (tree_file_path.empty())
            throw runtime_error("No tree file provided");
    }

    if (n_gamma_cats > 1 && use_error_model && error_model_file_path.empty())
    {
        throw runtime_error("Estimating an error model with a gamma distribution is not supported at this time");
    }
    //! Options -i and -f cannot be both specified. Either one or the other is used to specify the root eq freq distr'n.
    if (!input_file_path.empty() && !rootdist.empty()) {
        throw runtime_error("Options -i and -f are mutually exclusive.");
    }

    }
}

/* START: Reading in tree data */
//! Read tree from user-provided tree file
/*!
  This function is called by CAFExp's main function when "--tree"/"-t" is specified
*/
clade* read_tree(string tree_file_path, bool lambda_tree) {
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
    
    clade *p_tree = parse_newick(line, lambda_tree);
    
    if (p_tree->is_leaf())
        throw std::runtime_error(tree_file_path + " does not seem to be a valid tree");

    return p_tree;
}
/* END: Reading in tree data */

//! Read gene family data from user-provided tab-delimited file
/*!
  This function is called by execute::read_gene_family_data, which is itself called by CAFExp's main function when "--infile"/"-i" is specified  
*/
void read_gene_families(std::istream& input_file, clade *p_tree, std::vector<gene_family> &gene_families) {
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
                    for (size_t i = 0; i < tokens.size(); ++i) {
                        if (i == 0 || i == 1) {} // Ignoring description and ID columns
                        sp_col_map[i] = tokens[i];
                    }
                }
            }
        }
        
        // Header has ended, reading gene family counts
        else {
            gene_family genfam; 
            
            for (size_t i = 0; i < tokens.size(); ++i) {                
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
            
            gene_families.push_back(genfam);
        }
    }

    if (gene_families.empty())
        throw std::runtime_error("No families found");
}
/* END: Reading in gene family data */

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
            p_error_model->set_max_family_size(max_cnt);
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
                p_error_model->set_probabilities(sz, values);
            }
        }
    }
        // std::vector<std::string> tokens = tokenize_str(line, '\t'); 
}
/* END: Reading in error model data */

void write_error_model_file(std::ostream& ost, error_model& errormodel)
{
    ost << "maxcnt: " << errormodel.get_max_family_size()-1 << "\n";
    ost << "cntdiff:";
    for (int j : errormodel._deviations) {
        ost << " " << j;
    }
    ost << "\n";

    vector<double> last_probs;
    for (size_t j = 0; j < errormodel.get_max_family_size(); j++) {
        auto probs = errormodel.get_probs(j);
        if (probs == last_probs) continue;
        last_probs = probs;

        ost << j;
        for (auto p : probs)
            ost << " " << p;
        ost << endl;
    }
}

//! Split string into vector of strings given delimiter
std::vector<std::string> tokenize_str(std::string some_string, char some_delim) {
	std::istringstream ist(some_string);
	std::string token;
	std::vector<std::string> tokens;

	while (std::getline(ist, token, some_delim)) {
		tokens.push_back(token);
	}

	return tokens;
}

/// OS-specific. If mkdir succeeds it returns 0. If it returns -1
/// check to see if it failed because the directory already exists.
/// If it failed for some other reason, throw an exception.
void create_directory(std::string& dir)
{
    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if (errno != EEXIST)
            throw std::runtime_error("Failed to create directory");
    }
}
