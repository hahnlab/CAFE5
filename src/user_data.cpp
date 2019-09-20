#include "user_data.h"

#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "io.h"
#include "lambda.h"

using namespace std;

//! Read user provided gene family data (whose path is stored in input_parameters instance)
/// @param[in] my_input_parameters Parsed parameters passed to the application
/// @param[in] p_tree The tree to be used in calculations. Necessary for syncing tree data to gene family data
/// @param[out] p_gene_families Parsed data in the gene family file specified by my_input_parameters
/// @param[out] max_family_size Equal to the largest family size given in the file plus 20%, or plus 50 if the largest family size is more than 250
/// @param[out] max_root_family_size Equal to 5/4 the size of the largest family size given in the file (with a minimum of 30)
void user_data::read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree, std::vector<gene_family> *p_gene_families) {

    ifstream input_file(my_input_parameters.input_file_path);
    if (!input_file.is_open())
        throw std::runtime_error("Failed to open " + my_input_parameters.input_file_path + ". Exiting...");

    read_gene_families(input_file, p_tree, p_gene_families); // in io.cpp/io.h

                                                             // Iterating over gene families to get max gene family size
    for (std::vector<gene_family>::iterator it = p_gene_families->begin(); it != p_gene_families->end(); ++it) {
        int this_family_max_size = it->get_max_size();

        if (max_family_size < this_family_max_size)
            max_family_size = this_family_max_size;
    }

    max_root_family_size = std::max(30, static_cast<int>(std::rint(max_family_size*1.25)));
    max_family_size = max_family_size + std::max(50, max_family_size / 5);
    // cout << "Read input file " << my_input_parameters.input_file_path << "." << endl;
    // cout << "Max (parsed) family size is: " << max_family_size << endl;
    // cout << "Max root family size is: " << max_root_family_size << endl;
}

//! Read user provided error model file (whose path is stored in input_parameters instance)
void user_data::read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model) {

    ifstream error_model_file(my_input_parameters.error_model_file_path);
    if (!error_model_file.is_open()) {
        throw std::runtime_error("Failed to open " + my_input_parameters.error_model_file_path + ". Exiting...");
    }

    read_error_model_file(error_model_file, p_error_model);

} // GOTTA WRITE THIS!

  //! Read user provided phylogenetic tree (whose path is stored in input_parameters instance)
clade * user_data::read_input_tree(const input_parameters &my_input_parameters) {
    return read_tree(my_input_parameters.tree_file_path, false);
}

//! Read user provided lambda tree (lambda structure)
clade * user_data::read_lambda_tree(const input_parameters &my_input_parameters) {
    return read_tree(my_input_parameters.lambda_tree_file_path, true);
}

//! Read user provided single or multiple lambdas
lambda * user_data::read_lambda(const input_parameters &my_input_parameters, clade *p_lambda_tree) {

    lambda *p_lambda = NULL; // lambda is an abstract class, and so we can only instantiate it as single_lambda or multiple lambda -- therefore initializing it to NULL

                             // -l
    if (my_input_parameters.fixed_lambda > 0.0) {
        p_lambda = new single_lambda(my_input_parameters.fixed_lambda);
        // call_viterbi(max_family_size, max_root_family_size, 15, p_lambda, *p_gene_families, p_tree);
    }

    // -m
    if (!my_input_parameters.fixed_multiple_lambdas.empty()) {
        map<std::string, int> node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map(); // allows matching different lambda values to nodes in lambda tree
        vector<string> lambdastrings = tokenize_str(my_input_parameters.fixed_multiple_lambdas, ',');
        vector<double> lambdas(lambdastrings.size());

        // transform is like R's apply (vector lambdas takes the outputs, here we are making doubles from strings
        transform(lambdastrings.begin(), lambdastrings.end(), lambdas.begin(),
            [](string const& val) { return stod(val); } // this is the equivalent of a Python's lambda function
        );

        p_lambda = new multiple_lambda(node_name_to_lambda_index, lambdas);
    }

    return p_lambda;
}

//! Populate famdist_map with root family distribution read from famdist_file_path
void user_data::read_rootdist(string rootdist_file_path) {

    ifstream rootdist_file(rootdist_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    if (!rootdist_file.is_open())
        throw std::runtime_error("Failed to open file '" + rootdist_file_path + "'");
    string line;
    while (getline(rootdist_file, line)) {
        istringstream ist(line);
        int fam_size, fam_count;
        ist >> fam_size >> fam_count;
        rootdist[fam_size] = fam_count;
    }
}


void user_data::read_datafiles(const input_parameters& my_input_parameters)
{
    /* -t */
    if (!my_input_parameters.tree_file_path.empty()) {
        p_tree = read_input_tree(my_input_parameters); // populates p_tree (pointer to phylogenetic tree)
    }

    /* -i */
    if (!my_input_parameters.input_file_path.empty()) {
        // Populates (pointer to) vector of gene family instances, max_family_size and max_root_family_size (last two passed by reference)
        read_gene_family_data(my_input_parameters, max_family_size, max_root_family_size, p_tree, &gene_families);
    }

    /* -e */
    if (!my_input_parameters.error_model_file_path.empty()) {
        p_error_model = new error_model;
        read_error_model(my_input_parameters, p_error_model);
    }

    /* -y */
    if (!my_input_parameters.lambda_tree_file_path.empty()) {
        p_lambda_tree = read_lambda_tree(my_input_parameters);
		p_tree->validate_lambda_tree(p_lambda_tree);
    }

    /* -l/-m (in the absence of -l, estimate) */
    p_lambda = read_lambda(my_input_parameters, p_lambda_tree);

    if (!my_input_parameters.rootdist.empty())
        read_rootdist(my_input_parameters.rootdist);
}

