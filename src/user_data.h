#ifndef USER_DATA_H
#define USER_DATA_H

#include <stddef.h>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "gene_family.h"
#include "root_equilibrium_distribution.h"

class root_equilibrium_distribution;
class clade;
class lambda;
class error_model;
struct input_parameters;

/// Class holding data defined by the user, or derived from data defined by the user
class user_data {
public:
    user_data() : prior(100)
    {

    }
    int max_family_size = 120; //!<  The maximum family size for which probabilities will be calculated
    int max_root_family_size = 125; //!<  The maximum family size for which probabilities will be calculated at the root of the tree

    clade *p_tree = NULL; // instead of new clade(), o.w. mem leak
    lambda *p_lambda = NULL;
    clade *p_lambda_tree = NULL;
    error_model *p_error_model = NULL;
    root_equilibrium_distribution prior;

    std::vector<gene_family> gene_families;
    std::map<int, int> rootdist;

    void read_datafiles(const input_parameters& my_input_parameters);

    //! Read in gene family data
    void read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree, std::vector<gene_family> *p_gene_families);

    //! Read in error model file
    void read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model);

    //! Read in phylogenetic tree data
    clade * read_input_tree(const input_parameters &my_input_parameters);

    //! Read in lambda tree
    clade * read_lambda_tree(const input_parameters &my_input_parameters);

    //! Read in single or multiple lambda
    lambda * read_lambda(const input_parameters &my_input_parameters, clade *p_lambda_tree);

    void read_rootdist(std::string rootdist_file_path);
};

#endif
