#ifndef EXECUTE_H
#define EXECUTE_H

#include "io.h"

struct input_parameters;
class lambda;
class probability_calculator;
class model;
class root_equilibrium_distribution;

class execute {
public:
    //! Read in gene family data (-i)
    // std::vector<gene_family> * read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree);
    void read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree, std::vector<gene_family> *p_gene_families);
    
    //! Read in phylogenetic tree data (-t)
//    clade * read_input_tree(const input_parameters &my_input_parameters);
    clade * read_input_tree(const input_parameters &my_input_parameters);

    //! Read in lambda tree (-y)
    clade * read_lambda_tree(const input_parameters &my_input_parameters);

    //! Read in single or multiple lambda (-l or -m)
    lambda * read_lambda(const input_parameters &my_input_parameters, probability_calculator &my_calculator, clade *p_lambda_tree);

    void compute(std::vector<model *>& models, std::vector<gene_family> *p_gene_families, root_equilibrium_distribution *p_prior, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size);

    void simulate(std::vector<model *>& models, const input_parameters &my_input_parameters);

    void reconstruct(std::vector<model *>& models, const input_parameters &my_input_parameters, root_equilibrium_distribution *p_prior, probability_calculator& calculator);
    void estimate_lambda(model * p_model, const input_parameters &my_input_parameters, root_equilibrium_distribution *p_prior, clade *p_tree, clade *p_lambda_tree,
        std::vector<gene_family>* p_gene_families, int max_family_size, int max_root_family_size, probability_calculator& calculator);
};
#endif /* EXECUTE_H */
