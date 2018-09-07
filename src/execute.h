#ifndef EXECUTE_H
#define EXECUTE_H

#include "io.h"

struct input_parameters;
class lambda;
class matrix_cache;
class model;
class root_equilibrium_distribution;

class action
{
public:
    virtual void execute(std::vector<model *>& models) = 0;
    virtual ~action() {}
};

class simulator : public action
{
    const input_parameters &_user_input;
    void simulate(std::vector<model *>& models, const input_parameters &my_input_parameters);
public:
    simulator(const input_parameters &my_input_parameters) : _user_input(my_input_parameters)
    {

    }
    virtual void execute(std::vector<model *>& models);
};

class chisquare_compare : public action
{
    std::string _values;
public:
    chisquare_compare(std::string values) : _values(values)
    {
    }
    virtual void execute(std::vector<model *>& models);
};

action* get_executor(input_parameters& user_input);

class execute {
public:
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

    void compute(std::vector<model *>& models, std::vector<gene_family> *p_gene_families, root_equilibrium_distribution *p_prior, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size);

    void estimate_lambda(std::vector<model *>& models, std::vector<gene_family> &gene_families, error_model *p_error_model, clade *p_tree, clade *p_lambda_tree, root_equilibrium_distribution *p_prior);

    void reconstruct(std::vector<model *>& models, const input_parameters &my_input_parameters, int max_family_size, root_equilibrium_distribution *p_prior);

    void write_results(std::vector<model *>& models, const input_parameters &my_input_parameters, std::vector<double>& pvalues);
};
#endif /* EXECUTE_H */
