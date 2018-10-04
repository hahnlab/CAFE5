#ifndef EXECUTE_H
#define EXECUTE_H

#include "io.h"

struct input_parameters;
class lambda;
class matrix_cache;
class model;
class root_equilibrium_distribution;
class user_data;

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

class estimator : public action
{
    user_data& data;
    root_equilibrium_distribution *p_prior;
    vector<double> compute_pvalues(const user_data& data, int number_of_simulations);
    const input_parameters &_user_input;
public:
    estimator(user_data& d, root_equilibrium_distribution *dist, const input_parameters& ui) : data(d), p_prior(dist), _user_input(ui)
    {

    }
    virtual void execute(std::vector<model *>& models);

    void compute(std::vector<model *>& models, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size);

    void estimate_lambda(std::vector<model *>& models, std::vector<gene_family> &gene_families, error_model *p_error_model, clade *p_tree, clade *p_lambda_tree);
};

action* get_executor(input_parameters& user_input, user_data& data, root_equilibrium_distribution *dist);

#endif /* EXECUTE_H */
