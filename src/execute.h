#ifndef EXECUTE_H
#define EXECUTE_H

#include <vector>

#include "io.h"

struct input_parameters;
class lambda;
class matrix_cache;
class model;
class root_equilibrium_distribution;
class user_data;
class simulation_data;
    
/*! @brief All of the actions that the application can perform 

    Based on user input, CAFExp creates an action which drives the rest of the program.
    Action is a base class that defines which task is being performed.

    All actions take the input parameters that the user entered, and the 
    data that is described by those parameters (loaded files, clades, etc.)
*/
class action
{
protected:
    user_data& data;
    const input_parameters &_user_input;
public:
    /*! Performs an action on each of the models provided */
    virtual void execute(std::vector<model *>& models) = 0;

    /*! Construct an action and store references to the parameters the user asked for and any data already loaded */
    action(user_data& d, const input_parameters& ui) : data(d), _user_input(ui)
    {

    }
    virtual ~action() {}

    bool quiet = false;
};

class chisquare_compare : public action
{
    std::string _values;
public:
    chisquare_compare(user_data& d, const input_parameters& ui) : action(d, ui)
    {
        _values = ui.chisquare_compare;
    }
    virtual void execute(std::vector<model *>& models);
};

/*! @brief Estimator is used to guess at any missing values that the
    user did not provide. It will provide values for lambda, gamma, 
    and epsilon.
*/
class estimator : public action
{
public:
    estimator(user_data& d, const input_parameters& ui) : action(d, ui)
    {

    }

    virtual void execute(std::vector<model *>& models);

    void write_error_model_if_specified(const input_parameters& my_input_parameters, const model * p_model);

    void compute(std::vector<model *>& models, const input_parameters &my_input_parameters);

    void estimate_missing_variables(std::vector<model *>& models, user_data& data);
    void estimate_lambda_per_family(model *p_model, std::ostream& ost);

};

void initialization_failure_advice(std::ostream& ost, const std::vector<gene_family>& families);

action* get_executor(input_parameters& user_input, user_data& data);

#endif /* EXECUTE_H */
