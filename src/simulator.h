#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "execute.h"

class root_distribution;

/*! @brief Build simulated families based on the user's input
*/

class simulator : public action
{
    void simulate(std::vector<model *>& models, const input_parameters &my_input_parameters);
public:
    simulator(user_data& d, const input_parameters& ui) : action(d, ui)
    {

    }
    trial* create_trial(model *p_model, const root_distribution& rd, int family_number, const matrix_cache& cache);

    virtual void execute(std::vector<model *>& models);
    void print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<trial *>& results);

    //! Does the actual work of simulation. Calls the given model to load simulation parameters,
    //! and places the simulations in results. Every fifty simulations, the model's \ref model::perturb_lambda
    //! is called in order to provide a bit of additional randomness in the simulation.
    void simulate_processes(model *p_model, std::vector<trial *>& results);

};

int select_root_size(const user_data& data, const root_distribution& rd, int family_number);

#endif
