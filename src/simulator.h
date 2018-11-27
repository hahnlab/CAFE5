#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "execute.h"

class root_distribution;

class simulator : public action
{
    void simulate(std::vector<model *>& models, const input_parameters &my_input_parameters);
public:
    simulator(user_data& d, const input_parameters& ui) : action(d, ui)
    {

    }
    trial* create_trial(model *p_model, const root_distribution& rd, int family_number, const user_data& data, const matrix_cache& cache);

    virtual void execute(std::vector<model *>& models);
    void print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<trial *>& results);
    void simulate_processes(model *p_model, std::vector<trial *>& results);

};

int select_root_size(const user_data& data, const root_distribution& rd, int family_number);

#endif
