#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "execute.h"
#include "clade.h"

class root_distribution;

struct simulated_family
{
    clademap<int> values;
    double lambda;

    simulated_family()
    {

    }

    ~simulated_family()
    {

    }

    simulated_family(simulated_family&& other) 
    {
        *this = std::move(other);
    }

    // this is the move assignment operator
    simulated_family& operator=(simulated_family&& other) 
    { 
        values = std::move(other.values);
        lambda = other.lambda;
        return *this; 
    }
};

/*! @brief Build simulated families based on the user's input

The user asks for a simulation given a lambda. We generate a multiplier m from a normal distribution 
with a standard deviation of .3. If the user also specifies an alpha, A, we generate m from a 
gamma distribution described by alpha=A, beta = 1/A

If the user also specifies a number of clusters, K, we do the following:
•	Discretize the gamma value so we have a multiplier for each cluster.
•	Randomly select a multiplier and use that as the mean of a normal distribution. 
      The standard deviation is calculated as a percentage of the distance between the clusters.
•	We choose a final multiplier by selecting randomly from this normal distribution.

We generate 50 families by taking the user’s lambda and multiplying it by m. We repeat this process 
until we have the number of families requested by the user.


*/

/// The root distribution is specified by the user. If they provide a rootdist file, the number of families is taken directly
/// from this distribution and it matches the distribution exactly.
/// If they request more simulations than are in the distribution, random values are drawn from the distribution to match
/// the count. If they request fewer simulations, random values are deleted from the distribution.
/// If they do not specify a distribution, random family sizes are selected between 1 and 100.

class simulator : public action
{
    void simulate(std::vector<model *>& models, const input_parameters &my_input_parameters);
public:
    simulator(user_data& d, const input_parameters& ui);

    simulated_family create_trial(const lambda *p_lambda, int family_number, const matrix_cache& cache);

    virtual void execute(std::vector<model *>& models);
    void print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<simulated_family>& results);

    //! Does the actual work of simulation. Calls the given model to load simulation parameters,
    //! and places the simulations in results. Every fifty simulations, the model's \ref model::perturb_lambda
    //! is called in order to provide a bit of additional randomness in the simulation.
    void simulate_processes(model *p_model, std::vector<simulated_family>& results);

};

#endif
