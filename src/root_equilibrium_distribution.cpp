#include <numeric>
#include <iostream>

#include "root_equilibrium_distribution.h"
#include "poisson.h"
#include "io.h"


float uniform_distribution::compute(int val) const
{
    int sum = std::accumulate(_rootdist_vec.begin(), _rootdist_vec.end(), 0);
    return float(_rootdist_vec[val]) / float(sum);
}

::poisson_distribution::poisson_distribution(std::vector<gene_family> *p_gene_families)
{
    vector<double> root_poisson_lambda = find_poisson_lambda(*p_gene_families);
    _poisson_lambda = root_poisson_lambda[0];
    cout << "Estimated poisson lambda: " << _poisson_lambda << std::endl;
}

void ::poisson_distribution::initialize(std::vector<int> rootdist_vec)
{
    poisson = get_prior_rfsize_poisson_lambda(0, rootdist_vec.size(), _poisson_lambda);
}

root_equilibrium_distribution* root_eq_dist_factory(const input_parameters& my_input_parameters, std::vector<gene_family> *p_gene_families)
{
    root_equilibrium_distribution *p_prior = NULL;
    if (my_input_parameters.use_uniform_eq_freq)
    {
        p_prior = new uniform_distribution();
    }
    else
    {
        double pl = my_input_parameters.poisson_lambda;
        if (pl > 0)
            p_prior = new ::poisson_distribution(pl);
        else
            p_prior = new ::poisson_distribution(p_gene_families);
    }

    return p_prior;
}
