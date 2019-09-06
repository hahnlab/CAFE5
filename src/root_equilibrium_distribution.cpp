#include <numeric>
#include <iostream>

#include "root_distribution.h"
#include "root_equilibrium_distribution.h"
#include "poisson.h"
#include "io.h"
#include "optimizer.h"

uniform_distribution::uniform_distribution() : _p_root_distribution(new root_distribution())
{

}

uniform_distribution::~uniform_distribution()
{
    delete _p_root_distribution;
}

void uniform_distribution::initialize(const root_distribution* root_distribution)
{
    *_p_root_distribution = *root_distribution;
    _root_distribution_sum = _p_root_distribution->sum();
}

float uniform_distribution::compute(size_t val) const
{
    if (val >= _p_root_distribution->size())
        return 0;

    return float(_p_root_distribution->at(val)) / float(_root_distribution_sum);
}

::poisson_distribution::poisson_distribution(std::vector<gene_family> *p_gene_families)
{
    poisson_scorer scorer(*p_gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    cout << "\nEmpirical Prior Estimation Result : (" << result.num_iterations << " iterations)" << endl;
    cout << "Poisson lambda: " << result.values[0] << " &  Score: " << result.score << endl;

    _poisson_lambda = result.values[0];
}

void ::poisson_distribution::initialize(const root_distribution* root_distribution)
{
    poisson = get_prior_rfsize_poisson_lambda(0, root_distribution->size(), _poisson_lambda);
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
