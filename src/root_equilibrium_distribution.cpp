#include <numeric>
#include <iostream>

#include "root_distribution.h"
#include "root_equilibrium_distribution.h"
#include "poisson.h"
#include "io.h"
#include "lambda.h" // for definition of optimizer :/

float uniform_distribution::compute(int val) const
{
    if (val >= _root_distribution->size())
        return 0;

    return float(_root_distribution->at(val)) / float(_root_distribution->sum());
}

::poisson_distribution::poisson_distribution(std::vector<gene_family> *p_gene_families)
{
    poisson_scorer scorer(*p_gene_families);
    optimizer opt(&scorer);

    auto result = opt.optimize();

    cout << "Empirical Prior Estimation Result : (" << result.num_iterations << " iterations)" << endl;
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
