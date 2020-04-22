#include <numeric>
#include <iostream>
#include <random>
#include <algorithm>


#include "root_equilibrium_distribution.h"
#include "poisson.h"
#include "io.h"
#include "optimizer.h"

extern std::mt19937 randomizer_engine; // seeding random number engine

int uniform_distribution::select_root_size(int family_number) const
{
    // return a random root size
    std::uniform_int_distribution<> dis(0, _max_root_family_size+1);
    return dis(randomizer_engine); 
}

specified_distribution::specified_distribution(const map<int, int>& root_distribution)
{
    if (root_distribution.empty())
        throw std::runtime_error("No root distribution specified");

    for (auto it = root_distribution.begin(); it != root_distribution.end(); ++it) {
        for (int i = 0; i < it->second; ++i) {
            _vectorized_distribution.push_back(it->first);
        }
    }
}

float specified_distribution::compute(size_t val) const
{
    if (val >= _vectorized_distribution.size())
        return 0;

    return float(_vectorized_distribution.at(val)) / float(std::accumulate(_vectorized_distribution.begin(), _vectorized_distribution.end(), 0));
}

int specified_distribution::select_root_size(int family_number) const
{
    if ((size_t)family_number >= _vectorized_distribution.size())
        return 0;

    return _vectorized_distribution[family_number];
}

void specified_distribution::resize(size_t new_size)
{
    auto& v = _vectorized_distribution;
    if (new_size < v.size())
    {
        // pare back the distribution randomly
        shuffle(v.begin(), v.end(), randomizer_engine);
        v.erase(v.begin() + new_size, v.end());
    }
    else
    {
        std::uniform_int_distribution<> dis(0, v.size() - 1);
        for (size_t i = v.size(); i < new_size; ++i)
        {
            v.push_back(v.at(dis(randomizer_engine)));
        }
    }
    sort(v.begin(), v.end());
}


::poisson_distribution::poisson_distribution(double poisson_lambda, int num_values) : _poisson_lambda(poisson_lambda)
{
    poisson = get_prior_rfsize_poisson_lambda(0, num_values, _poisson_lambda);
}

::poisson_distribution::poisson_distribution(std::vector<gene_family> *p_gene_families, int num_values)
{
    poisson_scorer scorer(*p_gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    cout << "\nEmpirical Prior Estimation Result : (" << result.num_iterations << " iterations)" << endl;
    cout << "Poisson lambda: " << result.values[0] << " &  Score: " << result.score << endl;

    _poisson_lambda = result.values[0];

    poisson = get_prior_rfsize_poisson_lambda(0, num_values, _poisson_lambda);
}

int ::poisson_distribution::select_root_size(int family_number) const
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double my_random = distribution(randomizer_engine);
    double x = poisson[0];
    for (size_t i = 0; i < poisson.size()-1; ++i)
    {
        if (my_random < x)
            return i + 1;
        x += poisson[i + 1];
    }

    return 0;
}

root_equilibrium_distribution* root_eq_dist_factory(const input_parameters& params, std::vector<gene_family> *p_gene_families, const std::map<int, int>& root_distribution, int max_root_family_size)
{
    root_equilibrium_distribution *p_prior = NULL;
    if (params.use_uniform_eq_freq)
    {
        if (root_distribution.empty())
            p_prior = new uniform_distribution(max_root_family_size * 0.8);
        else
        {
            auto t = new specified_distribution(root_distribution);
            if (params.nsims > 0)
                t->resize(params.nsims);
            p_prior = t;
        }
    }
    else
    {
        int num_values = max_root_family_size * 0.8;
        if (!root_distribution.empty())
            num_values = accumulate(root_distribution.begin(), root_distribution.end(), 0,
                [](int acc, std::pair<int, int> p) { return (acc + p.second); });

        double pl = params.poisson_lambda;
        if (pl > 0)
            p_prior = new ::poisson_distribution(pl, num_values);
        else
            p_prior = new ::poisson_distribution(p_gene_families, num_values);
    }

    return p_prior;
}
