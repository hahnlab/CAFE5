#include <numeric>
#include <iostream>
#include <random>
#include <algorithm>

#include "easylogging++.h"

#include "root_equilibrium_distribution.h"
#include "poisson.h"
#include "io.h"
#include "optimizer.h"

extern std::mt19937 randomizer_engine; // seeding random number engine

root_equilibrium_distribution::root_equilibrium_distribution(const map<int, int>& root_distribution)
{
    if (root_distribution.empty())
        throw std::runtime_error("No root distribution specified");

    for (auto it = root_distribution.begin(); it != root_distribution.end(); ++it) {
        for (int i = 0; i < it->second; ++i) {
            _vectorized_distribution.push_back(it->first);
        }
    }

    build_percentages();

}

root_equilibrium_distribution::root_equilibrium_distribution(size_t max_size)
{
    _vectorized_distribution.resize(max_size);
    iota(_vectorized_distribution.begin(), _vectorized_distribution.end(), 1);
    build_percentages();
}

root_equilibrium_distribution::root_equilibrium_distribution(double poisson_lambda, int num_values)
{
    auto poisson = get_prior_rfsize_poisson_lambda(0, num_values, poisson_lambda);
    int n = poisson.size();
    for (size_t i = 0; i < poisson.size() - 1; ++i)
    {
        for (size_t j = 0; j < poisson[i] * n; ++j)
            _vectorized_distribution.push_back(i + 1);
    }
    build_percentages();

}

root_equilibrium_distribution::root_equilibrium_distribution(std::vector<gene_family> *p_gene_families, int num_values)
{
    poisson_scorer scorer(*p_gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    LOG(INFO) << "Empirical Prior Estimation Result : (" << result.num_iterations << " iterations)";
    LOG(INFO) << "Poisson lambda: " << result.values[0] << " &  Score: " << result.score;

    double poisson_lambda = result.values[0];

    auto poisson = get_prior_rfsize_poisson_lambda(0, num_values, poisson_lambda);
    int n = poisson.size();
    for (size_t i = 0; i < poisson.size() - 1; ++i)
    {
        for (size_t j = 0; j < poisson[0] * n; ++j)
            _vectorized_distribution.push_back(i + 1);
    }
    build_percentages();

}

void root_equilibrium_distribution::build_percentages()
{
    auto max = (size_t)*max_element(_vectorized_distribution.begin(), _vectorized_distribution.end()) + 1;
    _frequency_percentage.resize(max);
    for (size_t i = 0; i < max; ++i)
    {
        size_t c = count(_vectorized_distribution.begin(), _vectorized_distribution.end(), i);
        _frequency_percentage[i] = float(c) / float(_vectorized_distribution.size());
    }
}

float root_equilibrium_distribution::compute(size_t val) const
{
    if (val >= _frequency_percentage.size())
        return 0;

    return _frequency_percentage[val];
}

int root_equilibrium_distribution::select_root_size(int family_number) const
{
    if ((size_t)family_number >= _vectorized_distribution.size())
        return 0;

    return _vectorized_distribution[family_number];
}

void root_equilibrium_distribution::resize(size_t new_size)
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


/// Root distributions are affected by three parameters: -p, -i, -f
/// If a rootdist file is specified (-f), those values will be used for the root distribution and the other flags
/// are ignored. Otherwise, if a poisson distribution is specified (-p) with a value, the root
/// distribution will be based on that poisson distribution. If no poisson value
/// is specified, a family file must be given (-i) and those families will be used to calculate a
/// poisson distribution. If a Poisson distribution is used, values above a max family size
/// will be considered to be 0. The max family size defaults to 100, but is calculated from
/// family file if one is given.
root_equilibrium_distribution create_root_distribution(const input_parameters& params, std::vector<gene_family> *p_gene_families, const std::map<int, int>& root_distribution, int max_root_family_size)
{
    root_equilibrium_distribution result(1);

    if (!root_distribution.empty())
    {
        result = root_equilibrium_distribution(root_distribution);
    }
    else
    {
        if (params.use_uniform_eq_freq)
        {
            result = root_equilibrium_distribution(max_root_family_size);
        }
        else
        {
            double pl = params.poisson_lambda;
            if (pl > 0)
                result = root_equilibrium_distribution(pl, max_root_family_size);
            else
                result = root_equilibrium_distribution(p_gene_families, max_root_family_size);

        }

    }

    if (params.nsims > 0)
        result.resize(params.nsims);

    return result;
}
