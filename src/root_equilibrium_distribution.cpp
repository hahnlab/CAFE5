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

root_equilibrium_distribution::root_equilibrium_distribution(double poisson_lambda, size_t num_values)
{
    create_from_poisson(poisson_lambda, num_values);
    build_percentages();
}

root_equilibrium_distribution::root_equilibrium_distribution(const std::vector<gene_family>& gene_families, size_t num_values)
{
    poisson_scorer scorer(gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    LOG(INFO) << "\nEmpirical Prior Estimation Result : (" << result.num_iterations << " iterations)";
    LOG(INFO) << "Poisson lambda: " << result.values[0] << " &  Score: " << result.score << "\n";

    create_from_poisson(result.values[0], num_values);
    build_percentages();

}

void root_equilibrium_distribution::create_from_poisson(double poisson_lambda, size_t num_values)
{
    for (int i = 0; _vectorized_distribution.size() < num_values; ++i)
    {
        double pct = poisspdf(i, poisson_lambda);
        for (size_t j = 0; j < pct * num_values; ++j)
            _vectorized_distribution.push_back(i + 1);
    }

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
