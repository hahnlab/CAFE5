#include <numeric>
#include <set>
#include <cassert>
#include <algorithm>

#include "error_model.h"

using namespace std;

error_model::error_model()
{
    _deviations = { -1, 0, 1 };
}

void error_model::set_max_family_size(size_t max_cnt) {
    _max_family_size = max_cnt;
}

void error_model::set_deviations(std::vector<std::string> deviations) {
    _deviations.resize(deviations.size());
    transform(deviations.begin(), deviations.end(), _deviations.begin(), [](const string& s) {return std::stoi(s); });
}

inline bool is_nearly_equal(double x, double y)
{
    const double epsilon = 0.0000000000001;
    return std::abs(x - y) <= epsilon * std::abs(x);
}

void error_model::set_probabilities(size_t fam_size, std::vector<double> probs_deviation) {
    if ((fam_size == 0 || _error_dists.empty()) && !is_nearly_equal(probs_deviation[0], 0.0))
    {
        throw std::runtime_error("Cannot have a non-zero probability for family size 0 for negative deviation");
    }

    if (!is_nearly_equal(accumulate(probs_deviation.begin(), probs_deviation.end(), 0.0), 1.0))
    {
        throw std::runtime_error("Sum of probabilities must be equal to one");
    }

    if (_error_dists.empty())
        _error_dists.push_back(probs_deviation);

    if (_error_dists.size() <= fam_size)
    {
        _error_dists.resize(fam_size + 1, _error_dists.back());
    }
    _error_dists[fam_size] = probs_deviation; // fam_size starts at 0 at tips, so fam_size = index of vector
}

std::vector<double> error_model::get_probs(size_t fam_size) const {
    if (fam_size >= _error_dists.size() && fam_size <= _max_family_size)
        return _error_dists.back();

    return _error_dists[fam_size];
}

std::vector<double> error_model::get_epsilons() const {
    set<double> unique_values;
    for (auto& vec : _error_dists)
        unique_values.insert(vec.back());

    vector<double> result(unique_values.size());
    copy(unique_values.begin(), unique_values.end(), result.begin());
    return result;
}

// simple case where we have a single epsilon value in the tree
void error_model::update_single_epsilon(double new_epsilon)
{
    auto epsilons = get_epsilons();
    assert(epsilons.size() == 1);
    map<double, double> replacements;
    replacements[epsilons[0]] = new_epsilon;
    replace_epsilons(&replacements);
}

void error_model::replace_epsilons(std::map<double, double>* new_epsilons)
{
    vector<double> vec = _error_dists[0];
    assert(vec.size() == 3);
    for (auto kv : *new_epsilons)
    {
        if (is_nearly_equal(kv.first, vec.back()))
        {
            vec.back() = kv.second;
            vec[1] = 1 - kv.second;
            set_probabilities(0, vec);
        }
    }

    for (size_t i = 1; i < _error_dists.size(); ++i)
    {
        vector<double> vec = _error_dists[i];
        assert(vec.size() == 3);

        for (auto kv : *new_epsilons)
        {
            if (is_nearly_equal(kv.first, vec.back()))
            {
                vec.back() = kv.second;
                vec.front() = kv.second;
                vec[1] = 1 - (kv.second * 2);
                set_probabilities(i, vec);
            }
        }
    }
}

