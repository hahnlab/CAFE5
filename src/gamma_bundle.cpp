#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <functional>
#include <cmath>
#include <sstream>
#include <memory>

#include "gamma_bundle.h"
#include "gene_family_reconstructor.h"
#include "root_equilibrium_distribution.h"
#include "clade.h"
#include "gene_family.h"

using namespace std;

gamma_bundle::gamma_bundle(std::vector<double> lambda_multipliers, const clade *p_tree, const gene_family *p_gene_family,
    std::ostream & ost, const lambda* lambda, int max_family_size, int max_root_family_size) : 
    _p_tree(p_tree),
    _p_gene_family(p_gene_family),
    _lambda_multipliers(lambda_multipliers),
    _max_family_size(max_family_size),
    _max_root_family_size(max_root_family_size),
    _ost(ost),
    _p_lambda(lambda)
{
}

std::vector<const clade *> gamma_bundle::get_taxa()
{
    return _p_tree->find_internal_nodes();
}

string gamma_bundle::get_family_id() const
{ 
    return _p_gene_family->id();
}

bool gamma_bundle::prune(const vector<double>& _gamma_cat_probs, root_equilibrium_distribution *eq, matrix_cache& calc, const lambda *p_lambda) {
    //assert(_gamma_cat_probs.size() == _inf_processes.size());
    _category_likelihoods.clear();

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        auto partial_likelihood = inference_prune(*_p_gene_family, calc, p_lambda, _p_tree, _lambda_multipliers[k], _max_root_family_size, _max_family_size);
        if (accumulate(partial_likelihood.begin(), partial_likelihood.end(), 0.0) == 0.0)
            return false;   // saturation

        std::vector<double> full(partial_likelihood.size());
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = eq->compute(j);
            full[j] = partial_likelihood[j] * eq_freq;
        }

//        _category_likelihoods.push_back(accumulate(full.begin(), full.end(), 0.0) * _gamma_cat_probs[k]); // sum over all sizes (Felsenstein's approach)
        _category_likelihoods.push_back(*max_element(full.begin(), full.end()) * _gamma_cat_probs[k]); // get max (CAFE's approach)
    }

    return true;
}

clademap<double> get_weighted_averages(std::vector<clademap<int>>& m, const vector<double>& _gamma_cat_probs)
{
    cladevector nodes(m[0].size());
    std::transform(m[0].begin(), m[0].end(), nodes.begin(), [](std::pair<const clade *, int> v) { return v.first;  });

    clademap<double> result;
    for (auto node : nodes)
    {
        double val = 0.0;
        for (size_t i = 0; i<_gamma_cat_probs.size(); ++i)
        {
            val += _gamma_cat_probs[i] * double(m[i].at(node));
        }
        result[node] = val;
    }

    return result;
}


void gamma_bundle::reconstruct(const vector<double>& _gamma_cat_probs, matrix_cache *calc, root_equilibrium_distribution*prior)
{
    reconstructed_states.resize(_gamma_cat_probs.size());
    increase_decrease_map.resize(_gamma_cat_probs.size());
    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        unique_ptr<lambda> ml(_p_lambda->multiply(_lambda_multipliers[k]));
        reconstruct_gene_families(ml.get(), _p_tree, _max_family_size, _max_root_family_size, _p_gene_family, calc, prior, reconstructed_states[k]);
        compute_increase_decrease(reconstructed_states[k], increase_decrease_map[k]);
    }

    // multiply every reconstruction by gamma_cat_prob
    reconstruction = get_weighted_averages(reconstructed_states, _gamma_cat_probs);

    compute_increase_decrease(reconstruction, _increase_decrease_map);
}

string gamma_bundle::get_reconstructed_states(const clade *node) const
{
    std::ostringstream ost;

    if (node->is_leaf())
    {
        ost << _p_gene_family->get_species_size(node->get_taxon_name());
    }
    else
    {
        for (auto& r : reconstructed_states)
        {
            ost << r.at(node) << '_';
        }
        ost << std::round(reconstruction.at(node));
    }
    return ost.str();
}

void gamma_bundle::print_reconstruction(std::ostream & ost, cladevector& order)
{
    string family_id = _p_gene_family->id();

    auto g = [this](const clade *node) {
        return get_reconstructed_states(node);
    };

    auto f = [order, g, this](const clade *node) {
        return newick_node(node, order, g);
    };
    
    ost << "  TREE " << family_id << " = ";
    _p_tree->write_newick(ost, f);

    ost << ';' << endl;
}

increase_decrease gamma_bundle::get_increases_decreases(cladevector& order, double pvalue)
{
    increase_decrease result;
    result.change.resize(order.size());
    result.gene_family_id = _p_gene_family->id();
    result.pvalue = pvalue;

    transform(order.begin(), order.end(), result.change.begin(), [this](const clade *taxon)->family_size_change {
        if (taxon->is_leaf() || taxon->is_root())
            return Constant;
        else
            return _increase_decrease_map.at(taxon);
    });

    result.category_likelihoods = _category_likelihoods;
    return result;
}

double gamma_bundle::get_lambda_likelihood(int family_id)
{
    return _lambda_multipliers[family_id];
}
