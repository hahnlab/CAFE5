#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <functional>
#include <cmath>
#include <sstream>

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
    _max_root_family_size(max_root_family_size)
{
//    for (double m: lambda_multipliers)
//        _inf_processes.push_back(new inference_process(ost, lambda, m, p_tree, max_family_size, max_root_family_size, _p_gene_family, NULL)); 

    for (auto m : lambda_multipliers)
    {
        _rec_processes.push_back(new gene_family_reconstructor(ost, lambda, m, _p_tree, max_family_size, max_root_family_size, _p_gene_family,
            NULL, NULL));
    }
}

gamma_bundle::~gamma_bundle()
{
    for (auto r : _rec_processes)
        delete r;
}

std::vector<const clade *> gamma_bundle::get_taxa()
{
    return _p_tree->find_internal_nodes();
}

string gamma_bundle::get_family_id() const
{ 
    return _p_gene_family->id();
}

void gamma_bundle::set_values(matrix_cache *calc, root_equilibrium_distribution*prior)
{
    for (auto rec : _rec_processes)
    {
        rec->set_values(calc, prior);
    }
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

clademap<double> get_weighted_averages(std::vector<gene_family_reconstructor *> m, const vector<double>& _gamma_cat_probs)
{
    auto nodes = m[0]->get_nodes();

    clademap<double> result;
    for (auto node : nodes)
    {
        double val = 0.0;
        for (size_t i = 0; i<_gamma_cat_probs.size(); ++i)
        {
            val += _gamma_cat_probs[i] * double(m[i]->get_reconstructed_value(node));
        }
        result[node] = val;
    }

    return result;
}


void gamma_bundle::reconstruct(const vector<double>& _gamma_cat_probs)
{
    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        _rec_processes[k]->reconstruct();
    }

    // multiply every reconstruction by gamma_cat_prob
    reconstruction = get_weighted_averages(_rec_processes, _gamma_cat_probs);

    compute_increase_decrease(reconstruction, increase_decrease_map);
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
        for (auto r : _rec_processes)
            ost << r->get_reconstructed_value(node) << '_';

        ost << std::round(reconstruction.at(node));
    }
    return ost.str();
}

void gamma_bundle::print_reconstruction(std::ostream & ost, cladevector& order)
{
    string family_id = _p_gene_family->id();

    auto f = [order, this](const clade *node) { 
        return newick_node(node, order, this);
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
            return increase_decrease_map.at(taxon);
    });

    result.category_likelihoods = _category_likelihoods;
    return result;
}

double gamma_bundle::get_lambda_likelihood(int family_id)
{
    return _lambda_multipliers[family_id];
}
