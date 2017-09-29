#include <assert.h>
#include <numeric>

#include "gamma_core.h"
#include "gamma.h"
#include "process.h"

//! Simulation: gamma_core constructor when just alpha is provided.
gamma_core::gamma_core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
    int n_gamma_cats, double alpha) : core(ost, lambda, p_tree, max_family_size, total_n_families, rootdist_vec),
    _gamma_cat_probs(n_gamma_cats), _lambda_multipliers(n_gamma_cats) {

    if (!rootdist_vec.empty()) {
        _rootdist_bins.push_back(rootdist_vec); // just 1st element
    }

    else {
        _rootdist_vec = uniform_dist(total_n_families, 1, max_family_size); // the user did not specify one... using uniform from 1 to max_family_size!
        _rootdist_bins.push_back(_rootdist_vec); // just 1st element (but we could specify different root dists for each lambda bin)
    }

    if (n_gamma_cats > 1) {
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference
        auto cats = weighted_cat_draw(total_n_families, _gamma_cat_probs);
        _gamma_cats = *cats;
        delete cats;
    }

    for (auto i = _gamma_cat_probs.begin(); i != _gamma_cat_probs.end(); ++i) {
        cout << "Should be all the same probability: " << *i << endl;
    }

    for (auto i = _lambda_multipliers.begin(); i != _lambda_multipliers.end(); ++i) {
        cout << "Lambda multiplier (rho) is: " << *i << endl;
    }

    for (auto i = _gamma_cats.begin(); i != _gamma_cats.end(); ++i) {
        cout << "Gamma category is: " << *i << endl;
    }

}

//! Simulation: core constructor when alpha is not provided and membership is provided directly, in addition to the lambda multipliers (something the user decides to use)
gamma_core::gamma_core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
    vector<int>& cats, vector<double>&mul) : core(ost, lambda, p_tree, max_family_size, total_n_families, rootdist_vec),
    _gamma_cats(cats), _lambda_multipliers(mul) {}

gamma_core::~gamma_core()
{
    for (size_t i = 0; i < _inference_bundles.size(); ++i)
        _inference_bundles[i].clear();
    _inference_bundles.clear();
}

void gamma_core::print_results(std::ostream& ost)
{
    for (size_t i = 0; i < _inference_bundles.size(); ++i)
    {
        for (size_t k = 0; k < _inference_bundles[i].results.size(); ++k)
        {
            family_info_stash& r = _inference_bundles[i].results[k];
            ost << r.family_id << "\t" << r.lambda_multiplier << "\t" << r.category_likelihood << "\t" << r.family_likelihood << endl;
        }
    }
}


void gamma_core::initialize_with_alpha(int n_gamma_cats, int n_families, double alpha)
{
    adjust_n_gamma_cats(n_gamma_cats);
    adjust_family_gamma_membership(n_families);
    set_alpha(alpha, n_families);
}

void gamma_core::initialize_without_alpha(int n_gamma_cats, int n_families, vector<double> lambda_multipliers, std::vector<int> gamma_cats)
{
    assert(lambda_multipliers.size() == n_gamma_cats);
    assert(gamma_cats.size() == n_families);
    //    adjust_n_gamma_cats(n_gamma_cats);
    //    adjust_family_gamma_membership(n_families);
    //    set_alpha(alpha, n_families);

    set_lambda_multipliers(lambda_multipliers);
    set_gamma_cats(gamma_cats);
}

//! Resize all gamma-related vectors according to provided number (integer) of gamma categories
void gamma_core::adjust_n_gamma_cats(int n_gamma_cats) {
    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
}

//! Resize gamma_cats vector that assigns gamma class membership of families to be inferred/simulated
void gamma_core::adjust_family_gamma_membership(int n_families) {
    _gamma_cats.resize(n_families);
}

//! Set alpha for gamma distribution
void gamma_core::set_alpha(double alpha, int n_families) {
    _alpha = alpha;
    if (_gamma_cats.size() > 1)
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference

    vector<int>* cats = weighted_cat_draw(n_families, _gamma_cat_probs);
    _gamma_cats = *cats;
    delete cats;

    for (std::vector<double>::iterator it = _gamma_cat_probs.begin(); it != _gamma_cat_probs.end(); ++it) {
        cout << "Gamma cat prob is : " << *it << endl;
    }
}

//! Set lambda multipliers for each gamma category
void gamma_core::set_lambda_multipliers(std::vector<double> lambda_multipliers) {
    _lambda_multipliers = lambda_multipliers;
}

//! Set lambda bins (each int is a vector pointing to a gamma category)
void gamma_core::set_gamma_cats(std::vector<int> gamma_cats) {
    _gamma_cats = gamma_cats;
}

void gamma_core::start_inference_processes() {

    for (int i = 0; i < _p_gene_families->size(); ++i) {
        gamma_bundle bundle;

        cout << "Started inference bundle " << i + 1 << endl;

        for (int j = 0; j < _lambda_multipliers.size(); ++j) {
            //            double lambda_bin = _gamma_cat_probs[j];
            inference_process *p_new_process = new inference_process(_ost, _p_lambda, _lambda_multipliers[j], _p_tree, _max_family_size, _max_root_family_size, &_p_gene_families->at(i), _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
            bundle.add(p_new_process);

            cout << "  Started inference process " << j + 1 << endl;

        }

        _inference_bundles.push_back(bundle);
    }
}

//! Populate _processes (vector of processes)
simulation_process* gamma_core::create_simulation_process(int family_number) {
    double lambda_bin = _gamma_cats[family_number];
    return new simulation_process(_ost, _p_lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
}

//! Infer bundle
double gamma_core::infer_processes() {

    initialize_rootdist_if_necessary();

    equilibrium_frequency eq(_rootdist_vec);
    std::vector<double> all_bundles_likelihood(_inference_bundles.size());

    for (int i = 0; i < _inference_bundles.size(); ++i) {
        cout << endl << "About to prune a gamma bundle." << endl;
        _inference_bundles[i].prune(_gamma_cat_probs, &eq);

        all_bundles_likelihood[i] = 0;
        double family_likelihood = 0;
        for (size_t k = 0; k < _inference_bundles[i].results.size(); ++k)
        {
            _inference_bundles[i].results[k].family_likelihood += _inference_bundles[i].results[k].category_likelihood;
            _inference_bundles[i].results[k].family_id = i;
        }
        all_bundles_likelihood[i] = _inference_bundles[i].results[0].family_likelihood;

        cout << "Likelihood of family " << i << " = " << all_bundles_likelihood[i] << std::endl;
    }

    double multi = std::accumulate(all_bundles_likelihood.begin(), all_bundles_likelihood.end(), 1.0, std::multiplies<double>());

    cout << "Final answer: " << multi << std::endl;

    return multi;
}

