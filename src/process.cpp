#include <iostream>
#include <ostream>
#include <random>

#include "process.h"
#include "lambda.h"
#include "clade.h"
#include "probability.h"
#include "core.h"
#include "gene_family.h"
#include "root_distribution.h"

//using namespace std;

simulation_process::simulation_process(std::ostream &ost, lambda* lambda, double lambda_multiplier, const clade *p_tree,
    int max_root_family_size, int max_family_size_sim, int root_size, error_model* p_error_model ) : 
        process(ost, lambda, lambda_multiplier, p_tree, max_root_family_size), 
        _p_error_model(p_error_model),
        _max_family_size_sim(max_family_size_sim),
        _root_size(root_size)
{
}

//! Run process' simulation
trial* simulation_process::run_simulation(const matrix_cache& cache) {

    if (_lambda == NULL)
        throw std::runtime_error("No lambda specified for simulation");

    unique_ptr<lambda> multiplier(_lambda->multiply(_lambda_multiplier));
    if (_p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    auto *result = new trial();
    random_familysize_setter rfs(result, _max_family_size_sim, multiplier.get(), _p_error_model, cache);
    (*result)[_p_tree] = _root_size;
    _p_tree->apply_prefix_order(rfs); // this is where the () overload of random_familysize_setter is used

    return result;
}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by creating a likelihood_computer and calling it on all modes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_process::prune(matrix_cache& calc) {

    unique_ptr<lambda> multiplier(_lambda->multiply(_lambda_multiplier));

	likelihood_computer pruner(_max_root_family_size, _max_family_size, multiplier.get(), *_p_gene_family, calc, _p_error_model); // likelihood_computer has a pointer to a gene family as a member, that's why &(*p_gene_families)[0]
    pruner.initialize_memory(_p_tree);
	_p_tree->apply_reverse_level_order(pruner);

	return pruner.get_likelihoods(_p_tree); // likelihood of the whole tree = multiplication of likelihood of all nodes
}

std::string inference_process::family_id() const
{
    return _p_gene_family->id();
}
