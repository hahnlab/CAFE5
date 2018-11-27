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
#include "user_data.h"

//using namespace std;

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
