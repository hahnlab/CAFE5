#include <iostream>
#include <ostream>
#include <random>

#include "process.h"
#include "lambda.h"
#include "clade.h"
#include "probability.h"
#include "core.h"

/* START: Drawing random root size from uniform */
template<typename itr, typename random_generator>
itr select_randomly(itr start, itr end, random_generator& g) {
	std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
	std::advance(start, dis(g)); // advances iterator (start) by dis(g), where g is a seeded mt generator
	return start;
}

template<typename itr>
itr select_randomly(itr start, itr end) {
	static std::random_device rd; // randomly generates uniformly distributed ints (seed)
	static std::mt19937 gen(rd()); // seeding Mersenne Twister generator
	return select_randomly(start, end, gen); // plug in mt generator to advance our container and draw random element from it
}
/* END: Drawing random root size from uniform */


//using namespace std;

simulation_process::simulation_process(std::ostream &ost, lambda* lambda, double lambda_multiplier, const clade *p_tree, int max_family_size,
    int max_root_family_size, std::vector<int> rootdist, int family_number, error_model* p_error_model ) : process(ost, lambda, lambda_multiplier, p_tree,
        max_family_size, max_root_family_size, rootdist), _p_error_model(p_error_model) {

    // generating uniform root distribution when no distribution is provided 
    if (_rootdist_vec.empty()) {
        _max_family_size_sim = 100;

        cout << "Max family size to simulate: " << _max_family_size_sim << endl;
        _rootdist_vec.resize(_max_family_size_sim);

        for (size_t i = 0; i < _rootdist_vec.size(); ++i)
            _rootdist_vec[i] = i;

        _root_size = *select_randomly(_rootdist_vec.begin(), _rootdist_vec.end()); // getting a random root size from the provided (core's) root distribution
    }
    else {
        _max_family_size_sim = get_max_family_size_to_simulate();
        _root_size = _rootdist_vec[family_number];
    }

}

int simulation_process::get_max_family_size_to_simulate() const
{
    return 2 * *std::max_element(_rootdist_vec.begin(), _rootdist_vec.end());
}

//! Run process' simulation
void simulation_process::run_simulation() {
    if (_lambda == NULL)
        throw std::runtime_error("No lambda specified for simulation");

    lambda *multiplier = _lambda->multiply(_lambda_multiplier);
	_my_simulation = simulate_family_from_root_size(_p_tree, _root_size, _max_family_size_sim, multiplier, _p_error_model);
}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by creating a likelihood_computer and calling it on all modes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_process::prune(matrix_cache& calc) {

    unique_ptr<lambda> multiplier(_lambda->multiply(_lambda_multiplier));

	likelihood_computer pruner(_max_root_family_size, _max_family_size, multiplier.get(), *_p_gene_family, calc, _p_error_model); // likelihood_computer has a pointer to a gene family as a member, that's why &(*p_gene_families)[0]
	_p_tree->apply_reverse_level_order(pruner);

	return pruner.get_likelihoods(_p_tree); // likelihood of the whole tree = multiplication of likelihood of all nodes
}

std::string inference_process::family_id() const
{
    return _p_gene_family->id();
}

//! Printing process' simulation
void simulation_process::print_simulation(std::ostream & ost, int index) {

    // Printing gene counts
	for (trial::iterator it = _my_simulation->begin(); it != _my_simulation->end(); ++it) {
		ost << it->second << "\t";
	}

	ost << _lambda_multiplier << '\t' << index << endl;
}

//! Return simulation
trial * simulation_process::get_simulation() {
	return _my_simulation;
}

