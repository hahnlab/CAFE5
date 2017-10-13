#include <iostream>
#include <ostream>
#include <random>

#include "process.h"
#include "lambda.h"
#include "clade.h"
#include "probability.h"

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

simulation_process::simulation_process(std::ostream &ost, lambda* lambda, double lambda_multiplier, clade *p_tree, int max_family_size,
	int max_root_family_size, std::vector<int> rootdist, int family_number) : process(ost, lambda, lambda_multiplier, p_tree,
		max_family_size, max_root_family_size, rootdist) {

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
#if 0
            cout << "Using provided root distribution ";
            for (auto i : _rootdist_vec)
                cout << i << " ";
            cout << endl;
#endif
            _max_family_size_sim = *std::max_element(_rootdist_vec.begin(), _rootdist_vec.end());
            _root_size = _rootdist_vec[family_number];
    }

	cout << "_root_size is " << _root_size << endl;
}

//! Run process' simulation
void simulation_process::run_simulation() {
	single_lambda *sl = dynamic_cast<single_lambda*>(_lambda);	// we don't support multiple lambdas yet
	double lambda_m = sl->get_single_lambda() * _lambda_multiplier;
	_my_simulation = simulate_family_from_root_size(_p_tree, _root_size, _max_family_size_sim, lambda_m);
}

//! Prune process
std::vector<double> inference_process::prune() {
	single_lambda *sl = dynamic_cast<single_lambda*>(_lambda);	// we don't support multiple lambdas yet

	cout << endl << "Max root family size is: " << _max_root_family_size << endl;
	cout << "Max family size is: " << _max_family_size << endl;
	cout << "Lambda multiplier is: " << _lambda_multiplier << endl;

    single_lambda *multiplier = sl->multiply(_lambda_multiplier);
	likelihood_computer pruner(_max_root_family_size, _max_family_size, multiplier, _p_gene_family); // likelihood_computer has a pointer to a gene family as a member, that's why &(*p_gene_families)[0]

	cout << "  About to prune process." << endl;
	_p_tree->apply_reverse_level_order(pruner);

    delete multiplier;
	return pruner.get_likelihoods(_p_tree); // likelihood of the whole tree = multiplication of likelihood of all nodes

}

//! Printing process' simulation
void simulation_process::print_simulation(std::ostream & ost) {

    // Printing gene counts
	for (trial::iterator it = _my_simulation->begin(); it != _my_simulation->end(); ++it) {
		ost << it->second << "\t";
	}

	ost << _lambda_multiplier << endl;
}

//! Return simulation
trial * simulation_process::get_simulation() {
	return _my_simulation;
}

