#include <vector>
#include <map>
#include <random>

#include "clade.h"
#include "probability.h"
#include "family_generator.h"

/* TODO: Maybe change to a return value or something else, as global variables are to be avoided */
static map<clade *, int> family_sizes; // key: clades, value: family size

/* parameters to family size randomizer. TODO: Be more clever about passing these into set_node_familysize_random
 * Global variables are to be avoided */
int _max_family_size;
double _lambda;

/* Set the family size of a node to a random value */
void set_node_familysize_random(clade *node)
{
  if (node->get_parent() == NULL) { return; } // if node is root, get_parent() = NULL, and we do nothing

	std::default_random_engine gen(static_cast<long unsigned int>(time(0)));
	std::uniform_real_distribution<> dis(0, 1); // draw random number from uniform distribution
	double rnd = dis(gen);
	double cumul = 0;
	int parent_family_size = family_sizes[node->get_parent()];
	int c = 0;
	for (; c < _max_family_size - 1; c++)
	{
		cumul += the_probability_of_going_from_parent_fam_size_to_c(_lambda, node->get_branch_length(), parent_family_size, c);
		if (cumul >= rnd) break;
	}

	family_sizes[node] = c;
	node->apply_to_descendants(set_node_familysize_random);
}

map<clade *, int> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda)
{
	family_sizes.clear();
	_max_family_size = max_family_size;
	_lambda = lambda;
	for (int t = 0; t < num_trials; ++t)
	{
			family_sizes[tree] = root_family_size;
			tree->apply_to_descendants(set_node_familysize_random);
	}
	return family_sizes;
}

vector<trial> simulate_families_from_distribution(clade *tree, int num_trials, std::vector<int> root_dist, int max_family_size, double lambda)
{
	vector<trial> result;
	for (int i = 1; i <= root_dist.size(); i++)
	{
		result.push_back(simulate_families_from_root_size(tree, num_trials, root_dist[i], max_family_size, lambda));
	}
	return result;
}

void print_simulation(trial &sim, std::ostream& ost)
{
	trial::iterator it = sim.begin();
	for (; it != sim.end(); ++it)
	{
		ost << "#" << it->first->get_taxon_name() << endl;
	}
	for (it = sim.begin(); it != sim.end(); ++it)
	{
		ost << it->second << "\t";
	}
	ost << endl;

}

void print_simulation(std::vector<trial> &sim, std::ostream& ost)
{
	trial::iterator it = sim[0].begin();
	for (; it != sim[0].end(); ++it)
	{
		ost << "#" << it->first->get_taxon_name() << endl;
	}
	for (int i = 0; i < sim.size(); ++i)
	{
		for (it = sim[i].begin(); it != sim[i].end(); ++it)
		{
			ost << it->second << "\t";
		}
		ost << endl;
	}

}
