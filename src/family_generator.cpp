#include <vector>
#include <map>
#include "clade.h"

double unifrnd()
{
	return rand() / (RAND_MAX + 1.0);
}

double the_probability_of_going_from_parent_fam_size_to_c(int parent_size, int size)
{
	return 0.05;
}

map<clade *, int> family_sizes;
int _max_family_size;

/* Set the family size of a node to a random value */
void set_node_familysize_random(clade *node)
{
	if (!node->get_parent()) return;

	double rnd = unifrnd();
	double cumul = 0;
	int parent_family_size = family_sizes[node->get_parent()];
	int c = 0;
	for (; c < _max_family_size - 1; c++)
	{
		cumul += the_probability_of_going_from_parent_fam_size_to_c(parent_family_size, c);
		//cumul += square_matrix_get(pcnode->birthdeath_matrix, parent_family_size, c);
		if (cumul >= rnd) break;
	}
	family_sizes[node] = c;
	node->apply_to_descendants(set_node_familysize_random);
}

void simulate_families(clade *tree, int num_trials, std::vector<int> root_dist, int max_family_size)
{
	_max_family_size = max_family_size;
	for (int t = 0; t < num_trials; ++t)
	{
		for (int i = 1; i <= root_dist.size(); i++)
		{
			family_sizes[tree] = i;
			tree->apply_to_descendants(set_node_familysize_random);
		}
	}
}