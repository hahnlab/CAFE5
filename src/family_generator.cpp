#include <vector>
#include <map>
#include <random>
#include "clade.h"
#include "probability.h"
#include "family_generator.h"

static map<clade *, int> family_sizes; // key: clades, value: family size

std::default_random_engine gen(12);
std::uniform_real_distribution<> dis(0, 1); // draw random number from uniform distribution

/* Parameters to family size randomizer. TODO: Be more clever about passing these into set_node_familysize_random. Global variables are to be avoided */
int _max_family_size;
double _lambda;

/* Recursively set the family size of a node to a random value, using parent's family size */
void set_node_familysize_random(clade *node) {

  if (node->is_root()) { return; } // if node is root, we do nothing

  /* Drawing random number from uniform */
//  std::default_random_engine gen(static_cast<long unsigned int>(time(0)));
  double rnd = dis(gen);

  //cout << "Max family size: " << _max_family_size << " and rnd = " << rnd << endl;
  double cumul = 0;
  int parent_family_size = family_sizes[node->get_parent()];
  int c = 0; // c is the family size we will go to
  for (; c < _max_family_size - 1; c++) {
      cumul += the_probability_of_going_from_parent_fam_size_to_c(_lambda, node->get_branch_length(), parent_family_size, c);
	  if (cumul >= rnd)
	  {
		  break;
	  }
  }

  family_sizes[node] = c;
  node->apply_to_descendants(set_node_familysize_random); // recursion
}

map<clade *, int> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda) {

  family_sizes.clear(); // resets map {clade: family size}
  _max_family_size = max_family_size; // we could pass these variables values in a different way
  _lambda = lambda;
  for (int t = 0; t < num_trials; ++t) {
    family_sizes[tree] = root_family_size; // set root family size
    tree->apply_to_descendants(set_node_familysize_random); // setting the family size (random) of the descendants
  }

  return family_sizes;
}

vector<trial> simulate_families_from_distribution(clade *tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda) {

  vector<trial> result;

  std::map<int, int>::const_iterator it = root_dist.begin();
  while (it != root_dist.end())
  {
	  int root_dist = it->first;
	  cout << "Root dist: " << root_dist << endl;
	  int n_gene_families = it->second;
	  for (int i = 0; i<n_gene_families; ++i)
	    result.push_back(simulate_families_from_root_size(tree, num_trials, root_dist, max_family_size, lambda));
	  it++;
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
