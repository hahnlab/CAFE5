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

//! Set the family size of a node to a random value, using parent's family size */
/*!
  Starting from 0, the gene family size of the child (c) is increased until the cumulative probability of c (given the gene family size s of the parent) exceeds a random draw from a uniform distribution. When this happens, the last c becomes the child's gene family size.
 
  Note that the smaller the draw from the uniform, the higher the chance that c will be far away from s.
*/
void set_node_familysize_random(clade *node) {

    if (node->is_root()) { return; } // if node is root, we do nothing

    /* Drawing random number from uniform */
    //  std::default_random_engine gen(static_cast<long unsigned int>(time(0)));
    // double rnd = dis(gen);
    double rnd = unifrnd(); // Ben: the smaller rnd is, the further away from the parent family size will c be; the larger rnd is, the closer c will be of the parent family size (?)
    cout << "Max family size: " << _max_family_size << " and rnd = " << rnd << endl;
    cout << "Branch length: " << node->get_branch_length() << endl;
    double cumul = 0;
    int parent_family_size = family_sizes[node->get_parent()];
    int c = 0; // c is the family size we will go to
    
    if (parent_family_size > 0) {
        for (; c < _max_family_size - 1; c++) {
            double prob = the_probability_of_going_from_parent_fam_size_to_c(_lambda, node->get_branch_length(), parent_family_size, c);
            cumul += prob;
            cout << "Probability value for " << c << ": " << prob << " (cumulative " << cumul << ")" << endl;
        
            if (cumul >= rnd) {
                cout << "Stopping" << endl;
                break;
            }
        }

//    if (c == 349)
//    {
//      double to5 = the_probability_of_going_from_parent_fam_size_to_c(_lambda, node->get_branch_length(), parent_family_size, 5);
//      cout << "Setting " << node->get_taxon_name() << " to: " << c << endl;
//      cout << "Probability of going from " << parent_family_size << " to 5 is " << to5 << endl;
//    }
    }
    
    family_sizes[node] = c;
}

//! Simulate gene family evolution from a single root family size 
/*!
  Given a root gene family size and a lambda, simulates gene family counts for all nodes in the tree.
  Returns family_sizes map, key = pointer to node, value = gene family size
*/
map<clade *, int> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda) {

  family_sizes.clear(); // resets map {clade: family size}
  _max_family_size = max_family_size; // we could pass these variables values in a different way
  _lambda = lambda;
  for (int t = 0; t < num_trials; ++t) {
    family_sizes[tree] = root_family_size; // set root family size
    tree->apply_prefix_order(set_node_familysize_random);
    //tree->apply_to_descendants(set_node_familysize_random); // setting the family size (random) of the descendants
  }

  return family_sizes;
}

//! Simulate gene family evolution from a distribution of root family sizes
/*!
  Given a root family size distribution (provided as a map read using read_rootdist()), simulates gene family counts for all nodes in the tree, and the number of times specified by the distribution (e.g., 10 families with root size 2, 5 families with root size 3, etc.)
*/
vector<trial> simulate_families_from_distribution(clade *tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda) {

    vector<trial> result;
    std::map<int, int>::const_iterator it = root_dist.begin();
    // Ben: while loop visits all keys in root_dist map (?)
    while (it != root_dist.end()) {
        int root_size = it->first;
        cout << "Root size: " << root_size << endl;
        int n_fams_this_size = it->second;
	
        for (int i = 0; i<n_fams_this_size; ++i) {
            result.push_back(simulate_families_from_root_size(tree, num_trials, root_size, max_family_size, lambda));
            it++;
        }
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
