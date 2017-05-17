#include <vector>
#include <iosfwd>
#include <map>
#include "probability.h"

class clade;

typedef std::map<clade *, int> trial;

//vector<trial *> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda);
//Before: vector<trial> simulate_families_from_distribution(clade *tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda);
std::vector<std::vector<trial *> > simulate_families_from_distribution(clade *p_tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda);
//Before: void print_simulation(std::vector<trial> &sim, std::ostream& ost);


