
#include <vector>

class clade;

void simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda);
void simulate_families_from_distribution(clade *tree, int num_trials, std::vector<int> root_dist, int max_family_size, double lambda);