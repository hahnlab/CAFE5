#ifndef FAMILY_H
#define FAMILY_H

#include <vector>
#include <iosfwd>
#include <map>
//#include "probability.h"

class clade;
class lambda;

typedef std::map<clade *, int> trial;

trial * simulate_family_from_root_size(clade *tree, int root_family_size, int max_family_size, lambda * p_lambda);

std::vector<std::vector<trial *> > simulate_families_from_distribution(clade *p_tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda);

std::vector<trial *> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda);

#endif
