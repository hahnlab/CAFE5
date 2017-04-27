#ifndef POISSON_H
#define POISSON_H

#include <vector>

class clade;
class gene_family;

std::vector<double> find_poisson_lambda(clade* tree, std::vector<gene_family> gene_families);
std::vector<double> get_prior_rfsize_poisson_lambda(int min_family_size, int max_family_size, double poisson_lambda);

#endif
