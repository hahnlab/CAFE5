#ifndef POISSON_H
#define POISSON_H

#include <vector>

class clade;
class gene_family;

std::vector<double> find_poisson_lambda(clade* tree, std::vector<gene_family> gene_families);

#endif
