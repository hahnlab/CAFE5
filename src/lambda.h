#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>

class clade;
class gene_family;

struct lambda_search_params
{
	lambda_search_params(clade *pt, std::vector<gene_family> f, int mfs) : ptree(pt), families(f), max_family_size(mfs)
	{

	}
	clade *ptree;
	std::vector<gene_family> families;
	int max_family_size;
};


std::vector<double> get_posterior(std::vector<gene_family> gene_families, int max_family_size, double lambda, clade *p_tree);
double find_best_lambda(lambda_search_params *params);

#endif
