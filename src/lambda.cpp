#include <vector>
#include <map>
#include <cmath>
#include "lambda.h"
#include "clade.h"
#include "fminsearch.h"
#include "utils.h"
#include "probability.h"

using namespace std;


/// score of a lambda is the -log likelihood of the most likely resulting family size
double calculate_lambda_score(double* plambda, void* args)
{
	lambda_search_params *param = (lambda_search_params *)args;

	vector<double> posterior = get_posterior(param->families, param->max_family_size, *plambda, param->ptree);
	return -log(*max_element(posterior.begin(), posterior.end()));
}


double find_best_lambda(lambda_search_params *params)
{
	int lambda_len = 1;
	FMinSearch* pfm;
	pfm = fminsearch_new_with_eq(calculate_lambda_score, lambda_len, &params);
	pfm->tolx = 1e-6;
	pfm->tolf = 1e-6;
	double result;
	fminsearch_min(pfm, &result);
	double *re = fminsearch_get_minX(pfm);
	return *re;
}
