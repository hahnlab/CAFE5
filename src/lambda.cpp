#include <vector>
#include <map>
#include <cmath>
#include <iomanip>

#include "lambda.h"
#include "clade.h"
#include "fminsearch.h"
#include "utils.h"
#include "probability.h"
#include "core.h"

using namespace std;

/* START: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

std::vector<double> single_lambda::calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size)
{
	// cout << "Child node " << child->get_taxon_name() << " has " << probabilities.size() << " probabilities" << endl;
	std::vector<std::vector<double> > matrix = _p_calc->get_matrix(probabilities.size(), child->get_branch_length(), _lambda); // Ben: is _factors[child].size() the same as _max_root_family_size? If so, why not use _max_root_family_size instead?
	return matrix_multiply(matrix, probabilities, s_min_family_size, s_max_family_size, c_min_family_size, c_max_family_size);
}

std::string single_lambda::to_string()
{
    ostringstream ost;
    ost << setw(15) << setprecision(14) << _lambda;
    return ost.str();
}

std::vector<double> multiple_lambda::calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size)
{
	std::string nodename = child->get_taxon_name();
	int lambda_index = _node_name_to_lambda_index[nodename];
	double lambda = _lambdas[lambda_index];
	//cout << "Matrix for " << child->get_taxon_name() << endl;
	std::vector<std::vector<double> > matrix = _p_calc->get_matrix(probabilities.size(), child->get_branch_length(), lambda); // Ben: is _factors[child].size() the same as _max_root_family_size? If so, why not use _max_root_family_size instead?
	return matrix_multiply(matrix, probabilities, s_min_family_size, s_max_family_size, c_min_family_size, c_max_family_size);
}

void multiple_lambda::update(double* values)
{
    std::copy(values, values + _lambdas.size(), _lambdas.begin());
}

std::string multiple_lambda::to_string()
{
    ostringstream ost;
    ost << setw(15) << setprecision(14);
    for (size_t i = 0; i < _lambdas.size(); ++i)
    {
        ost << _lambdas[i];
        if (i != _lambdas.size() - 1) ost << ", ";
    }
    return ost.str();
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

/// score of a lambda is the -log likelihood of the most likely resulting family size
double calculate_lambda_score(double* p_lambda, void* args)
{
    std::pair<model *, root_equilibrium_distribution *>* vals = (std::pair<model *, root_equilibrium_distribution *>*)args;

    model *core = vals->first;
    root_equilibrium_distribution* dist = vals->second;

    core->set_current_guesses(p_lambda);
    core->start_inference_processes();

    return core->infer_processes(dist);
}

double* find_best_lambda(model * p_model, root_equilibrium_distribution *p_distribution)
{
    auto initial = p_model->initial_guesses();
	FMinSearch* pfm;
    std::pair<model *, root_equilibrium_distribution *> args(p_model, p_distribution);
	pfm = fminsearch_new_with_eq(calculate_lambda_score, initial.size(), &args);
	pfm->tolx = 1e-6;
	pfm->tolf = 1e-6;
    pfm->maxiters = 25;
	fminsearch_min(pfm, &initial[0]);
	double *re = fminsearch_get_minX(pfm);
    cout << "Best match" << (initial.size() == 1 ? " is: " : "es are: ") << setw(15) << setprecision(14);
    for (size_t i = 0; i<initial.size(); ++i)
        cout << re[i] << ',';
    cout << endl;

    return re;
}
