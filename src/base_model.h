#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

class gene_family_reconstructor;
class matrix_cache;

/*! @brief A Base model can simulate families or estimate lambdas and error models.

    The estimation method creates a best guess for the missing lambda and, optionally,
    epsilon values, by using an @optimizer to provide guesses for the missing values
    and calculating a score for each one.

    \defgroup base Base Model

 */
class base_model : public model {
    double simulation_lambda_multiplier = 1.0;

public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, const clade *p_tree, const vector<gene_family>* p_gene_families,
        int max_family_size, int max_root_family_size, error_model *p_error_model);

    virtual double infer_family_likelihoods(root_equilibrium_distribution *prior, const std::map<int, int>& root_distribution_map, const lambda *p_lambda);

    virtual std::string name() {
        return "Base";
    }

    virtual void write_family_likelihoods(std::ostream& ost);

    virtual inference_optimizer_scorer *get_lambda_optimizer(user_data& data);

    virtual reconstruction* reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior);

    virtual void prepare_matrices_for_simulation(matrix_cache& cache);

    //! Return the lambda provided by the user, but with a small randomizer
    /// provided by simulation_lambda_multiplier. That is modified by calling \ref perturb_lambda.
    virtual lambda* get_simulation_lambda(const user_data& data);

    //! Sets simulation_lambda_multiplier by a normal distribution mean=1, stddev = .3
    void perturb_lambda();

};

//! \ingroup base Base Model
class base_model_reconstruction : public reconstruction
{
public:

    base_model_reconstruction(size_t num_families) : families(num_families)
    {

    }

    void print_reconstructed_states(std::ostream& ost, const cladevector& order, const std::vector<gene_family>& gene_families, const clade *p_tree) override;

    std::vector<reconstructed_family<int>> families;

    void print_increases_decreases_by_family(std::ostream& ost, const cladevector& order, const std::vector<double>& pvalues) override;
    void print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order) override;
};

#endif
