#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

class gene_family_reconstructor;
class matrix_cache;

class base_model : public model {
    double simulation_lambda_multiplier = 1.0;

public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, const clade *p_tree, const vector<gene_family>* p_gene_families,
        int max_family_size, int max_root_family_size, error_model *p_error_model);

    virtual double infer_processes(root_equilibrium_distribution *prior, const std::map<int, int>& root_distribution_map, const lambda *p_lambda);

    virtual std::string name() {
        return "Base";
    }

    virtual void write_family_likelihoods(std::ostream& ost);

    virtual inference_optimizer_scorer *get_lambda_optimizer(user_data& data);

    virtual reconstruction* reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior);

    virtual void prepare_matrices_for_simulation(matrix_cache& cache);

    virtual lambda* get_simulation_lambda(const user_data& data);

    void perturb_lambda();

};


class base_model_reconstruction : public reconstruction
{
    void print_increases_decreases_by_family(std::ostream& ost, const std::vector<double>& pvalues);
    void print_increases_decreases_by_clade(std::ostream& ost);
public:
    base_model_reconstruction()
    {

    }

    ~base_model_reconstruction();

    vector<gene_family_reconstructor*> _rec_processes;

    void print_reconstructed_states(std::ostream& ost);
};

#endif
