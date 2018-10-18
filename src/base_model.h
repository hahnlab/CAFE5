#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

class gene_family_reconstructor;
class matrix_cache;

class base_model : public model {
    std::vector<inference_process *> processes;
    virtual simulation_process* create_simulation_process(int family_number);
public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, const clade *p_tree, const vector<gene_family>* p_gene_families,
        int max_family_size, int max_root_family_size, std::map<int, int> * p_rootdist_map, error_model *p_error_model);

    virtual void start_inference_processes();
    virtual double infer_processes(root_equilibrium_distribution *prior);

    virtual std::string name() {
        return "Base";
    }
    virtual ~base_model();

    virtual void write_family_likelihoods(std::ostream& ost);

    virtual optimizer_scorer *get_lambda_optimizer(user_data& data);

    virtual reconstruction* reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior);
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
