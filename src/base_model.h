#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

class reconstruction_process;

class base_model : public model {
    std::vector<inference_process *> processes;
    virtual simulation_process* create_simulation_process(int family_number);
public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, clade *p_tree, vector<gene_family> *p_gene_families,
        int max_family_size, int max_root_family_size, std::map<int, int> * p_rootdist_map);

    virtual void start_inference_processes();
    virtual double infer_processes(root_equilibrium_distribution *prior);

    virtual std::string name() {
        return "Base";
    }
    virtual ~base_model();

    virtual void print_results(std::ostream& ost);

    virtual std::vector<double> initial_guesses();
    void set_current_guesses(double * guesses);

    virtual void reconstruct_ancestral_states(probability_calculator *p_calc, root_equilibrium_distribution* p_prior);
    reconstruction_process* create_reconstruction_process(int family_number, probability_calculator *p_calc, root_equilibrium_distribution* p_prior);
};

#endif
