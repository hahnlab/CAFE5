#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

class base_model : public model {
    std::vector<inference_process *> processes;
    virtual simulation_process* create_simulation_process(int family_number);
public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, clade *p_tree, vector<gene_family> *p_gene_families, int max_family_size, int max_root_family_size) :
        model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size) {}

    virtual void start_inference_processes();
    virtual double infer_processes(root_equilibrium_distribution *prior);

    virtual std::string name() {
        return "Base";
    }
    virtual ~base_model();

    virtual void print_results(std::ostream& ost);

    virtual std::vector<double> initial_guesses();
    void set_current_guesses(double * guesses);
};

#endif
