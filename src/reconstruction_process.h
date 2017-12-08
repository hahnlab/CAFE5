#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H
#include "process.h"

class reconstruction_process : public process {

    gene_family *_gene_family;
    probability_calculator *_p_calc;

    /// Filled out when reconstruct() is called (NULL before then)
    root_equilibrium_distribution* _p_prior;

    void reconstruct_internal_node(clade * c, lambda * sl);
    void reconstruct_leaf_node(clade * c, lambda * sl);
    void reconstruct_root_node(clade * c);
    std::map<clade *, std::vector<int> > all_node_Cs;
    std::map<clade *, std::vector<double> > all_node_Ls;

    std::map<clade *, int> reconstructed_states;
public:
    void reconstruct();

    reconstruction_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree,
        int max_family_size,
        int max_root_family_size, std::vector<int> rootdist,
        gene_family *gf,
        probability_calculator *calc,
        root_equilibrium_distribution* p_prior);

    void set_values(probability_calculator *calc, root_equilibrium_distribution* prior)
    {
        _p_prior = prior;
        _p_calc = calc;
    }
    std::vector<clade *> get_taxa();

    void print_reconstruction(std::ostream & ost, std::vector<clade *>& order);
    std::map<clade *, int> get_reconstructed_states() const
    {
        return reconstructed_states;
    }
    void operator()(clade *c);

    static std::map<clade *, double> get_weighted_averages(std::vector<reconstruction_process *> m,
        const std::vector<double>& _gamma_cat_probs);

    std::string get_family_id() const;
};
#endif
