#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H

#include "process.h"
#include "core.h"

class matrix_cache;

class reconstruction_process : public process {

    gene_family *_gene_family;
    matrix_cache *_p_calc;

    /// Filled out when reconstruct() is called (NULL before then)
    root_equilibrium_distribution* _p_prior;

    void reconstruct_internal_node(clade * c, lambda * sl);
    void reconstruct_leaf_node(clade * c, lambda * sl);
    void reconstruct_root_node(clade * c);
    std::map<clade *, std::vector<int> > all_node_Cs;

    /// Ls hold a probability for each family size (values are probabilities of any given family size)
    std::map<clade *, std::vector<double> > all_node_Ls;

    std::map<clade *, int> reconstructed_states;

    std::map<clade *, family_size_change> increase_decrease_map;
public:
    void reconstruct();

    reconstruction_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree,
        int max_family_size,
        int max_root_family_size, std::vector<int> rootdist,
        gene_family *gf,
        matrix_cache *calc,
        root_equilibrium_distribution* p_prior);

    void set_values(matrix_cache *calc, root_equilibrium_distribution* prior)
    {
        _p_prior = prior;
        _p_calc = calc;
    }
    std::vector<clade *> get_taxa();

    void print_reconstruction(std::ostream & ost, std::vector<clade *>& order);

    increase_decrease get_increases_decreases(std::vector<clade *>& order);

    std::map<clade *, int> get_reconstructed_states() const
    {
        return reconstructed_states;
    }
    void operator()(clade *c);

    static std::map<clade *, double> get_weighted_averages(std::vector<reconstruction_process *> m,
        const std::vector<double>& _gamma_cat_probs);

    std::string get_family_id() const;

    std::vector<double> get_L(clade *c) const
    {
        return all_node_Ls.at(c);
    }
};

void compute_increase_decrease(std::map<clade *, int>& input, std::map<clade *, family_size_change>& output);
void compute_increase_decrease(std::map<clade *, double>& input, std::map<clade *, family_size_change>& output);
std::ostream& operator<<(std::ostream & ost, const increase_decrease& val);

#endif
