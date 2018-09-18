#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H

#include "process.h"
#include "core.h"
#include <map>

class matrix_cache;

class reconstruction_process : public process {

    gene_family *_gene_family;
    matrix_cache *_p_calc;

    /// Filled out when reconstruct() is called (NULL before then)
    root_equilibrium_distribution* _p_prior;

    void reconstruct_internal_node(const clade * c, lambda * sl);
    void reconstruct_leaf_node(const clade * c, lambda * sl);
    void reconstruct_root_node(const clade * c);
    std::map<const clade *, std::vector<int> > all_node_Cs;

    /// Ls hold a probability for each family size (values are probabilities of any given family size)
    clademap<std::vector<double>> all_node_Ls;
    clademap<int> reconstructed_states;
    clademap<family_size_change> increase_decrease_map;
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
    std::vector<const clade *> get_taxa();

    void print_reconstruction(std::ostream & ost, cladevector& order);

    increase_decrease get_increases_decreases(cladevector& order, double pvalue);

    clademap<int> get_reconstructed_states() const
    {
        return reconstructed_states;
    }
    void operator()(const clade *c);

    static clademap<double> get_weighted_averages(std::vector<reconstruction_process *> m,
        const std::vector<double>& _gamma_cat_probs);

    std::string get_family_id() const;

    std::vector<double> get_L(const clade *c) const
    {
        return all_node_Ls.at(c);
    }
};

void compute_increase_decrease(clademap<int>& input, clademap<family_size_change>& output);
void compute_increase_decrease(clademap<double>& input, clademap<family_size_change>& output);
std::ostream& operator<<(std::ostream & ost, const increase_decrease& val);

#endif
