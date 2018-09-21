#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H

#include "core.h"
#include <map>

class matrix_cache;

class gene_family_reconstructor {
    const lambda* _lambda;
    const clade *_p_tree;
    const gene_family *_gene_family;

    double _lambda_multiplier;
    int _max_family_size;
    int _max_root_family_size;
    matrix_cache *_p_calc;

    /// Filled out when reconstruct() is called (NULL before then)
    root_equilibrium_distribution* _p_prior;

    void reconstruct_internal_node(const clade * c, lambda * sl);
    void reconstruct_leaf_node(const clade * c, lambda * sl);
    void reconstruct_root_node(const clade * c);
    clademap<std::vector<int>> all_node_Cs;

    /// Ls hold a probability for each family size (values are probabilities of any given family size)
    clademap<std::vector<double>> all_node_Ls;
    clademap<int> reconstructed_states;
    clademap<family_size_change> increase_decrease_map;
public:
    void reconstruct();

    gene_family_reconstructor(std::ostream & ost, lambda* lambda, double lambda_multiplier, const clade *p_tree,
        int max_family_size,
        int max_root_family_size,
        const gene_family *gf,
        matrix_cache *p_calc,
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

    static clademap<double> get_weighted_averages(std::vector<gene_family_reconstructor *> m,
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
