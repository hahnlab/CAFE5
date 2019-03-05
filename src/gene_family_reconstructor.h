#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H

#include "core.h"
#include <map>

class matrix_cache;

/// Given a gene gamily and a tree, reconstructs the most likely values at each node on tree. Used in the base model to calculate values for each
/// gene family. Also used in a gamma bundle, one for each gamma category. Differences are represented by the lambda multiplier.
class gene_family_reconstructor {
    const gene_family *_gene_family;
    matrix_cache *_p_calc;
    /// Filled out when reconstruct() is called (NULL before then)
    root_equilibrium_distribution* _p_prior;

    const lambda* _lambda;
    const clade *_p_tree;

    int _max_family_size;
    int _max_root_family_size;

    double _lambda_multiplier;

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

    gene_family_reconstructor(std::ostream & ost, const lambda* lambda, double lambda_multiplier, const clade *p_tree,
        int max_family_size,
        int max_root_family_size,
        const gene_family *gf,
        matrix_cache *p_calc,
        root_equilibrium_distribution* p_prior);

    /// used to copy known values into the reconstructor
    gene_family_reconstructor(const gene_family *family, const clade *p_tree, clademap<int> states) {
        _gene_family = family;
        reconstructed_states = states;
        _p_tree = p_tree;
    }

    void set_values(matrix_cache *calc, root_equilibrium_distribution* prior)
    {
        _p_prior = prior;
        _p_calc = calc;
    }
    cladevector get_taxa() const;

    void print_reconstruction(std::ostream & ost, cladevector& order);

    increase_decrease get_increases_decreases(cladevector& order, double pvalue);

    cladevector get_nodes();
    int get_reconstructed_value(const clade *node)
    {
        return reconstructed_states[node];
    }
    void operator()(const clade *c);

    std::vector<double> get_L(const clade *c) const
    {
        return all_node_Ls.at(c);
    }
    std::string get_reconstructed_states(const clade *node) const;



};

void compute_increase_decrease(clademap<int>& input, clademap<family_size_change>& output);
void compute_increase_decrease(clademap<double>& input, clademap<family_size_change>& output);
std::ostream& operator<<(std::ostream & ost, const increase_decrease& val);

template<class T>
string newick_node(const clade *node, const cladevector& order, const T *reconstructor)
{
    auto node_id = distance(order.begin(), find(order.begin(), order.end(), node));
    ostringstream ost;
    if (node->is_leaf())
        ost << node->get_taxon_name();
    else
    {
        ost << node_id;
    }
    ost << "_" << reconstructor->get_reconstructed_states(node);
    if (!node->is_root())
        ost << ':' << node->get_branch_length();
    return ost.str();
}

#endif
