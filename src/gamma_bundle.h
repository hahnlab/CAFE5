#ifndef GAMMA_BUNDLE_H
#define GAMMA_BUNDLE_H

class inference_process;
class gene_family_reconstructor;
class root_equilibrium_distribution;
class matrix_cache;
class lambda;
class clade;
class gene_family;

#include <vector>
#include <iosfwd>
#include <map>

#include "core.h"

//! One gamma bundle per family
//! Should reconstruct values for all gamma category probabilities
class gamma_bundle {
//    std::vector<gene_family_reconstructor *> _rec_processes;
    clademap<family_size_change> _increase_decrease_map;
    clademap<double> reconstruction;
    const clade *_p_tree;
    const gene_family *_p_gene_family;

    std::vector<double> _lambda_multipliers;
    std::vector<double> _category_likelihoods;

    int _max_family_size;
    int _max_root_family_size;

    vector<clademap<int>> reconstructed_states;
    vector<clademap<family_size_change>> increase_decrease_map;
    std::ostream& _ost;
    const lambda *_p_lambda;
public:
    gamma_bundle(std::vector<double> lambda_multipliers, const clade *p_tree, const gene_family *p_gene_family,
        std::ostream & ost, const lambda* lambda, int max_family_size, int max_root_family_size);

    /// used to copy known values into the reconstructor
    gamma_bundle(vector<clademap<int>>& rs, clademap<double> rc, clademap<family_size_change> idm,
        const clade *p_tree, const gene_family *p_gene_family) :
        _increase_decrease_map(idm),
        reconstruction(rc),
        _p_tree(p_tree),
        _p_gene_family(p_gene_family),
        reconstructed_states(rs),
        _ost(cout)
        {}

    bool prune(const std::vector<double>& gamma_cat_probs, root_equilibrium_distribution *eq_freq,
        matrix_cache& calc, const lambda *p_lambda);

    void reconstruct(const std::vector<double>& _gamma_cat_probs, matrix_cache *calc, root_equilibrium_distribution*prior);

    double get_lambda_likelihood(int family_id);

    void print_reconstruction(std::ostream& ost, cladevector& order);

    increase_decrease get_increases_decreases(cladevector& order, double pvalue);

    std::vector<const clade *> get_taxa();

    std::vector<double> get_category_likelihoods() const
    {
        return _category_likelihoods;
    }

    std::string get_family_id() const;
    std::string get_reconstructed_states(const clade *node) const;
};

#endif
