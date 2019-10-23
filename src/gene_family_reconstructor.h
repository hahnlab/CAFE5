#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H

#include "core.h"
#include <map>

class matrix_cache;

void reconstruct_leaf_node(const clade * c, const lambda * _lambda, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int _max_family_size, const gene_family* _gene_family, const matrix_cache *_p_calc);
void reconstruct_at_node(const clade *c, const lambda *_lambda, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int max_family_size, int max_root_family_size, const matrix_cache* p_calc, const root_equilibrium_distribution* p_prior, const gene_family *p_family);
void reconstruct_internal_node(const clade * c, const lambda * _lambda, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int _max_family_size, const matrix_cache *_p_calc);

/// Given a gene gamily and a tree, reconstructs the most likely values at each node on tree. Used in the base model to calculate values for each
/// gene family. Also used in a gamma bundle, one for each gamma category. Differences are represented by the lambda multiplier.
void reconstruct_gene_family(const lambda* lambda, const clade *p_tree,
    int max_family_size,
    int max_root_family_size,
    const gene_family *gf,
    matrix_cache *p_calc,
    root_equilibrium_distribution* p_prior, clademap<int>& reconstructed_states);

string newick_node(const clade *node, const cladevector& order, std::function<std::string(const clade *c)> textwriter);

double compute_viterbi_sum(const clade* c, const gene_family& family, const reconstruction* rec, int max_family_size, const matrix_cache& cache, const lambda* p_lambda);

void print_branch_probabilities(std::ostream& ost, const cladevector& order, const vector<gene_family>& gene_families, const std::map<std::string, clademap<double>>& branch_probabilities);

#endif
