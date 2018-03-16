#ifndef GAMMA_BUNDLE_H
#define GAMMA_BUNDLE_H

class inference_process;
class reconstruction_process;
class root_equilibrium_distribution;
class matrix_cache;
class lambda;
class clade;
class gene_family;

#include <vector>
#include <iosfwd>
#include <map>

class inference_process_factory
{
    std::ostream & _ost;
    lambda* _lambda;
    clade *_p_tree;
    int _max_family_size;
    int _max_root_family_size;
    std::vector<int> _rootdist_vec; // distribution of relative values. probability can be found by dividing a single value by the total of all values
    int _root_size; // will be drawn from _rootdist_vec by process itself
    gene_family *_family;
public:

    inference_process_factory(std::ostream & ost, lambda* lambda, clade *p_tree, int max_family_size,
        int max_root_family_size, std::vector<int> rootdist);

    void set_gene_family(gene_family *family) {
        _family = family;
    }

    inference_process* operator()(double lambda_multiplier);

    reconstruction_process* create_reconstruction_process(double lambda_multiplier);
};

//! One gamma bundle per family
//! Should reconstruct values for all gamma category probabilities
class gamma_bundle {
    std::vector<inference_process *> _inf_processes;
    std::vector<reconstruction_process *> _rec_processes;

    std::map<clade *, double> reconstruction;

public:
    gamma_bundle(inference_process_factory& factory, std::vector<double> lambda_multipliers);
    ~gamma_bundle();

    void clear();

    bool prune(const std::vector<double>& gamma_cat_probs, root_equilibrium_distribution *eq_freq,
        matrix_cache& calc, std::vector<double>& cat_likelihoods);

    void reconstruct(const std::vector<double>& _gamma_cat_probs);

    double get_lambda_likelihood(int family_id);

    void set_values(matrix_cache *, root_equilibrium_distribution*);

    void print_reconstruction(std::ostream& ost, std::vector<clade *> order);
    void print_increases_decreases(std::ostream& ost, std::vector<clade *> order);

    std::vector<clade *> get_taxa();
};

#endif
