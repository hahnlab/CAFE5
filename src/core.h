#ifndef CORE_H
#define CORE_H

#include <set>

#include "clade.h"
#include "probability.h"

class simulation_process;
class inference_process;
class reconstruction_process;

struct family_info_stash {
    family_info_stash() : family_id(0), lambda_multiplier(0.0), category_likelihood(0.0), family_likelihood(0.0), 
        posterior_probability(0.0), significant(false) {}
    family_info_stash(int fam, double lam, double cat_lh, double fam_lh, double pp, bool signif) : 
        family_id(fam), lambda_multiplier(lam), category_likelihood(cat_lh), family_likelihood(fam_lh),
        posterior_probability(pp), significant(signif) {}
    int family_id;
    double lambda_multiplier;
    double category_likelihood;
    double family_likelihood;
    double posterior_probability;
    bool significant;
};

class branch_length_finder
{
    std::set<double> _result;
public:
    void operator()(clade *c);

    std::set<double> result() const
    {
        return _result;
    }

    double longest() const;
};

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r);

class model {
protected:
    std::ostream & _ost; 
    lambda *_p_lambda; // TODO: multiple lambdas for different branches
    clade *_p_tree;
    int _max_family_size;
    int _max_root_family_size;
    int _total_n_families_sim;
    vector<gene_family> * _p_gene_families;
    vector<int> _rootdist_vec; // in case the user wants to use a specific root size distribution for all simulations
    vector<vector<int> > _rootdist_bins; // holds the distribution for each lambda bin
    
    //! Simulations
    vector<simulation_process*> _sim_processes; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
    vector<reconstruction_process*> _rec_processes;

    void initialize_rootdist_if_necessary();

    std::vector<family_info_stash> results;

    error_model* _p_error_model;
public:
    //! Basic constructor
    model(): _ost(cout), _p_lambda(NULL), _p_tree(NULL), _p_gene_families(NULL), _total_n_families_sim(1) {}
    
    model(ostream & ost, 
        lambda* p_lambda, 
        clade *p_tree, 
        int max_family_size, 
        int total_n_families, 
        vector<int> rootdist_vec, 
        error_model *p_error_model) : 
        _ost(ost), _p_lambda(p_lambda), _p_tree(p_tree), _max_family_size(max_family_size), 
        _total_n_families_sim(total_n_families), _rootdist_vec(rootdist_vec), _p_error_model(p_error_model) {}
    
    model(lambda* p_lambda, 
        clade *p_tree, 
        vector<gene_family> *p_gene_families, 
        int max_family_size, 
        int max_root_family_size, 
        error_model *p_error_model) : 
        _ost(cout), _p_lambda(p_lambda), _p_tree(p_tree), _p_gene_families(p_gene_families), _max_family_size(max_family_size), 
        _max_root_family_size(max_root_family_size), _p_error_model(p_error_model) {}
    
    virtual ~model() {}
    //void estimate_processes(); 
    
    //! Setter methods
    void set_lambda(lambda *p_lambda);
    
    void set_tree(clade *p_tree);
    
    void set_gene_families(std::vector<gene_family> *p_gene_families);
    
    void set_max_sizes(int max_family_size, int max_root_family_size);
    
    void set_rootdist_vec(std::vector<int> rootdist_vec);
    
    void set_total_n_families_sim(int total_n_families_sim);
    
    //! Simulation methods
    void start_sim_processes();

    virtual simulation_process* create_simulation_process(int family_number) = 0;

    void simulate_processes();

    void print_processes(std::ostream& ost);

    //! Inference methods
    virtual void start_inference_processes() = 0;
    
    virtual double infer_processes(root_equilibrium_distribution *prior) = 0;  // return vector of likelihoods
    
    //! Printing methods
    void print_parameter_values();
    
    void adjust_family(ostream& ost);

    virtual std::string name() = 0;
    virtual void print_results(std::ostream& ost) = 0;

    virtual void reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior) = 0;
    virtual void print_reconstructed_states(std::ostream& ost) = 0;

    virtual optimizer *get_lambda_optimizer(root_equilibrium_distribution* p_distribution) = 0;
    void print_node_depths(std::ostream& ost);

    std::size_t get_gene_family_count() const {
        return _p_gene_families->size();
    }
};

std::vector<model *> build_models(const input_parameters& my_input_parameters, 
    clade *p_tree, 
    lambda *p_lambda, 
    std::vector<gene_family>* p_gene_families, 
    int max_family_size, 
    int max_root_family_size,
    error_model *p_error_model);
#endif /* CORE_H */

