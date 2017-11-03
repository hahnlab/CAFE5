#ifndef CORE_H
#define CORE_H

#include "clade.h"
#include "probability.h"

class clade; // core.cpp includes clade.h before including core.h

class process;
class inference_process;
class simulation_process;
class inference_process_factory;
class root_equilibrium_distribution;

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
    
    void initialize_rootdist_if_necessary();

    std::vector<family_info_stash> results;
public:
    //! Basic constructor
    model(): _ost(cout), _p_lambda(NULL), _p_tree(NULL), _p_gene_families(NULL), _total_n_families_sim(1) {}
    
    model(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec): _ost(ost), _p_lambda(p_lambda), _p_tree(p_tree), _max_family_size(max_family_size), _total_n_families_sim(total_n_families), _rootdist_vec(rootdist_vec) {}
    
    model(lambda* p_lambda, clade *p_tree, vector<gene_family> *p_gene_families, int max_family_size, int max_root_family_size): _ost(cout), _p_lambda(p_lambda), _p_tree(p_tree), _p_gene_families(p_gene_families), _max_family_size(max_family_size), _max_root_family_size(max_root_family_size) {}
    
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

    double initialize_lambda_guess();
};

class base_model : public model {
    std::vector<inference_process *> processes;
    virtual simulation_process* create_simulation_process(int family_number);
public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, clade *p_tree, vector<gene_family> *p_gene_families, int max_family_size, int max_root_family_size): 
    model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size) {}
    
    virtual void start_inference_processes();
    virtual double infer_processes(root_equilibrium_distribution *prior);

    virtual std::string name() {
        return "Base";
    }
    virtual ~base_model();

    virtual void print_results(std::ostream& ost);

};

#endif /* CORE_H */

