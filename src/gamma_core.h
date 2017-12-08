#include "core.h"

class inference_process_factory;

class gamma_bundle {
    std::vector<inference_process *> _inf_processes;
    std::vector<reconstruction_process *> _rec_processes;
public:
    gamma_bundle(inference_process_factory& factory, std::vector<double> lambda_multipliers);

    void clear();

    std::vector<double> prune(const vector<double>& gamma_cat_probs, root_equilibrium_distribution *eq_freq);

    void reconstruct(const vector<double>& _gamma_cat_probs);

    double get_lambda_likelihood(int family_id);

    void set_values(probability_calculator *, root_equilibrium_distribution*);

    void print_reconstruction(std::ostream& ost, std::vector<clade *> order);

    std::vector<clade *> get_taxa();
};


class gamma_model : public model {
private:
    //! Gamma
    std::vector<double> _lambda_multipliers;

    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    std::vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families

    double _alpha;

    vector<gamma_bundle> _family_bundles; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
                                             //! Basic constructor

    std::vector<double> get_posterior_probabilities(std::vector<double> cat_likelihoods);
public:
    
    //! Basic constructor
    gamma_model() { _alpha = 0; };

    //! Computation or estimation constructor
    gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
        int max_root_family_size, int n_gamma_cats, double fixed_alpha, std::map<int, int> *p_rootdist_map);
    
    //! Simulation constructors
    gamma_model(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, int n_gamma_cats, double alpha);
    
    gamma_model(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, vector<int>& cats, vector<double>&mul);

    ~gamma_model();
    //! Gamma methods
    void adjust_n_gamma_cats(int n_gamma_cats);

    void adjust_family_gamma_membership(int n_families);

    //! Setters
    void set_alpha(double alpha, int n_families);
    void initialize_with_alpha(int n_gamma_cats, int n_families, double alpha);
    void initialize_without_alpha(int n_gamma_cats, int n_families, vector<double> lambda_multipliers, std::vector<int> gamma_cats);

    void set_lambda_multipliers(std::vector<double> lambda_multipliers);

    //! Simulation methods
    virtual simulation_process* create_simulation_process(int family_number);

    //! Inference methods
    void start_inference_processes();

    double infer_processes(root_equilibrium_distribution *prior);

    std::vector<double> initial_guesses();

    void set_current_guesses(double * guesses);

    virtual std::string name() {
        return "Gamma";
    }

    virtual void print_results(std::ostream& ost);

    virtual void reconstruct_ancestral_states(probability_calculator *, root_equilibrium_distribution* p_prior);
    void print_reconstructed_states(std::ostream& ost);

};
