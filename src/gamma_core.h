#include "core.h"
#include "gamma_bundle.h"

class inference_process_factory;

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

    //! Computation or estimation constructor
    gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
        int max_root_family_size, int n_gamma_cats, double fixed_alpha, std::map<int, int> *p_rootdist_map);

    ~gamma_model();
    //! Gamma methods

    //! Setters
    void set_alpha(double alpha, int n_families);

    void write_probabilities(std::ostream& ost);

    //! Simulation methods
    virtual simulation_process* create_simulation_process(int family_number);

    //! Inference methods
    void start_inference_processes();

    double infer_processes(probability_calculator& calc, root_equilibrium_distribution *prior);

    std::vector<double> initial_guesses();

    void set_current_guesses(double * guesses);

    virtual std::string name() {
        return "Gamma";
    }

    virtual void print_results(std::ostream& ost);

    virtual void reconstruct_ancestral_states(probability_calculator *, root_equilibrium_distribution* p_prior);
    void print_reconstructed_states(std::ostream& ost);

};
