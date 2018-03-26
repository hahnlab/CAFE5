#include "core.h"

class inference_process_factory;
class gamma_bundle;

class gamma_model : public model {
private:
    //! Gamma
    std::vector<double> _lambda_multipliers;

    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    std::vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families

    double _alpha;

    vector<gamma_bundle *> _family_bundles; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
                                             //! Basic constructor

    std::vector<double> get_posterior_probabilities(std::vector<double> cat_likelihoods);
public:

    //! Computation or estimation constructor
    gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
        int max_root_family_size, int n_gamma_cats, double fixed_alpha, std::map<int, int> *p_rootdist_map, error_model *p_error_model);

    ~gamma_model();
    //! Gamma methods

    //! Setters
    void set_alpha(double alpha, int n_families);

    void write_probabilities(std::ostream& ost);

    //! Simulation methods
    virtual simulation_process* create_simulation_process(int family_number);

    //! Inference methods
    void start_inference_processes();

    double infer_processes(root_equilibrium_distribution *prior);

    virtual optimizer_scorer *get_lambda_optimizer(root_equilibrium_distribution* p_distribution);

    virtual std::string name() {
        return "Gamma";
    }

    virtual void write_family_likelihoods(std::ostream& ost);
    virtual void write_vital_statistics(std::ostream& ost, double final_likelihood);

    virtual void reconstruct_ancestral_states(matrix_cache *, root_equilibrium_distribution* p_prior);
    void print_reconstructed_states(std::ostream& ost);
    void print_increases_decreases_by_family(std::ostream& ost, const std::vector<double>& pvalues);
    void print_increases_decreases_by_clade(std::ostream& ost);

    std::size_t get_gamma_cat_probs_count() const {
        return _gamma_cat_probs.size();
    }

    std::size_t get_lambda_multiplier_count() const {
        return _lambda_multipliers.size();
    }
};

class gamma_lambda_optimizer : public optimizer_scorer
{
    gamma_model *_p_model;
    root_equilibrium_distribution *_p_distribution;
    clade *_p_tree;
    lambda *_p_lambda;
public:
    gamma_lambda_optimizer(clade *p_tree, lambda *p_lambda, gamma_model * p_model, root_equilibrium_distribution *p_distribution) :
        _p_tree(p_tree),
        _p_lambda(p_lambda),
        _p_model(p_model),
        _p_distribution(p_distribution)
    {

    }

    std::vector<double> initial_guesses();

    double calculate_score(double *values);

    /// results consists of the desired number of lambdas and one alpha value
    void finalize(double *results) {
        _p_lambda->update(results);
        double alpha = results[_p_lambda->count()];
        _p_model->set_alpha(alpha, _p_model->get_gene_family_count());
    }

};

