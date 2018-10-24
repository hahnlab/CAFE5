#include "core.h"

class inference_process_factory;
class gamma_bundle;

class gamma_model : public model {
private:
    //! Gamma
    std::vector<double> _lambda_multipliers;

    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    category_selector _gamma_cats;

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

    virtual optimizer_scorer *get_lambda_optimizer(user_data& data);

    virtual std::string name() {
        return "Gamma";
    }

    virtual void write_family_likelihoods(std::ostream& ost);
    virtual void write_vital_statistics(std::ostream& ost, double final_likelihood);

    virtual reconstruction* reconstruct_ancestral_states(matrix_cache *, root_equilibrium_distribution* p_prior);
    void print_reconstructed_states(std::ostream& ost);
    void print_increases_decreases_by_family(std::ostream& ost, const std::vector<double>& pvalues);
    void print_increases_decreases_by_clade(std::ostream& ost);

    std::size_t get_gamma_cat_probs_count() const {
        return _gamma_cat_probs.size();
    }

    std::size_t get_lambda_multiplier_count() const {
        return _lambda_multipliers.size();
    }

    void prepare_matrices_for_simulation(matrix_cache& cache);

    void initialize_simulations(size_t count);

};

class gamma_model_reconstruction : public reconstruction
{
    const std::vector<double>& _lambda_multipliers;
    const vector<gamma_bundle *>& _family_bundles;
    void print_increases_decreases_by_family(std::ostream& ost, const std::vector<double>& pvalues);
    void print_increases_decreases_by_clade(std::ostream& ost);
public:
    gamma_model_reconstruction(const std::vector<double>& lambda_multipliers, const vector<gamma_bundle *>& family_bundles) :
        _lambda_multipliers(lambda_multipliers), _family_bundles(family_bundles)
    {

    }

    void print_reconstructed_states(std::ostream& ost);
};

