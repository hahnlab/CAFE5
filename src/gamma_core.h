#include "core.h"

//! \defgroup gamma Gamma Model
//! @brief Extends the Base model by assuming lambda values belong to a gamma distribution
class inference_process_factory;
class gamma_bundle;


//! @brief Holds data for reconstructing a tree based on the Gamma model
//! \ingroup gamma
class gamma_model_reconstruction : public reconstruction
{
    const std::vector<double> _lambda_multipliers;
    virtual void write_nexus_extensions(std::ostream& ost) override;

public:
    gamma_model_reconstruction(const std::vector<double>& lambda_multipliers) :
        _lambda_multipliers(lambda_multipliers)
    {
    }

    void print_additional_data(const cladevector& order, familyvector& gene_families, std::string output_prefix) override;

    void print_category_likelihoods(std::ostream& ost, const cladevector& order, familyvector& gene_families);

    int get_node_count(const gene_family& gf, const clade* c) const override;

    struct gamma_reconstruction {
        std::vector<clademap<int>> category_reconstruction;
        clademap<double> reconstruction;
        std::vector<double> _category_likelihoods;
    };

    std::map<string, gamma_reconstruction> _reconstructions;
};

//! @brief Represents a model of species change in which lambda values are expected to belong to a gamma distribution
//! \ingroup gamma
class gamma_model : public model {
private:
    //! Gamma
    std::vector<double> _lambda_multipliers;

    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    std::vector<std::vector<double>> _category_likelihoods;

    double _alpha;

    std::vector<double> get_posterior_probabilities(std::vector<double> cat_likelihoods);
public:

    //! Calculate gamma categories and lambda multipliers based on category count and a fixed alpha
    gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
        int max_root_family_size, int n_gamma_cats, double fixed_alpha, error_model *p_error_model);

    //! Specify gamma categories and lambda multipliers explicitly
    gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
        int max_root_family_size, std::vector<double> gamma_categories, std::vector<double> multipliers, error_model *p_error_model);

    void set_alpha(double alpha);
    double get_alpha() const { return _alpha; }

    void write_probabilities(std::ostream& ost);

    //! Randomly select one of the multipliers to apply to the simulation
    virtual lambda* get_simulation_lambda() override;

    double infer_family_likelihoods(const root_equilibrium_distribution& prior, const lambda *p_lambda) override;

    virtual inference_optimizer_scorer *get_lambda_optimizer(const user_data& data) override;

    virtual std::string name() const override {
        return "Gamma";
    }

    virtual void write_family_likelihoods(std::ostream& ost) override;
    virtual void write_vital_statistics(std::ostream& ost, double final_likelihood) override;

    virtual reconstruction* reconstruct_ancestral_states(const vector<gene_family>& families, matrix_cache *, root_equilibrium_distribution* p_prior) override;

    std::size_t get_gamma_cat_probs_count() const {
        return _gamma_cat_probs.size();
    }

    std::size_t get_lambda_multiplier_count() const {
        return _lambda_multipliers.size();
    }

    std::vector<double> get_lambda_multipliers() const {
        return _lambda_multipliers;
    }

    void prepare_matrices_for_simulation(matrix_cache& cache) override;

    bool can_infer() const;

    bool prune(const gene_family& family, const root_equilibrium_distribution& eq, matrix_cache& calc, const lambda *p_lambda,
        std::vector<double>& category_likelihoods);
};

//! \ingroup gamma
clademap<double> get_weighted_averages(const std::vector<clademap<int>>& m, const vector<double>& probabilities);

