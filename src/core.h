#ifndef CORE_H
#define CORE_H

#include <set>

#include "clade.h"
#include "probability.h"
#include "root_distribution.h"

class simulation_data;
class inference_process;
class gene_family_reconstructor;
class reconstruction;
class user_data;
class root_distribution;
class inference_optimizer_scorer;

struct family_info_stash {
    family_info_stash() : lambda_multiplier(0.0), category_likelihood(0.0), family_likelihood(0.0), 
        posterior_probability(0.0), significant(false) {}
    family_info_stash(std::string fam, double lam, double cat_lh, double fam_lh, double pp, bool signif) : 
        family_id(fam), lambda_multiplier(lam), category_likelihood(cat_lh), family_likelihood(fam_lh),
        posterior_probability(pp), significant(signif) {}
    std::string family_id;
    double lambda_multiplier;
    double category_likelihood;
    double family_likelihood;
    double posterior_probability;
    bool significant;
};

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r);

//! The result of a model reconstruction. Should be able to (a) print reconstructed states with all available information;
/// (b) print increases and decreases by family; and (c) print increases and decreases by clade.
class reconstruction {
public:
    typedef const std::vector<gene_family> familyvector;

    void print_node_change(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree);

    void print_node_counts(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree);

    void print_reconstructed_states(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree, double test_pvalue, std::map<std::string, clademap<double>>& branch_probabilities);

    void print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order, familyvector& gene_families);

    void print_increases_decreases_by_family(std::ostream& ost, const cladevector& order, familyvector& gene_families, const std::vector<double>& pvalues, double test_pvalue);
        
    void print_family_clade_table(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree,
        std::function<string(int family_index, const clade* c)> get_family_clade_value);

    void write_results(std::string model_identifier, std::string output_prefix, const clade* p_tree, familyvector& families, std::vector<double>& pvalues, double test_pvalue, std::map<std::string, clademap<double>>& branch_probabilities);

    virtual int reconstructed_size(const gene_family& family, const clade* clade) const = 0;
    virtual ~reconstruction()
    {
    }

private:
    virtual void print_additional_data(const cladevector& order, familyvector& gene_families, std::string output_prefix) {};

    virtual int get_difference_from_parent(const gene_family* gf, const clade* c) = 0;
    virtual std::string get_reconstructed_state(const gene_family& gf, const clade* node) = 0;
    virtual void write_nexus_extensions(std::ostream& ost) {};
    virtual int get_node_count(const gene_family& gf, const clade* c) = 0;

};

class event_monitor
{
    std::map<string, int> failure_count;
    int attempts = 0;
    int rejects = 0;
public:
    void summarize(std::ostream& ost) const;

    void Event_InferenceAttempt_Started();
    void Event_InferenceAttempt_InvalidValues() { rejects++; }
    void Event_InferenceAttempt_Saturation(std::string family) { failure_count[family]++; }
    void Event_InferenceAttempt_Complete(double final_likelihood);

    void Event_Reconstruction_Started(std::string model);
    void Event_Reconstruction_Complete();
};

/*! @brief Describes the actions that are taken when estimating or simulating data

    A Model represents a way to calculate or simulate values in the data.
*/
class model {
protected:
    std::ostream & _ost; 
    lambda *_p_lambda;
    const clade *_p_tree;
    const std::vector<gene_family>* _p_gene_families;
    int _max_family_size;
    int _max_root_family_size;
    error_model* _p_error_model;
    vector<vector<int> > _rootdist_bins; // holds the distribution for each lambda bin

    /// Used to track gene families with identical species counts
    std::vector<size_t> references;

    std::vector<family_info_stash> results;

    event_monitor _monitor;

    //! Create a lambda based on the lambda tree model the user passed.
    /// Called when the user has provided no lambda value and one must
    /// be estimated. If the p_lambda_tree is NULL, uses a single
    /// lambda; otherwise uses the number of unique lambdas in the provided
    /// tree
    void initialize_lambda(clade *p_lambda_tree);
public:
    model(lambda* p_lambda,
        const clade *p_tree,
        const std::vector<gene_family>* p_gene_families,
        int max_family_size,
        int max_root_family_size,
        error_model *p_error_model);
    
    virtual ~model() {}
    
    /// Allows the replacement of the current set of families with a new set
    void set_families(const std::vector<gene_family>* p_gene_families)
    {
        _p_gene_families = p_gene_families;
    }

    lambda * get_lambda() const {
        return _p_lambda;
    }

    //! Returns a lambda suitable for creating a simulated family. Default case is simply to return the lambda provided by the user.
    virtual lambda* get_simulation_lambda();

    virtual void prepare_matrices_for_simulation(matrix_cache& cache) = 0;

    virtual double infer_family_likelihoods(root_equilibrium_distribution *prior, const std::map<int, int>& root_distribution_map, const lambda *p_lambda) = 0;  // return vector of likelihoods
    
    virtual std::string name() = 0;
    virtual void write_family_likelihoods(std::ostream& ost) = 0;
    virtual void write_vital_statistics(std::ostream& ost, double final_likelihood);
    void write_error_model(std::ostream& ost);

    //! Based on the model parameters, attempts to reconstruct the most likely counts of each family at each node
    virtual reconstruction* reconstruct_ancestral_states(const vector<const gene_family*>& families, matrix_cache *p_calc, root_equilibrium_distribution* p_prior) = 0;

    virtual inference_optimizer_scorer *get_lambda_optimizer(user_data& data) = 0;

    std::size_t get_gene_family_count() const;

    //! Tells the model to modify its lambdas slightly to provide a bit of extra randomness when simulating.
    //  Default is to do nothing.
    virtual void perturb_lambda() {}

    const event_monitor& get_monitor() { return _monitor;  }

    //! Returns true if the the model is interested in calculating pvalues for this family
    //! Default is to calculate pvalues for all families
    virtual bool should_calculate_pvalue(const gene_family& gf) const { return true; }

};

//! @brief Creates a list of families that are identical in all values
//!
//! With this information we can reduce the number of calculations required
//! and speed up the overall performance
std::vector<size_t> build_reference_list(const std::vector<gene_family>& families);

std::vector<model *> build_models(const input_parameters& my_input_parameters, user_data& user_data);

inline std::string filename(std::string base, std::string suffix, std::string extension)
{
    return (suffix.empty() ? std::string("results") : suffix) + "/" + base + "." + extension;
}

inline std::string filename(std::string base, std::string suffix)
{
    return filename(base, suffix, "txt");
}

std::vector<double> inference_prune(const gene_family& gf, matrix_cache& calc, const lambda *_lambda, const error_model *p_error_model, const clade *_p_tree, double _lambda_multiplier, int _max_root_family_size, int _max_family_size);

#endif /* CORE_H */

