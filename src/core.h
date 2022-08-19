#ifndef CORE_H
#define CORE_H

#include <set>

#include "easylogging++.h"

#include "clade.h"
#include "probability.h"

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

class branch_probabilities {
public:
    struct branch_probability {
        bool _is_valid;
        double _value;

        branch_probability(double value) : _is_valid(true), _value(value)
        {
            if (value < 0 || value > 1)
                throw std::runtime_error("Not a valid probability");
        }
        branch_probability() : _is_valid(false), _value(0.0)
        {

        }
    };


    bool contains(const gene_family& fam) const;
    branch_probability at(const gene_family& fam, const clade* c) const;
    void set(const gene_family& fam, const clade* c, branch_probability p);

    static branch_probability invalid() { return branch_probability(); }

private:
    std::map<std::string, clademap<branch_probability>> _probabilities;
};

//using probabilitymap = std::map<std::string, clademap<branch_probability>>;


//! The result of a model reconstruction. Should be able to (a) print reconstructed states with all available information;
/// (b) print increases and decreases by family; and (c) print increases and decreases by clade.
class reconstruction {
protected:
    const user_data& _data;
    const input_parameters& _user_input;

public:
    reconstruction(const user_data& d, const input_parameters& ui) : _data(d), _user_input(ui)
    {

    }

    void print_node_change(std::ostream& ost, const cladevector& order);

    void print_node_counts(std::ostream& ost, const cladevector& order);

    void print_reconstructed_states(std::ostream& ost, const cladevector& order, const branch_probabilities& branch_probabilities);

    void print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order);

    void print_increases_decreases_by_family(std::ostream& ost, const cladevector& order, const std::vector<double>& pvalues);
        
    void print_family_clade_table(std::ostream& ost, const cladevector& order,
        std::function<string(int family_index, const clade* c)> get_family_clade_value);

    void write_results(std::string model_identifier, std::vector<double>& pvalues, const branch_probabilities& branch_probabilities);

    virtual ~reconstruction()
    {
    }

    virtual int get_node_count(const gene_family& gf, const clade* c) const = 0;

    int get_difference_from_parent(const gene_family& gf, const clade* c);
private:
    virtual void print_additional_data(const cladevector& order) {};

    virtual void write_nexus_extensions(std::ostream& ost) {};

};

class event_monitor : public el::Loggable
{
    std::map<string, int> failure_count;
    int attempts = 0;
    int rejects = 0;
public:
    virtual void log(el::base::type::ostream_t& os) const;

    void Event_InferenceAttempt_Started();
    void Event_InferenceAttempt_InvalidValues() { rejects++; }
    void Event_InferenceAttempt_Saturation(std::string family) { failure_count[family]++; }
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

    virtual double infer_family_likelihoods(const root_equilibrium_distribution& prior, const lambda *p_lambda) = 0;  // return vector of likelihoods
    
    virtual std::string name() const = 0;
    virtual void write_family_likelihoods(std::ostream& ost) = 0;
    virtual void write_vital_statistics(std::ostream& ost, double final_likelihood);
    void write_error_model(std::ostream& ost) const;

    //! Based on the model parameters, attempts to reconstruct the most likely counts of each family at each node
    virtual reconstruction* reconstruct_ancestral_states(const user_data& data, const input_parameters& _user_input, matrix_cache *p_calc) = 0;

    virtual inference_optimizer_scorer *get_lambda_optimizer(const user_data& data) = 0;

    std::size_t get_gene_family_count() const;

    const event_monitor& get_monitor() { return _monitor;  }
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

void exclude_zero_root_families(const input_parameters& user_input, user_data& data);

#endif /* CORE_H */

