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

typedef clademap<int> trial;

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

class branch_length_finder
{
    std::set<double> _result;
public:
    void operator()(const clade *c);

    std::set<double> result() const
    {
        return _result;
    }

    double longest() const;
};

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r);

class reconstruction {
    virtual void print_reconstructed_states(std::ostream& ost) = 0;
    virtual void print_increases_decreases_by_family(std::ostream& ost, const std::vector<double>& pvalues) = 0;
    virtual void print_increases_decreases_by_clade(std::ostream& ost) = 0;

public:
    void write_results(model *p_model, std::string output_prefix, std::vector<double>& pvalues);
    virtual ~reconstruction()
    {
    }
};


class model {
protected:
    std::ostream & _ost; 
    lambda *_p_lambda; // TODO: multiple lambdas for different branches
    const clade *_p_tree;
    int _max_family_size;
    int _max_root_family_size;
    const std::vector<gene_family>* _p_gene_families;
    root_distribution _root_distribution; // in case the user wants to use a specific root size distribution for all simulations
    vector<vector<int> > _rootdist_bins; // holds the distribution for each lambda bin

    /// Used to track gene families with identical species counts
    std::vector<int> references;

    void initialize_rootdist_if_necessary();

    std::vector<family_info_stash> results;

    error_model* _p_error_model;
public:
    model(lambda* p_lambda,
        const clade *p_tree,
        const std::vector<gene_family>* p_gene_families,
        int max_family_size,
        int max_root_family_size,
        error_model *p_error_model);
    
    virtual ~model() {}
    
    void set_families(const std::vector<gene_family>* p_gene_families)
    {
        _p_gene_families = p_gene_families;
    }

    lambda * get_lambda() const {
        return _p_lambda;
    }
    void initialize_lambda(clade *p_lambda_tree);

    void set_max_sizes(int max_family_size, int max_root_family_size);
    
    //! Simulation methods
    virtual lambda* get_simulation_lambda(const user_data& data);

    virtual void prepare_matrices_for_simulation(matrix_cache& cache) = 0;

    //! Inference methods
    virtual void start_inference_processes(lambda *) = 0;
    
    virtual double infer_processes(root_equilibrium_distribution *prior) = 0;  // return vector of likelihoods
    
    virtual std::string name() = 0;
    virtual void write_family_likelihoods(std::ostream& ost) = 0;
    virtual void write_vital_statistics(std::ostream& ost, double final_likelihood);

    virtual reconstruction* reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior) = 0;

    virtual optimizer_scorer *get_lambda_optimizer(user_data& data) = 0;
    void print_node_depths(std::ostream& ost);

    std::size_t get_gene_family_count() const;

    int get_max_simulation_size() const;

    std::size_t get_rootdist_size() const {
        return _root_distribution.size();
    }
};

std::vector<int> build_reference_list(const std::vector<gene_family>& families);

enum family_size_change { Increase, Decrease, Constant };
std::ostream& operator<<(std::ostream& ost, family_size_change fsc);

struct increase_decrease
{
    std::string gene_family_id;
    double pvalue = 0.0;
    std::vector<family_size_change> change;
    std::vector<double> category_likelihoods;
};


std::vector<model *> build_models(const input_parameters& my_input_parameters, user_data& user_data);

template <class T>
void print_increases_decreases_by_family(std::ostream& ost, const std::vector<T>& printables, const std::vector<double>& pvalues)
{
    if (printables.size() != pvalues.size())
    {
        throw std::runtime_error("No pvalues found for family");
    }
    if (printables.empty())
    {
        ost << "No increases or decreases recorded\n";
        return;
    }
    auto rec = printables[0];
    auto order = rec->get_taxa();

    ost << "#FamilyID\tpvalue\t*\t";
    for (auto& it : order) {
        ost << it->get_taxon_name() << "\t";
    }
    ost << endl;

    for (size_t i = 0; i < printables.size(); ++i) {
        ost << printables[i]->get_increases_decreases(order, pvalues[i]);
    }
}

template<class T>
void print_increases_decreases_by_clade(std::ostream& ost, const std::vector<T>& printables) {
    if (printables.empty())
    {
        ost << "No increases or decreases recorded\n";
        return;
    }

    auto rec = printables[0];
    auto order = rec->get_taxa();

    clademap<pair<int, int>> increase_decrease_map;

    for (auto item : printables) {
        auto incdec = item->get_increases_decreases(order, 0.0);
        for (int i = 0; i < order.size(); ++i)
        {
            if (incdec.change[i] == Increase)
                increase_decrease_map[order[i]].first++;
            if (incdec.change[i] == Decrease)
                increase_decrease_map[order[i]].second++;
        }
    }

    ost << "#Taxon_ID\tIncrease/Decrease\n";
    for (auto& it : increase_decrease_map) {
        ost << it.first->get_taxon_name() << "\t";
        ost << it.second.first << "/" << it.second.second << endl;
    }
}

inline std::string filename(std::string base, std::string suffix)
{
    return base + (suffix.empty() ? "" : "_") + suffix + ".txt";
}


#endif /* CORE_H */

