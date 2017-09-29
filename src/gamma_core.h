#include "core.h"

class gamma_core : public core {
private:
    //! Gamma
    vector<double> _lambda_multipliers;

    vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families

    double _alpha;

    vector<gamma_bundle> _inference_bundles; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
                                             //! Basic constructor

public:
    gamma_core() { _alpha = 0; };

    //! Simulation constructor
    gamma_core(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, int n_gamma_cats, double alpha);

    //! Simulation constructor
    gamma_core(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, vector<int>& cats, vector<double>&mul);

    ~gamma_core();
    //! Gamma methods
    void adjust_n_gamma_cats(int n_gamma_cats);

    void adjust_family_gamma_membership(int n_families);

    //! Setters
    void set_alpha(double alpha, int n_families);
    void initialize_with_alpha(int n_gamma_cats, int n_families, double alpha);
    void initialize_without_alpha(int n_gamma_cats, int n_families, vector<double> lambda_multipliers, std::vector<int> gamma_cats);

    void set_lambda_multipliers(std::vector<double> lambda_multipliers);

    void set_gamma_cats(std::vector<int> gamma_cats);

    //! Simulation methods
    virtual simulation_process* create_simulation_process(int family_number);

    //! Inference methods
    void start_inference_processes();

    double infer_processes();

    virtual std::string name() {
        return "Gamma";
    }

    double get_lambda_multiplier(int family_id);

    virtual void print_results(std::ostream& ost);
};
