#ifndef io_h
#define io_h

#include "optimizer.h"

using namespace std;

extern struct option longopts[];

class clade;
class error_model;
class gene_family;

clade *read_tree(std::string tree_file_path, bool lambda_tree);

void read_gene_families(std::istream& input_file, clade *p_tree, std::vector<gene_family> *p_gene_families);

void read_error_model_file(std::istream& error_model_file, error_model *p_error_model);
void write_error_model_file(std::ostream& ost, error_model& errormodel);

struct input_parameters {
public:	
    std::string input_file_path;
    std::string error_model_file_path;
    std::string output_prefix;
    std::string tree_file_path;
    std::string lambda_tree_file_path;
    std::string fixed_multiple_lambdas;
    std::string chisquare_compare;
    std::string rootdist;
    double fixed_lambda = 0.0;
    double fixed_alpha = -1.0;
    double poisson_lambda = 0.0;
    double pvalue = 0.05;
    bool is_simulating = false;
    int nsims = 0;
    int n_gamma_cats = 1;
    bool do_log = false;
    bool use_uniform_eq_freq = true;
    bool exclude_zero_root_families = true;
    bool lambda_per_family = false;
    bool use_error_model = false;

    optimizer_parameters optimizer_params;
    bool help = false;

    //! Check calls
    void check_input();
};

std::vector<std::string> tokenize_str(std::string some_string, char some_delim);

#endif
