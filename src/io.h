#ifndef io_h
#define io_h

#include "utils.h"
#include "family_generator.h"
#include "clade.h"
#include <string>
#include <map>
#include <vector>
#include <map>
#include <iosfwd>

using namespace std;

extern struct option longopts[];

class clade;
class error_model;
class gene_family;

clade *read_tree(std::string tree_file_path, bool lambda_tree);

map<int, int> *read_rootdist(std::string famdist_file_path);

void read_gene_families(std::istream& input_file, clade *p_tree, std::vector<gene_family> *p_gene_families);

void read_error_model_file(std::istream& error_model_file, error_model *p_error_model);

/* START: Printing functions for simulation engine */
void print_simulation(std::vector<vector<trial *> >  &sim, std::ostream& ost);
/* END: Printing functions for simulation engine*/

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
    double fixed_alpha = 0.0;
    double poisson_lambda = 0.0;
    bool is_simulating = false;
    int nsims = 0;
    int n_gamma_cats = 1;
    bool do_log = false;
    bool use_uniform_eq_freq = true;

    //! Check calls
    void check_input();
};

/* START: Reading in error model file */
class error_model {
private:
  int _max_cnt;  
    
  std::vector<int> _deviations; //!< Deviations from the true gene family (e.g., -1 0 1)
  
  std::vector<std::vector<double> > _error_dists; //!< Each vector element will be a gene family size; the vector of doubles inside (e.g., 0.1 0.8 0.1) will be the probs of deviating from the true value
  
public:
  //! Set max family size for which deviations apply
  void set_max_cnt(int max_cnt);
  
  //! Set deviations
  void set_deviations(std::vector<std::string> deviations);
  
  //! Set deviation probability vector for a certain family size
  void set_probs(int fam_size, std::vector<double>);
    
  //! Get deviation probability vector for a certain family size
  std::vector<double> get_probs(int fam_size) const;

  size_t n_deviations() const {
      return _deviations.size();
  }

  size_t get_max_count() const {
      return _error_dists.size();
  }
  std::vector<double> get_epsilons() const;
  error_model* replace_epsilons(std::map<double,double> *new_epsilons);
  void update_single_epsilon(double new_epsilon);

};
#endif
