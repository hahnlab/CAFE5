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
    int nsims = 0;
    int n_gamma_cats = 1;
//    bool simulate = false;
    bool do_log = false;
    bool use_uniform_eq_freq = true;

    //! Check calls
    void check_input();
};

/* START: Reading in gene family data */
class gene_family {
private:
  std::string _id; //!< Gene family ID
  std::string _desc; //!< Gene family description
  int _max_family_size; //!< Gene family max size (largest observed gene count)
  int _parsed_max_family_size; //!< Gene family max as in CAFE (used for setting matrices dimensions)
  std::map<std::string, int> _species_size_map; //!< Map that stores each species gene family count: {sp1_name:count1, ...}

public:
  gene_family() { find_max_size(); } //!< Constructor for read user-provided input (sets _max_family_size and _parsed_max_family_size for this family)
  
  gene_family(trial *a_trial); //!< Constructor for simulations (trial is a typedef for a map = {clade *: int}
  
  void set_desc(std::string desc) { _desc = desc; }
  
  void set_id(std::string id) { _id = id; }
  
  void set_species_size(std::string species, int gene_count) {
      _species_size_map[species] = gene_count;
  }

  std::vector<std::string> get_species() const;

  //! Find and set _max_family_size and _parsed_max_family_size for this family
  void find_max_size();
  
  //! Getting max gene family size
  int get_max_size() { return _max_family_size; }

  //! Getting parsed max gene family size
  int get_parsed_max_size() { return _parsed_max_family_size; }

  std::string id() const { return _id; }
  
  //! Mainly for debugging: In case one want to grab the gene count for a given species
  int get_species_size(std::string species) const;
};
/* END: Reading in gene family data */

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
};
#endif
