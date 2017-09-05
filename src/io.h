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

clade *read_tree(std::string tree_file_path, bool lambda_tree);

map<int, int> *read_rootdist(std::string famdist_file_path);

std::vector<gene_family> * read_gene_families(std::string input_file_path);

/* START: Printing functions for simulation engine */
void print_simulation(std::vector<vector<trial *> >  &sim, std::ostream& ost);
/* END: Printing functions for simulation engine*/

/* START: Reading in gene family data */

struct input_parameters
{
	std::string input_file_path;
	std::string output_prefix;
	bool estimate = false;
	std::string tree_file_path;
	double fixed_lambda = 0.0;
	bool simulate = false;
	std::string fixed_multiple_lambdas;
	std::string rootdist;
	int nsims = 0;
	bool do_log = false;
	std::string lambda_tree_file_path;
};

class gene_family {
private:
  std::string _id; //!< Gene family ID
  std::string _desc; //!< Gene family description
  int _max_family_size; //!< Gene family max size
  int _parsed_max_family_size; //!< Gene family max as in CAFE (used for setting matrices dimensions)
  std::map<std::string, int> _species_size_map; //!< Map that stores each species gene family count: {sp1_name:count1, ...}

public:
  gene_family() { find_max_size(); } //!< Constructor
  gene_family(trial *a_trial);
  
  void set_desc(std::string desc) { _desc = desc; }
  
  void set_id(std::string id) { _id = id; }
  
  void set_species_size(std::string species, int gene_count) {
      _species_size_map[species] = gene_count;
  }

  std::vector<std::string> get_species() const;

  //! Find max gene family size
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
#endif
