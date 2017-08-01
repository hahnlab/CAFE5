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

//! This class holds data for ONE gene family
class gene_family {
private:
  std::string _id; //!< Gene family ID
  std::string _desc; //!< Gene family description
  std::map<std::string, int> _species_size_map; //!< Map that stores each species gene family count: {sp1_name:count1, ...}

public:
  gene_family() {} //!< Constructor
  
  void set_desc(std::string desc) {
      _desc = desc;
  }
  
  void set_id(std::string id) {
      _id = id;
  }
  
  void set_species_size(std::string species, int gene_count) {
      _species_size_map[species] = gene_count;
  }

  std::vector<std::string> get_species() const;

//  //! Getting max gene family size
//  int get_max_size()
  
  //! Mainly for debugging: In case one want to grab the gene count for a given species
  int get_species_size(std::string species) const {
      // First checks if species data has been entered (i.e., is key in map?)
      if (_species_size_map.find(species) == _species_size_map.end()) {
          throw std::runtime_error(species + " was not found in gene family " + _id);
      }
      
      return _species_size_map.at(species);
  }
  
  int max_family_size() const;

  std::string id() const { return _id; }
};

/* END: Reading in gene family data */
#endif
