#ifndef io_h
#define io_h

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

map<int, int> *read_rootdist(string famdist_file_path);

/* START: Printing functions for simulation engine */
void print_simulation(std::vector<vector<trial *> >  &sim, std::ostream& ost);
/* END: Printing functions for simulation engine*/

/* START: Reading in gene family data */

//! This class holds data for ONE gene family
class gene_family {
private:
  std::string _id; //!< Species name
  std::map<std::string, int> species_size_map; //!< Map that stores each species gene family count: {sp1_name:count1, ...}

public:
  gene_family(std::string id) : _id(id) {} //!< Constructor

  //! Mainly for debugging: In case one want to grab the gene count for a given species
  int get_species_size(std::string species) const {
      // First checks if species data has been entered (i.e., is key in map?)
      if (species_size_map.find(species) == species_size_map.end()) {
          throw std::runtime_error(species + " was not found in gene family " + _id);
      }
      
      return species_size_map.at(species);
  }
  
  void set_species_size(std::string species, int gene_count) {
      species_size_map[species] = gene_count;
  }

  std::vector<std::string> get_species() const;

  int max_family_size() const;

  std::string id() const { return _id; }
};

/* END: Reading in gene family data */
#endif
