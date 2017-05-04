#ifndef utils_h
#define utils_h

#include <regex>
#include <string>
#include <iostream>
#include <fstream>

class clade;

class gene_family
{
  std::string _id;
  std::map<std::string, int> species_size_map;
public:
  gene_family(std::string id) : _id(id) {}

  int get_species_size(std::string species) const
  {
    if (species_size_map.find(species) == species_size_map.end()) // checks if the key is in the map
    {
      throw std::runtime_error(species + " was not found in gene family " + _id);
    }
    return species_size_map.at(species);
  }
  void set_species_size(std::string species, int sz)
  {
    species_size_map[species] = sz;
  }

  std::vector<std::string> get_species() const;

  int max_family_size() const;

  std::string id() const { return _id; }
};

class newick_parser {

 private:
  int lp_count, rp_count;
  std::regex tokenizer;

  /* methods */
  clade *new_clade(clade *p_parent);
  
 public:
  std::string newick_string;

  /* methods */
  newick_parser(): tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;"), lp_count(0), rp_count(0) {} // constructor
  clade *parse_newick();
};

// these functions are intended to work with maps (key, value pairs)
template <typename T, typename U>
bool max_key(const std::pair<T, U> & p1, const std::pair<T, U> & p2)
{
  return p1.first < p2.first;
}

template <typename T, typename U>
bool max_value(const std::pair<T, U> & p1, const std::pair<T, U> & p2)
{
  return p1.second < p2.second;
}

#endif
