#ifndef GENE_FAMILY_H
#define GENE_FAMILY_H

#include <string>
#include <map>
#include <vector>

class gene_family {
private:
    std::string _id; //!< Gene family ID
    std::string _desc; //!< Gene family description
    std::map<std::string, int> _species_size_map; //!< Map that stores each species gene family count: {sp1_name:count1, ...}

public:
    gene_family() { }

    void set_desc(std::string desc) { _desc = desc; }

    void set_id(std::string id) { _id = id; }

    void set_species_size(std::string species, int gene_count) {
        _species_size_map[species] = gene_count;
    }

    std::vector<std::string> get_species() const;

    int get_max_size() const;

    std::string id() const { return _id; }

    //! Mainly for debugging: In case one want to grab the gene count for a given species
    int get_species_size(std::string species) const;

    //! Returns true if every species size for both gene families are identical
    std::map<std::string, int> get_species_map() const {
        return _species_size_map;
    }

    bool species_size_match(const gene_family& other) const
    {
        return _species_size_map == other._species_size_map;
    }
};
#endif
