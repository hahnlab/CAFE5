#include <algorithm>
#include <set>
#include <stdexcept>

#include "gene_family.h"
#include "clade.h"

using namespace std;

// these functions are intended to work with maps (key, value pairs)
template <typename T, typename U>
bool max_key(const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
	return p1.first < p2.first;
}

template <typename T, typename U>
bool max_value(const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
	return p1.second < p2.second;
}

//! Find and set _max_family_size and _parsed_max_family_size for this family
/*!
CAFE had a not_root_max (which we use = _parsed_max_family_size; see below) and a root_max = MAX(30, rint(max*1.25));
*/
int gene_family::get_max_size() const {
    // Max family size can only be found if there is data inside the object in the first place
    int max_family_size = 0;
    if (!_species_size_map.empty()) {
        max_family_size = (*max_element(_species_size_map.begin(), _species_size_map.end(), max_value<string, int>)).second;
        //_parsed_max_family_size = _max_family_size + std::max(50, _max_family_size/5);
    }
    return max_family_size;
}


//! Mainly for debugging: In case one want to grab the gene count for a given species
int gene_family::get_species_size(std::string species) const {
    // First checks if species data has been entered (i.e., is key in map?)
    if (_species_size_map.find(species) == _species_size_map.end()) {
        throw std::runtime_error(species + " was not found in gene family " + _id);
    }

    return _species_size_map.at(species);
}

//! Return first element of pair
template <typename type1, typename type2>
type1 do_get_species(const std::pair<type1, type2> & p1) {
    return p1.first;
}

//! Return vector of species names
vector<std::string> gene_family::get_species() const {
    vector<std::string> species_names(_species_size_map.size());
    transform(_species_size_map.begin(), _species_size_map.end(), species_names.begin(), do_get_species<string, int>); // Transform performs an operation on all elements of a container (here, the operation is the template)

    return species_names;
}

/// returns true if the family exists at the root, according to their parsimony reconstruction.
bool gene_family::exists_at_root(const clade *p_tree) const
{
    set<const clade *> exists;
    auto registered = [&exists](const clade *c) {
        return exists.find(c) != exists.end();
    };
    auto existence = [&](const clade *pc) {
        if (pc->is_leaf())
        {
            int sz = get_species_size(pc->get_taxon_name());
            if (sz > 0)
            {
                auto p = pc;
                do
                {
                    exists.insert(p);
                    //cout << "Registering node " << p->get_taxon_name() << endl;
                    p = p->get_parent();
                } while (p && !registered(p));
            }
        }
    };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), existence);

    bool exists_at_all_children = true;
    auto does_child_exist = [&](const clade *child) { exists_at_all_children &= registered(child); };
    p_tree->apply_to_descendants(does_child_exist);

    return exists_at_all_children;
}

int gene_family::species_size_differential() const
{
    auto compare = [](const std::pair<string, int>& a, const std::pair<string, int>& b) { return a.second < b.second; };
    int max_species_size = max_element(_species_size_map.begin(), _species_size_map.end(), compare)->second;
    int min_species_size = min_element(_species_size_map.begin(), _species_size_map.end(), compare)->second;
    return max_species_size - min_species_size;
}

void gene_family::init_from_clademap(const clademap<int>& values)
{
    for (auto& it : values) {
        if (it.first->is_leaf())
        {
            set_species_size(it.first->get_taxon_name(), it.second);
        }
    }

}