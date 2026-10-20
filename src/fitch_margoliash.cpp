#include <numeric>
#include <limits>

#include "doctest.h"
#include "easylogging++.h"

#include "clade.h"
#include "fitch_margoliash.h"

using namespace std;

fitch_margoliash::fitch_margoliash(distance_matrix distance_matrix, const std::vector<std::string>& species_names)
    : _distance_matrix(distance_matrix) {
    if (species_names.size() < 2)
        throw std::runtime_error("Fitch-Margoliash requires at least two species names");

    for (size_t i = 0; i < species_names.size(); ++i) {
        _clusters[species_names[i]] = new clade(species_names[i], -1);
    }
}   

double fitch_margoliash::distance(string sp1, string sp2) const
{
    auto it = _distance_matrix.find({sp1, sp2});
    if (it != _distance_matrix.end())
        return it->second;

    it = _distance_matrix.find({sp2, sp1});
    if (it != _distance_matrix.end())
        return it->second;

    throw std::runtime_error("Distance not found for " + sp1 + " and " + sp2);
}

set<string> get_tips(const clade *c) {
    if (c->is_leaf()) {
        return {c->get_taxon_name()}; // if it's a leaf, return its name
    }

    set<string> tips;
    // If the clade is not a leaf, we need to collect all tips in its descendants
    for (auto it = c->reverse_level_begin(); it != c->reverse_level_end(); ++it) {
        if ((*it)->is_leaf()) {
            tips.insert((*it)->get_taxon_name());
        }
    }
    return tips;
}

// Calculate the distance from a species to all other clusters
double fitch_margoliash::average_distance(string s) const {
    vector<double> distances;
    for (auto it = _clusters.begin(); it != _clusters.end(); ++it) {
        if (it->second->find_descendant(s)) 
            continue; // Skip if the species is a descendant of the cluster

        for(auto tip : get_tips(it->second)) {
            double d = distance(s, tip);
            distances.push_back(d);
        }
    }

    return accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
}

pair<string, string> fitch_margoliash::closest_clusters()
{
    double min_dist = std::numeric_limits<double>::max();
    pair<string, string> min_pair; // pair of species names and their distance
    for (auto& i : _clusters) {
        for (auto& j : _clusters) {
            if (i.first == j.first) continue; // Skip self-comparison
            double dist = distance(i.first, j.first);
            if (dist < min_dist) {
                min_dist = dist;
                min_pair = make_pair(i.first, j.first);
            }
        }
    }
    if (min_pair.first.empty() || min_pair.second.empty()) {
        throw std::runtime_error("No valid clusters found for merging");
    }
    return min_pair;
}

// Find the closest pair of clusters (A,B)
// create a new node to be their parent (C)
// compute branch lengths (AC, BC)
// compute distances from C to all other nodes
// Replace A and B in clusters with C
clade* fitch_margoliash::build_tree() 
{
    while (_clusters.size() > 1) {

        auto to_merge = closest_clusters();
        string sp1 = to_merge.first;
        string sp2 = to_merge.second;

        clade * new_cluster = new clade(sp1 + "_" + sp2, -1);
        _clusters[sp1]->_p_parent = new_cluster;
        _clusters[sp2]->_p_parent = new_cluster;
        new_cluster->add_descendant(_clusters[sp1]);
        new_cluster->add_descendant(_clusters[sp2]);

        compute_branch_lengths(sp1, sp2, new_cluster);

        // add distances to the new cluster
        for (auto cluster : _clusters) {
            if (cluster.first == sp1 || cluster.first == sp2)
                continue; // Skip the merged clusters

            double d = (distance(sp1, cluster.first) + distance(sp2, cluster.first)) / 2.0;
            _distance_matrix[{new_cluster->get_taxon_name(), cluster.first}] = d;
            auto c = cluster.second;
            for (auto it = c->reverse_level_begin(); it != c->reverse_level_end(); ++it) {
                string node_name = (*it)->get_taxon_name();
                double d = (distance(sp1, node_name) + distance(sp2, node_name)) / 2.0;
                _distance_matrix[{new_cluster->get_taxon_name(), node_name}] = d;
            }
        }

        // Update clusters
        _clusters.erase(sp1);
        _clusters.erase(sp2);
        _clusters[new_cluster->get_taxon_name()] = new_cluster;

    }

    return _clusters.begin()->second; // Return the root of the tree
}

void fitch_margoliash::compute_branch_lengths(std::string sp1, std::string sp2, clade *new_cluster)
{
    double dac = distance(sp1, sp2);
    double daw = average_distance(sp1);
    double dcw = average_distance(sp2);

    double b1 = (dac + daw - dcw) / 2.0;
    double b2 = (dac + dcw - daw) / 2.0;

    for (auto c = new_cluster->_descendants.begin(); c != new_cluster->_descendants.end(); ++c)
    {
        if ((*c)->get_taxon_name() == sp1)
        {
            (*c)->_branch_length = b1;
        }
        else if ((*c)->get_taxon_name() == sp2)
        {
            (*c)->_branch_length = b2;
        }
    }
}

TEST_CASE("fitch_margoliash")
{
    // https://www.cs.rice.edu/~nakhleh/COMP571/Slides/Phylogenetics-DistanceMethods-Full.pdf
    fitch_margoliash::distance_matrix distance_matrix = {
        {{"A", "B"}, 5.0},
        {{"A", "C"}, 4.0},
        {{"A", "D"}, 9.0},
        {{"A", "E"}, 8.0},
        {{"B", "C"}, 5.0},
        {{"B", "D"}, 10.0},
        {{"B", "E"}, 9.0},
        {{"C", "D"}, 7.0},
        {{"C", "E"}, 6.0},
        {{"D", "E"}, 7.0},
    };

    std::vector<std::string> species_names = {"A", "B", "C", "D", "E"};
    fitch_margoliash f(distance_matrix, species_names);
    auto tree = f.build_tree();

    ostringstream oss;
    tree->write_newick(oss, [](ostream& ost, const clade* c) {
        ost << c->get_taxon_name() << ":" << c->get_branch_length();
    });
    CHECK_EQ("(((A:2.5,C:1.5)AC:2.20833,B:2.79167)ABC:4.41667,(D:3.875,E:3.125)DE:4.08333)ABCDE:-1", oss.str());
}

TEST_CASE("compute_branch_lengths")
{
    fitch_margoliash::distance_matrix distance_matrix = {
        {{"A", "B"}, 5.0},
        {{"A", "C"}, 4.0},
        {{"B", "C"}, 3.0},
    };

    std::vector<std::string> species_names = {"A", "B", "C"};
    fitch_margoliash f(distance_matrix, species_names);
    clade *new_cluster = new clade("AB", -1);
    new_cluster->add_descendant(new clade("A", 0));
    new_cluster->add_descendant(new clade("B", 0));

    f.compute_branch_lengths("A", "B", new_cluster);

    CHECK_EQ(2.75, new_cluster->find_descendant("A")->get_branch_length());
    CHECK_EQ(2.25, new_cluster->find_descendant("B")->get_branch_length());
}

TEST_CASE("closest_clusters")
{
    fitch_margoliash::distance_matrix distance_matrix = {
        {{"A", "B"}, 5.0},
        {{"A", "C"}, 4.0},
        {{"B", "C"}, 3.0},
    };

    std::vector<std::string> species_names = {"A", "B", "C"};
    fitch_margoliash f(distance_matrix, species_names);
    
    auto closest = f.closest_clusters();
    CHECK_EQ("B", closest.first);
    CHECK_EQ("C", closest.second);
}

TEST_CASE("average_distance")
{
    fitch_margoliash::distance_matrix distance_matrix = {
        {{"A", "B"}, 5.0},
        {{"A", "C"}, 4.0},
        {{"B", "C"}, 3.0},
    };

    std::vector<std::string> species_names = {"A", "B", "C"};
    fitch_margoliash f(distance_matrix, species_names);
    
    double avg_dist = f.average_distance("A");
    CHECK_EQ(4.5, avg_dist); // Average distance from A to B and C
}

TEST_CASE("get_tips")
{
    string nwk = "((A:1,B:2)C:3,(D:4,E:5)F:6)G;"; // A tree with 5 tips
    unique_ptr<clade> p_tree(parse_newick(nwk));
    
    auto tips = get_tips(p_tree.get());
    CHECK_EQ(4, tips.size());
    CHECK(tips.count("A"));
    CHECK(tips.count("B"));
    CHECK(tips.count("D"));
    CHECK(tips.count("E"));
}

TEST_CASE("distance")
{
    fitch_margoliash::distance_matrix distance_matrix = {
        {{"A", "B"}, 5.0},
        {{"A", "C"}, 4.0},
        {{"B", "C"}, 3.0},
    };

    std::vector<std::string> species_names = {"A", "B", "C"};
    fitch_margoliash f(distance_matrix, species_names);
    
    CHECK_EQ(5.0, f.distance("A", "B")); // Direct distance from A to B
    CHECK_EQ(5.0, f.distance("B", "A")); // Direct distance from A to B
    CHECK_EQ(3.0, f.distance("B", "C")); // Direct distance from A to B
}
