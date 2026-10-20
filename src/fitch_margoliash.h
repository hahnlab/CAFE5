#ifndef FITCH_MARGOLIASH_H
#define FITCH_MARGOLIASH_H

#include <map>
#include <string>
#include <vector>

class clade;

/// @brief The fitch_margoliash algorithm builds a tree based on the distance matrix
/// It uses a distance matrix to iteratively merge the closest clusters until only one cluster remains.
/// The resulting tree is built from the clusters, with branch lengths computed based on the average distances
/// between the merged clusters.
class fitch_margoliash {
    // The fitch_margoliash algorithm builds a tree based on the distance matrix
public:
    using distance_matrix = std::map<std::pair<std::string, std::string>, double>;

    /// @brief Constructor that initializes the fitch_margoliash algorithm with a distance matrix and species names.
    fitch_margoliash(distance_matrix distance_matrix, const std::vector<std::string>& species_names);

    /// @brief Returns the average distance from a given species to all other species in the distance matrix.
    /// This is used to compute branch lengths in the tree.
    /// @param s The species name for which to compute the average distance.
    /// @return average distance from the species to all other species.
    double average_distance(std::string s) const;

    /// @brief Builds the tree using the fitch_margoliash algorithm.
    /// @return A pointer to the root of the tree.
    clade *build_tree();

    /// @brief Computes the branch lengths for the new cluster formed by merging two species.
    /// The branch lengths are computed based on the average distances from the merged species to all other
    /// species in the distance matrix.
    void compute_branch_lengths(std::string sp1, std::string sp2, clade *new_cluster);

    /// @brief Finds the distance between two species based on the distance matrix.
    /// @param sp1 The first species name.
    /// @param sp2 The second species name.
    /// @return The distance between the two species.
    double distance(std::string sp1, std::string sp2) const;

    /// @brief Finds the two closest clusters in the distance matrix.
    /// This is used to determine which clusters to merge next in the fitch_margoliash algorithm.
    /// @return A pair of species names representing the closest clusters.
    /// If there are multiple pairs with the same distance, the first one found is returned.
    /// If no valid clusters are found, an exception is thrown.
    std::pair<std::string, std::string> closest_clusters();

private:
    distance_matrix _distance_matrix;
    std::map<std::string, clade*> _clusters;
};


#endif