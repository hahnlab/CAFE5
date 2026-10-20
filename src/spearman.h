#ifndef SPEARMAN_H
#define SPEARMAN_H

#include <map>
#include <vector>
#include <string>

class gene_transcript;

/// @brief Computes the Spearman correlation for each pair of species based on the expression values of the transcripts.
/// @details The function iterates through the transcripts and computes the Spearman correlation for each pair of species.
/// It returns a map where the keys are pairs of species names and the values are the computed Spearman correlation coefficients.
/// @param transcripts A vector of gene_transcript objects containing the expression values for each species.
/// @return A map where the keys are pairs of species names and the values are the Spearman correlation coefficients.
std::map<std::pair<std::string, std::string>, double> spearman_correlation_by_species(const std::vector<gene_transcript>& transcripts);


#endif