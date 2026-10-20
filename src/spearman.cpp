#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

#include "doctest.h"
#include "easylogging++.h"

#include "gene_transcript.h"

using namespace std;

// Helper: compute ranks for a vector
std::vector<double> compute_ranks(const std::vector<double>& values) {
    std::vector<size_t> idx(values.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::vector<double> ranks(values.size());
    std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) { return values[i] < values[j]; });

    transform(idx.begin(), idx.end(), ranks.begin(), [&](size_t i) {
        return i+1;
    });
    return ranks;
}

// Main function: returns map from species to Spearman correlation
map<pair<string, string>, double> spearman_correlation_by_species(const std::vector<gene_transcript>& transcripts) {
    std::map<std::string, std::vector<double>> species_to_values;

    if (transcripts.size() < 2) {
        throw std::runtime_error("At least two transcripts are required to compute Spearman correlation.");
    }
    // Gather expression values for each species
    for (const auto& transcript : transcripts) {
        for (const auto& species : transcript.get_species()) {
            species_to_values[species].push_back(transcript.get_expression_value(species));
        }
    }

    std::vector<std::pair<std::string, std::string>> species_pairs;
    for (auto it1 = species_to_values.begin(); it1 != species_to_values.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != species_to_values.end(); ++it2) {
            species_pairs.emplace_back(it1->first, it2->first);
        }
    }

    map<pair<string, string>, double> spearman_correlations;
    for (auto species_pair : species_pairs) {
        const auto& sp1 = species_pair.first;
        const auto& sp2 = species_pair.second;

        auto ranks1 = compute_ranks(species_to_values[sp1]);
        auto ranks2 = compute_ranks(species_to_values[sp2]);
        if (ranks1.size() != ranks2.size()) continue;
        double sum_d_squared = 0.0;
        for (size_t i = 0; i < ranks1.size(); ++i) {
            double d = ranks1[i] - ranks2[i];
            sum_d_squared += d * d;
        }
        double n = static_cast<double>(ranks1.size());
        double spearman_corr = 1.0 - (6.0 * sum_d_squared) / (n * (n * n - 1));
        
        spearman_correlations[species_pair] = spearman_corr;
    }

    return spearman_correlations;
}


TEST_CASE("spearman correlation_by_species computes correct correlations")
{
    // Work the example of spearman correlation here: 
    // https://www.statisticshowto.com/probability-and-statistics/correlation-coefficient-formula/spearman-rank-correlation-definition-calculate/
    std::vector<double> sp1 = {35, 23, 47, 17, 10, 43, 9, 6, 28};
    std::vector<double> sp2 = {30, 33, 45, 23, 8, 49, 12, 4, 31};

    std::vector<gene_transcript> transcripts(sp1.size());
    for (size_t i = 0; i < sp1.size(); i++)
    {
        transcripts[i].set_expression_value("species1", sp1[i]);
        transcripts[i].set_expression_value("species2", sp2[i]);    
    }
     
    auto result = spearman_correlation_by_species(transcripts);
    CHECK_EQ(1, result.size());
    CHECK_EQ(doctest::Approx(0.766667), result[make_pair("species1","species2")]); 
}

TEST_CASE("spearman_correlation_by_species throws when less than two transcripts")
{
    string msg("At least two transcripts are required to compute Spearman correlation.");
    std::vector<gene_transcript> transcripts;
    CHECK_THROWS_WITH(spearman_correlation_by_species(transcripts), msg.c_str());

    transcripts.push_back(gene_transcript("transcript1", "", ""));
    CHECK_THROWS_WITH(spearman_correlation_by_species(transcripts), msg.c_str());
}

TEST_CASE("compute_ranks")
{
    std::vector<double> values = {10, 20, 30, 25, 15};
    auto ranks = compute_ranks(values);
    CHECK_EQ(ranks.size(), values.size());
    CHECK_EQ(ranks[0], 1); 
    CHECK_EQ(ranks[1], 5); 
    CHECK_EQ(ranks[2], 2); 
    CHECK_EQ(ranks[3], 4); 
    CHECK_EQ(ranks[4], 3); 
}

