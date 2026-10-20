#ifndef FREERATE_GLOBAL_MODEL_H
#define FREERATE_GLOBAL_MODEL_H

#include "core.h"
#include "clade.h"
#include "optimizer.h"

class gene_transcript;
typedef const std::vector<gene_transcript> transcript_vector;

class freerate_global_model : public model {

    clademap<std::pair<double, double>> _sigmas;
    bool _values_are_ratios;

    double optimize_sigmas(const user_data &ud, const clademap<prior> &priors);
    int _root_ape_index = -1;

    optimizer_parameters _optimizer_params;
public:
    //! Computation or estimation constructor
    freerate_global_model(bool values_are_ratios, optimizer_parameters params);
    
    virtual ~freerate_global_model() {}

    virtual double infer_transcript_likelihoods(const user_data &ud, const sigma_squared *p_sigma) override;

    virtual sigma_optimizer_scorer *get_sigma_optimizer(const user_data &data, const std::vector<std::string> &sample_groups) override;

    virtual reconstruction *reconstruct_ancestral_states(const user_data &ud, matrix_cache *p_calc);

    virtual void write_extra_vital_statistics(std::ostream& ost) override;

    virtual std::string get_name() const override { return "FreeRateGlobal"; }

    /// @brief Select initial values for sigmas at each node in the tree.
    /// @details This function initializes the sigma values for each node in the tree. It calculates
    /// a tree using the Fitch-Margoliash algorithm, basing the distance matrix on the Spearman
    /// correlation between the expression values of the transcripts across species.
    /// @param p_tree The tree for which to compute sigma values
    /// @param transcripts The transcripts used to compute the sigma values
    void initialize_sigmas(const clade* p_tree, const transcript_vector& transcripts);
};

#endif
