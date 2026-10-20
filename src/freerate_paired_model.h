#ifndef FREERATE_PAIRED_MODEL_H
#define FREERATE_PAIRED_MODEL_H

#include "core.h"
#include "clade.h"

class freerate_paired_model : public model {

    clademap<std::pair<double, double>> _sigmas;
    bool _values_are_ratios;
public:
    //! Computation or estimation constructor
    freerate_paired_model(bool values_are_ratios);
    
    virtual ~freerate_paired_model() {}

    virtual double infer_transcript_likelihoods(const user_data &ud, const sigma_squared *p_sigma) override;

    virtual sigma_optimizer_scorer *get_sigma_optimizer(const user_data &data, const std::vector<std::string> &sample_groups) override;

    virtual reconstruction *reconstruct_ancestral_states(const user_data &ud, matrix_cache *p_calc);

    virtual void write_extra_vital_statistics(std::ostream& ost) override;

    virtual std::string get_name() const override { return "FreeRatePaired"; }
};

#endif
