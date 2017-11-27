#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H
#include "process.h"

class reconstruction_process : public process {

    gene_family *_gene_family;
    probability_calculator *_p_calc;
    root_equilibrium_distribution* _p_prior;

    void reconstruct_internal_node(clade * c, lambda * sl);
    void reconstruct_leaf_node(clade * c, lambda * sl);
    void reconstruct_root_node(clade * c);
    std::map<clade *, std::vector<int> > all_node_Cs;
    std::map<clade *, std::vector<double> > all_node_Ls;
public:
    void reconstruct();

    reconstruction_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree,
        int max_family_size,
        int max_root_family_size, std::vector<int> rootdist,
        gene_family *gf,
        probability_calculator *calc,
        root_equilibrium_distribution* p_prior);


    void operator()(clade *c);
};
#endif
