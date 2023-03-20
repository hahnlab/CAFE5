#ifndef REPORT_H
#define REPORT_H

#include <iosfwd>

class clade;
class lambda;
class gene_family;
class reconstruction;
class branch_probabilities;

template<typename T>
using clademap = std::map<const clade*, T>;

struct family_line_item
{
    std::string node_id;
    std::string tree;
    double pvalue;
    std::string branch_pvalue_str;
};
family_line_item gene_family2report(const gene_family& gf, const clade* p_tree, reconstruction* r, double pvalue, const branch_probabilities& branch_probs);

struct Report {
    const clade* p_tree = nullptr;
    const clade* p_lambda_tree = nullptr;
    const lambda* p_lambda = nullptr;
    clademap<float> average_expansion;

    struct delta_counts { int expansion; int decrease; int remain; };
    clademap<delta_counts> delta_count;
    std::vector<family_line_item> families;

    void compute_expansion(const std::vector<gene_family>& gene_families, const reconstruction& reconstruct);
};

std::ostream& operator<<(std::ostream& ost, const Report& report);

#endif