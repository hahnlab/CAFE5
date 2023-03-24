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
    double pvalue = 0.0;
    std::string branch_pvalue_str;
};

struct clade_delta { 
    int expanded = 0; 
    int decreased = 0; 
    int same = 0; 
};

class Report {
    const clade* _p_tree = nullptr;
    const clade* _p_lambda_tree = nullptr;
    const lambda* _p_lambda = nullptr;
    clademap<float> average_expansion;

    clademap<clade_delta> delta_count;
    std::vector<family_line_item> families;

    void write_delta(std::ostream& ost, std::string header) const;

public:
    void compute_expansion(const std::vector<gene_family>& gene_families, const reconstruction& reconstruct);
    void add_line_item(const gene_family& gf, reconstruction* r, double pvalue, const branch_probabilities& branch_probs);

    Report(const clade *p_tree, const clade *p_lambda_tree, const lambda *p_lambda) : _p_tree(p_tree), _p_lambda_tree(p_lambda_tree), _p_lambda(p_lambda)
    {

    }

    float mean_expansion(const clade* c) const { return average_expansion.at(c); }
    clade_delta delta(const clade* c) const { return delta_count.at(c); }
    friend std::ostream& operator<<(std::ostream& ost, const Report& report);

};

std::ostream& operator<<(std::ostream& ost, const Report& report);

#endif