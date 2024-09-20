#include <ostream>
#include <iterator>
#include <numeric>

#include "easylogging++.h"
#include "doctest.h"

#include "clade.h"
#include "lambda.h"
#include "report.h"
#include "newick_ape_loader.h"
#include "gene_family.h"
#include "core.h"
#include "user_data.h"

using namespace std;

template<typename InputIt>
std::string join(InputIt first, InputIt last,
    const std::string& separator = ", ",
    const std::string& concluder = "")
{
    const std::string empty{};
    auto sep = &empty;

    std::stringstream ss;

    while (first != last)
    {
        ss << *sep << *first++;
        sep = &separator;
    }

    ss << concluder;

    return ss.str();
}

template<typename T, typename U>
std::ostream& operator<<(ostream& ost, const std::pair<T, U>& pair)
{
    ost << "(" << pair.first << "," << pair.second << ")";
    return ost;
}

std::ostream& operator<<(ostream& ost, const Report& report)
{
    ost << "Tree:";
    report._p_tree->write_newick(ost, [](const clade* c)
        {
            ostringstream ost;
            ost << (c->is_leaf() ? c->get_taxon_name() : "") << ":" << c->get_branch_length();
            return ost.str();
        });
    ost << "\n";

    if (report._p_lambda)
    {
        ost << "Lambda:\t";
        auto vals = get_lambda_values(report._p_lambda);
        copy(vals.begin(), vals.end(), ostream_iterator<double>(ost, "\t"));
    }
    ost << "\n";
    if (report._p_lambda_tree)
    {
        ost << "Lambda tree:\t";
        report._p_lambda_tree->write_newick(ost, [](const clade* c)
            {
                ostringstream ost;
                ost << c->get_lambda_index();
                return ost.str();
            });
        ost << "\n";
    }
    else
    {
        ost << "Lambda tree:\t";
        report._p_tree->write_newick(ost, [](const clade* c)
            {
                return "1";
            });
        ost << "\n";

    }

    auto t = get_ape_order(report._p_tree);
    map<const clade*, int> order;
    for (size_t i = 0; i < t.size(); ++i)
    {
        if (t[i]) order[t[i]] = i;
    }
    ost << "# IDs of nodes:";
    report._p_tree->write_newick(ost, [&order](const clade* c)
        {
            ostringstream ost;
            ost << (c->is_leaf() ? c->get_taxon_name() : "") << "<" << order[c] << ">";
            return ost.str();
        });
    ost << "\n";


    ost << "# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ";
    report._p_tree->apply_prefix_order([&ost, &order](const clade* c) {
        if (!c->is_leaf())
        {
            vector<int> v;
            transform(c->descendant_begin(), c->descendant_end(), back_inserter(v), [&order](const clade* c) { return order[c];  });
            ost << "(" << join(v.begin(), v.end(), ",") << ") ";
        }
    });

    ost << "\n";

    if (!report.average_expansion.empty())
    {
        ost << "Average Expansion:";
        report._p_tree->apply_prefix_order([&ost, &report](const clade* c) {
            if (!c->is_leaf())
            {
                vector<float> v;
                transform(c->descendant_begin(), c->descendant_end(), back_inserter(v), [&report](const clade* c) { return report.average_expansion.at(c);  });
                ost << "\t(" << join(v.begin(), v.end(), ",") << ")";
            }
            });
        ost << "\n";
    }

    if (!report.delta_count.empty())
    {
        report.write_delta(ost, "Expansion");
        report.write_delta(ost, "Remain");
        report.write_delta(ost, "Decrease");
    }

    ost << "'ID'\t'Newick'";

    for (auto item : report.families)
    {
        ost << item.node_id << "\t";
        ost << item.tree << "\t";
        ost << item.pvalue << "\t";
        report._p_tree->write_newick(ost, [&order](const clade* c)
            {
                ostringstream ost;
                ost << (c->is_leaf() ? c->get_taxon_name() : "") << "<" << order[c] << ">";
                return ost.str();
            });
        ost << endl;
#if 0
        for (size_t b = 0; b < item.pvalues.size(); b++)
        {
            std::pair<double, double> p = item.pvalues[b];
            if (p.first < 0)
            {
                ost << "(-,-)";
            }
            else
            {
                ost << "(" << p.first << "," << p.second << ")";
            }
            if (b < item.pvalues.size() - 1)
                ost << ",";
        }
        ost << ")\t";
#endif
    }
    return ost;
}

void Report::compute_expansion(const std::vector<gene_family>& gene_families, const reconstruction& reconstruct)
{
    _p_tree->apply_prefix_order([gene_families, &reconstruct, this](const clade* c) {
        if (!c->is_root())
        {
            vector<int> diffs(gene_families.size());
            transform(gene_families.begin(), gene_families.end(), diffs.begin(), [&reconstruct, c](const gene_family& gf) {
                return reconstruct.get_difference_from_parent(gf, c);
                });
            int total = accumulate(diffs.begin(), diffs.end(), 0);
            this->average_expansion[c] = float(total) / float(gene_families.size());
            this->delta_count[c].expanded = count_if(diffs.begin(), diffs.end(), [](int i) { return i > 0; });
            this->delta_count[c].decreased = count_if(diffs.begin(), diffs.end(), [](int i) { return i < 0; });
            this->delta_count[c].same = count_if(diffs.begin(), diffs.end(), [](int i) { return i == 0; });
        }
    });
}

family_line_item gene_family2report(const gene_family& gf, const clade* p_tree, reconstruction* r, double pvalue, const branch_probabilities& branch_probs)
{
    family_line_item fli;
    fli.node_id = gf.id();
    ostringstream ost;
    auto g = [gf, r](const clade* node) {
        ostringstream ost;
        if (node->is_leaf())
            ost << node->get_taxon_name();
        ost << "_" << r->get_node_count(gf, node);
        ost << ":" << node->get_branch_length();
        return ost.str();
    };

    p_tree->write_newick(ost, g);
    fli.tree = ost.str();

    fli.pvalue = pvalue;

    ostringstream ost2;
    p_tree->apply_prefix_order([&ost2, branch_probs, gf](const clade* c) {
        if (!c->is_leaf())
        {
            vector<double> v;
            transform(c->descendant_begin(), c->descendant_end(), back_inserter(v), [branch_probs, gf](const clade* c) {
                return branch_probs.at(gf, c)._value;  });

            ost2 << "(" << join(v.begin(), v.end(), ",") << ") ";
        }
        });
    fli.branch_pvalue_str = ost2.str();

    return fli;
}

void Report::add_line_item(const gene_family& gf, reconstruction* r, double pvalue, const branch_probabilities& branch_probs)
{
    if (!branch_probs.contains(gf))
        return;

    families.push_back(gene_family2report(gf, _p_tree, r, pvalue, branch_probs));
}

void Report::write_delta(ostream& ost, string header) const
{
    ost << header << " :";
    _p_tree->apply_prefix_order([&ost, this, header](const clade* c) {
        if (!c->is_leaf())
        {
            vector<float> v;
            transform(c->descendant_begin(), c->descendant_end(), back_inserter(v), [this, header](const clade* c)
                {
                    if (header == "Expansion") return delta_count.at(c).expanded;
                    if (header == "Remain") return delta_count.at(c).same;
                    if (header == "Decrease") return delta_count.at(c).decreased;
                    return 0;
                });
            ost << "\t(" << join(v.begin(), v.end(), ",") << ")";
        }
        });
    ost << "\n";
}


#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("Report writes tree")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r(p_tree.get(), nullptr, nullptr);

    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Tree:((A:1,B:1):1,(C:1,D:1):1):0");
}

TEST_CASE("Report writes lambda tree")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));
    unique_ptr<clade> p_tree2(parse_newick("((A:2,B:2):2,(C:3,D:3):1);", true));

    Report r(p_tree.get(), p_tree2.get(), nullptr);

    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Lambda tree:\t((2,2)2,(3,3)1)1");
}

TEST_CASE("Report writes dummy lambda tree if none specified")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:2,B:2):2,(C:3,D:3):1);"));

    Report r(p_tree.get(), nullptr, nullptr);

    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Lambda tree:\t((1,1)1,(1,1)1)1");
}

TEST_CASE("Report writes lambdas")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r(p_tree.get(), nullptr, new multiple_lambda(map<std::string, int>(), { 0.5, 0.3, 0.9 }));

    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Lambda:\t0.5\t0.3\t0.9");
}

TEST_CASE("Report writes IDs of nodes")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r(p_tree.get(), nullptr, new multiple_lambda(map<std::string, int>(), { 0.5, 0.3, 0.9 }));

    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "# IDs of nodes:((A<1>,B<2>)<6>,(C<3>,D<4>)<7>)<5>");
}

TEST_CASE("Report writes node pairs")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r(p_tree.get(), nullptr, new multiple_lambda(map<std::string, int>(), { 0.5, 0.3, 0.9 }));

    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, " (node ID, node ID): (6,7) (1,2) (3,4)");
}

class mock_reconstruction : public reconstruction   
{
    virtual int get_node_count(const gene_family& gf, const clade* c) const override 
    { 
        return values.at(gf.id())->find_descendant(c->get_taxon_name())->get_branch_length();
    }
    map<string, const clade*> values;
public:
    mock_reconstruction() : reconstruction(user_data(), input_parameters())
    {

    }

    void add(string id, string nwk)
    {
        values[id] = parse_newick(nwk);
    }
};

TEST_CASE("Report compute_expansion")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    mock_reconstruction r2;
    r2.add("f1", "((A:1,B:3):2,(C:2,D:1):2);");
    r2.add("f2", "((A:1,B:5):2,(C:7,D:1):4);");
    r2.add("f3", "((A:1,B:5):2,(C:3,D:1):3);");
    r2.add("f4", "((A:1,B:5):2,(C:5,D:1):3);");
    Report r(p_tree.get(), nullptr, nullptr);

    vector<gene_family> f(4);
    f[0].set_id("f1");
    f[1].set_id("f2");
    f[2].set_id("f3");
    f[3].set_id("f4");

    r.compute_expansion(f, r2);
    CHECK_EQ(-1, r.mean_expansion(p_tree->find_descendant("A")));
    CHECK_EQ(2.5, r.mean_expansion(p_tree->find_descendant("B")));
    CHECK_EQ(1.25, r.mean_expansion(p_tree->find_descendant("C")));
    CHECK_EQ(-2, r.mean_expansion(p_tree->find_descendant("D")));

    ostringstream ost;
    ost << r;
    CHECK_STREAM_CONTAINS(ost, "Average Expansion:\t(2,3)\t(-1,2.5)\t(1.25,-2)");
}

TEST_CASE("gene_family2report")
{
    vector<gene_family> f(4);
    f[0].set_id("f1");
    f[1].set_id("f2");
    f[2].set_id("f3");
    f[3].set_id("f4");

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    branch_probabilities bp;
    p_tree->apply_prefix_order([&bp, &f](const clade* c) { bp.set(f[0], c, 0.05); });

    mock_reconstruction r;
    r.add("f1", "((A:1,B:3):2,(C:2,D:1):2);");
    r.add("f2", "((A:1,B:5):2,(C:7,D:1):4);");
    r.add("f3", "((A:1,B:5):2,(C:3,D:1):3);");
    r.add("f4", "((A:1,B:5):2,(C:5,D:1):3);");
    // name_famsize:branchlength
    auto fli = gene_family2report(f[0], p_tree.get(), &r, 0.003, bp);
    CHECK_EQ("f1", fli.node_id);
    CHECK_EQ(0.003, fli.pvalue);
    CHECK_EQ("(0.05,0.05) (0.05,0.05) (0.05,0.05) ", fli.branch_pvalue_str);
    CHECK_EQ("((A_1:1,B_3:1)_2:1,(C_2:1,D_1:1)_2:1)_0:0", fli.tree);
}

TEST_CASE("Report compute_expansion expand count")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    mock_reconstruction r2;
    r2.add("f1", "((A:1,B:3):2,(C:2,D:1):2);");
    r2.add("f2", "((A:1,B:5):2,(C:7,D:1):4);");
    r2.add("f3", "((A:1,B:5):2,(C:3,D:1):3);");
    r2.add("f4", "((A:1,B:5):2,(C:5,D:1):3);");
    Report r(p_tree.get(), nullptr, nullptr);

    vector<gene_family> f(4);
    f[0].set_id("f1");
    f[1].set_id("f2");
    f[2].set_id("f3");
    f[3].set_id("f4");

    r.compute_expansion(f, r2);
    CHECK_EQ(0, r.delta(p_tree->find_descendant("A")).expanded);
    CHECK_EQ(4, r.delta(p_tree->find_descendant("B")).expanded);
    CHECK_EQ(2, r.delta(p_tree->find_descendant("C")).expanded);
    CHECK_EQ(0, r.delta(p_tree->find_descendant("D")).expanded);

    CHECK_EQ(4, r.delta(p_tree->find_descendant("A")).decreased);
    CHECK_EQ(0, r.delta(p_tree->find_descendant("B")).decreased);
    CHECK_EQ(0, r.delta(p_tree->find_descendant("C")).decreased);
    CHECK_EQ(4, r.delta(p_tree->find_descendant("D")).decreased);

    CHECK_EQ(0, r.delta(p_tree->find_descendant("A")).same);
    CHECK_EQ(0, r.delta(p_tree->find_descendant("B")).same);
    CHECK_EQ(2, r.delta(p_tree->find_descendant("C")).same);
    CHECK_EQ(0, r.delta(p_tree->find_descendant("D")).same);

    ostringstream ost;
    ost << r;
    CHECK_STREAM_CONTAINS(ost, "Expansion :\t(4,4)\t(0,4)\t(2,0)");
    CHECK_STREAM_CONTAINS(ost, "Remain :\t(0,0)\t(0,0)\t(2,0)");
    CHECK_STREAM_CONTAINS(ost, "Decrease :\t(0,0)\t(4,0)\t(0,4)");
}

