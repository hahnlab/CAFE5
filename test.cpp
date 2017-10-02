#include <numeric>
#include <cmath>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"

TEST_GROUP(GeneFamilies)
{
};

TEST_GROUP(Inference)
{
};


TEST(GeneFamilies, read_gene_families_reads_cafe_files)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    auto families = read_gene_families(ist, NULL);
    LONGS_EQUAL(5, families->at(0).get_species_size("A"));
    LONGS_EQUAL(10, families->at(0).get_species_size("B"));
    LONGS_EQUAL(2, families->at(0).get_species_size("C"));
    LONGS_EQUAL(6, families->at(0).get_species_size("D"));
    delete families;
}

TEST(GeneFamilies, read_gene_families_reads_simulation_files)
{
    std::string str = "#A\n#B\n#AB\n#CD\n#C\n#ABCD\n#D\n35\t36\t35\t35\t36\t34\t34\t1\n98\t96\t97\t98\t98\t98\t98\t1\n";
    std::istringstream ist(str);

    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    auto families = read_gene_families(ist, p_tree);
    LONGS_EQUAL(35, families->at(0).get_species_size("A"));
    LONGS_EQUAL(36, families->at(0).get_species_size("B"));
    LONGS_EQUAL(36, families->at(0).get_species_size("C"));
    LONGS_EQUAL(34, families->at(0).get_species_size("D"));
    delete families;
    delete p_tree;
}

TEST(Inference, infer_processes)
{
    vector<gene_family> families;
    gene_family fam;
    fam.set_species_size("A", 2);
    fam.set_species_size("B", 5);
    fam.set_species_size("C", 6);
    fam.set_species_size("D", 10);
    families.push_back(fam);
    base_core core;
    probability_calculator calc;
    single_lambda lambda(&calc, 0.01);
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    core.set_gene_families(&families);
    core.set_max_sizes(148, 122);
    core.set_lambda(&lambda);
    core.set_tree(p_tree);
    core.start_inference_processes();
    core.infer_processes();
    //core.get_likelihoods();
    //DOUBLES_EQUAL(0.05, likelihoods[1], 0.0001);
    delete p_tree;
}

TEST(Inference, equilibrium_frequency)
{
    base_core core;
    std::vector<int> rd(10);
    std::fill(rd.begin(), rd.end(), 1);
    equilibrium_frequency ef(rd);
    core.set_rootdist_vec(rd);
    DOUBLES_EQUAL(.1, ef.compute(5), 0.0001);
}

TEST(Inference, gamma_set_alpha)
{
    gamma_core core;
    core.adjust_n_gamma_cats(5);
    core.set_alpha(0.5, 3);
}

TEST(Inference, gamma_adjust_family_gamma_membership)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    auto families = read_gene_families(ist, NULL);
    delete families;

    gamma_core core;
    core.adjust_family_gamma_membership(5);
}

TEST(Inference, gamma)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6\n";
    std::istringstream ist(str);
    auto families = read_gene_families(ist, NULL);
    probability_calculator calc;
    single_lambda lambda(&calc, 0.5);
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    std::vector<int> rootdist;
//    gamma_core core(cout, &lambda, p_tree, 122, 4, rootdist, 2, 0.5);
    gamma_core core;
    core.set_gene_families(families);
    core.set_lambda(&lambda);
    core.set_tree(p_tree);
    core.initialize_with_alpha(2, 4, 0.5);
    core.set_max_sizes(148, 122);

    core.start_inference_processes();
    double actual = core.infer_processes();
    DOUBLES_EQUAL(-56.3469, std::log(actual), .0001);

    delete families;
    delete p_tree;
}

TEST(Inference, stash_stream)
{
    family_info_stash stash;
    stash.family_id = 1;
    stash.lambda_multiplier = 2.5;
    stash.family_likelihood = 3.7;
    stash.posterior_probability = 4.9;

    std::ostringstream ost;
    ost << stash;
    STRCMP_EQUAL("1\t2.5\t0\t3.7\t4.9\tN/S", ost.str().c_str());

}

int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}