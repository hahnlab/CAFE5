#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "src/io.h"

TEST_GROUP(GeneFamilies)
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

int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}