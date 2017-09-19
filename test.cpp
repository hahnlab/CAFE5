#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "src/io.h"

TEST_GROUP(GeneFamilies)
{
};


TEST(GeneFamilies, read_gene_families_reads_cafe_files)
{
    auto families = read_gene_families("examples/genefam_data.txt", NULL);
    LONGS_EQUAL(5, families->at(0).get_species_size("A"));
    delete families;
}

int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}