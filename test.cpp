
#include <numeric>
#include <cmath>
#include <getopt.h>
#include <sstream>
#include <random>

#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"
#include "src/root_equilibrium_distribution.h"
#include "src/base_model.h"
#include "src/process.h"
#include "src/gene_family_reconstructor.h"
#include "src/matrix_cache.h"
#include "src/gamma_bundle.h"
#include "src/probability.h"
#include "src/execute.h"
#include "src/user_data.h"
#include "src/optimizer_scorer.h"
#include "src/root_distribution.h"
#include "src/simulator.h"
#include "src/poisson.h"
#include "src/optimizer.h"

// these need to be at the end to stop weird STL errors
#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"

std::mt19937 randomizer_engine(10); // seeding random number engine

TEST_GROUP(GeneFamilies)
{
};

TEST_GROUP(Inference)
{
    user_data _user_data;
    single_lambda *_p_lambda;

    void setup()
    {
        _p_lambda = new single_lambda(0.05);
        newick_parser parser(false);
        parser.newick_string = "(A:1,B:1);";
        _user_data.p_tree = parser.parse_newick();
        _user_data.p_lambda = _p_lambda;
        _user_data.max_family_size = 10;
        _user_data.max_root_family_size = 8;
        _user_data.gene_families.resize(1);
        _user_data.gene_families[0].set_species_size("A", 1);
        _user_data.gene_families[0].set_species_size("B", 2);

        randomizer_engine.seed(10);
    }

    void teardown()
    {
        delete _user_data.p_tree;
        delete _p_lambda;
    }
};

TEST_GROUP(Simulation)
{
    void setup()
    {
        randomizer_engine.seed(10);
    }
};

TEST_GROUP(Probability)
{
};

TEST_GROUP(Clade)
{
};

TEST_GROUP(Optimizer)
{
};

TEST_GROUP(Options)
{
    char *values[100];
    size_t argc;

    void initialize(vector<string> arguments)
    {
        optind = 0;
        argc = arguments.size();
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            values[i] = strdup(arguments[i].c_str());
        }
    }

    void teardown()
    {
        for (size_t i = 0; i < argc; ++i)
        {
            free(values[i]);
        }
    }
};

input_parameters read_arguments(int argc, char *const argv[]);

TEST(Options, input_short)
{
    initialize({ "cafexp", "-ifile" });

    auto actual = read_arguments(argc, values);
    STRCMP_EQUAL("file", actual.input_file_path.c_str());
}

TEST(Options, input_long)
{
    initialize({ "cafexp", "--infile", "file" });

    auto actual = read_arguments(argc, values);
    STRCMP_EQUAL("file", actual.input_file_path.c_str());
}

TEST(Options, input_short_space_separated)
{
    initialize({ "cafexp", "-i", "file" });

    auto actual = read_arguments(argc, values);
    STRCMP_EQUAL("file", actual.input_file_path.c_str());
}

TEST(Options, simulate_long)
{
    initialize({ "cafexp", "--simulate=1000", "-l", "0.05" });

    auto actual = read_arguments(argc, values);
    CHECK_EQUAL(1000, actual.nsims);
}

TEST(Options, simulate_short)
{
    initialize({ "cafexp", "-s1000", "-l", "0.05"});

    auto actual = read_arguments(argc, values);
    CHECK_EQUAL(1000, actual.nsims);
}

TEST(Options, cannot_have_space_before_optional_parameter)
{
    try
    {
        initialize({ "cafexp", "-s", "1000" });

        auto actual = read_arguments(argc, values);
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Unrecognized parameter: '1000'", err.what());
    }
}

TEST(Options, must_specify_lambda_and_input_file_for_estimator)
{
    try
    {
        input_parameters params;
        params.fixed_lambda = 0.05;
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Options -l and -i must both be provided an argument.", err.what());
    }
}

TEST(Options, must_specify_lambda_for_simulation)
{
    try
    {
        input_parameters params;
        params.is_simulating = true;
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Cannot simulate without initial lambda values", err.what());
    }
}

TEST(Options, must_specify_alpha_for_gamma_simulation)
{
    try
    {
        input_parameters params;
        params.is_simulating = true;
        params.fixed_lambda = 0.05;
        params.n_gamma_cats = 3;
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Cannot simulate gamma clusters without an alpha value", err.what());
    }
}

TEST(Options, must_specify_alpha_and_k_for_gamma_inference)
{
    try
    {
        input_parameters params;
        params.fixed_alpha = 0.7;
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Alpha specified with 1 gamma category.", err.what());
    }
}

TEST(Options, can_specify_alpha_without_k_for_gamma_simulation)
{
    try
    {
        input_parameters params;
        params.fixed_alpha = 0.7;
        params.fixed_lambda = 0.01;
        params.is_simulating = true;
        params.check_input();
    }
    catch (runtime_error& err)
    {
        FAIL("Exception thrown checking input")
    }
}

TEST(Options, check_input_does_not_throw_when_simulating_with_multiple_lambdas)
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_multiple_lambdas = "0.01,0.05";
    params.lambda_tree_file_path = "./tree";
    params.check_input();
    CHECK(true);
}

TEST(Options, per_family_must_provide_families)
{
    try
    {
        input_parameters params;
        params.lambda_per_family = true;
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("No family file provided", err.what());
    }
}

TEST(Options, per_family_must_provide_tree)
{
    try
    {
        input_parameters params;
        params.lambda_per_family = true;
        params.input_file_path = "/tmp/test";
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("No tree file provided", err.what());
    }
}

TEST(GeneFamilies, read_gene_families_reads_cafe_files)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);
    LONGS_EQUAL(5, families.at(0).get_species_size("A"));
    LONGS_EQUAL(10, families.at(0).get_species_size("B"));
    LONGS_EQUAL(2, families.at(0).get_species_size("C"));
    LONGS_EQUAL(6, families.at(0).get_species_size("D"));
}

TEST(GeneFamilies, read_gene_families_reads_simulation_files)
{
    std::string str = "#A\n#B\n#AB\n#CD\n#C\n#ABCD\n#D\n35\t36\t35\t35\t36\t34\t34\t1\n98\t96\t97\t98\t98\t98\t98\t1\n";
    std::istringstream ist(str);

    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    std::vector<gene_family> families;
    read_gene_families(ist, p_tree, &families);
    LONGS_EQUAL(35, families.at(0).get_species_size("A"));
    LONGS_EQUAL(36, families.at(0).get_species_size("B"));
    LONGS_EQUAL(36, families.at(0).get_species_size("C"));
    LONGS_EQUAL(34, families.at(0).get_species_size("D"));
    delete p_tree;
}

TEST(Inference, infer_processes)
{
    vector<gene_family> families;
    gene_family fam;
    fam.set_species_size("A", 1);
    fam.set_species_size("B", 2);
    families.push_back(fam);
    gene_family fam2;
    fam2.set_species_size("A", 2);
    fam2.set_species_size("B", 1);
    families.push_back(fam2);
    gene_family fam3;
    fam3.set_species_size("A", 3);
    fam3.set_species_size("B", 6);
    families.push_back(fam3);
    gene_family fam4;
    fam4.set_species_size("A", 6);
    fam4.set_species_size("B", 3);
    families.push_back(fam4);

    single_lambda lambda(0.01);

    base_model core(&lambda, _user_data.p_tree, &families, 56, 30, NULL);
    core.start_inference_processes(&lambda);


    uniform_distribution frq;
    double multi = core.infer_processes(&frq, std::map<int, int>());
    //core.get_likelihoods();
    DOUBLES_EQUAL(41.7504, multi, 0.001);
}

TEST(Inference, uniform_distribution)
{
    root_distribution rd;
    rd.vectorize_uniform(10);
    uniform_distribution ef;
    ef.initialize(&rd);
    DOUBLES_EQUAL(.1, ef.compute(5), 0.0001);
}

TEST(Inference, gamma_set_alpha)
{
    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL);
    model.set_alpha(0.5);
}

TEST(Inference, gamma_adjust_family_gamma_membership)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);

    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL);
}

TEST(Inference, gamma_model_infers_processes_without_crashing)
{
    std::vector<int> rootdist;

    gamma_model core(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 0, 5, 1, 0, NULL);

    core.set_max_sizes(148, 122);

    core.start_inference_processes(_user_data.p_lambda);
    uniform_distribution frq;

    // TODO: make this return a non-infinite value and add a check for it
    core.infer_processes(&frq, std::map<int, int>());
    
}

TEST(Inference, stash_stream)
{
    family_info_stash stash;
    stash.family_id = "F01";
    stash.lambda_multiplier = 2.5;
    stash.family_likelihood = 3.7;
    stash.posterior_probability = 4.9;

    std::ostringstream ost;
    ost << stash;
    STRCMP_EQUAL("F01\t2.5\t0\t3.7\t4.9\tN/S", ost.str().c_str());

}

TEST(Probability, probability_of_some_values)
{
    matrix_cache calc(0);
    double lambda = 0.05;
    double branch_length = 5;
    DOUBLES_EQUAL(0.0152237, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 5, 9), 0.00001);

    DOUBLES_EQUAL(0.17573, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 9), 0.00001);

    DOUBLES_EQUAL(0.182728, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10), 0.00001);

    branch_length = 1;
    DOUBLES_EQUAL(0.465565, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10), 0.00001);
}

bool operator==(const matrix& m1, const matrix& m2)
{
    if (m1.size() != m2.size())
        return false;

    for (int i = 0; i < m1.size(); ++i)
    {
        for (int j = 0; j < m1.size(); ++j)
            if (abs(m1.get(i,j) - m2.get(i,j)) > 0.00001)
                return false;
    }

    return true;
}

TEST(Probability, matrices_take_fractional_branch_lengths_into_account)
{
    matrix_cache calc(141);
    single_lambda lambda(0.006335);
    std::set<double> branch_lengths{ 68, 68.7105 };
    calc.precalculate_matrices(get_lambda_values(&lambda), branch_lengths);
    DOUBLES_EQUAL(0.194661, calc.get_matrix(68.7105, 0.006335)->get(5,5), 0.00001); // a value 
    DOUBLES_EQUAL(0.195791, calc.get_matrix(68, 0.006335)->get(5, 5), 0.00001);
}

TEST(Probability, the_probability_of_going_from_parent_fam_size_to_c)
{
    DOUBLES_EQUAL(0.194661, the_probability_of_going_from_parent_fam_size_to_c(.006335, 68.7105, 5, 5), 0.00001);
}

TEST(Probability, probability_of_matrix)
{
    matrix_cache calc(5);
    single_lambda lambda(0.05);
    std::set<double> branch_lengths{ 5 };
    calc.precalculate_matrices(get_lambda_values(&lambda), branch_lengths);
    auto actual = calc.get_matrix(5, lambda.get_single_lambda());
    matrix expected(5);
    double values[5][5] = {
    {1,0,0,0,0},
    { 0.2,0.64,0.128,0.0256,0.00512 },
    { 0.04,0.256,0.4608,0.17408,0.0512 },
    { 0.008,0.0768,0.26112,0.36352,0.187392 },
    { 0.0016,0.02048,0.1024,0.249856,0.305562 } };
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            expected.set(i, j, values[i][j]);
    CHECK(*actual == expected);

    // a second call should get the same results as the first
    actual = calc.get_matrix(5, lambda.get_single_lambda());
    CHECK(*actual == expected);
}

TEST(Probability, get_random_probabilities)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lam(0.05);
    matrix_cache cache(15);
    cache.precalculate_matrices(vector<double>{0.05}, set<double>{1});
    auto probs = get_random_probabilities(p_tree.get(), 10, 3, 12, 8, &lam, cache, NULL);
    LONGS_EQUAL(10, probs.size());
    DOUBLES_EQUAL(0.00390567, probs[0], 0.0001);
}

TEST(Inference, base_optimizer_guesses_lambda_only)
{
    base_model model(_user_data.p_lambda, _user_data.p_tree, NULL, 0, 5, NULL);
    user_data data;
    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(1, guesses.size());
    DOUBLES_EQUAL(0.298761, guesses[0], 0.0001);
    delete model.get_lambda();
}

TEST(Inference, base_optimizer_guesses_lambda_and_unique_epsilons)
{
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    base_model model(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 10, 10, &err);

    _user_data.p_error_model = &err;
    _user_data.p_lambda = nullptr;
    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(_user_data));

    CHECK(opt);
    CHECK(dynamic_cast<lambda_epsilon_optimizer*>(opt.get()) != NULL);
#ifdef FALSE_MEMORY_LEAK_FIXED
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(3, guesses.size());
    DOUBLES_EQUAL(0.0808301, guesses[0], 0.0001);
    DOUBLES_EQUAL(0.3, guesses[1], 0.0001);
    DOUBLES_EQUAL(0.4, guesses[2], 0.0001);
    vector<double>().swap(guesses);
#endif
    delete model.get_lambda();
}


TEST(Inference, gamma_optimizer_guesses_lambda_and_alpha)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    gamma_model model(NULL, p_tree, NULL, 0, 5, 4, -1, NULL);
    user_data data;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));
    CHECK(opt);
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(2, guesses.size());
    DOUBLES_EQUAL(0.298761, guesses[0], 0.0001);
    DOUBLES_EQUAL(-0.172477, guesses[1], 0.0001);
    vector<double>().swap(guesses);
   delete p_tree;
   delete model.get_lambda();
}

TEST(Inference, gamma_optimizer_guesses_lambda_if_alpha_provided)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> p_tree(parser.parse_newick());

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, 0.25, NULL);

    user_data data;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));

    CHECK(opt);
    CHECK(dynamic_cast<lambda_optimizer *>(opt.get()));

    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(1, guesses.size());
    DOUBLES_EQUAL(0.298761, guesses[0], 0.0001);
    vector<double>().swap(guesses);
    delete model.get_lambda();
}

TEST(Inference, gamma_optimizer_guesses_alpha_if_lambda_provided)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> p_tree(parser.parse_newick());

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, -1, NULL);

    user_data data;
    single_lambda sl(0.05);
    data.p_lambda = &sl;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));

    CHECK(opt);
    CHECK(dynamic_cast<gamma_optimizer *>(opt.get()));

    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(1, guesses.size());
    DOUBLES_EQUAL(0.979494, guesses[0], 0.0001);
    vector<double>().swap(guesses);
    delete model.get_lambda();
}

TEST(Inference, gamma_optimizer_optimizes_nothing_if_lambda_and_alpha_provided)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> p_tree(parser.parse_newick());

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, .25, NULL);

    user_data data;
    single_lambda sl(0.05);
    data.p_lambda = &sl;

    CHECK(model.get_lambda_optimizer(data) == nullptr);
}

TEST(Inference, base_model_reconstruction)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1";
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda sl(0.05);

    std::vector<gene_family> families(1);
    families[0].set_species_size("A", 3);
    families[0].set_species_size("B", 4);

    base_model model(&sl, p_tree.get(), &families, 5, 5, NULL);

    matrix_cache calc(6);
    calc.precalculate_matrices(get_lambda_values(&sl), set<double>({ 1 }));
    root_distribution rd;
    rd.vectorize_increasing(6);
    uniform_distribution dist;
    dist.initialize(&rd);
    std::unique_ptr<reconstruction> rec(model.reconstruct_ancestral_states(&calc, &dist));

}

TEST(Inference, branch_length_finder)
{
    branch_length_finder finder;
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
    clade *p_tree = parser.parse_newick();
    p_tree->apply_prefix_order(finder);
    LONGS_EQUAL(finder.longest(), 23);
    auto expected = set<double>{ 1, 3, 7, 11, 17, 23 };
    CHECK(finder.result() == expected);
    delete p_tree;

}

TEST(Inference, increase_decrease)
{
    clademap<int> family_size;
    clademap<family_size_change> result;

    newick_parser parser(false);
    parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
    unique_ptr<clade> p_tree(parser.parse_newick());

     auto a = p_tree->find_descendant("A");
     auto b = p_tree->find_descendant("B");
     auto ab = p_tree->find_descendant("AB");
     auto abcd = p_tree->find_descendant("ABCD");

    family_size[a] = 2;
    family_size[b] = 4;
    family_size[ab] = 3;
    family_size[abcd] = 3;
    compute_increase_decrease(family_size, result);

    CHECK(result[a] == Decrease);
    CHECK(result[b] == Increase);
    CHECK(result[ab] == Constant);
}

TEST(Inference, increase_decrease_stream)
{
    ostringstream ost;
    increase_decrease id;
    id.gene_family_id = "1234";
    id.pvalue = 0.02;
    id.change = vector<family_size_change>{ Decrease, Constant, Increase, Increase};
    ost << id;
    STRCMP_EQUAL("1234\t0.02\ty\td\tc\ti\ti\t\n", ost.str().c_str());
}

TEST(Inference, gamma_increase_decrease_stream)
{
    ostringstream ost;
    increase_decrease id;
    id.gene_family_id = "1234";
    id.change = vector<family_size_change>{ Decrease, Constant, Increase, Increase };
    id.category_likelihoods = vector<double>{ .1,.2,.3,.4 };
    id.pvalue = 0.02;
    ost << id;
    STRCMP_EQUAL("1234\t0.02\ty\td\tc\ti\ti\t0.1\t0.2\t0.3\t0.4\t\n", ost.str().c_str());
}

TEST(Inference, precalculate_matrices_calculates_all_lambdas_all_branchlengths)
{
    matrix_cache calc(5);
    std::map<std::string, int> m;
    multiple_lambda lambda(m, vector<double>({ .1, .2, .3, .4 }));
    calc.precalculate_matrices(get_lambda_values(&lambda), set<double>({ 1,2,3 }));
    LONGS_EQUAL(12, calc.get_cache_size());
}

TEST_GROUP(Reconstruction)
{
    gene_family fam;
    unique_ptr<clade> p_tree;
    cladevector order;

    void setup()
    {
        newick_parser parser(false);
        parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
        p_tree.reset(parser.parse_newick());

        fam.set_id("Family5");
        fam.set_species_size("A", 11);
        fam.set_species_size("B", 2);
        fam.set_species_size("C", 5);
        fam.set_species_size("D", 6);

        vector<string> nodes{ "A", "B", "C", "D", "AB", "CD", "ABCD" };
        order.resize(nodes.size());
        const clade *t = p_tree.get();
        transform(nodes.begin(), nodes.end(), order.begin(), [t](string s) { return t->find_descendant(s); });

    }
};

TEST(Reconstruction, gene_family_reconstructor)
{
  single_lambda lambda(0.05);
  fam.set_species_size("Mouse", 3);

  clade leaf("Mouse", 7);

  matrix_cache calc(8);
  calc.precalculate_matrices({ .1 }, set<double>({ 7 }));
  gene_family_reconstructor process(cout, &lambda, 2.0, NULL, 7, 0, &fam, &calc, NULL);
  process(&leaf);
  auto L = process.get_L(&leaf);

  // L holds the probability of the leaf moving from size 3 to size n
  LONGS_EQUAL(8, L.size());
  DOUBLES_EQUAL(0.0, L[0], 0.0001);
  DOUBLES_EQUAL(0.0586679, L[1], 0.0001);
  DOUBLES_EQUAL(0.146916, L[2], 0.0001);
  DOUBLES_EQUAL(0.193072, L[3], 0.0001);
}

TEST(Reconstruction, gene_family_reconstructor_print_reconstruction)
{
    // auto a = p_tree->find_descendant("A");
    // cladevector order({ p_tree.get(), a });
    clademap<int> values;
    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    gene_family_reconstructor gfc(&fam, p_tree.get(), values);
    ostringstream ost;

    gfc.print_reconstruction(ost, order);
    STRCMP_EQUAL("  TREE Family5 = ((A_11:1,B_2:3)4_8:7,(C_5:11,D_6:17)5_6:23)6_7;\n", ost.str().c_str());
}

TEST(Reconstruction, gamma_bundle_print_reconstruction_prints_value_for_each_category_and_a_summation)
{
    // auto a = p_tree->find_descendant("A");
    clademap<int> values;
    values[p_tree.get()] = 7;

    clademap<double> values2;
    values2[p_tree.get()] = 7;
    values2[p_tree->find_descendant("AB")] = 8;
    values2[p_tree->find_descendant("CD")] = 6;

    auto gfc = new gene_family_reconstructor(&fam, p_tree.get(), values);
    vector<gene_family_reconstructor *> v{ gfc };

    ostringstream ost;
    gamma_bundle bundle(v, values2, clademap<family_size_change>(), p_tree.get(), &fam);
    bundle.print_reconstruction(ost, order);
    STRCMP_EQUAL("  TREE Family5 = ((A_11:1,B_2:3)4_0_8:7,(C_5:11,D_6:17)5_0_6:23)6_7_7;\n", ost.str().c_str());

}

TEST(Reconstruction, gamma_bundle_get_increases_decreases)
{
    auto a = p_tree->find_descendant("A");
    auto ab = p_tree->find_descendant("AB");
    order = { p_tree.get(), a, ab };

    clademap<int> values;
    values[p_tree.get()] = 7;
    values[a] = 11;
    values[ab] = 13;

    clademap<double> values2;
    values2[p_tree.get()] = 7;
    values2[a] = 11;
    values2[ab] = 13;

    auto gfc = new gene_family_reconstructor(&fam, p_tree.get(), values);
    vector<gene_family_reconstructor *> v{ gfc };
    clademap<family_size_change> c;
    c[ab] = Decrease;

    ostringstream ost;
    gamma_bundle bundle(v, values2, c, p_tree.get(), &fam);
    auto actual = bundle.get_increases_decreases(order, 0.05);
    DOUBLES_EQUAL(0.05, actual.pvalue, 0.000001);
    LONGS_EQUAL(3, actual.change.size());
    LONGS_EQUAL(Constant, actual.change[0]);
    LONGS_EQUAL(Constant, actual.change[1]);
    LONGS_EQUAL(Decrease, actual.change[2]);
}

TEST(Reconstruction, gamma_model_reconstruction)
{
    // auto a = p_tree->find_descendant("A");

    clademap<int> values;
    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    clademap<double> values2;
    values2[p_tree.get()] = 7;
    values2[p_tree->find_descendant("AB")] = 8;
    values2[p_tree->find_descendant("CD")] = 6;

    auto gfc = new gene_family_reconstructor(&fam, p_tree.get(), values);
    vector<gene_family_reconstructor *> v{ gfc };
    gamma_bundle bundle(v, values2, clademap<family_size_change>(), p_tree.get(), &fam);
    vector<gamma_bundle *> bundles{ &bundle };

    vector<double> multipliers{ 0.13, 1.4 };
    gamma_model_reconstruction gmr(multipliers, bundles);
    std::ostringstream ost;
    gmr.print_reconstructed_states(ost);

    STRCMP_CONTAINS("#NEXUS", ost.str().c_str());
    STRCMP_CONTAINS("BEGIN TREES;", ost.str().c_str());
    STRCMP_CONTAINS("  TREE Family5 = ((A_11:1,B_2:3)1_8_8:7,(C_5:11,D_6:17)2_6_6:23)0_7_7;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());

    STRCMP_CONTAINS("BEGIN LAMBDA_MULTIPLIERS;", ost.str().c_str());
    STRCMP_CONTAINS("  0.13;", ost.str().c_str());
    STRCMP_CONTAINS("  1.4;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST(Reconstruction, print_reconstructed_states_empty)
{
    base_model_reconstruction bmr;
    ostringstream ost;
    bmr.print_reconstructed_states(ost);
    STRCMP_EQUAL("", ost.str().c_str());
}

TEST(Reconstruction, print_reconstructed_states_no_print)
{
    // auto a = p_tree->find_descendant("A");

    clademap<int> values;
    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    base_model_reconstruction bmr;
    bmr._rec_processes.push_back(new gene_family_reconstructor(&fam, p_tree.get(), values));
    ostringstream ost;
    bmr.print_reconstructed_states(ost);
    STRCMP_CONTAINS("#NEXUS", ost.str().c_str());
    STRCMP_CONTAINS("BEGIN TREES;", ost.str().c_str());
    STRCMP_CONTAINS("  TREE Family5 = ((A_11:1,B_2:3)1_8:7,(C_5:11,D_6:17)2_6:23)0_7;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST(Reconstruction, reconstruction_process_internal_node)
{
    single_lambda s_lambda(0.05);
    double multiplier = 2.0;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);

    unique_ptr<lambda> m(s_lambda.multiply(multiplier));
    matrix_cache calc(25);
    calc.precalculate_matrices(get_lambda_values(m.get()), set<double>({ 1, 3, 7, 11, 17, 23 }));
    gene_family_reconstructor process(cout, &s_lambda, multiplier, NULL, 24, 24, &fam, &calc, NULL);

    process(p_tree->find_descendant("A"));
    process(p_tree->find_descendant("B"));

    auto internal_node = p_tree->find_descendant("AB");
    process(internal_node);
    auto L = process.get_L(internal_node);

    // L holds the probability of the leaf moving from size 3 to size n
    LONGS_EQUAL(25, L.size());
    DOUBLES_EQUAL(0.0, L[0], 0.0001);
    DOUBLES_EQUAL(0.00101688, L[1], 0.0001);
    DOUBLES_EQUAL(0.00254648, L[2], 0.0001);
    DOUBLES_EQUAL(0.0033465, L[3], 0.0001);
}

TEST(Inference, gamma_bundle_prune)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda lambda(0.005);
    inference_process_factory factory(cout, &lambda, p_tree.get(), 10, 8);
    factory.set_gene_family(&fam);
    gamma_bundle bundle(factory, { 0.1, 0.5 }, p_tree.get(), &fam);

    root_distribution rd;
    rd.vector({ 1,2,3,4,5,4,3,2,1 });
    uniform_distribution dist;
    dist.initialize(&rd);
    matrix_cache cache(11);
    multiple_lambda ml(map<string, int>(), {0.0005, 0.0025});
    cache.precalculate_matrices(get_lambda_values(&ml), set<double>{1, 3, 7});
    CHECK(bundle.prune({ 0.01, 0.05 }, &dist, cache)); 
    auto cat_likelihoods = bundle.get_category_likelihoods();

    LONGS_EQUAL(2, cat_likelihoods.size());
    DOUBLES_EQUAL(-23.3728, log(cat_likelihoods[0]), 0.0001);
    DOUBLES_EQUAL(-17.0086, log(cat_likelihoods[1]), 0.0001);

}

TEST(Inference, gamma_bundle_prune_returns_false_if_saturated)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda lambda(0.9);
    inference_process_factory factory(cout, &lambda, p_tree.get(), 10, 8);
    factory.set_gene_family(&fam);
    gamma_bundle bundle(factory, { 0.1, 0.5 }, p_tree.get(), &fam);

    root_distribution rd;
    rd.vector({ 1,2,3,4,5,4,3,2,1 });
    uniform_distribution dist;
    dist.initialize(&rd);
    matrix_cache cache(11);
    cache.precalculate_matrices({ 0.09, 0.45 }, set<double>{1, 3, 7});

    CHECK(!bundle.prune({ 0.01, 0.05 }, &dist, cache));
}

TEST(Inference, matrix_cache_key_handles_floating_point_imprecision)
{
    set<matrix_cache_key> keys;
    double t = 0.0;
    for (int i = 0; i < 31; i++)
    {
        t += 0.1;
        matrix_cache_key key(1, t, 0.3);
        keys.insert(key);
    }
    LONGS_EQUAL(31, keys.size());

    matrix_cache_key key(1, 3.0, 0.3);
    LONGS_EQUAL(1, keys.count(key));
}

TEST(Inference, birthdeath_rate_with_log_alpha)
{
    // alpha and coeff are derived values from lambda and t
    // (alpha = lambda*t / 1 + lambda*t, coeff = 1 - 2 * alpha);
    DOUBLES_EQUAL(-1.55455, log(birthdeath_rate_with_log_alpha(46, 45, -3.672556, 0.949177)), 0.00001);
    DOUBLES_EQUAL(-2.20436, log(birthdeath_rate_with_log_alpha(44, 46, -2.617970, 0.854098)), 0.00001);
    DOUBLES_EQUAL(-2.39974, log(birthdeath_rate_with_log_alpha(43, 43, -1.686354, 0.629613)), 0.00001);
    DOUBLES_EQUAL(-2.44301, log(birthdeath_rate_with_log_alpha(43, 44, -1.686354, 0.629613)), 0.00001);
    DOUBLES_EQUAL(-1.58253, log(birthdeath_rate_with_log_alpha(13, 14, -2.617970, 0.854098)), 0.00001);

    DOUBLES_EQUAL(0.107, birthdeath_rate_with_log_alpha(40, 42, -1.37, 0.5), .001);
    DOUBLES_EQUAL(0.006, birthdeath_rate_with_log_alpha(41, 34, -1.262, 0.4), .001);

    DOUBLES_EQUAL(0.194661, birthdeath_rate_with_log_alpha(5, 5, -1.1931291703283662, 0.39345841643135504), 0.0001);
}

TEST(Inference, create_one_model_if_lambda_is_null)
{
    input_parameters params;
    params.input_file_path = "foo";
    user_data data;
    auto models = build_models(params, data);
    LONGS_EQUAL(1, models.size());
    CHECK(dynamic_cast<base_model *>(models[0]));
    for (auto m : models)
        delete m;
}

TEST(Inference, create_gamma_model_if_alpha_provided)
{
    input_parameters params;
    params.input_file_path = "foo";
    params.fixed_alpha = 0.7;
    user_data data;
    auto models = build_models(params, data);
    LONGS_EQUAL(1, models.size());
    CHECK(dynamic_cast<gamma_model *>(models[0]));
    for (auto m : models)
        delete m;
}

TEST(Inference, create_gamma_model_if__n_gamma_cats__provided)
{
    input_parameters params;
    params.input_file_path = "foo";
    params.n_gamma_cats = 3;
    user_data data;
    auto models = build_models(params, data);
    LONGS_EQUAL(1, models.size());
    CHECK(dynamic_cast<gamma_model *>(models[0]));
    for (auto m : models)
        delete m;
}

void build_matrix(matrix& m)
{
    m.set(0, 0, 1);
    m.set(0, 1, 2);
    m.set(0, 2, 3);
    m.set(1, 0, 4);
    m.set(1, 1, 5);
    m.set(1, 2, 6);
    m.set(2, 0, 7);
    m.set(2, 1, 8);
    m.set(2, 2, 9);
}

TEST(Probability, matrix_multiply)
{
    matrix m1(3);
    build_matrix(m1);
    vector<double> m2({ 7, 9, 11 });
    auto result = m1.multiply(m2, 0, 2, 0, 2);
    LONGS_EQUAL(3, result.size());

    DOUBLES_EQUAL(58, result[0], .001);
    DOUBLES_EQUAL(139, result[1], .001);
    DOUBLES_EQUAL(220, result[2], .001);

    matrix m3(8);
    m3.set(3, 3, 1);
    m3.set(3, 4, 2);
    m3.set(3, 5, 3);
    m3.set(4, 3, 4);
    m3.set(4, 4, 5);
    m3.set(4, 5, 6);
    m3.set(5, 3, 7);
    m3.set(5, 4, 8);
    m3.set(5, 5, 9);

    result = m3.multiply(m2, 3, 5, 3, 5);
    LONGS_EQUAL(3, result.size());

    DOUBLES_EQUAL(58, result[0], .001);
    DOUBLES_EQUAL(139, result[1], .001);
    DOUBLES_EQUAL(220, result[2], .001);
}

TEST(Probability, error_model_set_probs)
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    auto vec = model.get_probs(0);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.0, vec[0], 0.00001);
    DOUBLES_EQUAL(0.7, vec[1], 0.00001);
    DOUBLES_EQUAL(0.3, vec[2], 0.00001);

    vec = model.get_probs(1);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.2, vec[0], 0.00001);
    DOUBLES_EQUAL(0.6, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);
}

TEST(Probability, error_model_get_epsilon)
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(2, { .1, .8, .1 });
    model.set_probabilities(3, { .2, .6, .2 });

    auto actual = model.get_epsilons();
    LONGS_EQUAL(3, actual.size());
    vector<double> expected{ .1, .2, .3 };
    CHECK(expected == actual);
}

TEST(Probability, error_model_get_epsilon_zero_zero_must_be_zero)
{
    error_model model;
    try
    {
        model.set_probabilities(0, { 0.4, 0.3, 0.3 });
        CHECK(false);
    }
    catch(runtime_error& err)
    {
        STRCMP_EQUAL("Cannot have a non-zero probability for family size 0 for negative deviation", err.what());
    }
}

TEST(Probability, error_model_rows_must_add_to_one)
{
    error_model model;
    try
    {
        model.set_probabilities(1, { 0.3, 0.3, 0.3 });
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Sum of probabilities must be equal to one", err.what());
    }
}

TEST(Probability, error_model_replace_epsilons)
{
    string input = "maxcnt: 10\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "2 0.2 0.6 0.2\n";

    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    map<double, double> replacements;
    replacements[.2] = .3;
    model.replace_epsilons(&replacements);

    auto actual = model.get_probs(0);
    DOUBLES_EQUAL(.7, actual[1], 0.0001);
    DOUBLES_EQUAL(.3, actual[2], 0.0001);

    actual = model.get_probs(1);
    DOUBLES_EQUAL(.3, actual[0], 0.0001);
    DOUBLES_EQUAL(.4, actual[1], 0.0001);
    DOUBLES_EQUAL(.3, actual[2], 0.0001);
}

TEST(Probability, read_error_model)
{
    string input = "maxcnt: 10\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "2 0.2 0.6 0.2\n"
        "3 0.2 0.6 0.2\n"
        "5 0.2 0.6 0.2\n";

    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);
    auto vec = model.get_probs(0);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.0, vec[0], 0.00001);
    DOUBLES_EQUAL(0.8, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);

    vec = model.get_probs(1);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.2, vec[0], 0.00001);
    DOUBLES_EQUAL(0.6, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);


    vec = model.get_probs(4);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.2, vec[0], 0.00001);
    DOUBLES_EQUAL(0.6, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);
}

TEST(Probability, matrix_cache_warns_on_saturation)
{
    matrix_cache m(10);
    m.precalculate_matrices({ 0.05, 0.01 }, { 25 });
    ostringstream ost;
    m.warn_on_saturation(ost);
    STRCMP_EQUAL("WARNING: Saturated branch using lambda 0.05 on branch length 25\n", ost.str().c_str());
}

TEST(Probability, matrix_is_saturated)
{
    matrix_cache c(10);
    CHECK(c.is_saturated(25, 0.05));
    CHECK_FALSE(c.is_saturated(25, 0.01));
}

TEST(Inference, build_reference_list)
{
    std::string str = "Desc\tFamily ID\tA\tB\n"
        "\t (null)1\t5\t10\n"
        "\t (null)2\t5\t7\n"
        "\t (null)3\t5\t10\n"
        "\t (null)4\t5\t7\n";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);
    auto actual = build_reference_list(families);
    vector<int> expected({ 0, 1, 0, 1 });
    LONGS_EQUAL(expected.size(), actual.size());

}

TEST(Inference, prune)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);
    inference_process process(ost, &lambda, 1.5, p_tree.get(), 20, 20, &fam, NULL);
    matrix_cache cache(21);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });
    auto actual = process.prune(cache);

    vector<double> log_expected{ -17.2771, -10.0323 , -5.0695 , -4.91426 , -5.86062 , -7.75163 , -10.7347 , -14.2334 , -18.0458 , 
        -22.073 , -26.2579 , -30.5639 , -34.9663 , -39.4472 , -43.9935 , -48.595 , -53.2439 , -57.9338 , -62.6597 , -67.4173 };

    LONGS_EQUAL(log_expected.size(), actual.size());
    for (size_t i = 0; i<log_expected.size(); ++i)
    {
        DOUBLES_EQUAL(log_expected[i], log(actual[i]), 0.0001);
    }
}

TEST(Inference, likelihood_computer_sets_leaf_nodes_correctly)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    likelihood_computer pruner(20, 20, &lambda, family, cache, NULL);
    pruner.initialize_memory(p_tree.get());
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });

    auto A = p_tree->find_descendant("A");
    pruner(A);
    auto actual = pruner.get_likelihoods(A);

    vector<double> expected{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i<expected.size(); ++i)
    {
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }

    auto B = p_tree->find_descendant("B");
    pruner(B);
    actual = pruner.get_likelihoods(B);

    expected = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i<expected.size(); ++i)
    {
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }
}

TEST(Inference, likelihood_computer_sets_root_nodes_correctly)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    likelihood_computer pruner(20, 20, &lambda, family, cache, NULL);
    pruner.initialize_memory(p_tree.get());

    cache.precalculate_matrices({ 0.03 }, { 1.0,3.0,7.0 });

    auto AB = p_tree->find_descendant("AB");
    pruner(p_tree->find_descendant("A"));
    pruner(p_tree->find_descendant("B"));
    pruner(AB);

    auto actual = pruner.get_likelihoods(AB);

    vector<double> log_expected{ -19.7743, -11.6688, -5.85672, -5.66748, -6.61256, -8.59725, -12.2301, -16.4424, -20.9882, -25.7574, 
        -30.6888, -35.7439, -40.8971, -46.1299, -51.4289, -56.7837, -62.1863, -67.6304, -73.1106, -78.6228
    };

    LONGS_EQUAL(log_expected.size(), actual.size());
    for (size_t i = 0; i<log_expected.size(); ++i)
    {
        DOUBLES_EQUAL(log_expected[i], log(actual[i]), 0.0001);
    }
}

TEST(Inference, likelihood_computer_sets_leaf_nodes_from_error_model_if_provided)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";

    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });

    string input = "maxcnt: 20\ncntdiff: -1 0 1\n"
        "1 0.2 0.6 0.2\n"
        "20 0.2 0.6 0.2\n";
    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    likelihood_computer pruner(20, 20, &lambda, family, cache, &model);
    pruner.initialize_memory(p_tree.get());

    auto A = p_tree->find_descendant("A");
    pruner(A);
    auto actual = pruner.get_likelihoods(A);

    vector<double> expected{ 0, 0, 0.2, 0.6, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        //cout << actual[i] << endl;
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }
}

TEST(Clade, get_lambda_index_throws_from_branch_length_tree)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    LONGS_EQUAL(7, p_tree->get_branch_length());
    try
    {
        p_tree->get_lambda_index();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Requested lambda index from branch length tree", err.what());
    }

}

TEST(Clade, get_branch_length_throws_from_lambda_tree)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    parser.parse_to_lambdas = true;
    unique_ptr<clade> p_tree(parser.parse_newick());
    LONGS_EQUAL(7, p_tree->get_lambda_index());

    try
    {
        p_tree->get_branch_length();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Requested branch length from lambda tree", err.what());
    }

}

TEST(Clade, exists_at_root_returns_false_if_not_all_children_exist)
{
    newick_parser parser(false);
    parser.newick_string= " ((((cat:68.710687,horse:68.710687):4.566771,cow:73.277458):20.722542,(((((chimp:4.444178,human:4.444178):6.682660,orang:11.126837):2.285866,gibbon:13.412704):7.211528,(macaque:4.567239,baboon:4.567239):16.056993):16.060691,marmoset:36.684923):57.315077)mancat:38.738115,(rat:36.302467,mouse:36.302467):96.435648)";
    unique_ptr<clade> p_tree(parser.parse_newick());

    istringstream ist(
    "Desc\tFamily ID\tcat\thorse\tcow\tchimp\thuman\torang\tgibbon\tmacaque\tbaboon\tmarmoset\trat\tmouse\n"
 "(null)\t1\t0\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0\n"
 "(null)\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\n");

    vector<gene_family> families;
    read_gene_families(ist, p_tree.get(), &families);
    CHECK_FALSE(families[0].exists_at_root(p_tree.get()));
    CHECK_FALSE(families[1].exists_at_root(p_tree.get()));
}

TEST(Clade, exists_at_root_returns_true_if_all_children_exist)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    CHECK(family.exists_at_root(p_tree.get()));
}

TEST(Inference, multiple_lambda_returns_correct_values)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    map<string, int> key;
    key["A"] = 5;
    key["B"] = 3;
    multiple_lambda ml(key, { .03, .05, .07, .011, .013, .017 });
    DOUBLES_EQUAL(.017, ml.get_value_for_clade(p_tree->find_descendant("A")), 0.0001);
    DOUBLES_EQUAL(.011, ml.get_value_for_clade(p_tree->find_descendant("B")), 0.0001);
}

class mock_optimizer : public inference_optimizer_scorer
{
    double _initial;
    // Inherited via optimizer
    virtual std::vector<double> initial_guesses() override
    {
        return std::vector<double>{_initial};
    }
    virtual void prepare_calculation(double * values) override
    {
    }
    virtual void report_precalculation() override
    {
    }

public:
    mock_optimizer(double initial) : 
        inference_optimizer_scorer(NULL, NULL, NULL, std::map<int, int>()),
        _initial(initial)
    {

    }

};

class mock_model : public model {
    // Inherited via model
    virtual void start_inference_processes(lambda *) override
    {
    }
    virtual std::string name() override
    {
        return "mockmodel";
    }
    virtual void write_family_likelihoods(std::ostream & ost) override
    {
    }
    virtual reconstruction* reconstruct_ancestral_states(matrix_cache * p_calc, root_equilibrium_distribution * p_prior) override
    {
        return nullptr;
    }
    virtual inference_optimizer_scorer * get_lambda_optimizer(user_data& data) override
    {
        branch_length_finder finder;
        _p_tree->apply_prefix_order(finder);

        initialize_lambda(data.p_lambda_tree);
        auto result = new lambda_optimizer(_p_lambda, this, data.p_prior.get(), finder.longest(), std::map<int, int>());
        result->quiet = true;
        return result;
    }
public:
    mock_model() : model(NULL, NULL, NULL, 0, 0, NULL)
    {

    }
    void set_lambda(lambda * lambda)
    {
        _p_lambda = lambda;
    }
    void set_tree(clade * tree)
    {
        _p_tree = tree;
    }

    // Inherited via model
    virtual void prepare_matrices_for_simulation(matrix_cache & cache) override
    {
        branch_length_finder lengths;
        _p_tree->apply_prefix_order(lengths);
        cache.precalculate_matrices(get_lambda_values(_p_lambda), lengths.result());
    }

    // Inherited via model
    virtual double infer_processes(root_equilibrium_distribution * prior, const std::map<int, int>& root_distribution_map) override
    {
        return 0.0;
    }
};

TEST(Simulation, select_root_size_returns_less_than_100_without_rootdist)
{
    root_distribution rd;
    rd.vectorize_increasing(100);

    user_data data;

    for (int i = 0; i<50; ++i)
    {
        int root_size = select_root_size(data, rd, 0);
        CHECK(root_size < 100);
    }
}

TEST(Simulation, print_process_prints_in_order)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    std::ostringstream ost;
    trial t;
    t[p_tree->find_descendant("B")] = 4;
    t[p_tree->find_descendant("A")] = 2;
    t[p_tree->find_descendant("AB")] = 6;

    vector<trial*> my_trials({ &t });

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, true, my_trials);

    STRCMP_CONTAINS("DESC\tFID\tB\tA\t2", ost.str().c_str());
    STRCMP_CONTAINS("NULL\tsimfam0\t4\t2\t6", ost.str().c_str());

}

TEST(Simulation, print_process_can_print_without_internal_nodes)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    std::ostringstream ost;
    trial t;
    t[p_tree->find_descendant("B")] = 4;
    t[p_tree->find_descendant("A")] = 2;
    t[p_tree->find_descendant("AB")] = 6;

    vector<trial*> my_trials({ &t });

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, false, my_trials);
    STRCMP_CONTAINS("DESC\tFID\tB\tA\n", ost.str().c_str());
    STRCMP_CONTAINS("NULL\tsimfam0\t4\t2\n", ost.str().c_str());

}

TEST(Simulation, gamma_model_get_simulation_lambda_selects_random_multiplier_based_on_alpha)
{
    gamma_model m(NULL, NULL, NULL, 0, 5, 3, 0.7, NULL);
    user_data data;
    single_lambda lam(0.05);
    data.p_lambda = &lam;
    unique_ptr<single_lambda> new_lam(dynamic_cast<single_lambda *>(m.get_simulation_lambda(data)));
    DOUBLES_EQUAL(0.00574028, new_lam->get_single_lambda(), 0.0000001);
}

TEST(Inference, model_vitals)
{
    mock_model model;
    single_lambda lambda(0.05);
    model.set_lambda(&lambda);
    std::ostringstream ost;
    model.write_vital_statistics(ost, 0.01);
    STRCMP_EQUAL("Model mockmodel Result: 0.01\nLambda:            0.05\n", ost.str().c_str());
}

class increasable {
    clade *_p_tree;
public:
    increasable(clade *p_tree) : _p_tree(p_tree)
    {

    }
    std::vector<const clade *> get_taxa()
    {
        vector<string> taxa{ "A", "B", "AB" };
        vector<const clade*> result(3);
        transform(taxa.begin(), taxa.end(), result.begin(), [this](string taxon)->const clade * {
            return _p_tree->find_descendant(taxon);});
        return result;
    }
    increase_decrease get_increases_decreases(std::vector<const clade *>& order, double pvalue)
    {
        increase_decrease result;
        result.gene_family_id = "myid";
        result.pvalue = pvalue;
        result.change = vector<family_size_change>{ Increase, Decrease, Constant };
        return result;
    }
};


TEST(Inference, print_increases_decreases_by_family)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    vector<increasable*> vec;
    vector<double> pvalues;
    ostringstream empty;
    ::print_increases_decreases_by_family(empty, vec, pvalues);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    unique_ptr<increasable> ic(new increasable(p_tree.get()));
    vec.push_back(ic.get());
    pvalues.push_back(0.07);
    ostringstream ost;
    ::print_increases_decreases_by_family(ost, vec, pvalues);
    STRCMP_CONTAINS("#FamilyID\tpvalue\t*\tA\tB\tAB", ost.str().c_str());
    STRCMP_CONTAINS("myid\t0.07\tn\ti\td\tc", ost.str().c_str());
}

TEST(Inference, print_increases_decreases_by_clade)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    vector<increasable*> vec;
    ostringstream empty;
    ::print_increases_decreases_by_clade(empty, vec);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    unique_ptr<increasable> ic(new increasable(p_tree.get()));
    vec.push_back(ic.get());
    ostringstream ost;
    ::print_increases_decreases_by_clade(ost, vec);
    STRCMP_CONTAINS("#Taxon_ID\tIncrease/Decrease", ost.str().c_str());
    STRCMP_CONTAINS("A\t1/0", ost.str().c_str());
    STRCMP_CONTAINS("B\t0/1", ost.str().c_str());
}

TEST(Inference, lambda_epsilon_optimizer)
{
    const double initial_epsilon = 0.01;
    error_model err;
    err.set_probabilities(0, { .0, .99, initial_epsilon });
    err.set_probabilities(1, { initial_epsilon, .98, initial_epsilon });

    mock_model model;

    single_lambda lambda(0.05);
    lambda_epsilon_optimizer optimizer(&model, &err, NULL, std::map<int, int>(), &lambda, 10);
    optimizer.initial_guesses();
    vector<double> values = { 0.05, 0.06 };
    optimizer.calculate_score(&values[0]);
    auto actual = err.get_probs(0);
    vector<double> expected{ 0, .94, .06 };
    CHECK(expected == actual);

    values[1] = 0.04;
    optimizer.calculate_score(&values[0]);
    actual = err.get_probs(0);
    expected = { 0, .96, .04 };
    CHECK(expected == actual);
}

TEST(Inference, gamma_optimizer)
{
    gamma_model m(NULL, NULL, NULL, 0, 0, 0, 0, NULL);
    gamma_optimizer optimizer(&m, NULL, std::map<int, int>());
    auto initial = optimizer.initial_guesses();
    DOUBLES_EQUAL(0.979494, initial[0], 0.00001);
}

TEST(GeneFamilies, model_set_families)
{
    mock_model m;

    vector<gene_family> fams;
    m.set_families(&fams);
    LONGS_EQUAL(0, m.get_gene_family_count());

    fams.resize(5);
    LONGS_EQUAL(5, m.get_gene_family_count());
}

TEST(Inference, lambda_per_family)
{
    user_data ud;
    
    ud.max_root_family_size = 10;
    ud.max_family_size = 10;
    ud.gene_families.resize(1);

    gene_family& family = ud.gene_families[0];
    family.set_id("test");
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);
    input_parameters params;
    params.lambda_per_family = true;

    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());
    ud.p_tree = p_tree.get();
    estimator v(ud, params);

    mock_model m;
    m.set_tree(ud.p_tree);
    ostringstream ost;
    v.estimate_lambda_per_family(&m, ost);
    STRCMP_EQUAL("test\t0.042680165523241\n", ost.str().c_str());
}

TEST(Inference, estimator_compute_pvalues)
{
    input_parameters params;

    estimator v(_user_data, params);
    auto values = v.compute_pvalues(_user_data, 3);
    LONGS_EQUAL(1, values.size());
    DOUBLES_EQUAL(0.666667, values[0], 0.00001);
}

TEST(Inference, gamma_lambda_optimizer)
{
    uniform_distribution frq;

    gamma_model m(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 10, 10, 4, 0.25, NULL);
    gamma_lambda_optimizer optimizer(_user_data.p_lambda, &m, &frq, std::map<int, int>(), 7);
    vector<double> values{ 0.01, 0.25 };
    DOUBLES_EQUAL(6.4168, optimizer.calculate_score(&values[0]), 0.0001);
}

TEST(Inference, poisson_scorer_optimizes_correct_value)
{
    poisson_scorer scorer(_user_data.gene_families);
    optimizer opt(&scorer);

    auto result = opt.optimize();

    DOUBLES_EQUAL(0.5, result.values[0], 0.0001)
}

class mock_scorer : public optimizer_scorer
{
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        return std::vector<double>{0.2};
    }
    virtual double calculate_score(double * values) override
    {
        if (force_scoring_error)
            return std::numeric_limits<double>::infinity();

        return 1000;
    }
public:
    bool force_scoring_error = false;
};

TEST(Inference, optimizer_gets_initial_guesses_from_scorer)
{
    mock_scorer scorer;
    optimizer opt(&scorer);
    auto guesses = opt.get_initial_guesses();
    LONGS_EQUAL(1, guesses.size())
    DOUBLES_EQUAL(0.2, guesses[0], 0.0001);
}

TEST(Inference, optimizer_disallows_bad_initializations)
{
    mock_scorer scorer;
    scorer.force_scoring_error = true;
    optimizer opt(&scorer);

    try
    {
        opt.get_initial_guesses();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Failed to find any reasonable values", err.what());
    }

}

TEST(Inference, root_distribution_copy)
{
    root_distribution rd;
    rd.vector({ 1,2,5 });
    LONGS_EQUAL(5, rd.max());

    root_distribution rd2 = rd;
    LONGS_EQUAL(3, rd2.size());
    LONGS_EQUAL(5, rd2.max());
}

TEST(Inference, uniform_root_distribution_does_not_hold_rootdist_pointer)
{
    uniform_distribution prior;
    {
        root_distribution rd;
        rd.vector({ 1,2,5 });
        prior.initialize(&rd);
    }
    DOUBLES_EQUAL(.25, prior.compute(1), 0.00001);
}

TEST(Inference, optimizer_result_stream)
{
    optimizer::result r;
    r.num_iterations = 10;
    r.score = 5;
    r.values = vector<double>{ .05, .03 };
    r.duration = chrono::seconds(5000);

    ostringstream ost;
    ost << r;
    STRCMP_CONTAINS("Completed 10 iterations", ost.str().c_str());
    STRCMP_CONTAINS("Time: 1H 23M 20S", ost.str().c_str());
    STRCMP_CONTAINS("Best matches are:            0.05,0.03", ost.str().c_str());
    STRCMP_CONTAINS("Final -lnL: 5", ost.str().c_str());
}

TEST(Simulation, base_prepare_matrices_for_simulation_creates_matrix_for_each_branch)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parser.parse_newick());
    base_model b(&lam, p_tree.get(), NULL, 0, 0, NULL);
    matrix_cache m(25);
    b.prepare_matrices_for_simulation(m);
    LONGS_EQUAL(3, m.get_cache_size());
}

TEST(Simulation, base_prepare_matrices_for_simulation_uses_perturbed_lambda)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    user_data data;
    single_lambda lam(0.05);
    data.p_lambda = &lam;
    unique_ptr<clade> p_tree(parser.parse_newick());
    base_model b(&lam, p_tree.get(), NULL, 0, 0, NULL);
    b.perturb_lambda();
    matrix_cache m(25);
    b.prepare_matrices_for_simulation(m);
    unique_ptr<single_lambda> sim_lambda(dynamic_cast<single_lambda *>(b.get_simulation_lambda(data)));
    try
    {
        m.get_matrix(7, sim_lambda->get_single_lambda());
    }
    catch (std::runtime_error& err)
    {
        FAIL("Failed to cache simulation lambda");
    }
}

TEST(Simulation, gamma_prepare_matrices_for_simulation_creates_matrix_for_each_branch_and_category)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parser.parse_newick());
    gamma_model g(&lam, p_tree.get(), NULL, 0, 0, 2, 0.5, NULL);
    matrix_cache m(25);
    g.prepare_matrices_for_simulation(m);
    LONGS_EQUAL(6, m.get_cache_size());
}

TEST(Simulation, set_random_node_size_without_error_model)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.05);
    trial t;
    matrix_cache cache(12);
    cache.precalculate_matrices({ 0.05 }, { 3 });
    auto b = p_tree->find_descendant("B");
    set_weighted_random_family_size(b, &t, &lambda, NULL, 10, cache);

    LONGS_EQUAL(0, t[b]);

    t[p_tree.get()] = 5;

    set_weighted_random_family_size(b, &t, &lambda, NULL, 10, cache);
    LONGS_EQUAL(4, t[b]);
}

TEST(Simulation, set_random_node_size_with_error_model)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.05);
    trial t;
    matrix_cache cache(20);
    cache.precalculate_matrices({ 0.05 }, { 3 });
    auto b = p_tree->find_descendant("B");

    error_model err;

    t[p_tree.get()] = 5;
    err.set_probabilities(5, { .9, .05, 0.05 });
    set_weighted_random_family_size(b, &t, &lambda, &err, 10, cache);
    LONGS_EQUAL(4, t[b]);

    err.set_probabilities(5, { .05, .05, 0.9 });
    set_weighted_random_family_size(b, &t, &lambda, &err, 10, cache);
    LONGS_EQUAL(6, t[b]);
}

TEST(Simulation, executor)
{
    input_parameters params;
    user_data ud;
    unique_ptr<action> act(get_executor(params, ud));
    CHECK(dynamic_cast<const estimator *>(act.get()))

    params.chisquare_compare = true;
    unique_ptr<action> act2(get_executor(params, ud));
    CHECK(dynamic_cast<const chisquare_compare *>(act2.get()))
}

TEST(Simulation, rootdist_vectorize_creates_matching_vector)
{
    root_distribution rd;
    std::map<int, int> m;
    m[2] = 3;
    rd.vectorize(m);
    CHECK(rd.size() == 3);
    CHECK(rd.at(0) == 2);
    CHECK(rd.at(1) == 2);
    CHECK(rd.at(2) == 2);
}

TEST(Simulation, rootdist_vectorize_uniform_creates_uniform_vector)
{
    root_distribution rd;
    rd.vectorize_uniform(5);
    CHECK(rd.size() == 5);
    CHECK(rd.at(0) == 1);
    CHECK(rd.at(1) == 1);
    CHECK(rd.at(2) == 1);
    CHECK(rd.at(3) == 1);
    CHECK(rd.at(4) == 1);
}

TEST(Simulation, rootdist_vectorize_increasing_creates_increasing_vector)
{
    root_distribution rd;
    rd.vectorize_increasing(5);
    CHECK(rd.size() == 5);
    CHECK(rd.at(0) == 0);
    CHECK(rd.at(1) == 1);
    CHECK(rd.at(2) == 2);
    CHECK(rd.at(3) == 3);
    CHECK(rd.at(4) == 4);
}

TEST(Simulation, rootdist_max_and_sum)
{
    root_distribution rd;
    std::map<int, int> m;

    m[2] = 3;
    m[4] = 1;
    m[8] = 1;
    rd.vectorize(m);
    CHECK(rd.max() == 8);
    CHECK(rd.sum() == 18);  // three twos, one four and one eight
}

TEST(Simulation, simulate_processes)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parser.parse_newick());
    mock_model m;
    m.set_tree(p_tree.get());
    m.set_lambda(&lam);

    user_data ud;
    ud.p_tree = p_tree.get();
    ud.p_lambda = &lam;
    input_parameters ip;
    ip.nsims = 100;
    simulator sim(ud, ip);
    vector<trial *> results;
    sim.simulate_processes(&m, results);
    LONGS_EQUAL(100, results.size());

    for (auto r : results) delete r;
}

TEST(Simulation, simulate_processes_uses_rootdist_if_available)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parser.parse_newick());
    mock_model m;
    m.set_tree(p_tree.get());
    m.set_lambda(&lam);

    user_data ud;
    ud.p_tree = p_tree.get();
    ud.p_lambda = &lam;
    ud.rootdist[5] = 50;
    ud.rootdist[10] = 50;
    input_parameters ip;
    ip.nsims = 100;
    simulator sim(ud, ip);
    vector<trial *> results;
    sim.simulate_processes(&m, results);
    LONGS_EQUAL(100, results.size());

    LONGS_EQUAL(5, results[0]->at(p_tree.get()));
    LONGS_EQUAL(10, results[75]->at(p_tree.get()));
    for (auto r : results) delete r;
}

TEST(Simulation, root_distribution_pare)
{
    root_distribution rd;
    rd.vectorize_increasing(100);
    LONGS_EQUAL(100, rd.size());

    rd.pare(10);
    LONGS_EQUAL(10, rd.size());

    for (int i = 1; i < 10; ++i)
    {
        CHECK(rd.at(i-1) <= rd.at(i));
    }

}

TEST(Simulation, gamma_model_perturb_lambda_with_clusters)
{
    gamma_model model(NULL, NULL, NULL, 0, 5, 3, 0.7, NULL);
    model.perturb_lambda();
    auto multipliers = model.get_lambda_multipliers();
    LONGS_EQUAL(3, multipliers.size());
    DOUBLES_EQUAL(0.112844, multipliers[0], 0.00001);
    DOUBLES_EQUAL(1.05493, multipliers[1], 0.00001);
    DOUBLES_EQUAL(2.06751, multipliers[2], 0.00001);
}

TEST(Simulation, gamma_model_perturb_lambda_without_clusters)
{
    gamma_model model(NULL, NULL, NULL, 0, 5, 1, 5, NULL);
    model.perturb_lambda();
    auto multipliers = model.get_lambda_multipliers();
    LONGS_EQUAL(1, multipliers.size());
    DOUBLES_EQUAL(0.911359, multipliers[0], 0.00001);
}

TEST(Optimizer, fminsearch_sort_sorts_scores_and_moves_values)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.variable_count_plus_one = 3;
    vector<double> scores({ 3.0, 5.0, 1.0 });
    fm.fv = &scores[0];
    vector<int> indices(3);
    fm.idx = &indices[0];
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.vsort = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.v[0][0] = 300;
    fm.v[1][0] = 200;
    fm.v[2][0] = 100;
    __fminsearch_sort(&fm);

    DOUBLES_EQUAL(1.0, fm.fv[0], 0.00001);
    DOUBLES_EQUAL(3.0, fm.fv[1], 0.00001);
    DOUBLES_EQUAL(5.0, fm.fv[2], 0.00001);

    DOUBLES_EQUAL(100.0, fm.v[0][0], 0.00001);
    DOUBLES_EQUAL(300.0, fm.v[1][0], 0.00001);
    DOUBLES_EQUAL(200.0, fm.v[2][0], 0.00001);
    free_2dim((void **)fm.v, 3, 2);
    free_2dim((void **)fm.vsort, 3, 2);
}

TEST(Optimizer, fminsearch_checkV_compares_value_difference_to_lx)
{
    FMinSearch fm;
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.variable_count = 2;
    fm.v[0][0] = 1;
    fm.v[1][0] = 2;
    fm.v[2][0] = 3;
    fm.v[0][1] = 3;
    fm.v[1][1] = 4;
    fm.v[2][1] = 5;
    fm.tolx = 3;
    CHECK(__fminsearch_checkV(&fm));
    fm.tolx = .5;
    CHECK_FALSE(__fminsearch_checkV(&fm));
    free_2dim((void **)fm.v, 3, 2);
}

TEST(Optimizer, fminsearch_checkF_compares_score_difference_to_lf)
{
    FMinSearch fm;
    fm.variable_count_plus_one = 3;
    vector<double> scores({ 1.0, 3.0, 5.0 });
    fm.fv = &scores[0];
    fm.tolf = 5;
    CHECK(__fminsearch_checkF(&fm));
    fm.tolf = 1;
    CHECK_FALSE(__fminsearch_checkF(&fm));
}

class multiplier_scorer : public optimizer_scorer
{
public:
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        return std::vector<double>{5, 3};
    }
    virtual double calculate_score(double * values) override
    {
        return values[0] * values[1];
    }
};

TEST(Optimizer, fminsearch_min_init)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.variable_count_plus_one = 3;
    fm.delta = 0.05;
    fm.zero_delta = 0.00025;
    vector<double> scores(3);
    fm.fv = &scores[0];
    vector<int> indices(3);
    fm.idx = &indices[0];
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.vsort = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.v[0][0] = 300;
    fm.v[1][0] = 200;
    fm.v[2][0] = 100;

    multiplier_scorer ms;
    fm.scorer = &ms;
    auto init = ms.initial_guesses();
    __fminsearch_min_init(&fm, &init[0]);

    DOUBLES_EQUAL(15, fm.fv[0], 0.0001);
    DOUBLES_EQUAL(15.75, fm.fv[1], 0.0001);
    DOUBLES_EQUAL(15.75, fm.fv[2], 0.0001);

    free_2dim((void **)fm.v, 3, 2);
    free_2dim((void **)fm.vsort, 3, 2);
}

TEST(Optimizer, __fminsearch_x_mean)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    vector<double> means(2);
    fm.x_mean = &means[0];
    fm.v[0][0] = 300;
    fm.v[1][0] = 200;
    fm.v[0][1] = 12;
    fm.v[1][1] = 44;

    __fminsearch_x_mean(&fm);

    DOUBLES_EQUAL(250, fm.x_mean[0], 0.0001);
    DOUBLES_EQUAL(28, fm.x_mean[1], 0.0001);

    free_2dim((void **)fm.v, 3, 2);
}

TEST(Optimizer, __fminsearch_x_reflection)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.rho = 1;
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    vector<double> means({ 250,28 });
    vector<double> reflections(2);
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.v[0][0] = 300;
    fm.v[1][0] = 200;
    fm.v[0][1] = 12;
    fm.v[1][1] = 44;

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_reflection(&fm);

    DOUBLES_EQUAL(500, fm.x_r[0], 0.0001);
    DOUBLES_EQUAL(56, fm.x_r[1], 0.0001);
    DOUBLES_EQUAL(28000, score, 0.0001);

    free_2dim((void **)fm.v, 3, 2);

}

TEST(Optimizer, __fminsearch_x_expansion)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.chi = 2;
    vector<double> means({ 250,28 });
    vector<double> expansions(2);
    vector<double> reflections({ 500, 56 });
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_expansion(&fm);

    DOUBLES_EQUAL(750, fm.x_tmp[0], 0.0001);
    DOUBLES_EQUAL(84, fm.x_tmp[1], 0.0001);
    DOUBLES_EQUAL(63000, score, 0.0001);
}

TEST(Optimizer, __fminsearch_x_contract_outside)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.psi = 0.5;
    vector<double> means({ 250,28 });
    vector<double> reflections({ 500, 56 });
    vector<double> expansions(2);
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_contract_outside(&fm);

    DOUBLES_EQUAL(375, fm.x_tmp[0], 0.0001);
    DOUBLES_EQUAL(42, fm.x_tmp[1], 0.0001);
    DOUBLES_EQUAL(15750, score, 0.0001);
}


TEST(Optimizer, __fminsearch_x_contract_inside)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.psi = 0.5;
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    vector<double> means({ 250,28 });
    vector<double> expansions(2);
    fm.v[2][0] = 26;
    fm.v[2][1] = 12;

    fm.x_mean = &means[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_contract_inside(&fm);

    DOUBLES_EQUAL(362, fm.x_tmp[0], 0.0001);
    DOUBLES_EQUAL(36, fm.x_tmp[1], 0.0001);
    DOUBLES_EQUAL(13032, score, 0.0001);
    free_2dim((void **)fm.v, 3, 2);


}

TEST(Optimizer, __fminsearch_x_shrink)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.variable_count_plus_one = 3;
    fm.sigma = 0.5;
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.vsort = (double**)calloc_2dim(3, 2, sizeof(double));
    vector<int> indices(3);
    fm.idx = &indices[0];
    vector<double> scores(3);
    fm.fv = &scores[0];

    fm.v[0][0] = 300;
    fm.v[0][1] = 200;
    fm.v[1][0] = 42;
    fm.v[1][1] = 64;
    fm.v[2][0] = 26;
    fm.v[2][1] = 12;

    multiplier_scorer ms;
    fm.scorer = &ms;
    __fminsearch_x_shrink(&fm);

    DOUBLES_EQUAL(300, fm.v[0][0], 0.0001);
    DOUBLES_EQUAL(200, fm.v[0][1], 0.0001);
    DOUBLES_EQUAL(163, fm.v[1][0], 0.0001);
    DOUBLES_EQUAL(106, fm.v[1][1], 0.0001);
    DOUBLES_EQUAL(171, fm.v[2][0], 0.0001);
    DOUBLES_EQUAL(132, fm.v[2][1], 0.0001);

    free_2dim((void **)fm.v, 3, 2);
    free_2dim((void **)fm.vsort, 3, 2);
}

TEST(Optimizer, __fminsearch_set_last_element)
{
    FMinSearch fm;
    fm.variable_count = 2;
    fm.variable_count_plus_one = 3;
    fm.v = (double**)calloc_2dim(3, 2, sizeof(double));
    fm.vsort = (double**)calloc_2dim(3, 2, sizeof(double));
    vector<int> indices(3);
    fm.idx = &indices[0];
    vector<double> scores({2,4,6});
    fm.fv = &scores[0];

    fm.v[0][0] = 300;
    fm.v[0][1] = 200;
    fm.v[1][0] = 42;
    fm.v[1][1] = 64;
    fm.v[2][0] = 26;
    fm.v[2][1] = 12;

    vector<double> new_vals({ 99, 14 });
    __fminsearch_set_last_element(&fm, &new_vals[0], 3);

    DOUBLES_EQUAL(2.0, fm.fv[0], 0.00001);
    DOUBLES_EQUAL(3.0, fm.fv[1], 0.00001);
    DOUBLES_EQUAL(4.0, fm.fv[2], 0.00001);

    free_2dim((void **)fm.v, 3, 2);
    free_2dim((void **)fm.vsort, 3, 2);

}


void init_lgamma_cache();

int main(int ac, char** av)
{
    init_lgamma_cache();
    return CommandLineTestRunner::RunAllTests(ac, av);
}