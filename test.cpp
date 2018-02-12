
#include <numeric>
#include <cmath>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"
#include "src/root_equilibrium_distribution.h"
#include "src/base_model.h"
#include "src/reconstruction_process.h"
#include "src/matrix_cache.h"
#include "src/gamma_bundle.h"
#include "src/probability.h"

TEST_GROUP(GeneFamilies)
{
};

TEST_GROUP(Inference)
{
    void setup()
    {
        srand(10);
    }
};

TEST_GROUP(Simulation)
{
};

TEST_GROUP(Probability)
{
};

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
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:1);";
    clade *p_tree = parser.parse_newick();

    base_model core(&lambda, p_tree, &families, 56, 30, NULL, NULL);
    core.start_inference_processes();


    uniform_distribution frq;
    double multi = core.infer_processes(&frq);
    //core.get_likelihoods();
    DOUBLES_EQUAL(39.9207, multi, 0.001);
    delete p_tree;
}

TEST(Inference, uniform_distribution)
{
    std::vector<int> rd(10);
    std::fill(rd.begin(), rd.end(), 1);
    uniform_distribution ef;
    ef.initialize(rd);
    DOUBLES_EQUAL(.1, ef.compute(5), 0.0001);
}

TEST(Inference, gamma_set_alpha)
{
    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL, NULL);
    model.set_alpha(0.5, 3);
}

TEST(Inference, gamma_adjust_family_gamma_membership)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);

    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL, NULL);
}

TEST(Inference, gamma)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6\n";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);

    single_lambda lambda(0.5);
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    std::vector<int> rootdist;
    gamma_model core(NULL, NULL, NULL, 0, 5, 0, 0, NULL, NULL);
    core.set_gene_families(&families);
    core.set_lambda(&lambda);
    core.set_tree(p_tree);
    core.set_max_sizes(148, 122);

    core.start_inference_processes();
    uniform_distribution frq;
    double actual = core.infer_processes(&frq);
    // DOUBLES_EQUAL(56.3469, actual, .0001);

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

TEST(Simulation, gamma_cats)
{
    std::vector<gene_family> families;
    std::ostringstream ost;
    vector<int> rootdist_vec;

    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL, NULL);
    model.start_sim_processes();
}

TEST(Probability, probability_of_some_values)
{
    matrix_cache calc;
    double lambda = 0.05;
    double branch_length = 5;
    DOUBLES_EQUAL(0.0152237, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 5, 9), 0.00001);

    DOUBLES_EQUAL(0.17573, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 9), 0.00001);

    DOUBLES_EQUAL(0.182728, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10), 0.00001);

    branch_length = 1;
    DOUBLES_EQUAL(0.465565, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10), 0.00001);
}

bool operator==(matrix& m1, matrix& m2)
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

TEST(Probability, probability_of_matrix2)
{
    matrix_cache calc;
    single_lambda lambda(0.006335);
    std::set<double> branch_lengths{ 68.7105 };
    calc.precalculate_matrices(141, get_lambda_values(&lambda), branch_lengths);
    DOUBLES_EQUAL(0.194661, calc.get_matrix(141, 68.7105, 0.006335).get(5,5), 0.00001);
}

TEST(Probability, the_probability_of_going_from_parent_fam_size_to_c)
{
    DOUBLES_EQUAL(0.194661, the_probability_of_going_from_parent_fam_size_to_c(.006335, 68.7105, 5, 5), 0.00001);
}

TEST(Probability, probability_of_matrix)
{
    matrix_cache calc;
    single_lambda lambda(0.05);
    std::set<double> branch_lengths{ 5 };
    calc.precalculate_matrices(5, get_lambda_values(&lambda), branch_lengths);
    matrix actual = calc.get_matrix(5, 5, lambda.get_single_lambda());
    matrix expected(5);
    double values[5][5] = {
    {0,0,0,0,0},
    { 0.2,0.64,0.128,0.0256,0.00512 },
    { 0.04,0.256,0.4608,0.17408,0.0512 },
    { 0.008,0.0768,0.26112,0.36352,0.187392 },
    { 0.0016,0.02048,0.1024,0.249856,0.305562 } };
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            expected.set(i, j, values[i][j]);
    CHECK(actual == expected);

    // a second call should get the same results as the first
    actual = calc.get_matrix(5, 5, lambda.get_single_lambda());
    CHECK(actual == expected);
}

TEST(Inference, gamma_model_initial_guesses)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();
    single_lambda sl(0.05);

    gamma_model model(&sl, p_tree, NULL, 0, 5, 4, 0.7, NULL, NULL);
    unique_ptr<optimizer> opt(model.get_lambda_optimizer(NULL));
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(2, guesses.size());
    DOUBLES_EQUAL(0.218321, guesses[0], 0.0001);
    DOUBLES_EQUAL(0.565811, guesses[1], 0.0001);

   delete p_tree;
}

TEST(Inference, base_model_initial_guesses)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> tree(parser.parse_newick());
    single_lambda sl(0.05);

    base_model model(&sl, tree.get(), NULL, 0, 5, NULL, NULL);
    unique_ptr<optimizer> opt(model.get_lambda_optimizer(NULL));
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(1, guesses.size());
    DOUBLES_EQUAL(0.565811, guesses[0], 0.0001);
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

    base_model model(&sl, p_tree.get(), &families, 5, 5, NULL, NULL);

    matrix_cache calc;
    calc.precalculate_matrices(6, get_lambda_values(&sl), set<double>({ 1 }));
    uniform_distribution dist;
    dist.initialize({ 1,2,3,4,5,6 });
    model.reconstruct_ancestral_states(&calc, &dist);

}

TEST(Inference, branch_length_finder)
{
    branch_length_finder finder;
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
    clade *p_tree = parser.parse_newick();
    p_tree->apply_prefix_order(finder);
    LONGS_EQUAL(finder.longest(), 23);
    auto expected = set<double>{ 0, 1, 3, 7, 11, 17, 23 };
    CHECK(finder.result() == expected);
    delete p_tree;

}

TEST(Inference, reconstruction_process)
{
    vector<int> rootdist;
    single_lambda lambda(0.05);
    gene_family fam;
    fam.set_species_size("Mouse", 3);

    clade leaf("Mouse",7);

    matrix_cache calc;
   reconstruction_process process(cout, &lambda, 2.0, NULL, 3, 0, rootdist, &fam, &calc, NULL);
   process(&leaf);
   auto L = process.get_L(&leaf);

   // L holds the probability of the leaf moving from size 3 to size n
   LONGS_EQUAL(4, L.size());
   DOUBLES_EQUAL(0.0, L[0], 0.0001);
   DOUBLES_EQUAL(0.0586679, L[1], 0.0001);
   DOUBLES_EQUAL(0.146916, L[2], 0.0001);
   DOUBLES_EQUAL(0.193072, L[3], 0.0001);
}

TEST(Inference, precalculate_matrices_calculates_all_lambdas_all_branchlengths)
{
    matrix_cache calc;
    std::map<std::string, int> m;
    multiple_lambda lambda(m, vector<double>({ .1, .2, .3, .4 }));
    calc.precalculate_matrices(5, get_lambda_values(&lambda), set<double>({ 1,2,3 }));
    LONGS_EQUAL(12, calc.get_cache_size());
}

TEST(Inference, reconstruction_process_internal_node)
{
    vector<int> rootdist;
    single_lambda s_lambda(0.05);
    double multiplier = 2.0;
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);

    newick_parser parser(false);
    parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
    unique_ptr<clade> p_tree(parser.parse_newick());
    unique_ptr<lambda> m(s_lambda.multiply(multiplier));
    matrix_cache calc;
    calc.precalculate_matrices(4, get_lambda_values(m.get()), set<double>({ 1, 3, 7, 11, 17, 23 }));
    reconstruction_process process(cout, &s_lambda, multiplier, NULL, 3, 0, rootdist, &fam, &calc, NULL);

    process(p_tree->find_descendant("A"));
    process(p_tree->find_descendant("B"));

    clade *internal_node = p_tree->find_descendant("AB");
    process(internal_node);
    auto L = process.get_L(internal_node);

    // L holds the probability of the leaf moving from size 3 to size n
    LONGS_EQUAL(4, L.size());
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
    inference_process_factory factory(cout, &lambda, p_tree.get(), 10, 8, {2,2,2});
    factory.set_gene_family(&fam);
    gamma_bundle bundle(factory, { 0.1, 0.5 });
    uniform_distribution dist;
    dist.initialize({ 1,2,3,4,5,4,3,2,1 });
    matrix_cache cache;
    multiple_lambda ml(map<string, int>(), {0.0005, 0.0025});
    cache.precalculate_matrices(11, get_lambda_values(&ml), set<double>{1, 3, 7});
    vector<double> cat_likelihoods;
    CHECK(bundle.prune({ 0.01, 0.05 }, &dist, cache, cat_likelihoods)); 

    LONGS_EQUAL(2, cat_likelihoods.size());
    DOUBLES_EQUAL(-22.5251, log(cat_likelihoods[0]), 0.0001);
    DOUBLES_EQUAL(-16.1766, log(cat_likelihoods[1]), 0.0001);

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
    inference_process_factory factory(cout, &lambda, p_tree.get(), 10, 8, { 2,2,2 });
    factory.set_gene_family(&fam);
    gamma_bundle bundle(factory, { 0.1, 0.5 });
    uniform_distribution dist;
    dist.initialize({ 1,2,3,4,5,4,3,2,1 });
    matrix_cache cache;
    cache.precalculate_matrices(11, { 0.09, 0.45 }, set<double>{1, 3, 7});
    vector<double> cat_likelihoods;

    CHECK(!bundle.prune({ 0.01, 0.05 }, &dist, cache, cat_likelihoods));
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
    clade c;
    single_lambda lambda(0.02);
    vector<gene_family> fams;
    auto models = build_models(params, &c, &lambda, &fams, 10, 10, NULL);
    LONGS_EQUAL(1, models.size());
    CHECK(dynamic_cast<base_model *>(models[0]));
    for (auto m : models)
        delete m;
}

TEST(Inference, build_models_sets_error_model)
{
    input_parameters params;
    params.input_file_path = "foo";
    clade c;
    single_lambda lambda(0.02);
    vector<gene_family> fams;
    auto models = build_models(params, &c, &lambda, &fams, 10, 10, NULL);
    LONGS_EQUAL(1, models.size());
    CHECK(dynamic_cast<base_model *>(models[0]));
    for (auto m : models)
        delete m;
}

TEST(Inference, base_model_computes_epsilon)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(1);
    error_model err;
    err.set_probs(0, { .0, .7, .3 });
    err.set_probs(1, { .4, .2, .4 });
    std::vector<gene_family> fams{ fam };
    uniform_distribution frq;

    base_model model(&lambda, p_tree.get(), &fams, 10, 10, NULL, &err);
    unique_ptr<optimizer> opt(model.get_epsilon_optimizer(&frq));

    opt->optimize();
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
    model.set_probs(0, { .0, .7, .3 });
    model.set_probs(1, { .2, .6, .2 });
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

TEST(Probability, base_model_initial_epsilon)
{
    error_model err;
    err.set_probs(0, { .0, .7, .3 });
    err.set_probs(1, { .4, .2, .4 });
    std::vector<gene_family> fams;
    base_model model(NULL, NULL, &fams, 10, 10, NULL, &err);
    optimizer* opt = model.get_epsilon_optimizer(NULL);
    auto actual = opt->initial_guesses();
    LONGS_EQUAL(2, actual.size());
    double a;
    opt->finalize(&a);
    delete opt;
}

TEST(Probability, error_model_get_epsilon)
{
    error_model model;
    model.set_probs(0, { .0, .7, .3 });
    model.set_probs(1, { .2, .6, .2 });
    model.set_probs(2, { .1, .8, .1 });
    model.set_probs(3, { .2, .6, .2 });

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
        model.set_probs(0, { 0.4, 0.3, 0.3 });
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
        model.set_probs(1, { 0.3, 0.3, 0.3 });
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
    inference_process process(ost, &lambda, 1.5, p_tree.get(), 20, 20, &fam, { 1,2,3 }, NULL);
    matrix_cache cache;
    cache.precalculate_matrices(21, { 0.045 }, { 1.0,3.0,7.0 });
    auto actual = process.prune(cache);
    vector<double> log_expected{ -22.5513, -14.7206, -8.85104, -4.67892, -4.57456, -5.52831, -7.37063, -10.152, 
        -13.4435, -17.0609, -20.9074, -24.9251, -29.076, -33.334, -37.6799, -42.0995, -46.5817, -51.1179, -55.7011, -60.3256 };

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
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache;
    likelihood_computer pruner(20, 20, &lambda, &fam, cache, NULL);
    cache.precalculate_matrices(21, { 0.045 }, { 1.0,3.0,7.0 });

    clade *A = p_tree->find_descendant("A");
    pruner(A);
    auto actual = pruner.get_likelihoods(A);

    vector<double> expected{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i<expected.size(); ++i)
    {
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }

    clade *B = p_tree->find_descendant("B");
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
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache;
    likelihood_computer pruner(20, 20, &lambda, &fam, cache, NULL);
    cache.precalculate_matrices(21, { 0.03 }, { 1.0,3.0,7.0 });

    clade *AB = p_tree->find_descendant("AB");
    try
    {
        pruner(AB);
        CHECK(false);   
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Child node probabilities not calculated", err.what());
    }

    pruner(p_tree->find_descendant("A"));
    pruner(p_tree->find_descendant("B"));
    pruner(AB);

    auto actual = pruner.get_likelihoods(AB);

    vector<double> log_expected{ -25.8046, -17.1121, -10.388, -5.37765, -5.24652, -6.21594, -8.18487, -11.5884,
        -15.5796, -19.9223, -24.5056, -29.2661, -34.1635, -39.1702, -44.2664, -49.4374, -54.6718, -59.9609, -65.2974, -70.6756 };

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
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache;
    cache.precalculate_matrices(21, { 0.045 }, { 1.0,3.0,7.0 });

    string input = "maxcnt: 20\ncntdiff: -1 0 1\n"
        "1 0.2 0.6 0.2\n"
        "20 0.2 0.6 0.2\n";
    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    likelihood_computer pruner(20, 20, &lambda, &fam, cache, &model);

    clade *A = p_tree->find_descendant("A");
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


int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}