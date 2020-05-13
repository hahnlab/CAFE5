
#include <numeric>
#include <cmath>
#include <getopt.h>
#include <sstream>
#include <random>
#include <algorithm>

#include <string.h>

#define ELPP_NO_CHECK_MACROS
#include "src/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

#define DOCTEST_CONFIG_IMPLEMENT
#define DOCTEST_CONFIG_IMPLEMENT
#include "src/doctest.h"

#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"
#include "src/root_equilibrium_distribution.h"
#include "src/base_model.h"
#include "src/gene_family_reconstructor.h"
#include "src/matrix_cache.h"
#include "src/probability.h"
#include "src/execute.h"
#include "src/user_data.h"
#include "src/optimizer_scorer.h"
#include "src/simulator.h"
#include "src/poisson.h"
#include "src/optimizer.h"
#include "src/error_model.h"
#include "src/likelihood_ratio.h"

std::mt19937 randomizer_engine(10); // seeding random number engine

class mock_model : public model {
    // Inherited via model
    virtual std::string name() const override
    {
        return "mockmodel";
    }
    virtual void write_family_likelihoods(std::ostream& ost) override
    {
    }
    virtual reconstruction* reconstruct_ancestral_states(const vector<gene_family>& families, matrix_cache* p_calc, root_equilibrium_distribution* p_prior) override
    {
        return nullptr;
    }
    virtual inference_optimizer_scorer* get_lambda_optimizer(const user_data& data) override
    {
        auto lengths = _p_tree->get_branch_lengths();
        auto longest_branch = *max_element(lengths.begin(), lengths.end());

        initialize_lambda(data.p_lambda_tree);
        auto result = new lambda_optimizer(_p_lambda, this, &data.prior, longest_branch);
        result->quiet = true;
        return result;
    }
    bool _invalid_likelihood = false;
public:
    mock_model() : model(NULL, NULL, NULL, 0, 0, NULL)
    {

    }
    void set_lambda(lambda* lambda)
    {
        _p_lambda = lambda;
    }
    void set_tree(clade* tree)
    {
        _p_tree = tree;
    }
    void set_invalid_likelihood() { _invalid_likelihood = true;  }
    // Inherited via model
    virtual void prepare_matrices_for_simulation(matrix_cache& cache) override
    {
        cache.precalculate_matrices(get_lambda_values(_p_lambda), _p_tree->get_branch_lengths());
    }

    // Inherited via model
    virtual double infer_family_likelihoods(const root_equilibrium_distribution& prior, const lambda* p_lambda) override
    {
        return _invalid_likelihood ? nan("") : 0.0;
    }
};

class Inference
{
public:
    user_data _user_data;
    single_lambda* _p_lambda;

    Inference()
    {
        _p_lambda = new single_lambda(0.05);
        _user_data.p_tree = parse_newick("(A:1,B:1);");
        _user_data.p_lambda = _p_lambda;
        _user_data.max_family_size = 10;
        _user_data.max_root_family_size = 8;
        _user_data.gene_families.resize(1);
        _user_data.gene_families[0].set_id("TestFamily1");
        _user_data.gene_families[0].set_species_size("A", 1);
        _user_data.gene_families[0].set_species_size("B", 2);

        randomizer_engine.seed(10);
    }

    ~Inference()
    {
        delete _user_data.p_tree;
        delete _p_lambda;
    }
};

class Optimizer
{
public:
    FMinSearch fm;
    vector<candidate> c;
    Optimizer()
    {
        for (int i = 0; i < 3; ++i)
            c.push_back(candidate(3));
        fm.variable_count = 2;
        fm.variable_count_plus_one = 3;
        fm.candidates.resize(3);
        transform(c.begin(), c.end(), fm.candidates.begin(), [](candidate& c) { return &c;  });
    }

};


input_parameters read_arguments(int argc, char *const argv[]);

struct option_test
{
    char* values[100];
    size_t argc;

    option_test(vector<string> arguments)
    {
        optind = 0;
        argc = arguments.size();
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            values[i] = strdup(arguments[i].c_str());
        }
    }

    ~option_test()
    {
        for (size_t i = 0; i < argc; ++i)
        {
            free(values[i]);
        }
    }
};

#define STRCMP_EQUAL(x, y) CHECK(strcmp(x,y) == 0)
#define STRCMP_CONTAINS(x, y) CHECK(strstr(y,x) != nullptr)

TEST_CASE("read_arguments translates short values ") {
    option_test c({ "cafexp", "-ifile" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}

TEST_CASE("read_arguments translates long values ") {
    option_test c({ "cafexp", "--infile", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}
TEST_CASE("Options, input_short_space_separated")
{
    option_test c({ "cafexp", "-i", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}
TEST_CASE("Options, simulate_long")
{
    option_test c({ "cafexp", "--simulate=1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, simulate_short")
{
    option_test c({ "cafexp", "-s1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, pvalue_long")
{
    option_test c({ "cafexp", "--pvalue=0.01" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.01, actual.pvalue);
}

TEST_CASE("Options, pvalue_short")
{
    option_test c({ "cafexp", "-P0.01" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.01, actual.pvalue);
}

TEST_CASE("Options, optimizer_long")
{
    option_test c({ "cafexp", "--optimizer_expansion=0.05", "--optimizer_reflection=3.2", "--optimizer_iterations=5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.neldermead_iterations);
}

TEST_CASE("Options, optimizer_short")
{
    option_test c({ "cafexp", "-E", "0.05", "-R", "3.2", "-I", "5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.neldermead_iterations);
}

TEST_CASE("Options, errormodel_accepts_argument")
{
    option_test c({ "cafexp", "-eerror.txt" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    STRCMP_EQUAL("error.txt", actual.error_model_file_path.c_str());
}

TEST_CASE("Options, errormodel_accepts_no_argument")
{
    option_test c({ "cafexp", "-e" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path.empty());
}

TEST_CASE("Options, zero_root_familes")
{
    input_parameters by_default;
    CHECK(by_default.exclude_zero_root_families);

    option_test c({ "cafexp", "-z" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_FALSE(actual.exclude_zero_root_families);
}

TEST_CASE("Options: cannot_have_space_before_optional_parameter")
{
    option_test c({ "cafexp", "-s", "1000" });

    CHECK_THROWS_WITH(read_arguments(c.argc, c.values), "Unrecognized parameter: '1000'");
}

TEST_CASE("Options: must_specify_lambda_and_input_file_for_estimator")
{
    input_parameters params;
    params.fixed_lambda = 0.05;
    CHECK_THROWS_WITH(params.check_input(), "Options -l and -i must both be provided an argument.");
}

TEST_CASE("Options: must_specify_lambda_for_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate without initial lambda values");
}

TEST_CASE("Options: must_specify_alpha_for_gamma_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_lambda = 0.05;
    params.n_gamma_cats = 3;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate gamma clusters without an alpha value");
}

TEST_CASE("Options: must_specify_alpha_and_k_for_gamma_inference")
{
    input_parameters params;
    params.fixed_alpha = 0.7;
    CHECK_THROWS_WITH(params.check_input(), "Alpha specified with 1 gamma category.");
}

TEST_CASE("Options: can_specify_alpha_without_k_for_gamma_simulation")
{
    input_parameters params;
    params.fixed_alpha = 0.7;
    params.fixed_lambda = 0.01;
    params.is_simulating = true;
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: check_input_does_not_throw_when_simulating_with_multiple_lambdas")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_multiple_lambdas = "0.01,0.05";
    params.lambda_tree_file_path = "./tree";
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: per_family_must_provide_families")
{
    input_parameters params;
    params.lambda_per_family = true;
    CHECK_THROWS_WITH_AS(params.check_input(), "No family file provided", runtime_error);
}

TEST_CASE("Options: per_family_must_provide_tree")
{
    input_parameters params;
    params.lambda_per_family = true;
    params.input_file_path = "/tmp/test";
    CHECK_THROWS_WITH_AS(params.check_input(), "No tree file provided", runtime_error);
}

TEST_CASE("Options: cannot_estimate_error_and_gamma_together")
{
    input_parameters params;
    params.n_gamma_cats = 3;
    params.use_error_model = true;
    params.error_model_file_path = "model.txt";
    params.check_input();
    CHECK(true);

    params.error_model_file_path.clear();
    CHECK_THROWS_WITH_AS(params.check_input(), "Estimating an error model with a gamma distribution is not supported at this time", runtime_error);

}

TEST_CASE("GeneFamilies: read_gene_families_reads_cafe_files")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, families);
    CHECK_EQ(5, families.at(0).get_species_size("A"));
    CHECK_EQ(10, families.at(0).get_species_size("B"));
    CHECK_EQ(2, families.at(0).get_species_size("C"));
    CHECK_EQ(6, families.at(0).get_species_size("D"));
}

TEST_CASE("GeneFamilies: init_from_clademap")
{
    clade* p_tree = parse_newick("((A:1,B:1):1,(C:1,D:1):1);");
    clademap<int> values;
    values[p_tree->find_descendant("A")] = 3;
    values[p_tree->find_descendant("B")] = 5;
    values[p_tree->find_descendant("C")] = 7;
    values[p_tree->find_descendant("D")] = 11;
    gene_family gf;
    gf.init_from_clademap(values);
    CHECK_EQ(3, gf.get_species_size("A"));
}

TEST_CASE("GeneFamilies: read_gene_families_reads_simulation_files")
{
    std::string str = "#A\n#B\n#AB\n#CD\n#C\n#ABCD\n#D\n35\t36\t35\t35\t36\t34\t34\t1\n98\t96\t97\t98\t98\t98\t98\t1\n";
    std::istringstream ist(str);

    clade* p_tree = parse_newick("((A:1,B:1):1,(C:1,D:1):1);");

    std::vector<gene_family> families;
    read_gene_families(ist, p_tree, families);
    CHECK_EQ(35, families.at(0).get_species_size("A"));
    CHECK_EQ(36, families.at(0).get_species_size("B"));
    CHECK_EQ(36, families.at(0).get_species_size("C"));
    CHECK_EQ(34, families.at(0).get_species_size("D"));
    delete p_tree;
}

TEST_CASE("GeneFamilies: read_gene_families_throws_if_no_families_found")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));
    std::vector<gene_family> families;
    CHECK_THROWS_WITH_AS(read_gene_families(ist, p_tree.get(), families), "No families found", runtime_error);
}

TEST_CASE("GeneFamilies: model_set_families")
{
    mock_model m;

    vector<gene_family> fams;
    m.set_families(&fams);
    CHECK_EQ(0, m.get_gene_family_count());

    fams.resize(5);
    CHECK_EQ(5, m.get_gene_family_count());
}

TEST_CASE("GeneFamilies: species_size_is_case_insensitive")
{
    gene_family gf;
    gf.set_species_size("Human", 5);
    CHECK_EQ(5, gf.get_species_size("human"));
    CHECK_EQ(5, gf.get_species_size("HUMAN"));
    CHECK_EQ(5, gf.get_species_size("hUmAn"));
}

TEST_CASE("GeneFamilies: species_size_differential")
{
    gene_family gf;
    gf.set_species_size("Cat", 5);
    gf.set_species_size("Horse", 3);
    gf.set_species_size("Cow", 1);

    CHECK_EQ(4, gf.species_size_differential());

    gf.set_species_size("Chicken", 12);
    CHECK_EQ(11, gf.species_size_differential());
}

TEST_CASE_FIXTURE(Inference, "infer_processes")
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

    base_model core(&lambda, _user_data.p_tree, &families, 56, _user_data.max_root_family_size, NULL);

    double multi = core.infer_family_likelihoods(_user_data.prior, &lambda);

    CHECK_EQ(doctest::Approx(46.56632), multi);
}

TEST_CASE("Inference: root_equilibrium_distribution__with_no_rootdist_is_uniform")
{
    root_equilibrium_distribution ef(10);
    CHECK_EQ(doctest::Approx(.1), ef.compute(5));
}

TEST_CASE_FIXTURE(Inference, "root_equilibrium_distribution__with_rootdist_uses_rootdist")
{
    _user_data.rootdist[1] = 3;
    _user_data.rootdist[2] = 5;

    root_equilibrium_distribution ef(_user_data.rootdist);

    CHECK_EQ(0.375, ef.compute(1));
}

TEST_CASE("Inference: gamma_set_alpha")
{
    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL);
    model.set_alpha(0.5);
    CHECK(true);
}

TEST_CASE( "Inference: gamma_adjust_family_gamma_membership")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, families);

    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL);
    CHECK(true);
}

TEST_CASE_FIXTURE(Inference, "gamma_model_infers_processes_without_crashing")
{
    gamma_model core(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 148, 122, 1, 0, NULL);

    // TODO: make this return a non-infinite value and add a check for it
    core.infer_family_likelihoods(_user_data.prior, _user_data.p_lambda);
    CHECK(true);

}

TEST_CASE("Inference: stash_stream")
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

TEST_CASE("Probability: probability_of_some_values")
{
    matrix_cache calc(0);
    double lambda = 0.05;
    double branch_length = 5;
    CHECK_EQ(doctest::Approx(0.0152237).scale(10000), calc.get_from_parent_fam_size_to_c(lambda, branch_length, 5, 9));

    CHECK_EQ(doctest::Approx(0.17573).scale(10000), calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 9));

    CHECK_EQ(doctest::Approx(0.182728).scale(10000), calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10));

    branch_length = 1;
    CHECK_EQ(doctest::Approx(0.465565).scale(10000), calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10));
}

TEST_CASE("Probability:matrices_take_fractional_branch_lengths_into_account")
{
    matrix_cache calc(141);
    single_lambda lambda(0.006335);
    std::set<double> branch_lengths{ 68, 68.7105 };
    calc.precalculate_matrices(get_lambda_values(&lambda), branch_lengths);
    CHECK_EQ(doctest::Approx(0.194661).epsilon(0.0001), calc.get_matrix(68.7105, 0.006335)->get(5, 5)); // a value 
    CHECK_EQ(doctest::Approx(0.195791).epsilon(0.0001), calc.get_matrix(68, 0.006335)->get(5, 5));
}

TEST_CASE("Probability: the_probability_of_going_from_parent_fam_size_to_c")
{
    CHECK_EQ(doctest::Approx(0.194661).epsilon(.00001), the_probability_of_going_from_parent_fam_size_to_c(.006335, 68.7105, 5, 5));
}

bool operator==(const matrix& m1, const matrix& m2)
{
    if (m1.size() != m2.size())
        return false;

    for (int i = 0; i < m1.size(); ++i)
    {
        for (int j = 0; j < m1.size(); ++j)
            if (abs(m1.get(i, j) - m2.get(i, j)) > 0.00001)
                return false;
    }

    return true;
}


TEST_CASE("Probability: probability_of_matrix")
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

TEST_CASE("Probability: get_random_probabilities")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    single_lambda lam(0.05);
    matrix_cache cache(15);
    cache.precalculate_matrices(vector<double>{0.05}, set<double>{1});
    auto probs = get_random_probabilities(p_tree.get(), 10, 3, 12, 8, &lam, cache, NULL);
    CHECK_EQ(10, probs.size());
    CHECK_EQ(doctest::Approx(0.001905924).scale(10000), probs[0]);
}

TEST_CASE_FIXTURE(Inference, "base_optimizer_guesses_lambda_only")
{
    _user_data.p_lambda = NULL;

    base_model model(_user_data.p_lambda, _user_data.p_tree, NULL, 0, 5, NULL);

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(_user_data));
    auto guesses = opt->initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(doctest::Approx(0.2498383).epsilon(0.00001), guesses[0]);
    delete model.get_lambda();
}

TEST_CASE_FIXTURE(Inference, "base_optimizer_guesses_lambda_and_unique_epsilons")
{
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    base_model model(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 10, 10, &err);

    _user_data.p_error_model = nullptr;
    _user_data.p_lambda = nullptr;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(_user_data));

    CHECK(opt);
    CHECK(dynamic_cast<lambda_epsilon_optimizer*>(opt.get()) != nullptr);
    auto guesses = opt->initial_guesses();
    CHECK_EQ(3, guesses.size());
    CHECK_EQ(doctest::Approx(0.2498383).epsilon(0.00001), guesses[0]);
    CHECK_EQ(0.3, guesses[1]);
    CHECK_EQ(0.4, guesses[2]);

    delete model.get_lambda();
}


TEST_CASE_FIXTURE(Inference, "gamma_model_creates__gamma_lambda_optimizer_if_nothing_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, -1, NULL);
    user_data data;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));
    CHECK(opt);
    CHECK(dynamic_cast<gamma_lambda_optimizer*>(opt.get()));

    delete model.get_lambda();
}

TEST_CASE("Inference: gamma_model__creates__lambda_optimizer__if_alpha_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, 0.25, NULL);

    user_data data;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));

    CHECK(opt);
    CHECK(dynamic_cast<lambda_optimizer*>(opt.get()));
    delete model.get_lambda();
}

TEST_CASE("Inference: gamma_model__creates__gamma_optimizer__if_lambda_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, -1, NULL);

    user_data data;

    single_lambda sl(0.05);
    data.p_lambda = &sl;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));

    CHECK(opt);
    CHECK(dynamic_cast<gamma_optimizer*>(opt.get()));

    delete model.get_lambda();
}

TEST_CASE("Inference: gamma_model_creates_nothing_if_lambda_and_alpha_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, p_tree.get(), NULL, 0, 5, 4, .25, NULL);

    user_data data;

    single_lambda sl(0.05);
    data.p_lambda = &sl;

    CHECK(model.get_lambda_optimizer(data) == nullptr);
}

TEST_CASE("Inference: gamma_lambda_optimizer__provides_two_guesses")
{
    single_lambda sl(0.05);
    gamma_model model(NULL, NULL, NULL, 0, 5, 4, .25, NULL);
    gamma_lambda_optimizer glo(&sl, &model, NULL, 5);
    auto guesses = glo.initial_guesses();
    CHECK_EQ(2, guesses.size());

    double lambda = guesses[0];
    CHECK_GT(lambda, 0);
    CHECK_LT(lambda, 1);

    double alpha = guesses[1];
    CHECK_GT(alpha, 0);
    CHECK_LT(alpha, 10);
}

TEST_CASE("Inference: gamma_optimizer__creates_single_initial_guess")
{
    gamma_model m(NULL, NULL, NULL, 0, 0, 0, 0, NULL);
    gamma_optimizer optimizer(&m, NULL);
    auto initial = optimizer.initial_guesses();
    CHECK_EQ(1, initial.size());
}

TEST_CASE_FIXTURE(Inference, "base_model_reconstruction")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1"));
    single_lambda sl(0.05);

    std::vector<gene_family> families(1);
    families[0].set_species_size("A", 3);
    families[0].set_species_size("B", 4);

    base_model model(&sl, p_tree.get(), &families, 5, 5, NULL);

    matrix_cache calc(6);
    calc.precalculate_matrices(get_lambda_values(&sl), set<double>({ 1 }));
    root_equilibrium_distribution dist(_user_data.max_root_family_size);

    std::unique_ptr<base_model_reconstruction> rec(dynamic_cast<base_model_reconstruction*>(model.reconstruct_ancestral_states(families, &calc, &dist)));

    CHECK_EQ(1, rec->_reconstructions.size());

}

TEST_CASE("Inference: branch_length_finder")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));
    auto actual = p_tree->get_branch_lengths();
    auto expected = set<double>{ 1, 3, 7, 11, 17, 23 };
    CHECK(actual == expected);
}

TEST_CASE("Inference: increase_decrease")
{
    base_model_reconstruction bmr;
    gene_family gf;

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

    auto a = p_tree->find_descendant("A");
    auto b = p_tree->find_descendant("B");
    auto ab = p_tree->find_descendant("AB");
    auto abcd = p_tree->find_descendant("ABCD");

    gf.set_id("myid");
    gf.set_species_size("A", 4);
    gf.set_species_size("B", 2);
    bmr._reconstructions["myid"][ab] = 3;
    bmr._reconstructions["myid"][abcd] = 3;

    CHECK_EQ(1, bmr.get_difference_from_parent(gf, a));
    CHECK_EQ(-1, bmr.get_difference_from_parent(gf, b));
    CHECK_EQ(0, bmr.get_difference_from_parent(gf, ab));
}

TEST_CASE( "Inference: precalculate_matrices_calculates_all_lambdas_all_branchlengths")
{
    matrix_cache calc(5);
    std::map<std::string, int> m;
    multiple_lambda lambda(m, vector<double>({ .1, .2, .3, .4 }));
    calc.precalculate_matrices(get_lambda_values(&lambda), set<double>({ 1,2,3 }));
    CHECK_EQ(12, calc.get_cache_size());
}

class Reconstruction
{
public:
    gene_family fam;
    unique_ptr<clade> p_tree;
    cladevector order;

    Reconstruction()
    {
        p_tree.reset(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

        fam.set_id("Family5");
        fam.set_species_size("A", 11);
        fam.set_species_size("B", 2);
        fam.set_species_size("C", 5);
        fam.set_species_size("D", 6);

        vector<string> nodes{ "A", "B", "C", "D", "AB", "CD", "ABCD" };
        order.resize(nodes.size());
        const clade* t = p_tree.get();
        transform(nodes.begin(), nodes.end(), order.begin(), [t](string s) { return t->find_descendant(s); });

    }
};

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_leaf_node")
{
    single_lambda lambda(0.1);
    fam.set_species_size("Mouse", 3);

    clade leaf("Mouse", 7);

    matrix_cache calc(8);
    calc.precalculate_matrices({ .1 }, set<double>({ 7 }));
    clademap<std::vector<int>> all_node_Cs;
    clademap<std::vector<double>> all_node_Ls;
    reconstruct_leaf_node(&leaf, &lambda, all_node_Cs, all_node_Ls, 7, &fam, &calc);

    // L holds the probability of the leaf moving from size 3 to size n
    auto L = all_node_Ls[&leaf];

    CHECK_EQ(8, L.size());
    CHECK_EQ(0.0, L[0]);
    CHECK_EQ(doctest::Approx(0.0586679), L[1]);
    CHECK_EQ(doctest::Approx(0.146916), L[2]);
    CHECK_EQ(doctest::Approx(0.193072), L[3]);
}

TEST_CASE_FIXTURE(Reconstruction, "print_reconstructed_states__prints_star_for_significant_values")
{
    gene_family gf;
    gf.set_id("Family5");
    base_model_reconstruction bmr;
    auto& values = bmr._reconstructions[gf.id()];

    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    branch_probabilities branch_probs;
    p_tree->apply_reverse_level_order([&branch_probs, &gf](const clade* c) {branch_probs.set(gf, c, branch_probabilities::branch_probability(.5)); });
    branch_probs.set(gf, p_tree->find_descendant("AB"), 0.02);
    branch_probs.set(gf, p_tree.get(), branch_probabilities::invalid());  /// root is never significant regardless of the value

    ostringstream sig;
    bmr.print_reconstructed_states(sig, order, { fam }, p_tree.get(), 0.05, branch_probs);
    STRCMP_CONTAINS("  TREE Family5 = ((A<0>_11:1,B<1>_2:3)<4>*_8:7,(C<2>_5:11,D<3>_6:17)<5>_6:23)<6>_7;", sig.str().c_str());

    ostringstream insig;
    bmr.print_reconstructed_states(insig, order, { fam }, p_tree.get(), 0.01, branch_probs);
    STRCMP_CONTAINS("  TREE Family5 = ((A<0>_11:1,B<1>_2:3)<4>_8:7,(C<2>_5:11,D<3>_6:17)<5>_6:23)<6>_7;", insig.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__print_reconstructed_states__prints_value_for_each_category_and_a_summation")
{
    gamma_model_reconstruction gmr(vector<double>({ 1.0 }));

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()] = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")] = 0;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")] = 0;

    rec.reconstruction[p_tree.get()] = 7;
    rec.reconstruction[p_tree->find_descendant("AB")] = 8;
    rec.reconstruction[p_tree->find_descendant("CD")] = 6;

    ostringstream ost;
    branch_probabilities branch_probs;
    gmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);
    STRCMP_CONTAINS("  TREE Family5 = ((A<0>_11:1,B<1>_2:3)<4>_8:7,(C<2>_5:11,D<3>_6:17)<5>_6:23)<6>_7;", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__print_additional_data__prints_likelihoods")
{
    gamma_model_reconstruction gmr(vector<double>({ 0.3, 0.9, 1.4, 2.0 }));
    gmr._reconstructions["Family5"]._category_likelihoods = { 0.01, 0.03, 0.09, 0.07 };
    ostringstream ost;
    gmr.print_category_likelihoods(ost, order, { fam });
    STRCMP_CONTAINS("Family ID\t0.3\t0.9\t1.4\t2\t\n", ost.str().c_str());
    STRCMP_CONTAINS("Family5\t0.01\t0.03\0.09\t0.07", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__prints_lambda_multipiers")
{
    vector<double> multipliers{ 0.13, 1.4 };
    gamma_model_reconstruction gmr(multipliers);

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()] = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")] = 8;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")] = 6;

    rec.reconstruction[p_tree.get()] = 7;
    rec.reconstruction[p_tree->find_descendant("AB")] = 8;
    rec.reconstruction[p_tree->find_descendant("CD")] = 6;

    branch_probabilities branch_probs;

    std::ostringstream ost;
    gmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);

    STRCMP_CONTAINS("BEGIN LAMBDA_MULTIPLIERS;", ost.str().c_str());
    STRCMP_CONTAINS("  0.13;", ost.str().c_str());
    STRCMP_CONTAINS("  1.4;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "base_model_reconstruction__print_reconstructed_states")
{
    base_model_reconstruction bmr;
    auto& values = bmr._reconstructions["Family5"];

    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    branch_probabilities branch_probs;

    ostringstream ost;

    bmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);
    STRCMP_CONTAINS("#nexus", ost.str().c_str());
    STRCMP_CONTAINS("BEGIN TREES;", ost.str().c_str());
    STRCMP_CONTAINS("  TREE Family5 = ((A<0>_11:1,B<1>_2:3)<4>_8:7,(C<2>_5:11,D<3>_6:17)<5>_6:23)<6>_7;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruction_process_internal_node")
{
    single_lambda s_lambda(0.1);
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);

    matrix_cache calc(25);
    calc.precalculate_matrices({ 0.1 }, set<double>({ 1, 3, 7, 11, 17, 23 }));

    clademap<std::vector<int>> all_node_Cs;
    clademap<std::vector<double>> all_node_Ls;

    reconstruct_leaf_node(p_tree->find_descendant("A"), &s_lambda, all_node_Cs, all_node_Ls, 24, &fam, &calc);
    reconstruct_leaf_node(p_tree->find_descendant("B"), &s_lambda, all_node_Cs, all_node_Ls, 24, &fam, &calc);

    auto internal_node = p_tree->find_descendant("AB");
    reconstruct_internal_node(internal_node, &s_lambda, all_node_Cs, all_node_Ls, 24, &calc);
    auto L = all_node_Ls[internal_node];

    // L holds the probability of the node moving from size 3 to size n
    CHECK_EQ(25, L.size());
    CHECK_EQ(0.0, L[0]);
    CHECK_EQ(doctest::Approx(0.00101688), L[1]);
    CHECK_EQ(doctest::Approx(0.00254648), L[2]);
    CHECK_EQ(doctest::Approx(0.0033465), L[3]);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_family")
{
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    single_lambda lambda(0.005);
    matrix_cache cache(11);

    cache.precalculate_matrices(get_lambda_values(&lambda), set<double>{1, 3, 7});

    user_data ud;
    vector<int> v({ 1,2,3,4,5,4,3,2,1 });
    for (size_t i = 0; i < v.size(); ++i)
        ud.rootdist[i] = v[i];
    ud.max_root_family_size = 8;
    root_equilibrium_distribution dist(ud.rootdist);

    clademap<int> result;
    reconstruct_gene_family(&lambda, p_tree.get(), 10, ud.max_root_family_size, &fam, &cache, &dist, result);
    auto AB = p_tree->find_descendant("AB");
    CHECK_EQ(4, result[AB]);
}

TEST_CASE_FIXTURE(Reconstruction, "get_weighted_averages")
{
    clade c1;
    clade c2;

    clademap<int> rc1;
    rc1[&c1] = 10;
    rc1[&c2] = 2;

    clademap<int> rc2;
    rc2[&c1] = 20;
    rc2[&c2] = 8;

    auto avg = get_weighted_averages({ rc1, rc2 }, { .25, .75 });
    CHECK_EQ(17.5, avg[&c1]);
    CHECK_EQ(6.5, avg[&c2]);
}

TEST_CASE_FIXTURE(Reconstruction, "print_node_counts")
{
    gamma_model_reconstruction gmr({ .5 });
    ostringstream ost;

    auto initializer = [&gmr](const clade* c) { gmr._reconstructions["Family5"].reconstruction[c] = 5;  };
    p_tree->apply_prefix_order(initializer);

    gmr.print_node_counts(ost, order, { fam }, p_tree.get());
    STRCMP_CONTAINS("FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>", ost.str().c_str());
    STRCMP_CONTAINS("Family5\t11\t2\t5\t6\t5\t5\t5", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_node_change")
{
    gamma_model_reconstruction gmr({ .5 });
    ostringstream ost;

    std::normal_distribution<float> dist(0, 10);
    clademap<int> size_deltas;
    p_tree->apply_prefix_order([&gmr, &dist](const clade* c) {
        if (!c->is_leaf())
            gmr._reconstructions["Family5"].reconstruction[c] = dist(randomizer_engine);
        });

    gmr.print_node_change(ost, order, { fam }, p_tree.get());
    CHECK_MESSAGE(ost.str().find("FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>") != string::npos, ost.str());
    CHECK_MESSAGE(ost.str().find("Family5\t+1\t-8\t+5\t+6\t+17\t+7\t+0") != string::npos, ost.str());
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_index_in_angle_brackets_for_non_leaf")
{
    STRCMP_EQUAL("<0>", clade_index_or_name(p_tree.get(), { p_tree.get() }).c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_name_plus_index_in_angle_brackets_for_leaf")
{
    auto a = p_tree->find_descendant("A");
    STRCMP_EQUAL("A<1>", clade_index_or_name(a, { p_tree.get(), a }).c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_branch_probabilities__shows_NA_for_invalids")
{
    gene_family gf;
    gf.set_id("Family5");
    std::ostringstream ost;
    branch_probabilities probs;
    for (auto c : order)
        probs.set(gf, c, 0.05);
    probs.set(gf, p_tree->find_descendant("B"), branch_probabilities::invalid());
    probs.set(gf, p_tree.get(), branch_probabilities::invalid());

    print_branch_probabilities(ost, order, { fam }, probs);
    STRCMP_CONTAINS("FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>", ost.str().c_str());
    STRCMP_CONTAINS("Family5\t0.05\tN/A\t0.05\t0.05\t0.05\t0.05\tN/A\n", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_branch_probabilities__skips_families_without_reconstructions")
{
    std::ostringstream ost;
    branch_probabilities probs;

    print_branch_probabilities(ost, order, { fam }, probs);
    CHECK(ost.str().find("Family5") == string::npos);
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities")
{
    matrix_cache cache(25);
    cache.precalculate_matrices({ 0.05 }, { 1,3,7 });
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree->find_descendant("AB")] = 10;
    single_lambda lm(0.05);
    CHECK_EQ(doctest::Approx(0.2182032), compute_viterbi_sum(p_tree->find_descendant("A"), fam, &rec, 24, cache, &lm)._value);
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities_returns_invalid_if_equal_parent_and_child_sizes")
{
    matrix_cache cache(25);
    cache.precalculate_matrices({ 0.05 }, { 1,3,7 });
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree->find_descendant("AB")] = 11;
    single_lambda lm(0.05);
    CHECK_FALSE(compute_viterbi_sum(p_tree->find_descendant("A"), fam, &rec, 24, cache, &lm)._is_valid);
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities_returns_invalid_if_root")
{
    matrix_cache cache(25);
    cache.precalculate_matrices({ 0.05 }, { 1,3,7 });
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree.get()] = 11;
    single_lambda lm(0.05);
    CHECK_FALSE(compute_viterbi_sum(p_tree.get(), fam, &rec, 24, cache, &lm)._is_valid);
}

TEST_CASE_FIXTURE(Reconstruction, "pvalues")
{
    vector<double> cd(10);
    double n = 0;
    std::generate(cd.begin(), cd.end(), [&n]() mutable { return n += 0.01; });
    CHECK_EQ(0.5, pvalue(0.05, cd));
    CHECK_EQ(0.0, pvalue(0.0001, cd));
    CHECK_EQ(0.9, pvalue(0.099, cd));
}

TEST_CASE_FIXTURE(Reconstruction, "tree_pvalues")
{
    vector<vector<double>> cd(10);
    for (auto& d : cd)
    {
        d.resize(10);
        double n = 0;
        std::generate(d.begin(), d.end(), [&n]() mutable { return n += 0.01; });
    }
    clademap<vector<double>> results;
    results[p_tree.get()] = { 0, 0, 0 };
    auto fn = [this, &results](const clade* c)
    {
        if (c == p_tree.get())
            results[c][1] = 0.05;
    };
    CHECK_EQ(0.5, compute_tree_pvalue(p_tree.get(), fn, 10, cd, results));
}

TEST_CASE_FIXTURE(Reconstruction, "tree_pvalues_clears_results_before_using")
{
    vector<vector<double>> cd(10);
    for (auto& d : cd)
    {
        d.resize(10);
    }
    clademap<vector<double>> results;
    results[p_tree.get()] = { 1, 1, 1 };
    auto fn = [this, &results](const clade* c)
    {
        if (c == p_tree.get())
            results[c][1] = 0.05;
    };
    compute_tree_pvalue(p_tree.get(), fn, 10, cd, results);

    CHECK(vector<double>({ 0, 0.05, 0 }) == results[p_tree.get()]);
}

TEST_CASE_FIXTURE(Inference, "gamma_model_prune")
{
    vector<gene_family> families(1);
    families[0].set_species_size("A", 3);
    families[0].set_species_size("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    single_lambda lambda(0.005);

    _user_data.rootdist[1] = 2;
    _user_data.rootdist[2] = 2;
    _user_data.rootdist[3] = 2;
    _user_data.rootdist[4] = 2;
    _user_data.rootdist[5] = 1;
    root_equilibrium_distribution dist(_user_data.rootdist);
    matrix_cache cache(11);
    cache.precalculate_matrices({ 0.0005, 0.0025 }, set<double>{1, 3, 7});

    gamma_model model(&lambda, p_tree.get(), &families, 10, 8, { 0.01, 0.05 }, { 0.1, 0.5 }, NULL);

    vector<double> cat_likelihoods;
    CHECK(model.prune(families[0], dist, cache, &lambda, cat_likelihoods));

    CHECK_EQ(2, cat_likelihoods.size());
    CHECK_EQ(doctest::Approx(-23.04433), log(cat_likelihoods[0]));
    CHECK_EQ(doctest::Approx(-16.68005), log(cat_likelihoods[1]));
}

TEST_CASE_FIXTURE(Inference, "gamma_model_prune_returns_false_if_saturated")
{
    vector<gene_family> families(1);
    families[0].set_species_size("A", 3);
    families[0].set_species_size("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    single_lambda lambda(0.9);

    matrix_cache cache(11);
    cache.precalculate_matrices({ 0.09, 0.45 }, set<double>{1, 3, 7});
    vector<double> cat_likelihoods;

    gamma_model model(&lambda, p_tree.get(), &families, 10, 8, { 1.0,1.0 }, { 0.1, 0.5 }, NULL);

    CHECK(!model.prune(families[0], _user_data.prior, cache, &lambda, cat_likelihoods));
}

TEST_CASE("Inference: matrix_cache_key_handles_floating_point_imprecision")
{
    set<matrix_cache_key> keys;
    double t = 0.0;
    for (int i = 0; i < 31; i++)
    {
        t += 0.1;
        matrix_cache_key key(1, t, 0.3);
        keys.insert(key);
    }
    CHECK_EQ(31, keys.size());

    matrix_cache_key key(1, 3.0, 0.3);
    CHECK_EQ(1, keys.count(key));
}

TEST_CASE("Inference: birthdeath_rate_with_log_alpha")
{
    // alpha and coeff are derived values from lambda and t
    // (alpha = lambda*t / 1 + lambda*t, coeff = 1 - 2 * alpha);
    CHECK_EQ(doctest::Approx(-1.55455), log(birthdeath_rate_with_log_alpha(46, 45, -3.672556, 0.949177)));
    CHECK_EQ(doctest::Approx(-2.20436), log(birthdeath_rate_with_log_alpha(44, 46, -2.617970, 0.854098)));
    CHECK_EQ(doctest::Approx(-2.39974), log(birthdeath_rate_with_log_alpha(43, 43, -1.686354, 0.629613)));
    CHECK_EQ(doctest::Approx(-2.44301), log(birthdeath_rate_with_log_alpha(43, 44, -1.686354, 0.629613)));
    CHECK_EQ(doctest::Approx(-1.58253), log(birthdeath_rate_with_log_alpha(13, 14, -2.617970, 0.854098)));

    CHECK_EQ(doctest::Approx(0.107933), birthdeath_rate_with_log_alpha(40, 42, -1.37, 0.5));
    CHECK_EQ(doctest::Approx(0.005714), birthdeath_rate_with_log_alpha(41, 34, -1.262, 0.4));

    CHECK_EQ(doctest::Approx(0.194661), birthdeath_rate_with_log_alpha(5, 5, -1.1931291703283662, 0.39345841643135504));
}

TEST_CASE("Inference: create_one_model_if_lambda_is_null")
{
    input_parameters params;
    params.input_file_path = "foo";
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<base_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("Inference: create_gamma_model_if_alpha_provided")
{
    input_parameters params;
    params.input_file_path = "foo";
    params.fixed_alpha = 0.7;
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<gamma_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("Inference: create_gamma_model_if__n_gamma_cats__provided")
{
    input_parameters params;
    params.input_file_path = "foo";
    params.n_gamma_cats = 3;
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<gamma_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("Inference: build_models__uses_error_model_if_provided")
{
    input_parameters params;
    params.use_error_model = true;
    error_model em;
    em.set_probabilities(0, { 0, 0.99, 0.01 });
    user_data data;
    data.p_error_model = &em;
    data.p_tree = new clade("A", 5);
    data.max_family_size = 5;
    single_lambda lambda(0.05);
    data.p_lambda = &lambda;
    auto model = build_models(params, data)[0];
    std::ostringstream ost;
    model->write_vital_statistics(ost, 0.07);
    STRCMP_CONTAINS("Epsilon: 0.01\n", ost.str().c_str());
    delete model;
}

TEST_CASE("Inference: build_models__creates_default_error_model_if_needed")
{
    input_parameters params;
    params.use_error_model = true;
    user_data data;
    data.p_tree = new clade("A", 5);
    data.max_family_size = 5;
    single_lambda lambda(0.05);
    data.p_lambda = &lambda;
    auto model = build_models(params, data)[0];
    std::ostringstream ost;
    model->write_vital_statistics(ost, 0.01);
    STRCMP_CONTAINS("Epsilon: 0.05\n", ost.str().c_str());
    delete model;
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

TEST_CASE("Probability, matrix_multiply")
{
    matrix m1(3);
    build_matrix(m1);
    vector<double> m2({ 7, 9, 11 });
    auto result = m1.multiply(m2, 0, 2, 0, 2);
    CHECK_EQ(3, result.size());

    CHECK_EQ(58, result[0]);
    CHECK_EQ(139, result[1]);
    CHECK_EQ(220, result[2]);

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
    CHECK_EQ(3, result.size());

    CHECK_EQ(58, result[0]);
    CHECK_EQ(139, result[1]);
    CHECK_EQ(220, result[2]);
}

TEST_CASE("Probability: error_model_set_probs")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    auto vec = model.get_probs(0);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.0, vec[0]);
    CHECK_EQ(0.7, vec[1]);
    CHECK_EQ(0.3, vec[2]);

    vec = model.get_probs(1);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);
}

TEST_CASE("Probability: error_model_get_epsilon")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(2, { .1, .8, .1 });
    model.set_probabilities(3, { .2, .6, .2 });

    auto actual = model.get_epsilons();
    CHECK_EQ(3, actual.size());
    vector<double> expected{ .1, .2, .3 };
    CHECK(expected == actual);
}

TEST_CASE("Probability: error_model_get_epsilon_zero_zero_must_be_zero")
{
    error_model model;
    CHECK_THROWS_WITH_AS(model.set_probabilities(0, { 0.4, 0.3, 0.3 }), "Cannot have a non-zero probability for family size 0 for negative deviation", runtime_error);
}

TEST_CASE("Probability: error_model__set_probabilities__cannot_set_higher_values_without_setting_zero")
{
    error_model model;
    CHECK_THROWS_WITH_AS(model.set_probabilities(5, { 0.4, 0.3, 0.3 }), "Cannot have a non-zero probability for family size 0 for negative deviation", runtime_error);
    model.set_probabilities(0, { 0, 0.7, 0.3 });
    model.set_probabilities(5, { 0.4, 0.3, 0.3 });
    CHECK(model.get_probs(5) == vector<double>({ 0.4, 0.3, 0.3 }));
}

TEST_CASE("Probability: error_model__set_probabilities__can_set_higher_values_if_valid_for_zero")
{
    error_model model;
    model.set_probabilities(5, { 0, 0.7, 0.3 });
    CHECK(model.get_probs(5) == vector<double>({ 0, 0.7, 0.3 }));
}

TEST_CASE("Probability: error_model_rows_must_add_to_one")
{
    error_model model;
    model.set_probabilities(0, { 0,1,0 });
    CHECK(true);
    CHECK_THROWS_WITH_AS(model.set_probabilities(1, { 0.3, 0.3, 0.3 }), "Sum of probabilities must be equal to one", runtime_error);
}

TEST_CASE("Probability: error_model_replace_epsilons")
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
    CHECK_EQ(.7, actual[1]);
    CHECK_EQ(.3, actual[2]);

    actual = model.get_probs(1);
    CHECK_EQ(.3, actual[0]);
    CHECK_EQ(.4, actual[1]);
    CHECK_EQ(.3, actual[2]);
}

TEST_CASE("Probability: read_error_model")
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
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.0, vec[0]);
    CHECK_EQ(0.8, vec[1]);
    CHECK_EQ(0.2, vec[2]);

    vec = model.get_probs(1);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);

    vec = model.get_probs(4);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);

    vec = model.get_probs(7);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);
}

TEST_CASE("Probability: write_error_model")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(2, { .1, .8, .1 });
    model.set_probabilities(3, { .2, .6, .2 });

    ostringstream ost;
    write_error_model_file(ost, model);

    const char* expected = "maxcnt: 3\ncntdiff: -1 0 1\n"
        "0 0 0.7 0.3\n"
        "1 0.2 0.6 0.2\n"
        "2 0.1 0.8 0.1\n"
        "3 0.2 0.6 0.2\n";

    STRCMP_EQUAL(expected, ost.str().c_str());
}

TEST_CASE("Probability: write_error_model_skips_unnecessary_lines")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(10, { .1, .8, .1 });
    model.set_probabilities(15, { .05, .9, .05 });
    model.set_probabilities(18, { .05, .9, .05 });

    ostringstream ost;
    write_error_model_file(ost, model);

    const char* expected = "maxcnt: 18\ncntdiff: -1 0 1\n"
        "0 0 0.7 0.3\n"
        "1 0.2 0.6 0.2\n"
        "10 0.1 0.8 0.1\n"
        "15 0.05 0.9 0.05\n";

    STRCMP_EQUAL(expected, ost.str().c_str());
}

TEST_CASE("Probability: matrix_cache_warns_on_saturation")
{
    matrix_cache m(10);
    m.precalculate_matrices({ 0.05, 0.01 }, { 25 });
    ostringstream ost;
    m.warn_on_saturation(ost);
    STRCMP_EQUAL("WARNING: Saturated branch using lambda 0.05 on branch length 25\n", ost.str().c_str());
}

TEST_CASE("Probability, matrix_is_saturated")
{
    matrix_cache c(10);
    CHECK(c.is_saturated(25, 0.05));
    CHECK_FALSE(c.is_saturated(25, 0.01));
}

TEST_CASE("Inference: build_reference_list")
{
    std::string str = "Desc\tFamily ID\tA\tB\n"
        "\t (null)1\t5\t10\n"
        "\t (null)2\t5\t7\n"
        "\t (null)3\t5\t10\n"
        "\t (null)4\t5\t7\n";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, families);
    auto actual = build_reference_list(families);
    vector<int> expected({ 0, 1, 0, 1 });
    CHECK_EQ(expected.size(), actual.size());

}

TEST_CASE("Inference: prune")
{
    ostringstream ost;
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);
    matrix_cache cache(21);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });
    auto actual = inference_prune(fam, cache, &lambda, nullptr, p_tree.get(), 1.5, 20, 20);

    vector<double> log_expected{ -17.2771, -10.0323 , -5.0695 , -4.91426 , -5.86062 , -7.75163 , -10.7347 , -14.2334 , -18.0458 ,
        -22.073 , -26.2579 , -30.5639 , -34.9663 , -39.4472 , -43.9935 , -48.595 , -53.2439 , -57.9338 , -62.6597 , -67.4173 };

    CHECK_EQ(log_expected.size(), actual.size());
    for (size_t i = 0; i < log_expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(log_expected[i]), log(actual[i]));
    }
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_correctly")
{
    ostringstream ost;
    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    std::map<const clade*, std::vector<double> > _probabilities;

    auto init_func = [&](const clade* node) { _probabilities[node].resize(node->is_root() ? 20 : 21); };
    p_tree->apply_reverse_level_order(init_func);

    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, NULL, _probabilities, 20, 20, &lambda, cache);
    auto& actual = _probabilities[A];

    vector<double> expected{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_EQ(expected[i], actual[i]);
    }

    auto B = p_tree->find_descendant("B");
    compute_node_probability(B, family, NULL, _probabilities, 20, 20, &lambda, cache);
    actual = _probabilities[B];

    expected = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_EQ(expected[i], actual[i]);
    }
}

TEST_CASE("Inference: likelihood_computer_sets_root_nodes_correctly")
{
    ostringstream ost;
    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    std::map<const clade*, std::vector<double> > _probabilities;
    auto init_func = [&](const clade* node) { _probabilities[node].resize(node->is_root() ? 20 : 21); };
    p_tree->apply_reverse_level_order(init_func);

    cache.precalculate_matrices({ 0.03 }, { 1.0,3.0,7.0 });

    auto AB = p_tree->find_descendant("AB");
    compute_node_probability(p_tree->find_descendant("A"), family, NULL, _probabilities, 20, 20, &lambda, cache);
    compute_node_probability(p_tree->find_descendant("B"), family, NULL, _probabilities, 20, 20, &lambda, cache);
    compute_node_probability(AB, family, NULL, _probabilities, 20, 20, &lambda, cache);

    auto& actual = _probabilities[AB];

    vector<double> log_expected{ -19.7743, -11.6688, -5.85672, -5.66748, -6.61256, -8.59725, -12.2301, -16.4424, -20.9882, -25.7574,
        -30.6888, -35.7439, -40.8971, -46.1299, -51.4289, -56.7837, -62.1863, -67.6304, -73.1106, -78.6228
    };

    CHECK_EQ(log_expected.size(), actual.size());
    for (size_t i = 0; i < log_expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(log_expected[i]), log(actual[i]));
    }
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_from_error_model_if_provided")
{
    ostringstream ost;

    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });

    string input = "maxcnt: 20\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "20 0.2 0.6 0.2\n";
    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    std::map<const clade*, std::vector<double> > _probabilities;
    auto init_func = [&](const clade* node) { _probabilities[node].resize(node->is_root() ? 20 : 21); };
    p_tree->apply_reverse_level_order(init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, &model, _probabilities, 20, 20, &lambda, cache);
    auto& actual = _probabilities[A];

    vector<double> expected{ 0, 0, 0.2, 0.6, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        //cout << actual[i] << endl;
        CHECK_EQ(expected[i], actual[i]);
    }
}

TEST_CASE("Clade: get_lambda_index_throws_from_branch_length_tree")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    CHECK_EQ(7, p_tree->get_branch_length());
    CHECK_THROWS_WITH_AS(p_tree->get_lambda_index(), "Requested lambda index from branch length tree", runtime_error);

}

TEST_CASE("Clade: get_branch_length_throws_from_lambda_tree")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7", true));
    CHECK_EQ(7, p_tree->get_lambda_index());
    CHECK_THROWS_WITH_AS(p_tree->get_branch_length(), "Requested branch length from lambda tree", runtime_error);

}

TEST_CASE("Clade: lambda_tree_root_index_is_1_if_not_specified")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:2)", true));
    CHECK_EQ(1, p_tree->get_lambda_index());
}

TEST_CASE("Clade: exists_at_root_returns_false_if_not_all_children_exist")
{
    unique_ptr<clade> p_tree(parse_newick(" ((((cat:68.710687,horse:68.710687):4.566771,cow:73.277458):20.722542,(((((chimp:4.444178,human:4.444178):6.682660,orang:11.126837):2.285866,gibbon:13.412704):7.211528,(macaque:4.567239,baboon:4.567239):16.056993):16.060691,marmoset:36.684923):57.315077)mancat:38.738115,(rat:36.302467,mouse:36.302467):96.435648)"));

    istringstream ist(
        "Desc\tFamily ID\tcat\thorse\tcow\tchimp\thuman\torang\tgibbon\tmacaque\tbaboon\tmarmoset\trat\tmouse\n"
        "(null)\t1\t0\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0\n"
        "(null)\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\n");

    vector<gene_family> families;
    read_gene_families(ist, p_tree.get(), families);
    CHECK_FALSE(families[0].exists_at_root(p_tree.get()));
    CHECK_FALSE(families[1].exists_at_root(p_tree.get()));
}

TEST_CASE("Clade: exists_at_root_returns_true_if_all_children_exist")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    CHECK(family.exists_at_root(p_tree.get()));
}

TEST_CASE("Clade: parse_newick_throws_exception_for_invalid_lambdas_in_tree")
{
    CHECK_THROWS_WITH_AS(parse_newick("(A:1,B:0):2", true), "Invalid lambda index set for B", runtime_error);
    CHECK_THROWS_WITH_AS(parse_newick("(A:-1,B:2)", true), "Invalid lambda index set for A", runtime_error);
}

TEST_CASE("Clade: parse_newick_throws_exception_for_invalid_branch_length_in_tree")
{
    CHECK_THROWS_WITH_AS(parse_newick("(A:1,B:0):2", false), "Invalid branch length set for B", runtime_error);
    CHECK_THROWS_WITH_AS(parse_newick("(A:-1,B:2)", false), "Invalid branch length set for A", runtime_error);
}

TEST_CASE("Clade: copy_constructor")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    unique_ptr<clade> copy(new clade(*p_tree.get()));
    CHECK_EQ(1, copy->find_descendant("A")->get_branch_length());
    CHECK_EQ(7, copy->find_descendant("AB")->get_branch_length());
}

TEST_CASE("Clade: copy_constructor_modifying_branch")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    unique_ptr<clade> copy(new clade(*p_tree.get(), nullptr, [](const clade& c) { return c.get_branch_length() * 2; }));
    CHECK_EQ(2, copy->find_descendant("A")->get_branch_length());
    CHECK_EQ(14, copy->find_descendant("AB")->get_branch_length());
}

TEST_CASE("Inference: multiple_lambda_returns_correct_values")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    map<string, int> key;
    key["A"] = 5;
    key["B"] = 3;
    multiple_lambda ml(key, { .03, .05, .07, .011, .013, .017 });
    CHECK_EQ(.017, ml.get_value_for_clade(p_tree->find_descendant("A")));
    CHECK_EQ(.011, ml.get_value_for_clade(p_tree->find_descendant("B")));
}

TEST_CASE("Simulation: uniform_distribution__select_root_size__returns_sequential_values")
{
    root_equilibrium_distribution ud(20);
    CHECK_EQ(1, ud.select_root_size(0));
    CHECK_EQ(2, ud.select_root_size(1));
    CHECK_EQ(3, ud.select_root_size(2));
    CHECK_EQ(4, ud.select_root_size(3));
    CHECK_EQ(5, ud.select_root_size(4));
    CHECK_EQ(0, ud.select_root_size(20));
}

TEST_CASE("Simulation: specified_distribution__select_root_size__returns_exact_selection")
{
    user_data data;
    for (int i = 0; i < 20; ++i)
        data.rootdist[i] = 1;

    root_equilibrium_distribution sd(data.rootdist);
    for (size_t i = 0; i < 20; ++i)
        CHECK_EQ(i, sd.select_root_size(i));
}

TEST_CASE("Simulation: print_process_prints_in_order")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    std::ostringstream ost;
    clademap<int> t;
    t[p_tree->find_descendant("B")] = 4;
    t[p_tree->find_descendant("A")] = 2;
    t[p_tree->find_descendant("AB")] = 6;

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, true, my_trials);

    STRCMP_CONTAINS("DESC\tFID\tB\tA\t2", ost.str().c_str());
    STRCMP_CONTAINS("L0\tsimfam0\t4\t2\t6", ost.str().c_str());

}

TEST_CASE("Simulation: print_process_can_print_without_internal_nodes")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    std::ostringstream ost;
    clademap<int> t;
    t[p_tree->find_descendant("B")] = 4;
    t[p_tree->find_descendant("A")] = 2;
    t[p_tree->find_descendant("AB")] = 6;

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, false, my_trials);
    STRCMP_CONTAINS("DESC\tFID\tB\tA\n", ost.str().c_str());
    STRCMP_CONTAINS("L0\tsimfam0\t4\t2\n", ost.str().c_str());

}

TEST_CASE("Simulation: gamma_model_get_simulation_lambda_selects_random_multiplier_based_on_alpha")
{
    single_lambda lam(0.05);
    gamma_model m(&lam, NULL, NULL, 0, 5, 3, 0.7, NULL);
    unique_ptr<single_lambda> new_lam(dynamic_cast<single_lambda*>(m.get_simulation_lambda()));
    CHECK_EQ(doctest::Approx(0.00574028), new_lam->get_single_lambda());
}

TEST_CASE("Simulation: create_trial")
{
    randomizer_engine.seed(10);

    single_lambda lam(0.25);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    user_data data;
    data.p_tree = p_tree.get();
    data.rootdist[1] = 1;
    data.rootdist[2] = 1;
    data.rootdist[5] = 1;
    data.max_family_size = 10;
    data.max_root_family_size = 10;

    data.prior = root_equilibrium_distribution(data.rootdist);
    input_parameters params;
    simulator sim(data, params);

    matrix_cache cache(100);
    cache.precalculate_matrices(get_lambda_values(&lam), { 1,3,7 });

    simulated_family actual = sim.create_trial(&lam, 2, cache);

    CHECK_EQ(5, actual.values.at(p_tree.get()));
    CHECK_EQ(4, actual.values.at(p_tree->find_descendant("A")));
    CHECK_EQ(4, actual.values.at(p_tree->find_descendant("B")));
}

TEST_CASE("Inference: model_vitals")
{
    mock_model model;
    model.set_tree(new clade("A", 5));
    single_lambda lambda(0.05);
    model.set_lambda(&lambda);
    std::ostringstream ost;
    model.write_vital_statistics(ost, 0.01);
    STRCMP_CONTAINS("Model mockmodel Final Likelihood (-lnL): 0.01", ost.str().c_str());
    STRCMP_CONTAINS("Lambda:            0.05", ost.str().c_str());
    STRCMP_CONTAINS("Maximum possible lambda for this topology: 0.2", ost.str().c_str());
    STRCMP_CONTAINS("No attempts made", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gene_family_reconstrctor__print_increases_decreases_by_family__adds_flag_for_significance")
{
    ostringstream insignificant;
    base_model_reconstruction bmr;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 5;
    gene_family gf;
    gf.set_id("myid");
    gf.set_species_size("A", 7);
    order.clear();
    bmr.print_increases_decreases_by_family(insignificant, order, { gf }, { 0.03 }, 0.01);
    STRCMP_CONTAINS("myid\t0.03\tn", insignificant.str().c_str());

    ostringstream significant;
    bmr.print_increases_decreases_by_family(insignificant, order, { gf }, { 0.03 }, 0.05);
    STRCMP_CONTAINS("myid\t0.03\ty", insignificant.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_increases_decreases_by_family__prints_significance_level_in_header")
{
    ostringstream ost;
    base_model_reconstruction bmr;
    gene_family gf;

    bmr.print_increases_decreases_by_family(ost, order, { gf }, { 0.07 }, 0.00001);
    STRCMP_CONTAINS("#FamilyID\tpvalue\tSignificant at 1e-05\n", ost.str().c_str());
}

TEST_CASE("Reconstruction: base_model_print_increases_decreases_by_family")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;
    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    base_model_reconstruction bmr;
    bmr.print_increases_decreases_by_family(empty, order, {}, {}, 0.05);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 5;

    gene_family gf;
    gf.set_id("myid");
    gf.set_species_size("A", 7);
    gf.set_species_size("B", 2);

    ostringstream ost;
    bmr.print_increases_decreases_by_family(ost, order, { gf }, { 0.07 }, 0.05);
    STRCMP_CONTAINS("#FamilyID\tpvalue\tSignificant at 0.05\n", ost.str().c_str());
    STRCMP_CONTAINS("myid\t0.07\tn\n", ost.str().c_str());
}

TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_family")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;
    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    vector<double> multipliers({ .2, .75 });
    vector<gamma_bundle*> bundles; //  ({ &bundle });
    vector<double> em;
    gamma_model_reconstruction gmr(em);
    gmr.print_increases_decreases_by_family(empty, order, {}, {}, 0.05);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    gmr._reconstructions["myid"].reconstruction[p_tree->find_descendant("AB")] = 5;

    gene_family gf;
    gf.set_id("myid");
    gf.set_species_size("A", 7);
    gf.set_species_size("B", 2);

    ostringstream ost;
    gmr.print_increases_decreases_by_family(ost, order, { gf }, { 0.07 }, 0.05);
    STRCMP_CONTAINS("#FamilyID\tpvalue\tSignificant at 0.05", ost.str().c_str());
    STRCMP_CONTAINS("myid\t0.07\tn", ost.str().c_str());
}

TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    ostringstream empty;

    vector<double> multipliers({ .2, .75 });
    vector<gamma_bundle*> bundles; //  ({ &bundle });
    vector<double> em;
    gamma_model_reconstruction gmr(em);

    gmr.print_increases_decreases_by_clade(empty, order, {});
    STRCMP_EQUAL("#Taxon_ID\tIncrease\tDecrease\n", empty.str().c_str());

    gmr._reconstructions["myid"].reconstruction[p_tree->find_descendant("AB")] = 5;

    gene_family gf;
    gf.set_id("myid");
    gf.set_species_size("A", 7);
    gf.set_species_size("B", 2);

    ostringstream ost;
    gmr.print_increases_decreases_by_clade(ost, order, { gf });
    STRCMP_CONTAINS("#Taxon_ID\tIncrease\tDecrease", ost.str().c_str());
    STRCMP_CONTAINS("A<0>\t1\t0", ost.str().c_str());
    STRCMP_CONTAINS("B<1>\t0\t1", ost.str().c_str());
}

TEST_CASE("Reconstruction: base_model_print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    clade invalid;
    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    ostringstream empty;

    base_model_reconstruction bmr;

    bmr.print_increases_decreases_by_clade(empty, order, {});
    STRCMP_EQUAL("#Taxon_ID\tIncrease\tDecrease\n", empty.str().c_str());

    //    bmr._reconstructions["myid"][p_tree->find_descendant("A")] = 4;
    //    bmr._reconstructions["myid"][p_tree->find_descendant("B")] = -3;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 5;

    gene_family gf;
    gf.set_id("myid");
    gf.set_species_size("A", 7);
    gf.set_species_size("B", 2);

    ostringstream ost;
    bmr.print_increases_decreases_by_clade(ost, order, { gf });
    STRCMP_CONTAINS("#Taxon_ID\tIncrease\tDecrease", ost.str().c_str());
    STRCMP_CONTAINS("A<0>\t1\t0", ost.str().c_str());
    STRCMP_CONTAINS("B<1>\t0\t1", ost.str().c_str());
}

TEST_CASE("Inference: lambda_epsilon_optimizer")
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

TEST_CASE("Inference: lambda_per_family")
{
    randomizer_engine.seed(10);
    user_data ud;

    ud.max_root_family_size = 10;
    ud.max_family_size = 10;
    ud.gene_families.resize(1);
    ud.prior = root_equilibrium_distribution(ud.max_root_family_size);

    gene_family& family = ud.gene_families[0];
    family.set_id("test");
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);
    input_parameters params;
    params.lambda_per_family = true;

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    ud.p_tree = p_tree.get();
    estimator v(ud, params);

    mock_model m;
    m.set_tree(ud.p_tree);
    ostringstream ost;
    v.estimate_lambda_per_family(&m, ost);
    CHECK_EQ(std::string("test\t0.00053526736161992\n"), ost.str());
}

TEST_CASE_FIXTURE(Inference, "estimator_compute_pvalues")
{
    input_parameters params;
    matrix_cache cache(max(_user_data.max_family_size, _user_data.max_root_family_size) + 1);
    cache.precalculate_matrices(get_lambda_values(_user_data.p_lambda), _user_data.p_tree->get_branch_lengths());

    auto values = compute_pvalues(_user_data.p_tree, _user_data.gene_families, _user_data.p_lambda, cache, 3, _user_data.max_family_size, _user_data.max_root_family_size);
    CHECK_EQ(1, values.size());
    CHECK_EQ(doctest::Approx(0.666667), values[0]);
}

TEST_CASE_FIXTURE(Inference, "gamma_lambda_optimizer")
{
    _user_data.max_root_family_size = 10;
    _user_data.prior = root_equilibrium_distribution(_user_data.max_root_family_size);

    gamma_model m(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 10, _user_data.max_root_family_size, 4, 0.25, NULL);
    gamma_lambda_optimizer optimizer(_user_data.p_lambda, &m, &_user_data.prior, 7);
    vector<double> values{ 0.01, 0.25 };
    CHECK_EQ(doctest::Approx(6.4168), optimizer.calculate_score(&values[0]));
}

TEST_CASE("Inference: inference_optimizer_scorer__calculate_score__translates_nan_to_inf")
{
    single_lambda lam(0.05);
    mock_model m;
    m.set_invalid_likelihood();
    double val;
    lambda_optimizer opt(&lam, &m, NULL, 0);
    CHECK(std::isinf(opt.calculate_score(&val)));
}

TEST_CASE_FIXTURE(Inference, "poisson_scorer_optimizes_correct_value")
{
    poisson_scorer scorer(_user_data.gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    // DOUBLES_EQUAL(0.5, result.values[0], 0.0001)
}

TEST_CASE_FIXTURE(Inference, "poisson_scorer__lnlPoisson")
{
    poisson_scorer scorer(_user_data.gene_families);
    double lambda = 0.05;
    CHECK_EQ(doctest::Approx(3.095732), scorer.lnLPoisson(&lambda));
}

TEST_CASE_FIXTURE(Inference, "poisson_scorer__lnlPoisson_skips_incalculable_family_sizes")
{
    _user_data.gene_families.resize(2);
    _user_data.gene_families[1].set_id("TestFamily2");
    _user_data.gene_families[1].set_species_size("A", 3);
    _user_data.gene_families[1].set_species_size("B", 175);

    poisson_scorer scorer(_user_data.gene_families);
    double lambda = 0.05;
    CHECK_EQ(doctest::Approx(9.830344), scorer.lnLPoisson(&lambda));
}

class mock_scorer : public optimizer_scorer
{
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        return std::vector<double>{0.2};
    }
    virtual double calculate_score(const double* values) override
    {
        if (force_scoring_error)
            return std::numeric_limits<double>::infinity();

        return 1000;
    }
public:
    bool force_scoring_error = false;
};

TEST_CASE("Inference, optimizer_gets_initial_guesses_from_scorer")
{
    mock_scorer scorer;
    optimizer opt(&scorer);
    auto guesses = opt.get_initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(0.2, guesses[0]);
}

TEST_CASE("Inference, optimizer_disallows_bad_initializations")
{
    mock_scorer scorer;
    scorer.force_scoring_error = true;
    optimizer opt(&scorer);

    CHECK_THROWS_WITH_AS(opt.get_initial_guesses(), "Failed to initialize any reasonable values", runtime_error);
}

TEST_CASE("Inference, poisson_distribution__compute")
{
    root_equilibrium_distribution pd(0.75, 100);

    CHECK_EQ(doctest::Approx(0.2436548).scale(1000), pd.compute(1));
    CHECK_EQ(doctest::Approx(0.071).scale(1000), pd.compute(3));
    CHECK_EQ(doctest::Approx(0.005).scale(1000), pd.compute(5));
    CHECK_EQ(0.0, pd.compute(100));
}

TEST_CASE("Inference, poisson_distribution__select_root_size")
{
    root_equilibrium_distribution pd(0.75, 9);

    CHECK_EQ(1, pd.select_root_size(1));
    CHECK_EQ(1, pd.select_root_size(3));
    CHECK_EQ(2, pd.select_root_size(5));
    CHECK_EQ(2, pd.select_root_size(7));
    CHECK_EQ(0, pd.select_root_size(100));
}

TEST_CASE("Inference, optimizer_result_stream")
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

TEST_CASE("Inference, event_monitor_shows_no_attempts")
{
    event_monitor evm;

    ostringstream ost;
    evm.log(ost);
    STRCMP_EQUAL("No attempts made\n", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_one_attempt")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    ostringstream ost;

    evm.log(ost);
    STRCMP_EQUAL("1 values were attempted (0% rejected)\n", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_rejected_attempts")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_InvalidValues();
    ostringstream ost;

    evm.log(ost);
    STRCMP_EQUAL("2 values were attempted (50% rejected)\n", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_poor_performing_families")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Saturation("test");
    ostringstream ost;

    evm.log(ost);
    STRCMP_CONTAINS("2 values were attempted (0% rejected)", ost.str().c_str());
    STRCMP_CONTAINS("The following families had failure rates >20% of the time:", ost.str().c_str());
    STRCMP_CONTAINS("test had 1 failures", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_does_not_show_decent_performing_families")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Saturation("test");
    ostringstream ost;

    evm.log(ost);
    STRCMP_EQUAL("5 values were attempted (0% rejected)\n", ost.str().c_str());
}

TEST_CASE_FIXTURE(Inference, "initialization_failure_advice_shows_20_families_with_largest_differentials")
{
    std::ostringstream ost;
    _user_data.gene_families.resize(2);
    _user_data.gene_families[1].set_id("TestFamily2");
    _user_data.gene_families[1].set_species_size("A", 34);
    _user_data.gene_families[1].set_species_size("B", 86);

    initialization_failure_advice(ost, _user_data.gene_families);
    STRCMP_CONTAINS("Families with largest size differentials:", ost.str().c_str());
    STRCMP_CONTAINS("\nYou may want to try removing the top few families with the largest difference\nbetween the max and min counts and then re-run the analysis.\n", ost.str().c_str());
    STRCMP_CONTAINS("TestFamily2: 52\nTestFamily1: 1", ost.str().c_str());
}

TEST_CASE("Simulation, base_prepare_matrices_for_simulation_creates_matrix_for_each_branch")
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    base_model b(&lam, p_tree.get(), NULL, 0, 0, NULL);
    matrix_cache m(25);
    b.prepare_matrices_for_simulation(m);
    CHECK_EQ(3, m.get_cache_size());
}

TEST_CASE("Simulation, base_prepare_matrices_for_simulation_uses_perturbed_lambda")
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    base_model b(&lam, p_tree.get(), NULL, 0, 0, NULL);
    b.perturb_lambda();
    matrix_cache m(25);
    b.prepare_matrices_for_simulation(m);
    unique_ptr<single_lambda> sim_lambda(dynamic_cast<single_lambda*>(b.get_simulation_lambda()));
    try
    {
        m.get_matrix(7, sim_lambda->get_single_lambda());
    }
    catch (std::runtime_error & err)
    {
        FAIL("Failed to cache simulation lambda");
    }
}

TEST_CASE("Simulation, gamma_prepare_matrices_for_simulation_creates_matrix_for_each_branch_and_category")
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    gamma_model g(&lam, p_tree.get(), NULL, 0, 0, 2, 0.5, NULL);
    matrix_cache m(25);
    g.prepare_matrices_for_simulation(m);
    CHECK_EQ(6, m.get_cache_size());
}

TEST_CASE("Simulation, set_random_node_size_without_error_model")
{
    randomizer_engine.seed(10);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.05);
    clademap<int> t;
    matrix_cache cache(12);
    cache.precalculate_matrices({ 0.05 }, { 3 });
    auto b = p_tree->find_descendant("B");
    set_weighted_random_family_size(b, &t, &lambda, NULL, 10, cache);

    CHECK_EQ(0, t[b]);

    t[p_tree.get()] = 5;

    set_weighted_random_family_size(b, &t, &lambda, NULL, 10, cache);
    CHECK_EQ(4, t[b]);
}

TEST_CASE("Simulation, set_random_node_size_with_error_model")
{
    randomizer_engine.seed(10);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.05);
    clademap<int> t;
    matrix_cache cache(20);
    cache.precalculate_matrices({ 0.05 }, { 3 });
    auto b = p_tree->find_descendant("B");

    error_model err;
    err.set_probabilities(0, { 0, .95, 0.05 });

    t[p_tree.get()] = 5;
    err.set_probabilities(5, { .9, .05, 0.05 });
    set_weighted_random_family_size(b, &t, &lambda, &err, 10, cache);
    CHECK_EQ(4, t[b]);

    err.set_probabilities(5, { .05, .05, 0.9 });
    set_weighted_random_family_size(b, &t, &lambda, &err, 10, cache);
    CHECK_EQ(6, t[b]);
}

TEST_CASE("Simulation, executor")
{
    randomizer_engine.seed(10);

    input_parameters params;
    user_data ud;
    unique_ptr<action> act(get_executor(params, ud));
    CHECK(dynamic_cast<const estimator*>(act.get()));

    params.chisquare_compare = true;
    unique_ptr<action> act2(get_executor(params, ud));
    CHECK(dynamic_cast<const chisquare_compare*>(act2.get()));
}

TEST_CASE("Simulation, specified_distribution__with_rootdist_creates_matching_vector")
{
    std::map<int, int> m;
    m[2] = 3;
    m[4] = 1;
    m[8] = 1;
    root_equilibrium_distribution rd(m);
    CHECK_EQ(rd.select_root_size(0), 2);
    CHECK_EQ(rd.select_root_size(1), 2);
    CHECK_EQ(rd.select_root_size(2), 2);
    CHECK_EQ(rd.select_root_size(3), 4);
    CHECK_EQ(rd.select_root_size(4), 8);
    CHECK_EQ(rd.select_root_size(5), 0);
}

TEST_CASE("Simulation, specified_distribution__pare")
{
    randomizer_engine.seed(10);

    std::map<int, int> m;
    m[2] = 5;
    m[4] = 3;
    m[8] = 3;
    root_equilibrium_distribution rd(m);
    rd.resize(5);
    CHECK_EQ(rd.select_root_size(0), 2);
    CHECK_EQ(rd.select_root_size(1), 2);
    CHECK_EQ(rd.select_root_size(2), 2);
    CHECK_EQ(rd.select_root_size(3), 4);
    CHECK_EQ(rd.select_root_size(4), 8);
    CHECK_EQ(rd.select_root_size(5), 0);
}

TEST_CASE("Simulation, specified_distribution__expand")
{
    std::map<int, int> m;
    m[2] = 5;
    m[4] = 3;
    m[8] = 3;
    root_equilibrium_distribution rd(m);
    rd.resize(15);
    CHECK_EQ(rd.select_root_size(14), 8);
    CHECK_EQ(rd.select_root_size(15), 0);
}


TEST_CASE("Simulation, simulate_processes")
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    mock_model m;
    m.set_tree(p_tree.get());
    m.set_lambda(&lam);

    user_data ud;
    ud.p_tree = p_tree.get();
    ud.p_lambda = &lam;
    ud.prior = root_equilibrium_distribution(100);
    ud.max_family_size = 101;
    ud.max_root_family_size = 101;

    input_parameters ip;
    ip.nsims = 100;
    simulator sim(ud, ip);
    vector<simulated_family> results(1);
    sim.simulate_processes(&m, results);
    CHECK_EQ(100, results.size());
}

TEST_CASE("Simulation, simulate_processes_uses_rootdist_if_available")
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    mock_model m;
    m.set_tree(p_tree.get());
    m.set_lambda(&lam);

    user_data ud;
    ud.p_tree = p_tree.get();
    ud.p_lambda = &lam;
    ud.rootdist[5] = 50;
    ud.rootdist[10] = 50;
    ud.max_family_size = 60;
    ud.max_root_family_size = 60;
    ud.prior = root_equilibrium_distribution(ud.rootdist);
    input_parameters ip;
    ip.nsims = 100;
    simulator sim(ud, ip);
    vector<simulated_family> results;
    sim.simulate_processes(&m, results);
    CHECK_EQ(100, results.size());

    CHECK_EQ(5, results[0].values.at(p_tree.get()));
    CHECK_EQ(10, results[75].values.at(p_tree.get()));
}

TEST_CASE("Simulation, gamma_model_perturb_lambda_with_clusters")
{
    randomizer_engine.seed(10);

    gamma_model model(NULL, NULL, NULL, 0, 5, 3, 0.7, NULL);
    model.perturb_lambda();
    auto multipliers = model.get_lambda_multipliers();
    CHECK_EQ(3, multipliers.size());
    CHECK_EQ(doctest::Approx(0.1136284), multipliers[0]);
    CHECK_EQ(doctest::Approx(0.8763229), multipliers[1]);
    CHECK_EQ(doctest::Approx(2.151219), multipliers[2]);
}

TEST_CASE("Simulation, gamma_model_perturb_lambda_without_clusters")
{
    randomizer_engine.seed(10);

    gamma_model model(NULL, NULL, NULL, 0, 5, 1, 5, NULL);
    model.perturb_lambda();
    auto multipliers = model.get_lambda_multipliers();
    CHECK_EQ(1, multipliers.size());
    CHECK_EQ(doctest::Approx(0.911359), multipliers[0]);
}

TEST_CASE("create_root_distribution__creates__uniform_distribution")
{
    input_parameters params;
    map<int, int> rootdist;
    auto dist = create_root_distribution(params, nullptr, rootdist, 10);
    CHECK_EQ(doctest::Approx(.1), dist.compute(1));
    CHECK_EQ(doctest::Approx(.1), dist.compute(9));
    CHECK_EQ(0, dist.compute(11));
}

TEST_CASE("Simulation, create_root_distribution__creates__poisson_distribution_if_given")
{
    input_parameters params;
    params.use_uniform_eq_freq = false;
    params.poisson_lambda = 0.75;
    map<int, int> rootdist;
    auto dist = create_root_distribution(params, nullptr, rootdist, 100);
    CHECK_EQ(doctest::Approx(0.2436548), dist.compute(1));
    CHECK_EQ(0, dist.compute(101));
}

TEST_CASE("Simulation, create_root_distribution__creates__specifed_distribution_if_given")
{
    input_parameters params;
    map<int, int> rootdist;
    rootdist[2] = 11;
    rootdist[3] = 5;
    rootdist[4] = 7;
    rootdist[6] = 2;
    auto dist = create_root_distribution(params, nullptr, rootdist, 10);
    CHECK_EQ(doctest::Approx(0.44), dist.compute(2));
    CHECK_EQ(doctest::Approx(0.2), dist.compute(3));
    CHECK_EQ(doctest::Approx(0.28), dist.compute(4));
    CHECK_EQ(0, dist.compute(5));
    CHECK_EQ(doctest::Approx(.08), dist.compute(6));
}


TEST_CASE("Simulation, create_root_distribution__creates__specifed_distribution_if_given_distribution_and_poisson")
{
    input_parameters params;
    params.use_uniform_eq_freq = false;
    params.poisson_lambda = 0.75;
    auto dist = create_root_distribution(params, nullptr, map<int, int>(), 10);
    CHECK_EQ(1, dist.select_root_size(0));
    CHECK_EQ(1, dist.select_root_size(1));
    CHECK_EQ(1, dist.select_root_size(2));
    CHECK_EQ(1, dist.select_root_size(3));
    CHECK_EQ(1, dist.select_root_size(4));
    CHECK_EQ(2, dist.select_root_size(5));
    CHECK_EQ(2, dist.select_root_size(6));
    CHECK_EQ(2, dist.select_root_size(7));
    CHECK_EQ(2, dist.select_root_size(8));
}

TEST_CASE("Simulation, create_root_distribution__resizes_distribution_if_nsims_specified")
{
    randomizer_engine.seed(10);

    input_parameters params;
    params.nsims = 10;
    auto dist = create_root_distribution(params, nullptr, map<int, int>(), 100);
    /// GCC's implementation of shuffle changed so the numbers that are
    /// returned are in a slightly different order, even with the same seed
#if __GNUC__ >= 7
    CHECK_EQ(86, dist.select_root_size(9));
#else
    CHECK_EQ(81, dist.select_root_size(9));
#endif
    CHECK_EQ(0, dist.select_root_size(10));
}

TEST_CASE_FIXTURE(Optimizer, "fminsearch_sort_sorts_scores_and_moves_values")
{
    c[0].score = 3;
    c[1].score = 5;
    c[2].score = 1;

    __fminsearch_sort(&fm);

    CHECK_EQ(&c[2], fm.candidates[0]);
    CHECK_EQ(&c[0], fm.candidates[1]);
    CHECK_EQ(&c[1], fm.candidates[2]);
}

TEST_CASE_FIXTURE(Optimizer, "fminsearch_checkV_compares_value_difference_to_lx")
{
    c[0].values[0] = 1;
    c[1].values[0] = 2;
    c[2].values[0] = 3;
    c[0].values[1] = 3;
    c[1].values[1] = 4;
    c[2].values[1] = 5;
    fm.tolx = 3;
    CHECK(__fminsearch_checkV(&fm));
    fm.tolx = .5;
    CHECK_FALSE(__fminsearch_checkV(&fm));
}

TEST_CASE_FIXTURE(Optimizer, "fminsearch_checkF_compares_score_difference_to_lf")
{
    c[0].score = 1.0;
    c[1].score = 3.0;
    c[2].score = 5.0;

    fm.tolf = 5;
    CHECK(__fminsearch_checkF(&fm));
    fm.tolf = 1;
    CHECK_FALSE(__fminsearch_checkF(&fm));
}

class multiplier_scorer : public optimizer_scorer
{
    int num_scores;
public:
    multiplier_scorer() : multiplier_scorer(2) {}
    multiplier_scorer(int sz) : num_scores(sz)
    {

    }
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        vector<double> result(num_scores);
        result[0] = 5;
        result[1] = 3;
        return result;
    }
    virtual double calculate_score(const double* values) override
    {
        return std::accumulate(values, values + num_scores, 1.0, std::multiplies<double>());
    }
};

TEST_CASE_FIXTURE(Optimizer, "fminsearch_min_init")
{
    fm.delta = 0.05;
    fm.zero_delta = 0.00025;
    vector<int> indices(3);
    fm.idx = &indices[0];
    c[0].values[0] = 300;
    c[1].values[0] = 200;
    c[2].values[0] = 100;

    multiplier_scorer ms;
    fm.scorer = &ms;
    auto init = ms.initial_guesses();
    __fminsearch_min_init(&fm, &init[0]);

    CHECK_EQ(15, fm.candidates[0]->score);
    CHECK_EQ(15.75, fm.candidates[1]->score);
    CHECK_EQ(doctest::Approx(15.75), fm.candidates[2]->score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_mean")
{
    fm.variable_count = 2;
    vector<double> means(2);
    fm.x_mean = &means[0];
    fm.candidates[0]->values[0] = 300;
    fm.candidates[1]->values[0] = 200;
    fm.candidates[0]->values[1] = 12;
    fm.candidates[1]->values[1] = 44;

    __fminsearch_x_mean(&fm);

    CHECK_EQ(250, fm.x_mean[0]);
    CHECK_EQ(28, fm.x_mean[1]);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_reflection")
{
    fm.rho = 1;
    vector<double> means({ 250,28 });
    vector<double> reflections(2);
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.candidates[0]->values[0] = 300;
    fm.candidates[1]->values[0] = 200;
    fm.candidates[0]->values[1] = 12;
    fm.candidates[1]->values[1] = 44;

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_reflection(&fm);

    CHECK_EQ(500, fm.x_r[0]);
    CHECK_EQ(56, fm.x_r[1]);
    CHECK_EQ(28000, score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_expansion")
{
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

    CHECK_EQ(750, fm.x_tmp[0]);
    CHECK_EQ(84, fm.x_tmp[1]);
    CHECK_EQ(63000, score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_contract_outside")
{
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

    CHECK_EQ(375, fm.x_tmp[0]);
    CHECK_EQ(42, fm.x_tmp[1]);
    CHECK_EQ(15750, score);
}


TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_contract_inside")
{
    fm.psi = 0.5;
    vector<double> means({ 250,28 });
    vector<double> expansions(2);
    fm.candidates[2]->values[0] = 26;
    fm.candidates[2]->values[1] = 12;

    fm.x_mean = &means[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_contract_inside(&fm);

    CHECK_EQ(362, fm.x_tmp[0]);
    CHECK_EQ(36, fm.x_tmp[1]);
    CHECK_EQ(13032, score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_shrink")
{
    fm.variable_count_plus_one = 3;
    fm.sigma = 0.5;
    vector<int> indices(3);
    fm.idx = &indices[0];
    vector<double> scores(3);

    fm.candidates[0]->values[0] = 300;
    fm.candidates[0]->values[1] = 200;
    fm.candidates[1]->values[0] = 42;
    fm.candidates[1]->values[1] = 64;
    fm.candidates[2]->values[0] = 26;
    fm.candidates[2]->values[1] = 12;

    multiplier_scorer ms;
    fm.scorer = &ms;
    __fminsearch_x_shrink(&fm);

    CHECK_EQ(300, fm.candidates[0]->values[0]);
    CHECK_EQ(200, fm.candidates[0]->values[1]);
    CHECK_EQ(163, fm.candidates[1]->values[0]);
    CHECK_EQ(106, fm.candidates[1]->values[1]);
    CHECK_EQ(171, fm.candidates[2]->values[0]);
    CHECK_EQ(132, fm.candidates[2]->values[1]);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_set_last_element")
{
    fm.variable_count_plus_one = 3;
    vector<int> indices(3);
    fm.idx = &indices[0];

    fm.candidates[0]->values[0] = 300;
    fm.candidates[0]->values[1] = 200;
    fm.candidates[0]->score = 2;
    fm.candidates[1]->values[0] = 42;
    fm.candidates[1]->values[1] = 64;
    fm.candidates[1]->score = 4;
    fm.candidates[2]->values[0] = 26;
    fm.candidates[2]->values[1] = 12;
    fm.candidates[2]->score = 6;

    vector<double> new_vals({ 99, 14 });
    __fminsearch_set_last_element(&fm, &new_vals[0], 3);

    CHECK_EQ(2.0, fm.candidates[0]->score);
    CHECK_EQ(3.0, fm.candidates[1]->score);
    CHECK_EQ(4.0, fm.candidates[2]->score);

}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_false_on_first_nine_attempts")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < 9; ++i)
        CHECK_FALSE(strat.threshold_achieved_checking_similarity(&fm));

}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_true_if_all_attempts_match")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < OPTIMIZER_SIMILARITY_CUTOFF_SIZE - 1; ++i)
        strat.threshold_achieved_checking_similarity(&fm);

    CHECK(strat.threshold_achieved_checking_similarity(&fm));
}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_false_if_attempt_varies")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < 9; ++i)
        strat.threshold_achieved_checking_similarity(&fm);

    fm.candidates[0]->score = 100.1;
    CHECK_FALSE(strat.threshold_achieved_checking_similarity(&fm));
}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_true_if_attempt_varies_minimally")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < OPTIMIZER_SIMILARITY_CUTOFF_SIZE - 1; ++i)
        strat.threshold_achieved_checking_similarity(&fm);

    fm.candidates[0]->score = 100.0001;
    CHECK(strat.threshold_achieved_checking_similarity(&fm));
}

TEST_CASE("LikelihoodRatioTest, update_branchlength")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    unique_ptr<clade> actual(LikelihoodRatioTest::update_branchlength(p_tree.get(), .5, 3));

    CHECK_EQ(15.5, actual->find_descendant("AB")->get_branch_length());
    CHECK_EQ(3.5, actual->find_descendant("A")->get_branch_length());
    CHECK_EQ(7.5, actual->find_descendant("B")->get_branch_length());
}

TEST_CASE("LikelihoodRatioTest, get_likelihood_for_diff_lambdas")
{
    mock_scorer s;
    optimizer opt(&s);
    gene_family gf;
    gf.set_species_size("A", 5);
    gf.set_species_size("B", 9);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    std::vector<lambda*> cache(100);
    CHECK_EQ(0.0, LikelihoodRatioTest::get_likelihood_for_diff_lambdas(gf, p_tree.get(), 0, 0, cache, &opt, 12, 12));
}

TEST_CASE("LikelihoodRatioTest, compute_for_diff_lambdas")
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    user_data data;
    data.p_lambda = &lam;
    data.p_tree = p_tree.get();
    data.gene_families.resize(1);
    data.gene_families[0].set_species_size("A", 5);
    data.gene_families[0].set_species_size("B", 9);
    data.max_root_family_size = 12;
    data.max_family_size = 12;
    vector<int> lambda_index(data.gene_families.size(), -1);
    vector<double> pvalues(data.gene_families.size());
    vector<lambda*> lambdas(100);
    mock_scorer scorer;
    optimizer opt(&scorer);
    LikelihoodRatioTest::compute_for_diff_lambdas_i(data, lambda_index, pvalues, lambdas, &opt);
    CHECK_EQ(0, lambda_index[0]);
    CHECK(isinf(pvalues[0]));
}

void init_lgamma_cache();

int main(int argc, char** argv)
{
    init_lgamma_cache();

    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Global, el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", defaultConf);

    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int res = context.run(); // run

    if (context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests

    return res;
}
