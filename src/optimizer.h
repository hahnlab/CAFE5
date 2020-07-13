#ifndef FMINSEARCH_H
#define FMINSEARCH_H

#include <vector>
#include <map>
#include <chrono>
#include <iosfwd>
#include <functional>
#include <deque>
#include <stdexcept>

//! \defgroup optimizer Optimization
//! @brief Classes and functions designed to calculate optimal values for various parameters
//!
//! One of CAFE's main functions is to calculate the most
//! likely lambda value for a given tree and group of 
//! families. Generally the infer_processes method of 
//! the model is called to calculate the score, but
//! the optimizer_scorer is generic enough for other
//! purposes.
class optimizer_scorer;

//! \ingroup optimizer
enum strategies { RangeWidely, InitialVar, Perturb, Standard, SimilarityCutoff, NLOpt, LBFGS };

struct optimizer_parameters {
    double neldermead_expansion;
    double neldermead_reflection;
    int neldermead_iterations = 300;
    strategies strategy;
    optimizer_parameters();
};

//! @brief Values and score for a potential optimization
//! \ingroup optimizer
struct candidate {
    std::vector<double> values;
    double score;
    candidate(int size);
};

//! @brief Options for an optimization. In general,
//! the Nelder-Mead optimization strategy is used
//! but other possibilities are available.
//! \ingroup optimizer
struct FMinSearch
{
  int maxiters;
  int bymax;
  double rho, chi, psi, sigma;
  double tolx, tolf;
  double delta, zero_delta;

  int 	variable_count, variable_count_plus_one;
  int 	iters;
  std::vector<candidate *> candidates;
  double* x_mean;
  double* x_r;
  double* x_tmp;
  int*    idx;

  optimizer_scorer* scorer;
};

FMinSearch* fminsearch_new();
void fminsearch_free(FMinSearch* pfm);
FMinSearch* fminsearch_new_with_eq(optimizer_scorer* eq, int Xsize);
void fminsearch_set_equation(FMinSearch* pfm, optimizer_scorer* eq, int Xsize);
candidate *get_best_result(FMinSearch* pfm);
void __fminsearch_sort(FMinSearch* pfm);
void** calloc_2dim(int row, int col, int size);
void free_2dim(void** data, int row, int col);
int __fminsearch_checkV(FMinSearch* pfm);
int __fminsearch_checkF(FMinSearch* pfm);
void __fminsearch_min_init(FMinSearch* pfm, double* X0);
void __fminsearch_x_mean(FMinSearch* pfm);
double __fminsearch_x_reflection(FMinSearch* pfm);
double __fminsearch_x_expansion(FMinSearch* pfm);
double __fminsearch_x_contract_outside(FMinSearch* pfm);
double __fminsearch_x_contract_inside(FMinSearch* pfm);
void __fminsearch_x_shrink(FMinSearch* pfm);
void __fminsearch_set_last_element(FMinSearch* pfm, double* x, double f);
bool threshold_achieved(FMinSearch* pfm);
int fminsearch_min(FMinSearch* pfm, double* X0, std::function<bool(FMinSearch*)> threshold_func = threshold_achieved);

class OptimizerStrategy;

//! @brief Provides routines allowing the optimization of some function
//! \ingroup optimizer
//!
//! Comprises a optimizer_scorer and a OptimizerStrategy. The
//! optimizer guesses at a value, uses the scorer to get a score for that
//! value, and then uses the strategy and the score to guess at the 
//! a new value.
class optimizer {
    FMinSearch* pfm;
    optimizer_scorer *_p_scorer;
public:
    optimizer(optimizer_scorer *scorer);
    ~optimizer();

    //! @brief information on the results of an optimization run
    struct result {

        //! the optimized values
        std::vector<double> values;

        //! The score that the scorer returned for the optimized values
        double score;

        //! The number of iterations it took for the optimizer to complete
        int num_iterations;

        //! Length of time (in seconds) that the optimizer took to run
        std::chrono::seconds duration;
    };

    //! The main function that returns optimized values
    result optimize(const optimizer_parameters& params);

    bool quiet = false;

    //! Asks the scorer for initial values to try.
    //! Calculates the score for initial values to 
    //! determine if they are acceptable. Keeps
    //! asking until scorable values are found
    std::vector<double> get_initial_guesses();


    OptimizerStrategy* get_strategy(const optimizer_parameters& params);
};

//! @brief Base class for the optimizer strategy to be used
//! \ingroup optimizer
class OptimizerStrategy
{
public:
    virtual void Run(FMinSearch* pfm, optimizer::result& r, std::vector<double>& initial) = 0;

    virtual std::string Description() const = 0;

    virtual ~OptimizerStrategy() {}
};

class NelderMeadSimilarityCutoff : public OptimizerStrategy
{
    std::deque<double> scores;
public:
    void Run(FMinSearch* pfm, optimizer::result& r, std::vector<double>& initial) override;

    bool threshold_achieved_checking_similarity(FMinSearch* pfm);

    virtual std::string Description() const override { return "Nelder-Mead with similarity cutoff"; };
};

std::ostream& operator<<(std::ostream& ost, const optimizer::result& r);

class OptimizerInitializationFailure : public std::runtime_error {
public:
    OptimizerInitializationFailure() : std::runtime_error("Failed to initialize any reasonable values")
    {

    }
};

#define OPTIMIZER_SIMILARITY_CUTOFF_SIZE    15

#endif
