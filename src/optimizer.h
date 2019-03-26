#ifndef FMINSEARCH_H
#define FMINSEARCH_H

#include <vector>
#include <map>
#include <chrono>
#include <iosfwd>

#include "../config.h"

class optimizer_scorer;

struct candidate {
    std::vector<double> values;
    double score;
    candidate(int size);
    ~candidate();
};
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
int fminsearch_min(FMinSearch* pfm, double* X0);
bool threshold_achieved(FMinSearch* pfm);

class OptimizerStrategy;

class optimizer {
    FMinSearch* pfm;
    optimizer_scorer *_p_scorer;
    mutable int phase = 1;
public:
    optimizer(optimizer_scorer *scorer);
    ~optimizer();

    struct result {
        std::vector<double> values;
        double score;
        int num_iterations;
        std::chrono::seconds duration;
    };

    result optimize();

    bool quiet = false;
//    bool explode = false;

    std::vector<double> get_initial_guesses();

    OptimizerStrategy* get_strategy();
};

std::ostream& operator<<(std::ostream& ost, const optimizer::result& r);

enum strategies { RangeWidely, InitialVar, Perturb, Standard, NLOpt };
#ifdef HAVE_NLOPT_HPP
const strategies strategy = NLOpt;
#elif defined(OPTIMIZER_STRATEGY_RANGE_WIDELY_THEN_HOME_IN)
const strategies strategy = RangeWidely;
#elif defined(OPTIMIZER_STRATEGY_INITIAL_VARIANTS)
const strategies strategy = InitialVar;
#elif defined(OPTIMIZER_STRATEGY_PERTURB_WHEN_CLOSE)
const strategies strategy = Perturb;
#else
const strategies strategy = Standard;
#endif


#endif
