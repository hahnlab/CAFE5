#include <limits>
#include <stdlib.h>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iomanip>
#include <chrono>
#include <memory>

#include "../config.h"

#include "optimizer.h"
#include "io.h"
#include "optimizer_scorer.h"

#define PHASED_OPTIMIZER_PHASE2_PRECISION 1e-6

using namespace std;

const double MAX_DOUBLE = std::numeric_limits<double>::max();

optimizer_parameters::optimizer_parameters() : neldermead_expansion(2.0), neldermead_reflection(1.0)
{
#if defined(OPTIMIZER_STRATEGY_RANGE_WIDELY_THEN_HOME_IN)
    strategy = RangeWidely;
#elif defined(OPTIMIZER_STRATEGY_INITIAL_VARIANTS)
    strategy = InitialVar;
#elif defined(OPTIMIZER_STRATEGY_PERTURB_WHEN_CLOSE)
    strategy = Perturb;
#elif defined(OPTIMIZER_STRATEGY_SIMILARITY_CUTOFF)
    strategy = SimilarityCutoff;
#else
    strategy = Standard;
#endif
}

// TODO: If fminsearch is slow, replacing this with a single array might be more efficient
void** calloc_2dim(int row, int col, int size)
{
  void** data = (void**)calloc(row, sizeof(void*));
  int i;
  for (i = 0; i < row; i++)
  {
    data[i] = calloc(col, size);
  }
  return data;
}

void free_2dim(void** data, int row, int col)
{
  for (int r = 0; r < row; r++)
  {
    free(data[r]);
    data[r] = NULL;
  }
  free(data);
}

FMinSearch* fminsearch_new()
{
	FMinSearch* pfm = (FMinSearch*)calloc(1,sizeof(FMinSearch));	
	pfm->rho = 1;				// reflection
	pfm->chi = 2;				// expansion
	pfm->psi = 0.5;				// contraction
	pfm->sigma = 0.5;			// shrink
	pfm->tolx = 1e-6;
	pfm->tolf = 1e-6;
	pfm->delta = 0.05;
	pfm->zero_delta = 0.00025;
	pfm->maxiters = 250;
	return pfm;
}

FMinSearch* fminsearch_new_with_eq(optimizer_scorer* eq, int Xsize)
{
	FMinSearch* pfm = fminsearch_new();
	fminsearch_set_equation(pfm,eq,Xsize);
	return pfm;
}

void fminsearch_clear_memory(FMinSearch* pfm)
{
    for (auto c : pfm->candidates)
    {
        delete c;
    }
	free(pfm->x_mean);
	pfm->x_mean = NULL;
	free(pfm->x_r);
	pfm->x_r = NULL;
	free(pfm->x_tmp);
	pfm->x_tmp = NULL;
	free(pfm->idx);
	pfm->idx = NULL;
}

void fminsearch_free(FMinSearch* pfm)
{
	fminsearch_clear_memory(pfm);	
	free(pfm);
	pfm = NULL;
}

void fminsearch_set_equation(FMinSearch* pfm, optimizer_scorer* eq, int Xsize)
{
	if ( pfm->variable_count != Xsize )
	{
        if ( pfm->scorer ) fminsearch_clear_memory(pfm);
        pfm->candidates.resize(Xsize + 1);
        for (auto& c : pfm->candidates)
            c = new candidate(Xsize);
        pfm->x_mean = (double*)calloc(Xsize, sizeof(double));
		pfm->x_r = (double*)calloc(Xsize, sizeof(double));
		pfm->x_tmp = (double*)calloc(Xsize, sizeof(double));
		pfm->idx = (int*)calloc(Xsize+1,sizeof(int));
	}
	pfm->scorer = eq;
	pfm->variable_count = Xsize;
	pfm->variable_count_plus_one = Xsize + 1;
}

void __fminsearch_sort(FMinSearch* pfm)
{
    sort(pfm->candidates.begin(), pfm->candidates.end(), [](candidate * a, candidate *b) {
        return a->score < b->score;
    });
}


int __fminsearch_checkV(FMinSearch* pfm)
{
	int i,j;
	double t;
	double max = -MAX_DOUBLE;
	
	for ( i = 0 ; i < pfm->variable_count  ; i++ )
	{
		for ( j = 0 ; j < pfm->variable_count ; j++ )
		{
			t = fabs(pfm->candidates[i+1]->values[j] - pfm->candidates[i]->values[j] );
			if ( t > max ) max = t;
		}
	}
    return max <= pfm->tolx;
}

int __fminsearch_checkF(FMinSearch* pfm)
{
    using namespace std;
	int i;
	double t;
	double max = -MAX_DOUBLE;
	for ( i = 1 ; i < pfm->variable_count_plus_one ; i++ )
	{
		t = fabs( pfm->candidates[i]->score - pfm->candidates[0]->score );
		if ( t > max ) max = t;
	}

	return max <= pfm->tolf;
}

void __fminsearch_min_init(FMinSearch* pfm, double* X0)
{
    // run the optimizer a few times, tweaking the initial values to get an idea of what direction we should move
	int i,j;
	for ( i = 0 ; i < pfm->variable_count_plus_one ; i++ )
	{
		for ( j = 0 ; j < pfm->variable_count ; j++ )
		{
            if ( i > 1 && std::isinf(pfm->candidates[i-1]->score)) {
                if ( (i - 1)  == j )
                {
                    pfm->candidates[i]->values[j] = X0[j] ? ( 1 + pfm->delta*100 ) * X0[j] : pfm->zero_delta;
                }
                else
                {
                    pfm->candidates[i]->values[j] = X0[j];
                }                
            }
            else {
                if ( (i - 1)  == j )
                {
                    pfm->candidates[i]->values[j] = X0[j] ? ( 1 + pfm->delta ) * X0[j] : pfm->zero_delta;
                }
                else
                {
                    pfm->candidates[i]->values[j] = X0[j];
                }
            }
		}
		pfm->candidates[i]->score = pfm->scorer->calculate_score(&pfm->candidates[i]->values[0]);
	}
	__fminsearch_sort(pfm);
}

void __fminsearch_x_mean(FMinSearch* pfm)
{
    // compute the mean of the stored values of each of the variables being optimized
	int i,j;
	for ( i = 0 ; i < pfm->variable_count ; i++ )
	{
		pfm->x_mean[i] = 0;
		for ( j = 0 ; j < pfm->variable_count ; j++ )
		{
			pfm->x_mean[i] += pfm->candidates[j]->values[i];
		}
		pfm->x_mean[i] /= pfm->variable_count;
	}
}

double __fminsearch_x_reflection(FMinSearch* pfm)
{   
	int i;
	for ( i = 0 ; i < pfm->variable_count ; i++ )
	{
		pfm->x_r[i] = pfm->x_mean[i] + pfm->rho * ( pfm->x_mean[i] - pfm->candidates[pfm->variable_count]->values[i] );
	}
	return pfm->scorer->calculate_score(pfm->x_r);
}


double __fminsearch_x_expansion(FMinSearch* pfm)
{
	int i;
	for ( i = 0 ; i < pfm->variable_count ; i++ )
	{
		pfm->x_tmp[i] = pfm->x_mean[i] + pfm->chi * ( pfm->x_r[i] - pfm->x_mean[i] );
	}
	return pfm->scorer->calculate_score(pfm->x_tmp);
}

double __fminsearch_x_contract_outside(FMinSearch* pfm)
{
	int i;
	for ( i = 0 ; i < pfm->variable_count; i++ )
	{
		pfm->x_tmp[i] = pfm->x_mean[i] + pfm->psi * ( pfm->x_r[i] - pfm->x_mean[i] );
	}
	return pfm->scorer->calculate_score(pfm->x_tmp);
}

double __fminsearch_x_contract_inside(FMinSearch* pfm)
{
	int i;
	for ( i = 0 ; i < pfm->variable_count ; i++ )
	{
		pfm->x_tmp[i] = pfm->x_mean[i] + pfm->psi * ( pfm->x_mean[i] - pfm->candidates[pfm->variable_count]->values[i] );
	}
	return pfm->scorer->calculate_score(pfm->x_tmp);
}

void __fminsearch_x_shrink(FMinSearch* pfm)
{
	int i, j;
	for ( i = 1 ; i < pfm->variable_count_plus_one ; i++ )
	{
		for ( j = 0 ; j < pfm->variable_count ; j++ )
		{
			pfm->candidates[i]->values[j] = pfm->candidates[0]->values[j] + pfm->sigma * ( pfm->candidates[i]->values[j] - pfm->candidates[0]->values[j] );
		}
		pfm->candidates[i]->score = pfm->scorer->calculate_score(&pfm->candidates[i]->values[0]);
	}
	__fminsearch_sort(pfm);
}

/// Replace the worst performing row with a new set of values
void __fminsearch_set_last_element(FMinSearch* pfm, double* x, double f)
{
    candidate* c = pfm->candidates.back();
    copy(x, x + pfm->variable_count, c->values.begin());
	c->score = f;
	__fminsearch_sort(pfm);
}

int fminsearch_min(FMinSearch* pfm, double* X0, std::function<bool(FMinSearch *)> threshold_func)
{
	int i;
	__fminsearch_min_init(pfm, X0);
	for ( i = 0 ; i < pfm->maxiters; i++ )
	{
		if (threshold_func(pfm)) 
            break;

        __fminsearch_x_mean(pfm);
		double reflection_score = __fminsearch_x_reflection(pfm);
		if ( reflection_score < pfm->candidates[0]->score )
		{
			double expansion_score = __fminsearch_x_expansion(pfm);
			if (expansion_score < reflection_score ) 
                __fminsearch_set_last_element(pfm,pfm->x_tmp, expansion_score);
			else 
                __fminsearch_set_last_element(pfm,pfm->x_r, reflection_score);
		}
		else if ( reflection_score >= pfm->candidates[pfm->variable_count]->score )
		{
			if ( reflection_score > pfm->candidates[pfm->variable_count]->score )
			{
				double contract_inside_score = __fminsearch_x_contract_inside(pfm);
				if (contract_inside_score < pfm->candidates[pfm->variable_count]->score ) 
                    __fminsearch_set_last_element(pfm,pfm->x_tmp, contract_inside_score);
				else 
                    __fminsearch_x_shrink(pfm);
			}
			else
			{
				double contract_outside_score = __fminsearch_x_contract_outside(pfm);
				if (contract_outside_score <= reflection_score ) 
                    __fminsearch_set_last_element(pfm,pfm->x_tmp, contract_outside_score);
				else 
                    __fminsearch_x_shrink(pfm);
			}
		}
		else
		{
			__fminsearch_set_last_element(pfm,pfm->x_r, reflection_score);
		}
	}
	pfm->bymax = i == pfm->maxiters;
	pfm->iters = i;
	return pfm->bymax;
}

candidate *get_best_result(FMinSearch* pfm)
{
    return pfm->candidates[0];
}

candidate::candidate(int size) : values(size), score(0)
{
}

optimizer::optimizer(optimizer_scorer *p_scorer) : _p_scorer(p_scorer)
{
#ifdef SILENT
    quiet = true;
#endif
    pfm = fminsearch_new();
}

optimizer::~optimizer()
{
    fminsearch_free(pfm);
}


std::vector<double> optimizer::get_initial_guesses()
{
    auto initial = _p_scorer->initial_guesses();

    int i = 0;
    double first_run = _p_scorer->calculate_score(&initial[0]);
    while (std::isinf(first_run) && i < NUM_OPTIMIZER_INITIALIZATION_ATTEMPTS)
    {
        initial = _p_scorer->initial_guesses();
        first_run = _p_scorer->calculate_score(&initial[0]);
        i++;
    }
    if (std::isinf(first_run))
    {
        throw OptimizerInitializationFailure();
    }

    return initial;
}

class StandardNelderMead : public OptimizerStrategy
{
    bool explode = false;

public:
    void Run(FMinSearch *pfm, optimizer::result& r, std::vector<double>& initial)
    {
        if (explode)
        {
            pfm->rho *= 1.5;				// reflection
            pfm->chi *= 25;				// expansion
            pfm->delta = 0.4;
        }
        pfm->tolx = 1e-6;
        pfm->tolf = 1e-6;
        pfm->maxiters = 25;
        fminsearch_min(pfm, &initial[0]);
        auto result = get_best_result(pfm);

        r.score = result->score;
        r.values = result->values;
        r.num_iterations = pfm->iters;
    }

    virtual std::string Description() const override { return "Standard Nelder-Mead"; };
};

void NelderMeadSimilarityCutoff::Run(FMinSearch* pfm, optimizer::result& r, std::vector<double>& initial)
{
    pfm->tolx = 1e-6;
    pfm->tolf = 1e-6;
    pfm->maxiters = 300;
    fminsearch_min(pfm, &initial[0], [this](FMinSearch* pfm) { return threshold_achieved_checking_similarity(pfm); });
    auto result = get_best_result(pfm);

    r.score = result->score;
    r.values = result->values;
    r.num_iterations = pfm->iters;
}

bool NelderMeadSimilarityCutoff::threshold_achieved_checking_similarity(FMinSearch* pfm)
{
    if (threshold_achieved(pfm))
        return true;

    double current = get_best_result(pfm)->score;
    scores.push_back(current);
    if (scores.size() < OPTIMIZER_SIMILARITY_CUTOFF_SIZE)
        return false;

    if (scores.size() > OPTIMIZER_SIMILARITY_CUTOFF_SIZE)
        scores.pop_front();

    double mx = *std::max_element(scores.begin(), scores.end());
    double mn = *std::min_element(scores.begin(), scores.end());
    return mx-mn < OPTIMIZER_LOW_PRECISION;
}

class PerturbWhenClose : public OptimizerStrategy
{
    bool explode = false;
public:
    void Run(FMinSearch *pfm, optimizer::result& r, std::vector<double>& initial)
    {
        if (explode)
        {
            pfm->rho *= 1.5;				// reflection
            pfm->chi *= 25;				// expansion
            pfm->delta = 0.4;
        }
        pfm->tolf = OPTIMIZER_LOW_PRECISION;
        pfm->tolx = OPTIMIZER_LOW_PRECISION;

        fminsearch_min(pfm, &initial[0]);

        cout << "\n*****Threshold achieved, move to Phase 2*****\n\n";
        int phase1_iters = pfm->iters;
        pfm->rho *= 1.3;				// reflection
        pfm->chi *= 15;				// expansion
        pfm->delta = 0.4;
        pfm->tolf = OPTIMIZER_HIGH_PRECISION;
        pfm->tolx = OPTIMIZER_HIGH_PRECISION;

        auto phase1_result = get_best_result(pfm);
        r.score = phase1_result->score;
        r.values = phase1_result->values;

        fminsearch_min(pfm, &initial[0]);
        r.num_iterations = phase1_iters + pfm->iters;

        auto phase2_result = get_best_result(pfm);
        r.score = phase2_result->score;
        r.values = phase2_result->values;
    }
    virtual std::string Description() const override { return "Standard Nelder-Mead"; };
};

class InitialVariants : public OptimizerStrategy
{
    bool explode = false;
    optimizer& _opt;
public:
    InitialVariants(optimizer& opt) : _opt(opt)
    {

    }
    void Run(FMinSearch *pfm, optimizer::result& r, std::vector<double>& initial)
    {
        vector<optimizer::result> results(PHASED_OPTIMIZER_PHASE1_ATTEMPTS);

        for (auto& r : results)
        {
            pfm->tolf = OPTIMIZER_LOW_PRECISION;
            pfm->tolx = OPTIMIZER_LOW_PRECISION;

            initial = _opt.get_initial_guesses();
            fminsearch_min(pfm, &initial[0]);
            auto phase1_result = get_best_result(pfm);
            r.score = phase1_result->score;
            r.values = phase1_result->values;
            r.num_iterations = pfm->iters;
            //        cout << "Threshold achieved, move to Phase 2";

        }

        int phase1_iters = accumulate(results.begin(), results.end(), 0, [](int prev, const optimizer::result& r) { return prev + r.num_iterations;  });
        auto best = min_element(results.begin(), results.end(), [](const optimizer::result& r1, const optimizer::result& r2) { return r1.score < r2.score;  });

        pfm->tolf = OPTIMIZER_HIGH_PRECISION;
        pfm->tolx = OPTIMIZER_HIGH_PRECISION;

        fminsearch_min(pfm, &(best->values)[0]);
        auto phase2_result = get_best_result(pfm);
        r.score = phase2_result->score;
        r.values = phase2_result->values;
        r.num_iterations = pfm->iters + phase1_iters;
    }

    virtual std::string Description() const override { return "Vary initial conditions"; };
};

class RangeWidelyThenHomeIn : public OptimizerStrategy
{
public:
    void Run(FMinSearch *pfm, optimizer::result& r, std::vector<double>& initial)
    {
        pfm->rho *= 1.5;				// reflection
        pfm->chi *= 25;				// expansion
        pfm->delta = 0.4;

        pfm->tolf = OPTIMIZER_LOW_PRECISION;
        pfm->tolx = OPTIMIZER_LOW_PRECISION;

        fminsearch_min(pfm, &initial[0]);

        cout << "\n*****Threshold achieved, move to Phase 2*****\n\n";

        pfm->rho /= 1.5;				// reflection
        pfm->chi /= 25;				// expansion
        pfm->delta = 0.05;
        pfm->tolf = OPTIMIZER_HIGH_PRECISION;
        pfm->tolx = OPTIMIZER_HIGH_PRECISION;
        int phase1_iters = pfm->iters;
        auto phase1_result = get_best_result(pfm);
        initial = phase1_result->values;

        fminsearch_min(pfm, &initial[0]);

        r.num_iterations = phase1_iters + pfm->iters;
        auto phase2_result = get_best_result(pfm);
        r.score = phase2_result->score;
        r.values = phase2_result->values;
    }
    virtual std::string Description() const override { return "Search a wider area when close to a solution"; };
};    


optimizer::result optimizer::optimize(const optimizer_parameters& params)
{
    unique_ptr<OptimizerStrategy> strat(get_strategy(params));

    if (!quiet)
    {
        cout << "\nOptimizer strategy: " << strat->Description() << "\n" << endl;
    }

    using clock = std::chrono::system_clock;

    const auto before = clock::now();
    result r;

    auto initial = get_initial_guesses();
    fminsearch_set_equation(pfm, _p_scorer, initial.size());

    strat->Run(pfm, r, initial);
    r.duration = chrono::duration_cast<chrono::seconds>(clock::now() - before);

    if (!quiet)
    {
        cout << r;
    }
    return r;
}

std::ostream& operator<<(std::ostream& ost, const optimizer::result& r)
{
    if (r.score == -log(0))
    {
        ost << "\nFailed to find reasonable starting values to initialize search." << endl;
    }
    else
    {
        ost << "\nCompleted " << r.num_iterations << " iterations" << endl;
        ost << "Time: " << chrono::duration_cast<chrono::hours>(r.duration).count() << "H";
        ost << " " << chrono::duration_cast<chrono::minutes>(r.duration).count() % 60 << "M";
        ost << " " << chrono::duration_cast<chrono::seconds>(r.duration).count() % 60 << "S" << endl;
        ost << "Best match" << (r.values.size() == 1 ? " is: " : "es are: ") << setw(15) << setprecision(14);
        for (size_t i = 0; i < r.values.size() - 1; ++i)
            ost << r.values[i] << ',';
        ost << r.values[r.values.size() - 1] << endl; 
        ost << "Final -lnL: " << r.score << endl;
    }

    return ost;
}

bool threshold_achieved(FMinSearch* pfm)
{
    return __fminsearch_checkV(pfm) && __fminsearch_checkF(pfm);
}

OptimizerStrategy *optimizer::get_strategy(const optimizer_parameters& params)
{
    pfm->chi = params.neldermead_expansion;
    pfm->rho = params.neldermead_reflection;

    switch (params.strategy)
    {
    case RangeWidely:
        return new RangeWidelyThenHomeIn();
    case InitialVar:
        return new InitialVariants(*this);
    case Perturb:
        return new PerturbWhenClose();
    case SimilarityCutoff:
        return new NelderMeadSimilarityCutoff();
    case Standard:
        return new StandardNelderMead();
    case NLOpt:
    case LBFGS:
        throw std::runtime_error("Optimizer strategy not supported");
    }

    return new StandardNelderMead();
}

