#pragma once

#include <iostream>
#include <limits>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp: public gsl_function {
  const F func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->func(x);
  }
  public:
  gsl_function_pp(const F& f) : func(f) {
    function = &gsl_function_pp::invoke; //inherited from gsl_function
    params   = this;                     //inherited from gsl_function
  }
  operator gsl_function*(){return this;}
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
  return gsl_function_pp<F>(func);
}

// gsl_fsolve operating on lambda functions for template construction
template <typename F>
double gsl_fsolve(const F & func, double x_lo, double x_hi, const double tolerance=1e-15, const size_t max_iter = 100) {
  class gsl_root_fsolver_pp {
    gsl_root_fsolver * solver;
    public:
    gsl_root_fsolver_pp(const gsl_root_fsolver_type * s): solver(gsl_root_fsolver_alloc(s)) {}
    ~gsl_root_fsolver_pp() {gsl_root_fsolver_free(solver);}
    operator gsl_root_fsolver*() {return solver;}
  };
  gsl_root_fsolver_pp solver(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(solver, make_gsl_function(func), x_lo, x_hi);
  double result = std::numeric_limits<double>::quiet_NaN();
  int status = GSL_CONTINUE;
  for (size_t i = 0; i < max_iter && status == GSL_CONTINUE; ++i) {
    /* iterate one step of the solver */
    status = gsl_root_fsolver_iterate(solver);
    if (status != GSL_SUCCESS) break;

    /* get the solver's current best solution and bounds */
    result = gsl_root_fsolver_root(solver);
    x_lo   = gsl_root_fsolver_x_lower(solver);
    x_hi   = gsl_root_fsolver_x_upper(solver);

    /* Check to see if the solution is within tolerance */
    status = gsl_root_test_interval(x_lo, x_hi, 0, tolerance);
  }

  /* Check for errors */
  if (status == GSL_CONTINUE) {
      std::cout << "not enough iterations\n";
  } else if (status != GSL_SUCCESS) {
      std::cout << "error: " << gsl_strerror(status) << '\n';
  }

  return result;
}


// gsl_fmin operating on lambda functions for template construction
template <typename F>
double gsl_fmin(const F & func, double x_lo, double x_hi, const double tolerance=1e-6, const size_t max_iter = 100) {
  class gsl_min_fminimizer_pp {
    gsl_min_fminimizer * min;
    public:
    gsl_min_fminimizer_pp(const gsl_min_fminimizer_type * s): min(gsl_min_fminimizer_alloc(s)) {}
    ~gsl_min_fminimizer_pp() {gsl_min_fminimizer_free(min);}
    operator gsl_min_fminimizer*() {return min;}
  };
  gsl_min_fminimizer_pp solver(gsl_min_fminimizer_brent);
  gsl_min_fminimizer_set(solver, make_gsl_function(func), 0.5*(x_lo+x_hi), x_lo, x_hi);
  double result = std::numeric_limits<double>::quiet_NaN();
  int status = GSL_CONTINUE;
  for (size_t i = 0; i < max_iter && status == GSL_CONTINUE; ++i) {
    /* iterate one step of the solver */
    status = gsl_min_fminimizer_iterate(solver);
    if (status != GSL_SUCCESS) break;

    /* get the solver's current best solution and bounds */
    result = gsl_min_fminimizer_x_minimum(solver);
    x_lo   = gsl_min_fminimizer_x_lower(solver);
    x_hi   = gsl_min_fminimizer_x_upper(solver);

    /* Check to see if the solution is within tolerance */
    status = gsl_min_test_interval(x_lo, x_hi, tolerance, 0);
  }

  /* Check for errors */
  if (status == GSL_CONTINUE) {
      std::cout << "not enough iterations\n";
  } else if (status != GSL_SUCCESS) {
      std::cout << "error: " << gsl_strerror(status) << '\n';
  }

  return result;
}


