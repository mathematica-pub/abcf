#ifndef GUARD_param_draws_H
#define GUARD_funs_h

#include "tree.h"
#include "info.h"
#include "rng.h"

void drct_test_func(void);

void log_status(Rcpp::NumericVector& z_, std::vector<double>& y, double* allfit, double* ri, double mscale, double bscale0, double bscale1, Logger& logger);

void print_trees(std::string step, tree& t, xinfo& xi, bool verbose, Logger& logger);

void update_trees(std::string context, std::vector<tree>& t, 
                    xinfo& xi, dinfo& di, pinfo& pi, 
                    int ntrt, double* weight, Rcpp::NumericVector& z_, std::vector<double>& y,
                    double* allfit, double* allfit_spec, double* ri, 
                    double mscale, double bscale0, double bscale1,
                    RNG& gen, Logger& logger, bool verbose);

#endif