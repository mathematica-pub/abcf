#ifndef GUARD_param_draws_H
#define GUARD_funs_h

#include "tree.h"
#include "info.h"
#include "rng.h"

// Define two wrappers to simplify function inputs

// General bits, not specific to con or mod
struct ginfo {
   // Basic descriptives
   size_t n;
   int ntrt;
   // Input data
   Rcpp::NumericVector& z_;
   std::vector<double>& y;
   double* w;
   // Helper bits
   RNG& gen;
   Logger& logger;
};

// Wrapper around other pieces passed to many functions
struct winfo {
   int ntree;
   std::vector<tree>& t;
   double sd;
   double scale_prec;
   xinfo& xi;
   dinfo& di;
   pinfo& pi;
   double* ri;
   double* weight;
   double& delta;
};

void log_trees(std::string step, tree& t, xinfo& xi, bool verbose, Logger& logger);

void update_trees(std::string context,
                  double* allfit, double* allfit_spec, 
                  double mscale, double bscale0, double bscale1,
                  ginfo& gi, winfo& wi, bool verbose);

void calculate_rwww(int start, int stop, double s2, double scale, double* allfit_spec, double* allfit_alt, std::vector<double>& y, double* w, double& ww, double& rw);

void draw_scale(double& scale, double scale_prec, double ww, double rw, RNG& gen, Logger& logger, bool verbose);

void draw_delta(std::vector<tree>& t, pinfo& pi, double& delta, RNG& gen);

void update_pi(pinfo& pi, double spec_sd, double delta, int ntree, Logger& logger, bool verbose);

void update_mscale(double& mscale, double sigma,
                    double* allfit_con, double* allfit_mod,
                    ginfo& ginfo, winfo& wi, bool verbose);

void update_bscale(double& bscale0, double& bscale1,
                    bool b_half_normal, double sigma,
                    double* allfit_con, double* allfit_mod,
                    ginfo& gi, winfo& wi, bool verbose);

void update_sigma(double* allfit, double& sigma, double nu, double lambda, double mscale, pinfo& pi_con, pinfo& pi_mod, ginfo& ginfo);

void save_values(size_t& save_ctr, int n, int ntrt,
                Rcpp::NumericVector& msd_post, Rcpp::NumericVector& bsd_post, 
                Rcpp::NumericVector& b0_post, Rcpp::NumericVector& b1_post, Rcpp::NumericVector& sigma_post,
                double mscale, double bscale1, double bscale0, double sigma,
                Rcpp::NumericMatrix& m_post, Rcpp::NumericMatrix& yhat_post, Rcpp::NumericMatrix& b_post,
                double* allfit, double* allfit_con, double* allfit_mod);

#endif