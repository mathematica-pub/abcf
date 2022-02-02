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
   double* u;
   double* v;
   double* sigma2_i;
   double& sigma_y;
   double& sigma_u;
   double& sigma_v;
   double& rho;
   double ls_sigma_y;
   double ls_sigma_u;
   double ls_sigma_v;
   double ls_rho;
   int ac_sigma_y;
   int ac_sigma_u;
   int ac_sigma_v;
   int ac_rho;
   // Helper bits
   arma::vec& xform_sigma_v;
   arma::vec& xform_rho;
   arma::mat& xcov_sigma_v_rho;
   double* ftemp;
   double* prop_sig2;
   RNG& gen;
   Logger& logger;
};

// Wrapper around other pieces passed to many functions
struct winfo {
   int ntree;
   std::vector<tree>& t;
   double sd;
   double scale_prec;
   std::vector<std::reference_wrapper<double>> scale_idx;
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

void calculate_rwww(int start, int stop, double* sigma2_i, double scale, double* allfit_spec, double* allfit_alt, std::vector<double>& y, double* w, double& ww, double& rw);

void draw_scale(double& scale, double scale_prec, double ww, double rw, RNG& gen, Logger& logger, bool verbose);

void draw_delta(std::vector<tree>& t, pinfo& pi, double& delta, RNG& gen);

void update_pi(pinfo& pi, double spec_sd, double delta, int ntree, Logger& logger, bool verbose);

void update_mscale(double& mscale,
                    double* allfit_con, double* allfit_mod,
                    ginfo& gi, winfo& wi, bool verbose);

void update_bscale(double& bscale0, double& bscale1,
                    bool b_half_normal,
                    double* allfit_con, double* allfit_mod,
                    ginfo& gi, winfo& wi, bool verbose);

void update_bscale_block(double& bscale0, double& bscale1,
                        bool b_half_normal,
                        double* allfit_con, double* allfit_mod,
                        ginfo& gi, winfo& wi, bool verbose);

void initialize_sigmas(double& sigma_y, double& sigma_u, double& sigma_v, double& rho, double sigu_hyperprior, double sigv_hyperprior, RNG& gen);

double propose_sigma(double sigma_current, double ls_proposal, RNG& gen);

double propose_rho(double rho_current, double ls_proposal, RNG& gen);

arma::vec propose_sigma_v_rho(double sigma_v_current, double rho_current, arma::mat& xcov_sigma_v_rho, RNG& gen);

void update_sigma_y_conj(double* allfit, double& sigma, double nu, double lambda, double mscale, pinfo& pi_con, pinfo& pi_mod, ginfo& gi);

void update_sigma_y(ginfo& gi, double* allfit, double nu, double lambda);

void update_sigma_u(ginfo& gi, double* allfit, double hyperprior);

void update_sigma_v(ginfo& gi, double* allfit, double hyperprior);

void update_rho(ginfo& gi, double* allfit);

void update_sigma_v_rho(ginfo& gi, double* allfit, double hyperprior);

double calculate_lp_diff(ginfo& gi, double* allfit, double log_prior_current, double log_prior_proposed);

void calculate_sigma2_i(ginfo& gi, double sigma_y, double sigma_u, double sigma_v, double rho, double* return_loc);

void draw_uv(double* u, double* v, double* allfit, ginfo& gi);

void update_mh_cov(arma::mat& cov_loc, arma::vec par1, arma::vec par2);

void update_adaptive_ls(ginfo& gi, size_t iter, int batch_size, double ac_target=0.44);

double calculate_adaptive_ls(int accepted, double target, double log_sigma, double increment);

void save_values(size_t& save_ctr, int n, int ntrt,
                Rcpp::NumericVector& msd_post, Rcpp::NumericVector& bsd_post, 
                Rcpp::NumericVector& b0_post, Rcpp::NumericVector& b1_post, Rcpp::NumericVector& sigma_y_post,
                Rcpp::NumericVector& sigma_u_post, Rcpp::NumericVector& sigma_v_post, 
                Rcpp::NumericVector& rho_post, Rcpp::NumericMatrix& sigma_i_post,
                Rcpp::NumericMatrix& m_post, Rcpp::NumericMatrix& yhat_post, Rcpp::NumericMatrix& b_post,
                Rcpp::NumericMatrix& u_post, Rcpp::NumericMatrix& v_post, Rcpp::NumericVector& delta_con_post, 
                double mscale, double bscale1, double bscale0, ginfo& gi,
                double* allfit, double* allfit_con, double* allfit_mod, double delta_con);

#endif