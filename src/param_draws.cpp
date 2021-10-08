// Todo: what is the best order/location (cpp vs h) of includes? Can I drop any of these?
#include "funs.h"
#include "bd.h"
#include "logging.h"
#include <iostream>
#include <RcppArmadillo.h>
#include "param_draws.h"

void log_trees(std::string step, tree& t, xinfo& xi, bool verbose, Logger& logger) {
    logger.log("Attempting to print tree " + step);
    if(verbose){
        t.pr(xi);
        Rcpp::Rcout << "\n\n";
      }
}

// Desired behavior of verbose is unclear to me, since it is hardcoded originally.
void update_trees(std::string context,
                  double* allfit, double* allfit_spec, 
                  double mscale, double bscale0, double bscale1,
                  ginfo& gi, winfo& wi, bool verbose) {
    char logBuff[100];
    double* ftemp  = new double[gi.n];

    for(size_t iTree=0;iTree<wi.ntree;iTree++) {

      gi.logger.log("==================================");
      sprintf(logBuff, "%d of %d", iTree + 1, wi.ntree);
      gi.logger.log("Updating " + context + " tree: " + logBuff);
      gi.logger.log("==================================");
      gi.logger.startContext();

      log_trees("pre update", wi.t[iTree], wi.xi, verbose, gi.logger);
      
      fit(wi.t[iTree], // tree& t
          wi.xi, // xinfo& xi
          wi.di, // dinfo& di
          ftemp);
      
      log_trees("post first call to fit", wi.t[iTree], wi.xi, verbose, gi.logger);

      for(size_t k=0;k<gi.n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcpp::Rcout << context << " tree " << iTree <<" obs "<< k <<" "<< std::endl;
          Rcpp::Rcout << wi.t[iTree] << std::endl;
          Rcpp::stop("nan in ftemp");
        }

        allfit[k]      = allfit[k]      -wi.scale_idx[k]*ftemp[k];
        allfit_spec[k] = allfit_spec[k] -wi.scale_idx[k]*ftemp[k];
        
        wi.ri[k] = (gi.y[k]-allfit[k])/wi.scale_idx[k];

        if(wi.ri[k] != wi.ri[k]) {
          Rcpp::Rcout << (gi.y[k]-allfit[k]) << std::endl;
          Rcpp::Rcout << wi.scale_idx[k] << std::endl;
          Rcpp::Rcout << wi.ri[k] << std::endl;
          Rcpp::stop("NaN in resid");
        }
      }

      if(verbose){
        gi.logger.getVectorHead(wi.weight, logBuff);
        Rcpp::Rcout << "\n weight: " <<  logBuff << "\n\n";
      } 
      gi.logger.log("Starting birth / death processing");
      gi.logger.startContext();
      bd(wi.t[iTree], // tree& x
         wi.xi, // xinfo& xi
         wi.di, // dinfo& di
         wi.weight, // phi
         wi.pi, // pinfo& pi
         gi.gen,
         gi.logger); // RNG& gen
      gi.logger.stopContext();

      log_trees("post bd", wi.t[iTree], wi.xi, verbose, gi.logger);
      if(verbose){
        log_status(gi.z_, gi.y, allfit, wi.ri, mscale, bscale0, bscale1, gi.logger);
      }

      gi.logger.log("Starting to draw mu");
      gi.logger.startContext();

      drmu(wi.t[iTree],  // tree& x
           wi.xi, // xinfo& xi
           wi.di, // dinfo& di
           wi.pi, // pinfo& pi,
           wi.weight,
           gi.gen); // RNG& gen

      gi.logger.stopContext();

      log_trees("post drmu", wi.t[iTree], wi.xi, verbose, gi.logger);

      fit(wi.t[iTree],
          wi.xi,
          wi.di,
          ftemp);

      for(size_t k=0;k<gi.n;k++) {
        allfit[k] += wi.scale_idx[k]*ftemp[k];
        allfit_spec[k] += wi.scale_idx[k]*ftemp[k];
      }

      log_trees("post second call to fit", wi.t[iTree], wi.xi, verbose, gi.logger);
      
      gi.logger.stopContext();
    }
}

void calculate_rwww(int start, int stop, double* sigma_i, double scale, double* allfit_spec, double* allfit_alt, std::vector<double>& y, double& ww, double& rw) {
    for(size_t k=start; k<stop; ++k) {
        double scale_factor = (allfit_spec[k]*allfit_spec[k])/(sigma_i[k]*sigma_i[k]*scale*scale);
        
        if(scale_factor!=scale_factor) {
          Rcpp::Rcout << " scale_factor " << scale_factor << endl;
          Rcpp::stop("NaN in scale factor");
        }

        // numerator is what's unexplained by the other factor
        double r = (y[k] - allfit_alt[k])*scale/allfit_spec[k];
        
        if(r!=r) {
          Rcpp::Rcout << " individual " << k << " r " << r << endl;
          Rcpp::stop("NaN in r");
        }

        ww += scale_factor;
        rw += r*scale_factor;
    }
}

void draw_scale(double& scale, double scale_prec, double ww, double rw, RNG& gen, Logger& logger, bool verbose) {
    logger.startContext();

    double scale_old = scale;
    double scale_fc_var = 1/(ww + scale_prec);
    scale = scale_fc_var*rw + gen.normal(0., 1.)*sqrt(scale_fc_var);
    if(verbose){
        Rcpp::Rcout << "Original : " << scale_old << "\n";
        Rcpp::Rcout << "scale_prec : " << scale_prec << ", ww : " << ww << ", rw : " << rw << "\n";
        Rcpp::Rcout << "New : " << scale << "\n\n";
    }
    logger.stopContext();
}

void draw_delta(std::vector<tree>& t, pinfo& pi, double& delta, RNG& gen) {
  int ntree = t.size();
  double ssq = 0.0;
  tree::npv bnv;
  typedef tree::npv::size_type bvsz;
  double endnode_count = 0.0;

  for(size_t iTree=0;iTree<ntree;iTree++) {
    bnv.clear();
    t[iTree].getbots(bnv);
    bvsz nb = bnv.size();
    for(bvsz ii = 0; ii<nb; ++ii) {
      double mm = bnv[ii]->getm(); //node parameter
      ssq += mm*mm/(pi.tau*pi.tau);
      endnode_count += 1.0;
    }
  }
  delta = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
}

void update_pi(winfo& wi, Logger& logger, bool verbose) {
    if(verbose){
        logger.log("Updating pi.tau");
        Rcpp::Rcout << "Original pi.tau : " <<  wi.pi.tau << "\n";
    }
      
    wi.pi.tau = wi.sd/(sqrt(wi.delta)*sqrt((double) wi.ntree));
      
    if(verbose){
        Rcpp::Rcout << "New pi.tau : " <<  wi.pi.tau << "\n\n";
    }
}

void update_mscale(double& mscale,
                    double* allfit_con, double* allfit_mod,
                    ginfo& gi, winfo& wi, bool verbose) {
    double ww = 0.0;
    double rw = 0.0;

    calculate_rwww(0, gi.n, gi.sigma_i, mscale, allfit_con, allfit_mod, gi.y, ww, rw);

    double mscale_old = mscale;
    gi.logger.log("Drawing mscale");
    draw_scale(mscale, wi.scale_prec, ww, rw, gi.gen, gi.logger, verbose);

    for(size_t k=0; k<gi.n; ++k) {
        allfit_con[k] = allfit_con[k] * mscale / mscale_old;
    }

    draw_delta(wi.t, wi.pi, wi.delta, gi.gen) ;
    
    update_pi(wi, gi.logger, verbose);
}

void update_bscale(double& bscale0, double& bscale1,
                    bool b_half_normal,
                    double* allfit_con, double* allfit_mod,
                    ginfo& gi, winfo& wi, bool verbose) {
    double ww0 = 0.0, ww1 = 0.0;
    double rw0 = 0.0, rw1 = 0.0;

    calculate_rwww(0,       gi.ntrt, gi.sigma_i, bscale1, allfit_mod, allfit_con, gi.y, ww1, rw1);
    calculate_rwww(gi.ntrt, gi.n,    gi.sigma_i, bscale0, allfit_mod, allfit_con, gi.y, ww0, rw0);

    double bscale0_old = bscale0;
    double bscale1_old = bscale1;
    gi.logger.log("Drawing bscale1");
    draw_scale(bscale1, wi.scale_prec, ww1, rw1, gi.gen, gi.logger, verbose);
    gi.logger.log("Drawing bscale0");
    draw_scale(bscale0, wi.scale_prec, ww0, rw0, gi.gen, gi.logger, verbose);

    for(size_t k=0; k<gi.n; ++k) {
        double scale_ratio = (k<gi.ntrt) ? bscale1/bscale1_old : bscale0/bscale0_old;
        allfit_mod[k] = allfit_mod[k]*scale_ratio;
    }

    // Not currently accesible from R - hardcoded to true
    if(!b_half_normal) {
       draw_delta(wi.t, wi.pi, wi.delta, gi.gen) ;
    }

    update_pi(wi, gi.logger, verbose);
}

void initialize_sigmas(double& sigma_y, double& sigma_u, double& sigma_v, double& rho, RNG& gen) {
  // sigma_y is not changed
  sigma_u = abs(gen.normal(0., 1.));
  sigma_v = abs(gen.normal(0., 1.));
  rho = rc_invcdf(gen.uniform(), 0., 1.);
}

double propose_sigma(double sigma_current, double ls_proposal, RNG& gen) {
  double delta = sqrt(exp(2*ls_proposal));
  double log_proposal = log(sigma_current) + gen.normal(0., 1.) * delta;
  double proposal = exp(log_proposal);
  return(proposal);
}

double propose_rho(double rho_current, double ls_proposal, RNG& gen) {
  double delta = sqrt(exp(2*ls_proposal));
  double xformed_proposal = log((rho_current + 1) / (1 - rho_current)) + gen.normal(0., 1.) * delta;
  double proposal = (exp(xformed_proposal) - 1) / (exp(xformed_proposal) + 1);
  return(proposal);
}

void update_sigma_y_conj(double* allfit, double& sigma, double nu, double lambda, double mscale, pinfo& pi_con, pinfo& pi_mod, ginfo& gi) {
  gi.logger.log("Draw sigma");
  double rss = 0.0;
  double restemp = 0.0;
  for(size_t k=0;k<gi.n;k++) {
    restemp = gi.y[k]-allfit[k];
    rss += gi.w[k]*restemp*restemp;
  }

  sigma = sqrt((nu*lambda + rss)/gi.gen.chi_square(nu+gi.n));
  pi_con.sigma = sigma/fabs(mscale);
  pi_mod.sigma = sigma;
}

void update_sigma_y(ginfo& gi, double* allfit, double nu, double lambda) {
  // Proposal is an adaptive MH draw, scaled by ls_sigma_y
  double proposal = propose_sigma(gi.sigma_y, gi.ls_sigma_y, gi.gen);
  double* sigma_i_proposed = calculate_sigma_i(gi, proposal, gi.sigma_u, gi.sigma_v, gi.rho);

  double log_prior_current  = - (nu/2 + 1) * log(gi.sigma_y*gi.sigma_y) - nu*lambda / (2*gi.sigma_y*gi.sigma_y);
  double log_prior_proposed = - (nu/2 + 1) * log(proposal  *proposal)   - nu*lambda / (2*proposal  *proposal);
  double lp_diff = calculate_lp_diff(gi, allfit, log_prior_current, log_prior_proposed, gi.sigma_i, sigma_i_proposed);
  double log_ratio = lp_diff + log(proposal) - log(gi.sigma_y);

  //Accept or reject
  double cut = gi.gen.uniform();
  if (log(cut) < log_ratio) {
    gi.logger.log("Accepting proposed sigma_y " + std::to_string(proposal));
    gi.sigma_y = proposal;
    gi.ac_sigma_y += 1;
    gi.sigma_i = sigma_i_proposed;
  } else {
    gi.logger.log("Rejecting proposed sigma_y " + std::to_string(proposal));
  }
}

void update_sigma_u(ginfo& gi, double* allfit) {
  // Proposal is an adaptive MH draw, scaled by ls_sigma_u
  double proposal = propose_sigma(gi.sigma_u, gi.ls_sigma_u, gi.gen);
  double* sigma_i_proposed = calculate_sigma_i(gi, gi.sigma_y, proposal, gi.sigma_v, gi.rho);

  double log_prior_current =  -gi.sigma_u*gi.sigma_u/2;
  double log_prior_proposed = -proposal*proposal/2;
  double lp_diff = calculate_lp_diff(gi, allfit, log_prior_current, log_prior_proposed, gi.sigma_i, sigma_i_proposed);
  double log_ratio = lp_diff + log(proposal) - log(gi.sigma_u);

  //Accept or reject
  double cut = gi.gen.uniform();
  if (log(cut) < log_ratio) {
    gi.logger.log("Accepting proposed sigma_u " + std::to_string(proposal));
    gi.sigma_u = proposal;
    gi.ac_sigma_u += 1;
    gi.sigma_i = sigma_i_proposed;
  } else {
    gi.logger.log("Rejecting proposed sigma_u " + std::to_string(proposal));
  }
}

void update_sigma_v(ginfo& gi, double* allfit) {
  // Proposal is an adaptive MH draw, scaled by ls_sigma_v
  double proposal = propose_sigma(gi.sigma_v, gi.ls_sigma_v, gi.gen);
  double* sigma_i_proposed = calculate_sigma_i(gi, gi.sigma_y, gi.sigma_u, proposal, gi.rho);

  double log_prior_current =  -gi.sigma_v*gi.sigma_v/2;
  double log_prior_proposed = -proposal*proposal/2;
  double lp_diff = calculate_lp_diff(gi, allfit, log_prior_current, log_prior_proposed, gi.sigma_i, sigma_i_proposed);
  double log_ratio = lp_diff + log(proposal) - log(gi.sigma_v);

  //Accept or reject
  double cut = gi.gen.uniform();
  if (log(cut) < log_ratio) {
    gi.logger.log("Accepting proposed sigma_v " + std::to_string(proposal));
    gi.sigma_v = proposal;
    gi.ac_sigma_v += 1;
    gi.sigma_i = sigma_i_proposed;
  } else {
    gi.logger.log("Rejecting proposed sigma_v " + std::to_string(proposal));
  }
}

void update_rho(ginfo& gi, double* allfit) {
  // Proposal is an adaptive MH draw, scaled by ls_rho
  double proposal = propose_rho(gi.rho, gi.ls_rho, gi.gen);
  double* sigma_i_proposed = calculate_sigma_i(gi, gi.sigma_y, gi.sigma_u, gi.sigma_v, proposal);

  double log_prior_current  = log(1 + cos(M_PI*gi.rho));
  double log_prior_proposed = log(1 + cos(M_PI*proposal));
  double lp_diff = calculate_lp_diff(gi, allfit, log_prior_current, log_prior_proposed, gi.sigma_i, sigma_i_proposed);
  double log_ratio = lp_diff + log((proposal + 1) * (1 - proposal) / ((gi.rho + 1) * (1 - gi.rho)));

  //Accept or reject
  double cut = gi.gen.uniform();
  if (log(cut) < log_ratio) {
    gi.logger.log("Accepting proposed rho " + std::to_string(proposal));
    gi.rho = proposal;
    gi.ac_rho += 1;
    gi.sigma_i = sigma_i_proposed;
  } else {
    gi.logger.log("Rejecting proposed rho " + std::to_string(proposal));
  }
}

// program returns the difference in the log conditional posterior betweeen the propsal and the current value
double calculate_lp_diff(ginfo& gi, double* allfit, double log_prior_current, double log_prior_proposed, double* sigma_i_current, double* sigma_i_proposed) {
  // Log likelihood requires two different sums: sum of the log of sigma_i^2, and sum of resid/sigma_i^2
  double sum_log_sig2_i_current     = 0;
  double sum_r_over_sig2_i_current  = 0;
  double sum_log_sig2_i_proposed    = 0;
  double sum_r_over_sig2_i_proposed = 0;

  for (size_t i=0;i<gi.n;i++) {
    double r = gi.y[i] - allfit[i];
    double r2 = r*r;
    double v_current  = sigma_i_current[i]  * sigma_i_current[i];
    double v_proposed = sigma_i_proposed[i] * sigma_i_proposed[i];
    
    sum_log_sig2_i_current  += log(v_current);
    sum_log_sig2_i_proposed += log(v_proposed);

    sum_r_over_sig2_i_current  += r2/v_current;
    sum_r_over_sig2_i_proposed += r2/v_proposed;
  }
  // Now compose the log posteriors: log prior + log likelihood
  double lp_current  = log_prior_current  -0.5 * sum_log_sig2_i_current  - 0.5 * sum_r_over_sig2_i_current;
  double lp_proposed = log_prior_proposed -0.5 * sum_log_sig2_i_proposed - 0.5 * sum_r_over_sig2_i_proposed;

  double lp_diff = lp_proposed - lp_current;
  return(lp_diff);
}

double* calculate_sigma_i(ginfo& gi, double sigma_y, double sigma_u, double sigma_v, double rho) {
  // precalculate squares rather than calculating inside loop
  double v_y = sigma_y*sigma_y;
  double v_u = sigma_u*sigma_u;
  double v_v = sigma_v*sigma_v;
  double twocov_uv = 2*rho*sigma_u*sigma_v;
  double* sigma_i = new double[gi.n];
  for (size_t i=0;i<gi.n;i++) {
    // since we're working in variances, sigma_v^2 should be multiplied by z^2, but z is binary so no need
    sigma_i[i] = sqrt(v_y/gi.w[i] + v_u + v_v*gi.z_[i] + twocov_uv*gi.z_[i]);
  }

  return(sigma_i);
}

void draw_uv(double* u, double* v, ginfo& gi){
  for(size_t i=0;i<gi.n;i++) {
    u[i] = 0;
    v[i] = 0;
  }
}

void update_adaptive_ls(ginfo& gi, size_t iter, int batch_size, double ac_target) {
  // Convert from target acceptance % to target # accepted iters
  double ac_count = ac_target * (iter+1);
  // We will incrememnt the log_sigma terms by 1/# batches completed
  double ls_incr = batch_size / (iter + 1.);
  
  gi.ls_sigma_y = calculate_adaptive_ls(gi.ac_sigma_y, ac_count, gi.ls_sigma_y, ls_incr);
  gi.ls_sigma_u = calculate_adaptive_ls(gi.ac_sigma_u, ac_count, gi.ls_sigma_u, ls_incr);
  gi.ls_sigma_v = calculate_adaptive_ls(gi.ac_sigma_v, ac_count, gi.ls_sigma_v, ls_incr);
  gi.ls_rho     = calculate_adaptive_ls(gi.ac_rho,     ac_count, gi.ls_rho,     ls_incr);
}

double calculate_adaptive_ls(int accepted, double target, double log_sigma, double increment) {
  if (accepted < target) {
    log_sigma += -increment;
  } else if (accepted > target) {
    log_sigma += increment;
  }
  return(log_sigma);
}

void save_values(size_t& save_ctr, int n, int ntrt,
                Rcpp::NumericVector& msd_post, Rcpp::NumericVector& bsd_post, 
                Rcpp::NumericVector& b0_post, Rcpp::NumericVector& b1_post, Rcpp::NumericVector& sigma_y_post,
                Rcpp::NumericVector& sigma_u_post, Rcpp::NumericVector& sigma_v_post, 
                Rcpp::NumericVector& rho_post, Rcpp::NumericMatrix& sigma_i_post,
                double mscale, double bscale1, double bscale0, ginfo& gi,
                Rcpp::NumericMatrix& m_post, Rcpp::NumericMatrix& yhat_post, Rcpp::NumericMatrix& b_post,
                Rcpp::NumericMatrix& u_post, Rcpp::NumericMatrix& v_post,
                double* allfit, double* allfit_con, double* allfit_mod) {

  msd_post(save_ctr) = mscale;
  bsd_post(save_ctr) = bscale1-bscale0;
  b0_post(save_ctr)  = bscale0;
  b1_post(save_ctr)  = bscale1;
  sigma_y_post(save_ctr) = gi.sigma_y;
  sigma_u_post(save_ctr) = gi.sigma_u;
  sigma_v_post(save_ctr) = gi.sigma_v;
  rho_post(save_ctr) = gi.rho;

  for(size_t k=0;k<n;k++) {
    m_post(save_ctr, k) = allfit_con[k];
    yhat_post(save_ctr, k) = allfit[k];
    double bscale = (k<ntrt) ? bscale1 : bscale0;
    b_post(save_ctr, k) = (bscale1-bscale0)*allfit_mod[k]/bscale;
    u_post(save_ctr, k) = gi.u[k];
    v_post(save_ctr, k) = gi.v[k];
    sigma_i_post(save_ctr,k) = gi.sigma_i[k];
  }

  save_ctr += 1;
}
