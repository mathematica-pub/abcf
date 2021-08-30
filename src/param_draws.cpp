// Todo: what is the best order/location (cpp vs h) of includes? Can I drop any of these?
#include "funs.h"
#include "bd.h"
#include "logging.h"
#include <cmath>
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
void update_trees(update_tree_args& args) {
    std::string context = args.context;
    std::vector<tree>& t = args.t;
    xinfo& xi = args.xi;
    dinfo& di = args.di;
    pinfo& pi = args.pi;
    int ntrt = args.ntrt;
    double* weight = args.weight;
    Rcpp::NumericVector& z_ = args.z_;
    std::vector<double>& y = args.y;
    double* allfit = args.allfit;
    double* allfit_spec = args.allfit_spec;
    double* ri = args.ri;
    double mscale = args.mscale;
    double bscale0 = args.bscale0;
    double bscale1 = args.bscale1;
    RNG& gen = args.gen;
    Logger& logger = args.logger;
    bool verbose = args.verbose;
    // No need to pass these from parent
    int n = y.size();
    int ntree = t.size();
    char logBuff[100];
    double* ftemp  = new double[n];

    if (context!="control" && context!="moderate") {
        Rcpp::stop("context must be control or moderate");
    }

    for(size_t iTree=0;iTree<ntree;iTree++) {

      logger.log("==================================");
      sprintf(logBuff, "%d of %d", iTree + 1, ntree);
      logger.log("Updating " + context + " tree: " + logBuff);
      logger.log("==================================");
      logger.startContext();

      log_trees("pre update", t[iTree], xi, verbose, logger);
      
      fit(t[iTree], // tree& t
          xi, // xinfo& xi
          di, // dinfo& di
          ftemp);
      
      log_trees("post first call to fit", t[iTree], xi, verbose, logger);

      for(size_t k=0;k<n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcpp::Rcout << context << " tree " << iTree <<" obs "<< k<<" "<< std::endl;
          Rcpp::Rcout << t[iTree] << std::endl;
          Rcpp::stop("nan in ftemp");
        }

        // If we're updating control trees, use mscale; otherwise use bscale1 for treats and bscale0 for controls
        double scale = (context=="control") ? mscale : (k<ntrt) ? bscale1 : bscale0;

        allfit[k]     = allfit[k]     -scale*ftemp[k];
        allfit_spec[k] = allfit_spec[k] -scale*ftemp[k];
        
        ri[k] = (y[k]-allfit[k])/scale;

        if(ri[k] != ri[k]) {
          Rcpp::Rcout << (y[k]-allfit[k]) << std::endl;
          Rcpp::Rcout << scale << std::endl;
          Rcpp::Rcout << ri[k] << std::endl;
          Rcpp::stop("NaN in resid");
        }
      }

      //Before we didn't print this for mod, now we will be, but it will be weight_het
      if(verbose){
        logger.getVectorHead(weight, logBuff);
        Rcpp::Rcout << "\n weight: " <<  logBuff << "\n\n";
      } 
      logger.log("Starting birth / death processing");
      logger.startContext();
      bd(t[iTree], // tree& x
         xi, // xinfo& xi
         di, // dinfo& di
         weight, // phi
         pi, // pinfo& pi
         gen,
         logger); // RNG& gen
      logger.stopContext();

      log_trees("post bd", t[iTree], xi, verbose, logger);
      if(verbose){
        log_status(z_, y, allfit, ri, mscale, bscale0, bscale1, logger);
      }

      logger.log("Starting to draw mu");
      logger.startContext();

      drmu(t[iTree],  // tree& x
           xi, // xinfo& xi
           di, // dinfo& di
           pi, // pinfo& pi,
           weight,
           gen); // RNG& gen

      logger.stopContext();

      log_trees("post drmu", t[iTree], xi, verbose, logger);

      fit(t[iTree],
          xi,
          di,
          ftemp);

      for(size_t k=0;k<n;k++) {
        // If we;re updating control trees, use mscale; otherwise use bscale1 for treats and bscale0 for controls
        double scale = (context=="control") ? mscale : (k<ntrt) ? bscale1 : bscale0;
        allfit[k] += scale*ftemp[k];
        allfit_spec[k] += scale*ftemp[k];
      }

      log_trees("post second call to fit", t[iTree], xi, verbose, logger);
      
      logger.stopContext();
    }
}

void calculate_rwww(int start, int stop, double s2, double scale, double* allfit_spec, double* allfit_alt, std::vector<double>& y, double* w, double& ww, double& rw) {
    for(size_t k=start; k<stop; ++k) {
        double scale_factor = (w[k]*allfit_spec[k]*allfit_spec[k])/(s2*scale*scale);
        
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

void update_pi(pinfo& pi, double spec_sd, double delta, int ntree, Logger& logger, bool verbose) {
    if(verbose){
        logger.log("Updating pi.tau");
        Rcpp::Rcout << "Original pi.tau : " <<  pi.tau << "\n";
    }
      
    pi.tau = spec_sd/(sqrt(delta)*sqrt((double) ntree));
      
    if(verbose){
        Rcpp::Rcout << "New pi.tau : " <<  pi.tau << "\n\n";
    }
}

void update_mscale(double& mscale, std::vector<tree>& t, 
                    double scale_prec, double spec_sd, double sigma,
                    double* allfit_con, double* allfit_mod, pinfo& pi, double& delta,
                    std::vector<double>& y, double* w, 
                    RNG& gen, Logger& logger, bool verbose) {
    int n = y.size();
    int ntree = t.size();
    char logBuff[100];

    double ww = 0.0;
    double rw = 0.0;
    double s2 = sigma*sigma;

    calculate_rwww(0, n, s2, mscale, allfit_con, allfit_mod, y, w, ww, rw);

    double mscale_old = mscale;
    logger.log("Drawing mscale");
    draw_scale(mscale, scale_prec, ww, rw, gen, logger, verbose);

    for(size_t k=0; k<n; ++k) {
        allfit_con[k] = allfit_con[k] * mscale / mscale_old;
    }

    draw_delta(t, pi, delta, gen) ;
    
    update_pi(pi, spec_sd, delta, ntree, logger, verbose);
}

void update_bscale(double& bscale0, double& bscale1, std::vector<tree>& t, 
                    double scale_prec, double spec_sd, bool b_half_normal, double sigma,
                    double* allfit_con, double* allfit_mod, pinfo& pi, double& delta,
                    int ntrt, std::vector<double>& y, double* w, 
                    RNG& gen, Logger& logger, bool verbose) {
    int n = y.size();
    int ntree = t.size();
    char logBuff[100];

    double ww0 = 0.0, ww1 = 0.0;
    double rw0 = 0.0, rw1 = 0.0;
    double s2 = sigma*sigma;

    calculate_rwww(0, ntrt, s2, bscale1, allfit_mod, allfit_con, y, w, ww1, rw1);
    calculate_rwww(ntrt, n, s2, bscale0, allfit_mod, allfit_con, y, w, ww0, rw0);

    double bscale0_old = bscale0;
    double bscale1_old = bscale1;
    logger.log("Drawing bscale1");
    draw_scale(bscale1, scale_prec, ww1, rw1, gen, logger, verbose);
    logger.log("Drawing bscale0");
    draw_scale(bscale0, scale_prec, ww0, rw0, gen, logger, verbose);

    for(size_t k=0; k<n; ++k) {
        double scale_ratio = (k<ntrt) ? bscale1/bscale1_old : bscale0/bscale0_old;
        allfit_mod[k] = allfit_mod[k]*scale_ratio;
    }

    // Not currently accesible from R - hardcoded to true
    if(!b_half_normal) {
       draw_delta(t, pi, delta, gen) ;
    }

    update_pi(pi, spec_sd, delta, ntree, logger, verbose);
}

void update_sigma(std::vector<double>& y, double* w, double* allfit, double& sigma, double nu, double lambda, double mscale, pinfo& pi_con, pinfo& pi_mod, RNG& gen, Logger& logger) {
  logger.log("Draw sigma");
  int n = y.size();
  double rss = 0.0;
  double restemp = 0.0;
  for(size_t k=0;k<n;k++) {
    restemp = y[k]-allfit[k];
    rss += w[k]*restemp*restemp;
  }

  sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
  pi_con.sigma = sigma/fabs(mscale);
  pi_mod.sigma = sigma;
}

void save_values(size_t& save_ctr, int n, int ntrt,
                Rcpp::NumericVector& msd_post, Rcpp::NumericVector& bsd_post, 
                Rcpp::NumericVector& b0_post, Rcpp::NumericVector& b1_post, Rcpp::NumericVector& sigma_post,
                double mscale, double bscale1, double bscale0, double sigma,
                Rcpp::NumericMatrix& m_post, Rcpp::NumericMatrix& yhat_post, Rcpp::NumericMatrix& b_post,
                double* allfit, double* allfit_con, double* allfit_mod) {

  msd_post(save_ctr) = mscale;
  bsd_post(save_ctr) = bscale1-bscale0;
  b0_post(save_ctr)  = bscale0;
  b1_post(save_ctr)  = bscale1;
  sigma_post(save_ctr) = sigma;

  for(size_t k=0;k<n;k++) {
    m_post(save_ctr, k) = allfit_con[k];
    yhat_post(save_ctr, k) = allfit[k];
    double bscale = (k<ntrt) ? bscale1 : bscale0;
    b_post(save_ctr, k) = (bscale1-bscale0)*allfit_mod[k]/bscale;
  }

  save_ctr += 1;
}