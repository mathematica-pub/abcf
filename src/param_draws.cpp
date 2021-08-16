// Todo: what is the best order/location (cpp vs h) of includes? Can I drop any of these?
#include "funs.h"
#include "bd.h"
#include "logging.h"
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>
#include "param_draws.h"

void drct_test_func(void){
    Rcpp::Rcout << "I am the DRCT test function" << std::endl;
}

void log_iter(std::string context, int current, int final, double sigma, double mscale, double bscale0, double bscale1, Logger& logger) {
    char logBuff[100];
    logger.log("==============================================");
    sprintf(logBuff, "MCMC iteration: %d of %d ", current, final);
    logger.log(logBuff + context);
    sprintf(logBuff, "sigma %f, mscale %f, bscale0 %f, bscale1 %f",sigma, mscale, bscale0, bscale1);
    logger.log(logBuff);
    logger.log("==============================================");
}

void log_fit(std::vector<double>& y, double* allfit, double* allfit_con, double* allfit_mod, Logger& logger, bool verbose) {
    char logBuff[100];
    if (verbose) {
        logger.getVectorHead(y, logBuff);
        Rcpp::Rcout << "           y: " <<  logBuff << "\n";

        logger.getVectorHead(allfit, logBuff);
        Rcpp::Rcout << "Current Fit : " <<  logBuff << "\n";

        logger.getVectorHead(allfit_con, logBuff);
        Rcpp::Rcout << "allfit_con  : " <<  logBuff << "\n";

        logger.getVectorHead(allfit_mod, logBuff);
        Rcpp::Rcout << "allfit_mod  : " <<  logBuff << "\n";
    }
}

void log_status(Rcpp::NumericVector& z_, std::vector<double>& y, double* allfit, double* ri, double mscale, double bscale0, double bscale1, Logger& logger) {
    char logBuff[100];

    logger.log("Printing current status of fit");

    logger.getVectorHead(z_, logBuff);
    Rcpp::Rcout << "\n          z : " <<  logBuff << "\n";

    logger.getVectorHead(y, logBuff);
    Rcpp::Rcout << "          y : " <<  logBuff << "\n";

    logger.getVectorHead(allfit, logBuff);
    Rcpp::Rcout << "Fit - Tree  : " <<  logBuff << "\n";

    logger.getVectorHead(ri, logBuff);
    Rcpp::Rcout << "     resid  : " <<  logBuff << "\n\n";

    Rcpp::Rcout <<" MScale: " << mscale << "\n";

    Rcpp::Rcout <<" bscale0 : " << bscale0 << "\n";

    Rcpp::Rcout <<" bscale1 : " << bscale1 << "\n\n";
}

void print_trees(std::string step, tree& t, xinfo& xi, bool verbose, Logger& logger) {
    logger.log("Attempting to print tree " + step);
    if(verbose){
        t.pr(xi);
        Rcpp::Rcout << "\n\n";
      }
}

// Desired behavior of verbose is unclear to me, since it is hardcoded originally.
void update_trees(std::string context, std::vector<tree>& t, 
                    xinfo& xi, dinfo& di, pinfo& pi,
                    int ntrt, double* weight, Rcpp::NumericVector& z_, std::vector<double>& y,
                    double* allfit, double* allfit_spec, double* ri, 
                    double mscale, double bscale0, double bscale1,
                    RNG& gen, Logger& logger, bool verbose) {
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

      print_trees("pre update", t[iTree], xi, verbose, logger);
      
      fit(t[iTree], // tree& t
          xi, // xinfo& xi
          di, // dinfo& di
          ftemp);
      
      print_trees("post first call to fit", t[iTree], xi, verbose, logger);

      for(size_t k=0;k<n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcpp::Rcout << context << " tree " << iTree <<" obs "<< k<<" "<< std::endl;
          Rcpp::Rcout << t[iTree] << std::endl;
          Rcpp::stop("nan in ftemp");
        }

        // Scale to use depends on context and k
        double scale = (context=="control") ? mscale : (k<ntrt) ? bscale1 : bscale0;

        allfit[k]     = allfit[k]     -scale*ftemp[k];
        allfit_spec[k] = allfit_spec[k] -scale*ftemp[k];
        
        ri[k] = (y[k]-allfit[k])/scale;

        // DRCT: we originally only do this for con, not mod. Any reason?
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

      print_trees("post bd", t[iTree], xi, verbose, logger);
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

      print_trees("post drmu", t[iTree], xi, verbose, logger);

      fit(t[iTree],
          xi,
          di,
          ftemp);

      for(size_t k=0;k<n;k++) {
        double scale = (context=="control") ? mscale : (k<ntrt) ? bscale1 : bscale0;
        allfit[k] += scale*ftemp[k];
        allfit_spec[k] += scale*ftemp[k];
      }

      print_trees("post second call to fit", t[iTree], xi, verbose, logger);
      
      logger.stopContext();
    }
}

void draw_scale(std::string context, double& scale, double scale_prec, double ww, double rw, RNG& gen, Logger& logger, bool verbose) {
    logger.log("Drawing " + context);
    logger.startContext();

    double scale_old = scale;
    double scale_fc_var = 1/(ww + scale_prec);
    scale = scale_fc_var*rw + gen.normal(0., 1.)*sqrt(scale_fc_var);
    if(verbose){
        Rcpp::Rcout << "Original " << context << " : " << scale_old << "\n";
        Rcpp::Rcout << "scale_prec : " << scale_prec << ", ww : " << ww << ", rw : " << rw << "\n";
        Rcpp::Rcout << "New  " << context << " : " << scale << "\n\n";
    }
    logger.stopContext();
}

// Annoying to have to put in all 3 scales, even if only working on one of them.
// For allfit - "spec" = applicable allfit (con for m, mod for b), alt = opposite allfit (mod for m, con for b)
void update_scale(std::string context, std::vector<tree>& t, 
                    double scale_prec, double spec_sd, bool b_half_normal,
                    double sigma, double& mscale, double& bscale0, double& bscale1,
                    double* allfit_spec, double* allfit_alt, pinfo& pi, double& delta,
                    int ntrt, std::vector<double>& y, double* w, 
                    RNG& gen, Logger& logger, bool verbose) {
    
    // Basics
    int n = y.size();
    int ntree = t.size();
    char logBuff[100];

    if (context!="mscale" && context!="bscale") {
        Rcpp::stop("context must be mscale or bscale");
    }

    double ww = 0.0, ww0 = 0.0, ww1 = 0.;
    double rw = 0.0, rw0 = 0.0, rw1 = 0.;
    double s2 = sigma*sigma;

    for(size_t k=0; k<n; ++k) {
        double scale = (context=="mscale") ? mscale : (k<ntrt) ? bscale1 : bscale0;
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

        if(context=="mscale") {
          ww += scale_factor;
          rw += r*scale_factor;
        } else if(k<ntrt) {
          ww1 += scale_factor;
          rw1 += r*scale_factor;
        } else {
          ww0 += scale_factor;
          rw0 += r*scale_factor;
        }
    }

    double mscale_old, bscale0_old, bscale1_old;
    if(context=="mscale") {
        mscale_old = mscale;
        draw_scale("mscale", mscale, scale_prec, ww, rw, gen, logger, verbose);
    } else {
        bscale0_old = bscale0;
        bscale1_old = bscale1;
        draw_scale("bscale1", bscale1, scale_prec, ww1, rw1, gen, logger, verbose);
        draw_scale("bscale0", bscale0, scale_prec, ww0, rw0, gen, logger, verbose);
    }

    for(size_t k=0; k<n; ++k) {
        double scale_ratio = (context=="mscale") ? mscale/mscale_old : (k<ntrt) ? bscale1/bscale1_old : bscale0/bscale0_old;
        allfit_spec[k] = allfit_spec[k]*scale_ratio;
    }

    if(context=="mscale" || !b_half_normal) {
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

    if(verbose){
        logger.log("Updating pi.tau");
        Rcpp::Rcout << "Original pi.tau : " <<  pi.tau << "\n";
    }
      
    pi.tau = spec_sd/(sqrt(delta)*sqrt((double) ntree));
      
    if(verbose){
        Rcpp::Rcout << "New pi.tau : " <<  pi.tau << "\n\n";
    }
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