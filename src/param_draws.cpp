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
