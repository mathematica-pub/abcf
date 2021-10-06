#ifndef GUARD_logging_h
#define GUARD_logging_h
#include <RcppArmadillo.h>

class Logger
{
    private:

    int level;
    int depth;

    public:
    Logger();

    void log(std::string text);
    void setLevel(int inLevel);
    void startContext();
    void stopContext();
    void getVectorHead(Rcpp::NumericVector x, char s[100]);
    void getVectorHead(std::vector<double> x, char s[100]);
    void getVectorHead(double* x, char s[100]);


};

void log_iter(std::string context, int current, int final, double sigma_y, double sigma_u, double sigma_v, double rho, double mscale, double bscale0, double bscale1, Logger& logger);
void log_fit(std::vector<double>& y, double* allfit, double* allfit_con, double* allfit_mod, double* sigma_i, Logger& logger, bool verbose);
void log_status(Rcpp::NumericVector& z_, std::vector<double>& y, double* allfit, double* ri, double mscale, double bscale0, double bscale1, Logger& logger);

#endif