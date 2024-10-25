#ifndef GUARD_funs_h
#define GUARD_funs_h

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>
#include "tree.h"
#include "info.h"
#include "rng.h"

inline double logsumexp(const double &a, const double &b){
  return a < b ? b + log(1.0 + exp(a - b)) : a + log(1.0 + exp(b - a));
}

using std::cout;
using std::endl;

//pi and log(2*pi)
//#define PI 3.1415926535897931
#define LTPI 1.83787706640934536

//--------------------------------------------------
//my functions
//void fit(tree& t, xinfo& xi, dinfo& di, double* fv);
typedef std::vector<std::vector<int> > lookup_t;

lookup_t make_lookup(Rcpp::IntegerMatrix lookup_table, Rcpp::IntegerVector cx);
void impute_x(int v, //variable index
              std::vector<int>& mask,
              int n, xinfo& xi, std::vector<double>& x, std::vector<vector<int> >& x_cat,
              std::vector<int>& cx, std::vector<int>& offsets, std::vector<int>& x_type,
              std::vector<tree>& t, std::vector<double>& y, double& sigma, RNG& rng);

//--------------------------------------------------
//normal density
double pn(
   double x,    //variate
   double m,    //mean
   double v     //variance
);
//--------------------------------------------------
//draw from a discrete distribution
int rdisc(
   double *p,   //vector of probabilities
   RNG& gen     //random number generator
);
//--------------------------------------------------
//evaluate tree tr on grid xi, write to os
void grm(tree& tr, xinfo& xi, std::ostream& os);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars);
//calibart
void getinternalvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var);

//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi);


//From calibart
//--------------------------------------------------
// Get L,U for variable that pertnode splits on given the entire tree
void getLU(tree::tree_p pertnode, xinfo& xi, int* L, int* U);
void getpertLU(tree::tree_p pertnode, size_t pertvar, xinfo& xi, int* L, int* U);



//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv);
//--------------------------------------------------
//get counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di);
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv);
//--------------------------------------------------
//update counts (inc or dec) to reflect observation i
// deprecated:
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, int sign);
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, int sign);

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp,size_t>& bnmap, int sign);
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp,size_t>& bnmap, int sign, tree::tree_cp &tbn);

//--------------------------------------------------
//check minimum leaf size
bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di);
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuffBirth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr);
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuffDeath(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr);

//--------------------------------------------------
//log of the integreted likelihood
double lil(double n, double sy, double tau);
//sy isn't needed, but convenient to maintain fcn sig
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv);
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, double* fv);


/*
template<class T>
double fit_i(T i, tree& t, xinfo& xi, dinfo& di)
{
  double *xx;
  double fv = 0.0;
  tree::tree_cp bn;
	//for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		bn = t.bn(xx,xi);
		fv = bn->getm();
	//}
  return fv;
}
*/


template<class T>
double fit_i(T i, tree& t, xinfo& xi, dinfo& di)
{
  double *xx;
  double fv = 0.0;
  tree::tree_cp bn;
  xx = di.x + i*di.p;
  //for (size_t j=0; j<t.size(); ++j) {
  bn = t.bn(xx,xi);
	fv = bn->getm();
  //}
  return fv;
}
template<class T>
double fit_i(T i, std::vector<tree>& t, xinfo& xi, dinfo& di)
{
  double *xx;
  double fv = 0.0;
  tree::tree_cp bn;
  xx = di.x + i*di.p;
  for (size_t j=0; j<t.size(); ++j) {
		bn = t[j].bn(xx,xi);
		fv += bn->getm();
  }
  return fv;
}
template<class T>
double fit_i_mult(T i, std::vector<tree>& t, xinfo& xi, dinfo& di)
{
  double *xx;
  double fv = 1.0;
  tree::tree_cp bn;
  xx = di.x + i*di.p;
  for (size_t j=0; j<t.size(); ++j) {
  	bn = t[j].bn(xx,xi);
		fv *= bn->getm();
  }
  return fv;
}
//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv);
//--------------------------------------------------
// draw all the bottom node mu's
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double* weight, RNG& gen);
void drphi(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//get min/max for p predictors needed to make cutpoints.
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx);
//make xinfo = cutpoints given minx/maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx);

#endif
