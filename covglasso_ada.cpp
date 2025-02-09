/*  covglasso: Sparse covariance matrix estimation via the covariance graphical lasso
 *  Copyright (C) 2021  Michael Fop <michael.fop@ucd.ie>
 *  Partly based on R code by Hao Wang (former ass. prof. in  University of South Carolina)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List profileloglik(arma::mat sigma, arma::mat S, int n)
{
  int V = sigma.n_rows;

  arma::mat inv = inv_sympd(sigma);
  double val;
  double sign;
  log_det(val, sign, sigma);
  val = -n*0.5*val -n*0.5*trace( inv*S ) ;

  return Rcpp::List::create( Named("loglik") = val,
                             Named("inv") =  inv);
}


double calcu_pen(const arma::mat &sig,  // ## current sigma
                    double lam,         // ## lambda
                    int pen_type,       // ## 1 denotes L1, 2 denotes SCAD, 3 denotes MCP
                    double ga           // ## tuning parameter for SCAD or MCP
                    
){
  
  int pp = sig.n_cols;
  int j, k;
  arma::mat pen = sig * lam;
  
  if(pen_type == 2){
    for(j=0;j<pp;j++){
      for(k=0;k<pp;k++){
        if(sig(j,k) > ga*lam){
          pen(j,k) = lam * lam * (ga + 1) / 2;
        }
        else if(sig(j,k) > lam){
          pen(j,k) = (- sig(j,k) * sig(j,k) + 2 * ga * lam * sig(j,k) - lam * lam) / (2*ga - 1);
        }
      }
    }
  }
  
  if(pen_type == 3){
    for(j=0;j<pp;j++){
      for(k=0;k<pp;k++){
        if(sig(j,k) > ga*lam){
          pen(j,k) = ga * lam * lam / 2;
        }
        else{
          pen(j,k) = lam * sig(j,k) - sig(j,k) * sig(j,k) / ga / 2;
        }
      }
    }
  }
  
  
  return(arma::accu(pen));
  
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List covglasso(arma::mat S, arma::mat lambda, arma::mat start, int n,
                     double tolout, double tolin, double iterout, double iterin,
                     double lam, int pen_type, double ga)
{
  int V = S.n_rows;
  arma::uvec ind = linspace<uvec>(0, V-1, V);

  arma::mat beta;
  arma::mat sigma = start;
  Rcpp::List out = profileloglik(S, sigma, n);
  arma::mat inv = out["inv"];
  double llkprev = out["loglik"];
  double pen = calcu_pen(sigma, lam, pen_type, ga); // ## added for compute different type of penalties.
  llkprev = llkprev - pen;  // ## modified for computing llkpen
  // double pen;

  double llk;
  bool crit = true;
  double err;
  int it = 0;
  double sign;

  while ( crit ) {  // ## outer iteration
    
    
    for ( int i = 0; i < V; i++ ) {

      arma::uvec sub = find(ind != i);    // set elements [-i]
      arma::uvec h = find(ind == i);      // element [i]

      arma::mat inv11 = inv(sub,sub) - inv.submat(sub,h)*inv.submat(sub,h).t() / as_scalar( inv(h,h) );
      arma::mat inv11s12 = inv11 * S.submat(sub,h);
      arma::mat s11inv11 = S(sub,sub) * inv11;
      arma::mat W = inv11 * s11inv11;

      // update for gamma
      // arma::mat beta = sigma.submat(sub,h);
      beta = sigma.submat(sub,h);
      double a = 0.5;          // Nk in case of mixture models
      double b = as_scalar(lambda(h,h))*0.5;
      double c = as_scalar( beta.t()*W*beta*0.5 - beta.t()*inv11s12 + S(h,h)*0.5 );
      double gam;
      if ( b == 0.0 ) {
        gam = c/a;
      } else {
        gam = ( -a + sqrt(a*a + 4*b*c) )/(2*b);   // ( -a + sqrt(a^2 + 8*Nk*b*c) ) / (4*b)   in case of mixture
      }

      // update for beta
      arma::mat Vmat = W/gam + b*2*inv11;    // 0.5*Nk * W1/gam + b*2*inv11    in case of mixture (also below)
      arma::mat u = inv11s12/gam;            // 0.5*Nk * inv11s12/gam
      arma::mat Vbeta = Vmat * beta;
      bool critin = true;
      double errin;
      int itin = 0;
      arma::mat betaprev = beta;
      while ( critin ) {
        for ( int j = 0; j < (V-1); j++ ) {
          double betajprev = as_scalar(beta(j));
          double tmp = u(j) - Vbeta(j) + Vmat(j,j)*beta(j);
          double sgn;
          if (tmp < 0) sgn = -1/Vmat(j,j); else if (tmp > 0) sgn = 1/Vmat(j,j); else sgn = 0;
          beta(j) = max( 0.0, abs(tmp) - lambda(i,sub(j)) ) * sgn;
          double betadiff = beta(j) - betajprev;
          if ( betadiff != 0 ) Vbeta = Vbeta + betadiff * Vmat.col(j);
        }
        errin = as_scalar(max( abs(beta - betaprev) ));
        itin++;
        critin = ( (errin > tolin) & (itin < iterin) );
        betaprev = beta;
      }

      sigma(sub,h) = beta;
      sigma(h,sub) = beta.t();
      arma::mat inv11beta = inv11*beta;
      sigma(h,h) = gam + beta.t()*inv11beta;
      inv(sub,sub) = inv11 + ( inv11beta*inv11beta.t() )/gam;
      inv(sub,h) = -inv11beta/gam;
      inv(h,sub) = (-inv11beta.t())/gam;
      arma::mat tmp(1,1);
      tmp(0,0) = 1/gam;
      inv(h,h) = tmp;
    }

    out = profileloglik(sigma, S, n);
    //inv = out["inv"];
    llk = out["loglik"];
    pen = calcu_pen(sigma, lam, pen_type, ga); // ## need to be modified
    llk = llk - pen;
    err = std::abs(llk - llkprev) / (1 + std::abs(llk));
    llkprev = llk;
    it++;
    crit = ( (err > tolout) & (it < iterout) );

  }

  return Rcpp::List::create( Named("omega") = inv,
                             Named("sigma") = sigma,
                             Named("loglik") = llk + pen,
                             Named("loglikpen") = llk,
                             Named("pen") = pen,
                             Named("it") = it,
                             Named("err") = err );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List covglassopath_bic_ada(arma::mat S, arma::cube lambda, arma::mat start, int L, int n,
                             double tolout, double tolin, double iterout, double iterin,
                             double lam, int pen_type, double ga)
{
  Rcpp::List out(L);
  Rcpp::List fit;
  arma::uvec par;
  arma::uvec npar(L);
  arma::vec loglik(L);
  int V = S.n_rows;


  for ( int l = 0; l < L; l++ ) {
    fit = covglasso(S, lambda.slice(l), start, n, tolout, tolin, iterout, iterin, lam, pen_type, ga);
    out[l] = fit;
    arma::mat sigma = fit["sigma"];
    par = find(sigma);
    npar(l) = (par.n_elem - V)/2.0 + V;
    loglik(l) = fit["loglik"];
  }

  return Rcpp::List::create(
    Named("out") = out,
    Named("loglik") = loglik,
    Named("npar") = npar
  );
}
