#ifndef ADAPTUTILS_mats_utils_H
#define ADAPTUTILS_mats_utils_H

#include <RcppArmadillo.h>
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;





// Index of functions
// 1. List mazurekd  Estimate All Elements of Stage- and Function-based Population Projection Matrix
// 2. List mdabrowskiego  Estimate All Elements of Function-based Leslie Population Projection Matrix
// 3. List thenewpizzle  Create Element Index for Matrix Estimation with Trait Variants



namespace AdaptMats {
  
  //' Estimate All Elements of Stage- and Function-based Population Projection Matrix
  //' 
  //' Function \code{mazurekd()} swiftly calculates matrix elements in
  //' function-based population projection matrices involving stages.
  //' 
  //' @name mazurekd
  //' 
  //' @param AllStages A large data frame giving all required inputs for vital
  //' rate estimation other than the vital rate model coefficients themselves.
  //' Contains a row for each ultimate matrix element.
  //' @param stageframe The modified stageframe used in matrix calculations.
  //' @param matrixformat An integer representing the style of matrix to develop.
  //' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
  //' and age-by-stage MPM (4).
  //' @param survproxy List of coefficients estimated in model of survival.
  //' @param obsproxy List of coefficients estimated in model of observation.
  //' @param sizeproxy List of coefficients estimated in model of size.
  //' @param repstproxy List of coefficients estimated in model of reproductive 
  //' status.
  //' @param fecproxy List of coefficients estimated in model of fecundity.
  //' @param jsurvproxy List of coefficients estimated in model of juvenile
  //' survival.
  //' @param jobsproxy List of coefficients estimated in model of juvenile
  //' observation.
  //' @param jsizeproxy List of coefficients estimated in model of juvenile size.
  //' @param jrepstproxy List of coefficients estimated in model of juvenile
  //' reproductive status.
  //' @param jmatstproxy List of coefficients estimated in model of juvenile
  //' maturity probability.
  //' @param f2_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param r2_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{a} at each time \emph{t} to be used in analysis.
  //' @param r1_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{a} at each time \emph{t}-1 to be used in analysis.
  //' @param r2_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{b} at each time \emph{t} to be used in analysis.
  //' @param r1_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{b} at each time \emph{t}-1 to be used in analysis.
  //' @param r2_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{c} at each time \emph{t} to be used in analysis.
  //' @param r1_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{c} at each time \emph{t}-1 to be used in analysis.
  //' @param dev_terms A numeric vector containing the deviations to the linear
  //' models input by the user. The order is: survival, observation status, size,
  //' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
  //' observation status, juvenile size, juvenile size_b, juvenile size_c,
  //' juvenile reproductive status, and juvenile maturity status.
  //' @param dens_vr A logical value indicating whether any vital rates are
  //' density dependent.
  //' @param dvr_yn A logical vector indicating whether each vital rate is density
  //' dependent.
  //' @param dvr_style An integer vector indicating the style of density
  //' dependence for each vital rate.
  //' @param dvr_alpha A numeric vector indicating the value of alpha to use in
  //' density dependence for each vital rate.
  //' @param dvr_beta A numeric vector indicating the value of beta to use in
  //' density dependence for each vital rate.
  //' @param dens_n A numeric vector corresponding to the population size to use
  //' in vital rate density dependence calculations.
  //' @param dens A numeric value equal to the spatial density to be used in
  //' calculations.
  //' @param fecmod A scalar multiplier for fecundity.
  //' @param maxsize The maximum primary size to be used in element estimation.
  //' @param maxsizeb The maximum secondary size to be used in element estimation.
  //' @param maxsizec The maximum tertiary size to be used in element estimation.
  //' @param firstage The first age to be included in age-by-stage MPM estimation.
  //' @param finalage The final age to be included in age-by-stage MPM estimation.
  //' @param negfec A logical value denoting whether to change negative estimated
  //' fecundity to 0.
  //' @param yearnumber An integer specifying which time at time \emph{t} to
  //' develop matrices for. Must be in reference to the \code{listofyears} object
  //' developed in the \code{R} matrix estimator function.
  //' @param patchnumber An integer specifying which patch to develop matrices
  //' for. Must be in reference to the \code{listofyears} object developed in the
  //' \code{R} matrix estimator function.
  //' @param exp_tol A numeric value indicating the maximum limit for the
  //' \code{exp()} function to be used in vital rate calculations. Defaults to
  //' \code{700.0}.
  //' @param theta_tol A numeric value indicating a maximum value for theta in
  //' negative binomial probability density estimation. Defaults to
  //' \code{100000000.0}.
  //' @param ipm_cdf A logical value indicating whether to estimate size
  //' transitions using the cumulative density function in cases with continuous
  //' distributions. Defaults to \code{TRUE}, with the midpoint method used if
  //' \code{FALSE}.
  //' @param err_check A logical value indicating whether to export a matrix of
  //' conditional probabilities used to develop the \code{U} matrix. Defaults to
  //' \code{FALSE}.
  //' @param sparse A logical value indicating whether to output matrices in
  //' standard Armadillo format or in sparse matrix format. Defaults to
  //' \code{FALSE}.
  //' @param A_only A Boolean value indicating whether to export U and F matrices
  //' for alteration, or only A matrices. Defaults to \code{TRUE}.
  //' 
  //' @return A list with 2, 3, or 4 elements. If \code{A_only} is set to
  //' \code{FALSE}, then the first 3 elements are matrices, including the main MPM
  //' (A), the survival-transition matrix (U), and a fecundity matrix (F). If
  //' \code{A_only} is set to \code{TRUE}, then only the main MPM matrix (A) is
  //' output. If \code{err_check} is set to \code{TRUE}, then another element is
  //' added, which is a 7 column matrix showing survival probability, observation
  //' probability, reproduction probability, sizea transition probability, sizeb
  //' transition probability, sizec transition probability, and juvenile
  //' transition probability to maturity for each element of the final MPM. It is
  //' possible that, due to evolving development strategy, further columns are
  //' output, as well. Note that if \code{sparse = TRUE}, then output matrices are
  //' in sparse format.
  //' 
  //' @section Notes:
  //' The data frame AllStages introduces variables used in size and fecundity
  //' calculations. This DataFrame is broken up into long vectors composed of
  //' input sizes and related variables for these calculations. The "model" Lists
  //' bring in the vital rate models, and include random coefficients where
  //' needed. We also have a number of extra variables, that include such info as
  //' whether to use the Poisson, negative binomial, Gamma, or Gaussian
  //' distributions for size and fecundity calculations. If \code{sizedist},
  //' \code{sizebdist}, \code{sizecdist}, or \code{fecdist} equals 0, 1, 2, or 3,
  //' then the Poisson, negative binomial, Gaussian, or Gamma is used,
  //' respectively.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List mazurekd(const DataFrame& AllStages, const DataFrame& stageframe,
    int matrixformat, const List& survproxy, const List& obsproxy,
    const List& sizeproxy, const List& sizebproxy, const List& sizecproxy,
    const List& repstproxy, const List& fecproxy, const List& jsurvproxy,
    const List& jobsproxy, const List& jsizeproxy, const List& jsizebproxy,
    const List& jsizecproxy, const List& jrepstproxy, const List& jmatstproxy,
    NumericVector f2_inda, NumericVector f1_inda,
    NumericVector f2_indb, NumericVector f1_indb, NumericVector f2_indc,
    NumericVector f1_indc, StringVector r2_inda, StringVector r1_inda,
    StringVector r2_indb, StringVector r1_indb, StringVector r2_indc,
    StringVector r1_indc, const NumericVector& dev_terms, bool dens_vr,
    LogicalVector dvr_yn, IntegerVector dvr_style, NumericVector dvr_alpha,
    NumericVector dvr_beta, NumericVector dens_n, double dens, double fecmod,
    double maxsize, double maxsizeb, double maxsizec,
    unsigned int firstage, unsigned int finalage, bool negfec, int yearnumber,
    int patchnumber, double exp_tol = 700.0, double theta_tol = 100000000.0,
    bool ipm_cdf = true, bool err_check = false, bool sparse = false,
    bool A_only = true) {
    
    NumericMatrix out;
    
    StringVector stagenames = stageframe["stage"];
    int nostages = stagenames.length();
    unsigned long matrixdim {0};
    
    //Rcout << "mazurekd a ";
    
    int nostages_counter = nostages;
    for (int i = 0; i < nostages_counter; i++) {
      if (LefkoUtils::stringcompare_hard(as<std::string>(stagenames(i)), "AlmostBorn")) nostages -= 1;  
      if (LefkoUtils::stringcompare_hard(as<std::string>(stagenames(i)), "Dead")) nostages -= 1;
    }
    
    //Rcout << "mazurekd - matrixformat: " << matrixformat << endl;
    
    
    if (matrixformat == 1) { // Ehrlen-format hMPM
      matrixdim = nostages * nostages;
    } else if (matrixformat == 2) { // deVries-format hMPM
      matrixdim = nostages * (nostages + 1);
    } else if (matrixformat == 3) { // ahMPM
      matrixdim = nostages;
    } else if (matrixformat == 4) { // age-by-stage MPM
      matrixdim = nostages * (finalage - firstage + 1);
    }
    
    //Rcout << "matrixformat: " << matrixformat << endl;
    //Rcout << "matrixdim: " << matrixdim << endl;
    //Rcout << "mazurekd b ";
    
    // Proxy model imports and settings
    bool sizezero = as<bool>(sizeproxy["zero_inflated"]);
    bool sizebzero = as<bool>(sizebproxy["zero_inflated"]);
    bool sizeczero = as<bool>(sizecproxy["zero_inflated"]);
    bool feczero = as<bool>(fecproxy["zero_inflated"]);
    bool jsizezero = as<bool>(jsizeproxy["zero_inflated"]);
    bool jsizebzero = as<bool>(jsizebproxy["zero_inflated"]);
    bool jsizeczero = as<bool>(jsizecproxy["zero_inflated"]);
    
    bool sizetrunc = as<bool>(sizeproxy["zero_truncated"]);
    bool sizebtrunc = as<bool>(sizebproxy["zero_truncated"]);
    bool sizectrunc = as<bool>(sizecproxy["zero_truncated"]);
    bool fectrunc = as<bool>(fecproxy["zero_truncated"]);
    bool jsizetrunc = as<bool>(jsizeproxy["zero_truncated"]);
    bool jsizebtrunc = as<bool>(jsizebproxy["zero_truncated"]);
    bool jsizectrunc = as<bool>(jsizecproxy["zero_truncated"]);
    
    NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
    NumericVector obscoefs = as<NumericVector>(obsproxy["coefficients"]);
    NumericVector sizecoefs = as<NumericVector>(sizeproxy["coefficients"]);
    NumericVector sizebcoefs = as<NumericVector>(sizebproxy["coefficients"]);
    NumericVector sizeccoefs = as<NumericVector>(sizecproxy["coefficients"]);
    NumericVector repstcoefs = as<NumericVector>(repstproxy["coefficients"]);
    NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
    NumericVector jsurvcoefs = as<NumericVector>(jsurvproxy["coefficients"]);
    NumericVector jobscoefs = as<NumericVector>(jobsproxy["coefficients"]);
    NumericVector jsizecoefs = as<NumericVector>(jsizeproxy["coefficients"]);
    NumericVector jsizebcoefs = as<NumericVector>(jsizebproxy["coefficients"]);
    NumericVector jsizeccoefs = as<NumericVector>(jsizecproxy["coefficients"]);
    NumericVector jrepstcoefs = as<NumericVector>(jrepstproxy["coefficients"]);
    NumericVector jmatstcoefs = as<NumericVector>(jmatstproxy["coefficients"]);
    
    double survsigma = as<double>(survproxy["sigma"]);
    double obssigma = as<double>(obsproxy["sigma"]);
    double sizesigma = as<double>(sizeproxy["sigma"]);
    double sizebsigma = as<double>(sizebproxy["sigma"]);
    double sizecsigma = as<double>(sizecproxy["sigma"]);
    double repstsigma = as<double>(repstproxy["sigma"]);
    double fecsigma = as<double>(fecproxy["sigma"]);
    double jsurvsigma = as<double>(jsurvproxy["sigma"]);
    double jobssigma = as<double>(jobsproxy["sigma"]);
    double jsizesigma = as<double>(jsizeproxy["sigma"]);
    double jsizebsigma = as<double>(jsizebproxy["sigma"]);
    double jsizecsigma = as<double>(jsizecproxy["sigma"]);
    double jrepstsigma = as<double>(jrepstproxy["sigma"]);
    double jmatstsigma = as<double>(jmatstproxy["sigma"]);
    
    int survdist = as<int>(survproxy["dist"]);
    int obsdist = as<int>(obsproxy["dist"]);
    int sizedist = as<int>(sizeproxy["dist"]);
    int sizebdist = as<int>(sizebproxy["dist"]);
    int sizecdist = as<int>(sizecproxy["dist"]);
    int repstdist = as<int>(repstproxy["dist"]);
    int fecdist = as<int>(fecproxy["dist"]);
    int jsurvdist = as<int>(jsurvproxy["dist"]);
    int jobsdist = as<int>(jobsproxy["dist"]);
    int jsizedist = as<int>(jsizeproxy["dist"]);
    int jsizebdist = as<int>(jsizebproxy["dist"]);
    int jsizecdist = as<int>(jsizecproxy["dist"]);
    int jrepstdist = as<int>(jrepstproxy["dist"]);
    int jmatstdist = as<int>(jmatstproxy["dist"]);
    
    //Rcout << "mazurekd c ";
    
    if (NumericVector::is_na(sizesigma)) {
      if (sizedist == 1) {
        sizesigma = 1.0;
      } else {
        sizesigma = 0.0;
      }
    }
    if (NumericVector::is_na(sizebsigma)) {
      if (sizebdist == 1) {
        sizebsigma = 1.0;
      } else {
        sizebsigma = 0.0;
      }
    }
    if (NumericVector::is_na(sizecsigma)) {
      if (sizecdist == 1) {
        sizecsigma = 1.0;
      } else {
        sizecsigma = 0.0;
      }
    }
    if (NumericVector::is_na(fecsigma)) {
      if (fecdist == 1) {
        fecsigma = 1.0;
      } else {
        fecsigma = 0.0;
      }
    }
    if (NumericVector::is_na(jsizesigma)) {
      if (sizedist == 1) {
        jsizesigma = 1.0;
      } else {
        jsizesigma = 0.0;
      }
    }
    if (NumericVector::is_na(jsizebsigma)) {
      if (sizebdist == 1) {
        jsizebsigma = 1.0;
      } else {
        jsizebsigma = 0.0;
      }
    }
    if (NumericVector::is_na(jsizecsigma)) {
      if (sizecdist == 1) {
        jsizecsigma = 1.0;
      } else {
        jsizecsigma = 0.0;
      }
    }
    
    //Rcout << "mazurekd d ";
    
    NumericMatrix vital_year = LefkoUtils::revelations(survproxy, obsproxy, sizeproxy,
      sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
      jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 1);
    
    NumericMatrix vital_patch = LefkoUtils::revelations(survproxy, obsproxy, sizeproxy,
      sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
      jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 2);
    
    arma::imat rand_index = LefkoUtils::foi_index(survproxy, obsproxy, sizeproxy,
      sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
      jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy);
    
    //Rcout << "mazurekd e ";
    
    // NumericVector imports from model_proxy objects
    NumericVector sizeyearzi = as<NumericVector>(sizeproxy["zeroyear"]);
    NumericVector sizebyearzi = as<NumericVector>(sizebproxy["zeroyear"]);
    NumericVector sizecyearzi = as<NumericVector>(sizecproxy["zeroyear"]);
    NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
    NumericVector jsizeyearzi = as<NumericVector>(jsizeproxy["zeroyear"]);
    NumericVector jsizebyearzi = as<NumericVector>(jsizebproxy["zeroyear"]);
    NumericVector jsizecyearzi = as<NumericVector>(jsizecproxy["zeroyear"]);
    
    NumericVector dud_yearzi(sizeyearzi.length());
    
    NumericVector unisyzi = unique(sizeyearzi);
    NumericVector unisyzbi = unique(sizebyearzi);
    NumericVector unisyzci = unique(sizecyearzi);
    NumericVector unijsyzi = unique(jsizeyearzi);
    NumericVector unijsyzbi = unique(jsizebyearzi);
    NumericVector unijsyzci = unique(jsizecyearzi);
    NumericVector unifeci = unique(fecyearzi);
    
    NumericVector sizepatchzi = as<NumericVector>(sizeproxy["zeropatch"]);
    NumericVector sizebpatchzi = as<NumericVector>(sizebproxy["zeropatch"]);
    NumericVector sizecpatchzi = as<NumericVector>(sizecproxy["zeropatch"]);
    NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
    NumericVector jsizepatchzi = as<NumericVector>(jsizeproxy["zeropatch"]);
    NumericVector jsizebpatchzi = as<NumericVector>(jsizebproxy["zeropatch"]);
    NumericVector jsizecpatchzi = as<NumericVector>(jsizecproxy["zeropatch"]);
    
    NumericVector dud_patchzi(sizepatchzi.length());
    
    NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
    NumericVector obsgroups2 = as<NumericVector>(obsproxy["groups2"]);
    NumericVector sizegroups2 = as<NumericVector>(sizeproxy["groups2"]);
    NumericVector sizebgroups2 = as<NumericVector>(sizebproxy["groups2"]);
    NumericVector sizecgroups2 = as<NumericVector>(sizecproxy["groups2"]);
    NumericVector repstgroups2 = as<NumericVector>(repstproxy["groups2"]);
    NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
    NumericVector jsurvgroups2 = as<NumericVector>(jsurvproxy["groups2"]);
    NumericVector jobsgroups2 = as<NumericVector>(jobsproxy["groups2"]);
    NumericVector jsizegroups2 = as<NumericVector>(jsizeproxy["groups2"]);
    NumericVector jsizebgroups2 = as<NumericVector>(jsizebproxy["groups2"]);
    NumericVector jsizecgroups2 = as<NumericVector>(jsizecproxy["groups2"]);
    NumericVector jrepstgroups2 = as<NumericVector>(jrepstproxy["groups2"]);
    NumericVector jmatstgroups2 = as<NumericVector>(jmatstproxy["groups2"]);
    
    NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
    NumericVector obsgroups1 = as<NumericVector>(obsproxy["groups1"]);
    NumericVector sizegroups1 = as<NumericVector>(sizeproxy["groups1"]);
    NumericVector sizebgroups1 = as<NumericVector>(sizebproxy["groups1"]);
    NumericVector sizecgroups1 = as<NumericVector>(sizecproxy["groups1"]);
    NumericVector repstgroups1 = as<NumericVector>(repstproxy["groups1"]);
    NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
    NumericVector jsurvgroups1 = as<NumericVector>(jsurvproxy["groups1"]);
    NumericVector jobsgroups1 = as<NumericVector>(jobsproxy["groups1"]);
    NumericVector jsizegroups1 = as<NumericVector>(jsizeproxy["groups1"]);
    NumericVector jsizebgroups1 = as<NumericVector>(jsizebproxy["groups1"]);
    NumericVector jsizecgroups1 = as<NumericVector>(jsizecproxy["groups1"]);
    NumericVector jrepstgroups1 = as<NumericVector>(jrepstproxy["groups1"]);
    NumericVector jmatstgroups1 = as<NumericVector>(jmatstproxy["groups1"]);
    
    NumericVector sizegroups2zi = as<NumericVector>(sizeproxy["zerogroups2"]);
    NumericVector sizebgroups2zi = as<NumericVector>(sizebproxy["zerogroups2"]);
    NumericVector sizecgroups2zi = as<NumericVector>(sizecproxy["zerogroups2"]);
    NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
    NumericVector jsizegroups2zi = as<NumericVector>(jsizeproxy["zerogroups2"]);
    NumericVector jsizebgroups2zi = as<NumericVector>(jsizebproxy["zerogroups2"]);
    NumericVector jsizecgroups2zi = as<NumericVector>(jsizecproxy["zerogroups2"]);
    
    NumericVector dud_groups2zi(jsizecyearzi.length());
    
    NumericVector sizegroups1zi = as<NumericVector>(sizeproxy["zerogroups1"]);
    NumericVector sizebgroups1zi = as<NumericVector>(sizebproxy["zerogroups1"]);
    NumericVector sizecgroups1zi = as<NumericVector>(sizecproxy["zerogroups1"]);
    NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
    NumericVector jsizegroups1zi = as<NumericVector>(jsizeproxy["zerogroups1"]);
    NumericVector jsizebgroups1zi = as<NumericVector>(jsizebproxy["zerogroups1"]);
    NumericVector jsizecgroups1zi = as<NumericVector>(jsizecproxy["zerogroups1"]);
    
    NumericVector dud_groups1zi(jsizecyearzi.length());
    
    //Rcout << "mazurekd f ";
    
    NumericVector survind = LefkoUtils::flightoficarus(survproxy);
    NumericVector obsind = LefkoUtils::flightoficarus(obsproxy);
    NumericVector sizeind = LefkoUtils::flightoficarus(sizeproxy);
    NumericVector sizebind = LefkoUtils::flightoficarus(sizebproxy);
    NumericVector sizecind = LefkoUtils::flightoficarus(sizecproxy);
    NumericVector repstind = LefkoUtils::flightoficarus(repstproxy);
    NumericVector fecind = LefkoUtils::flightoficarus(fecproxy);
    NumericVector jsurvind = LefkoUtils::flightoficarus(jsurvproxy);
    NumericVector jobsind = LefkoUtils::flightoficarus(jobsproxy);
    NumericVector jsizeind = LefkoUtils::flightoficarus(jsizeproxy);
    NumericVector jsizebind = LefkoUtils::flightoficarus(jsizebproxy);
    NumericVector jsizecind = LefkoUtils::flightoficarus(jsizecproxy);
    NumericVector jrepstind = LefkoUtils::flightoficarus(jrepstproxy);
    NumericVector jmatstind = LefkoUtils::flightoficarus(jmatstproxy);
    
    //Rcout << "mazurekd g ";
    
    NumericVector sizeindzi = zero_flightoficarus(sizeproxy);
    NumericVector sizebindzi = zero_flightoficarus(sizebproxy);
    NumericVector sizecindzi = zero_flightoficarus(sizecproxy);
    NumericVector fecindzi = zero_flightoficarus(fecproxy);
    NumericVector jsizeindzi = zero_flightoficarus(jsizeproxy);
    NumericVector jsizebindzi = zero_flightoficarus(jsizebproxy);
    NumericVector jsizecindzi = zero_flightoficarus(jsizecproxy);
    
    //Rcout << "mazurekd h ";
    
    StringVector survind_rownames = bootson(survproxy);
    StringVector obsind_rownames = bootson(obsproxy);
    StringVector sizeind_rownames = bootson(sizeproxy);
    StringVector sizebind_rownames = bootson(sizebproxy);
    StringVector sizecind_rownames = bootson(sizecproxy);
    StringVector repstind_rownames = bootson(repstproxy);
    StringVector fecind_rownames = bootson(fecproxy);
    StringVector jsurvind_rownames = bootson(jsurvproxy);
    StringVector jobsind_rownames = bootson(jobsproxy);
    StringVector jsizeind_rownames = bootson(jsizeproxy);
    StringVector jsizebind_rownames = bootson(jsizebproxy);
    StringVector jsizecind_rownames = bootson(jsizecproxy);
    StringVector jrepstind_rownames = bootson(jrepstproxy);
    StringVector jmatstind_rownames = bootson(jmatstproxy);
    
    //Rcout << "mazurekd i ";
    
    StringVector sizeind_rownames_zi = zero_bootson(sizeproxy);
    StringVector sizebind_rownames_zi = zero_bootson(sizebproxy);
    StringVector sizecind_rownames_zi = zero_bootson(sizecproxy);
    StringVector fecind_rownames_zi = zero_bootson(fecproxy);
    StringVector jsizeind_rownames_zi = zero_bootson(jsizeproxy);
    StringVector jsizebind_rownames_zi = zero_bootson(jsizebproxy);
    StringVector jsizecind_rownames_zi = zero_bootson(jsizecproxy);
    
    //Rcout << "mazurekd j ";
    
    // AllStages import and settings
    Rcpp::NumericVector stage3_num = as<NumericVector>(AllStages["stage3"]);
    Rcpp::NumericVector stage2n_num = as<NumericVector>(AllStages["stage2n"]);
    Rcpp::NumericVector stage2o_num = as<NumericVector>(AllStages["stage2o"]);
    arma::vec stage3 = as<arma::vec>(stage3_num);
    arma::vec stage2n = as<arma::vec>(stage2n_num);
    arma::vec stage2o = as<arma::vec>(stage2o_num);
    
    Rcpp::NumericVector sz3 = as<NumericVector>(AllStages["size3"]);
    Rcpp::NumericVector sz2n = as<NumericVector>(AllStages["size2n"]);
    Rcpp::NumericVector sz2o = as<NumericVector>(AllStages["size2o"]);
    Rcpp::NumericVector sz1 = as<NumericVector>(AllStages["size1"]);
    Rcpp::NumericVector szb3 = as<NumericVector>(AllStages["sizeb3"]);
    Rcpp::NumericVector szb2n = as<NumericVector>(AllStages["sizeb2n"]);
    Rcpp::NumericVector szb2o = as<NumericVector>(AllStages["sizeb2o"]);
    Rcpp::NumericVector szb1 = as<NumericVector>(AllStages["sizeb1"]);
    Rcpp::NumericVector szc3 = as<NumericVector>(AllStages["sizec3"]);
    Rcpp::NumericVector szc2n = as<NumericVector>(AllStages["sizec2n"]);
    Rcpp::NumericVector szc2o = as<NumericVector>(AllStages["sizec2o"]);
    Rcpp::NumericVector szc1 = as<NumericVector>(AllStages["sizec1"]);
    Rcpp::NumericVector ob3 = as<NumericVector>(AllStages["obs3"]);
    Rcpp::NumericVector fl3 = as<NumericVector>(AllStages["rep3"]);
    Rcpp::NumericVector fl2n = as<NumericVector>(AllStages["rep2n"]);
    Rcpp::NumericVector fl2o = as<NumericVector>(AllStages["rep2o"]);
    Rcpp::NumericVector fl1 = as<NumericVector>(AllStages["rep1"]);
    Rcpp::NumericVector mat3 = as<NumericVector>(AllStages["mat3"]);
    Rcpp::NumericVector mat2n = as<NumericVector>(AllStages["mat2n"]);
    Rcpp::NumericVector mat2o = as<NumericVector>(AllStages["mat2o"]);
    Rcpp::NumericVector mat1 = as<NumericVector>(AllStages["mat1"]);
    Rcpp::NumericVector immat2n = as<NumericVector>(AllStages["imm2n"]);
    Rcpp::NumericVector immat2o = as<NumericVector>(AllStages["imm2o"]);
    Rcpp::NumericVector immat1 = as<NumericVector>(AllStages["imm1"]);
    
    Rcpp::NumericVector repentry = as<NumericVector>(AllStages["repentry3"]);
    Rcpp::NumericVector indata2n = as<NumericVector>(AllStages["indata2n"]);
    Rcpp::NumericVector indata2o = as<NumericVector>(AllStages["indata2o"]);
    Rcpp::NumericVector binwidth3 = as<NumericVector>(AllStages["binwidth"]);
    Rcpp::NumericVector binbwidth3 = as<NumericVector>(AllStages["binbwidth"]);
    Rcpp::NumericVector bincwidth3 = as<NumericVector>(AllStages["bincwidth"]);
    Rcpp::NumericVector actualage2 = as<NumericVector>(AllStages["actualage"]);
    
    Rcpp::NumericVector grp3 = as<NumericVector>(AllStages["group3"]);
    Rcpp::NumericVector grp2n = as<NumericVector>(AllStages["group2n"]);
    Rcpp::NumericVector grp2o = as<NumericVector>(AllStages["group2o"]);
    Rcpp::NumericVector grp1 = as<NumericVector>(AllStages["group1"]);
    
    Rcpp::NumericVector ovestt_num = as<NumericVector>(AllStages["ovest_t"]);
    arma::vec ovestt = as<arma::vec>(ovestt_num);
    
    Rcpp::NumericVector ovestf_num = as<NumericVector>(AllStages["ovest_f"]);
    arma::vec ovestf = as<arma::vec>(ovestf_num);
    
    Rcpp::NumericVector indata = as<NumericVector>(AllStages["indata"]);
    Rcpp::NumericVector ovgivent = as<NumericVector>(AllStages["ovgiven_t"]);
    Rcpp::NumericVector ovgivenf = as<NumericVector>(AllStages["ovgiven_f"]);
    
    Rcpp::NumericVector ovsurvmult = as<NumericVector>(AllStages["ovsurvmult"]);
    Rcpp::NumericVector ovfecmult = as<NumericVector>(AllStages["ovfecmult"]);
    arma::vec ovsurvmult_arma = as<arma::vec>(ovsurvmult);
    arma::vec ovfecmult_arma = as<arma::vec>(ovfecmult);
    
    Rcpp::IntegerVector index321_int = as<IntegerVector>(AllStages["index321"]);
    arma::uvec index321 = as<arma::uvec>(index321_int);
    
    Rcpp::IntegerVector aliveandequal = as<IntegerVector>(AllStages["aliveandequal"]);
    
    int n = static_cast<int>(stage3.n_elem);
    
    //Rcout << "n (stage3.n_elem): " << n << endl;
    //Rcout << "mazurekd l ";
    
    arma::uvec replacetvec = find(ovestt != -1.0);
    arma::uvec replacefvec = find(ovestf != -1.0);
    int replacementst = static_cast<int>(replacetvec.n_elem);
    int replacementsf = static_cast<int>(replacefvec.n_elem);
    
    arma::uvec tmults = find(ovsurvmult_arma != -1.0);
    arma::uvec fmults = find(ovfecmult_arma != -1.0);
    arma::uvec noreplacetvec = find(ovestt == -1.0);
    arma::uvec noreplacefvec = find(ovestf == -1.0);
    arma::uvec tmults_only = intersect(tmults, noreplacetvec);
    arma::uvec fmults_only = intersect(fmults, noreplacefvec);
    int tmults_only_st = static_cast<int>(tmults_only.n_elem);
    int fmults_only_st = static_cast<int>(fmults_only.n_elem);
    
    int repindex {0};
    int properindex {0};
    int proxyindex {0};
    
    //Rcout << "mazurekd m ";
    
    // Determination of choices of fixed and random individual covariates
    double inda1 = f1_inda(0);
    double indb1 = f1_indb(0);
    double indc1 = f1_indc(0);
    double inda2 = f2_inda(0);
    double indb2 = f2_indb(0);
    double indc2 = f2_indc(0);
    
    String chosen_f2inda_cat("none");
    String chosen_f1inda_cat("none");
    String chosen_f2indb_cat("none");
    String chosen_f1indb_cat("none");
    String chosen_f2indc_cat("none");
    String chosen_f1indc_cat("none");
    
    String chosen_r2inda(r2_inda(0));
    String chosen_r1inda(r1_inda(0));
    String chosen_r2indb(r2_indb(0));
    String chosen_r1indb(r1_indb(0));
    String chosen_r2indc(r2_indc(0));
    String chosen_r1indc(r1_indc(0));
    
    // Density corrections under vital rate density dependence
    double vr1_dcorr {1.0};
    double vr2_dcorr {1.0};
    double vr3_dcorr {1.0};
    double vr4_dcorr {1.0};
    double vr5_dcorr {1.0};
    double vr6_dcorr {1.0};
    double vr7_dcorr {1.0};
    double vr8_dcorr {1.0};
    double vr9_dcorr {1.0};
    double vr10_dcorr {1.0};
    double vr11_dcorr {1.0};
    double vr12_dcorr {1.0};
    double vr13_dcorr {1.0};
    double vr14_dcorr {1.0};
    
    //Rcout << "mazurekd n ";
    
    if (dens_vr) {
      
      //Rcout << "mazurekd o ";
      
      // Adult survival
      if (dvr_yn(0)) {
        if (dvr_style(0) == 1) {
          vr1_dcorr = dvr_alpha(0) * exp(-1 * dvr_beta(0) * dens_n(0));
        } else if (dvr_style(0) == 2) {
          vr1_dcorr = dvr_alpha(0) / (1 + dvr_beta(0) * dens_n(0));
        } else if (dvr_style(0) == 3) {
          vr1_dcorr = 1 / (1 + exp(dvr_alpha(0) * dens_n(0) + dvr_beta(0)));
        } else if (dvr_style(0) == 4) {
          vr1_dcorr = 1 - (dens_n(0) / dvr_alpha(0));
          if (dvr_beta(0) != 0 && vr1_dcorr < -1.0) vr1_dcorr = 0.0; 
        }
      }
      
      // Adult observation
      if (dvr_yn(1)) {
        if (dvr_style(1) == 1) {
          vr2_dcorr = dvr_alpha(1) * exp(-1 * dvr_beta(1) * dens_n(1));
        } else if (dvr_style(1) == 2) {
          vr2_dcorr = dvr_alpha(1) / (1 + dvr_beta(1) * dens_n(1));
        } else if (dvr_style(1) == 3) {
          vr2_dcorr = 1 / (1 + exp(dvr_alpha(1) * dens_n(1) + dvr_beta(1)));
        } else if (dvr_style(1) == 4) {
          vr2_dcorr = 1 - (dens_n(1) / dvr_alpha(1));
          if (dvr_beta(1) != 0 && vr2_dcorr < -1.0) vr2_dcorr = 0.0; 
        }
      }
      
      // Adult sizea
      if (dvr_yn(2)) {
        if (dvr_style(2) == 1) {
          vr3_dcorr = dvr_alpha(2) * exp(-1 * dvr_beta(2) * dens_n(2));
        } else if (dvr_style(2) == 2) {
          vr3_dcorr = dvr_alpha(2) / (1 + dvr_beta(2) * dens_n(2));
        } else if (dvr_style(2) == 3) {
          vr3_dcorr = 1 / (1 + exp(dvr_alpha(2) * dens_n(2) + dvr_beta(2)));
        } else if (dvr_style(2) == 4) {
          vr3_dcorr = 1 - (dens_n(2) / dvr_alpha(2));
          if (dvr_beta(2) != 0 && vr3_dcorr < -1.0) vr3_dcorr = 0.0; 
        }
      }
      
      // Adult sizeb
      if (dvr_yn(3)) {
        if (dvr_style(3) == 1) {
          vr4_dcorr = dvr_alpha(3) * exp(-1 * dvr_beta(3) * dens_n(3));
        } else if (dvr_style(3) == 2) {
          vr4_dcorr = dvr_alpha(3) / (1 + dvr_beta(3) * dens_n(3));
        } else if (dvr_style(3) == 3) {
          vr4_dcorr = 1 / (1 + exp(dvr_alpha(3) * dens_n(3) + dvr_beta(3)));
        } else if (dvr_style(3) == 4) {
          vr4_dcorr = 1 - (dens_n(3) / dvr_alpha(3));
          if (dvr_beta(3) != 0 && vr4_dcorr < -1.0) vr4_dcorr = 0.0; 
        }
      }
      
      // Adult sizec
      if (dvr_yn(4)) {
        if (dvr_style(4) == 1) {
          vr5_dcorr = dvr_alpha(4) * exp(-1 * dvr_beta(4) * dens_n(4));
        } else if (dvr_style(4) == 2) {
          vr5_dcorr = dvr_alpha(4) / (1 + dvr_beta(4) * dens_n(4));
        } else if (dvr_style(4) == 3) {
          vr5_dcorr = 1 / (1 + exp(dvr_alpha(4) * dens_n(4) + dvr_beta(4)));
        } else if (dvr_style(4) == 4) {
          vr5_dcorr = 1 - (dens_n(4) / dvr_alpha(4));
          if (dvr_beta(4) != 0 && vr5_dcorr < -1.0) vr5_dcorr = 0.0; 
        }
      }
      
      // Adult reproduction
      if (dvr_yn(5)) {
        if (dvr_style(5) == 1) {
          vr6_dcorr = dvr_alpha(5) * exp(-1 * dvr_beta(5) * dens_n(5));
        } else if (dvr_style(5) == 2) {
          vr6_dcorr = dvr_alpha(5) / (1 + dvr_beta(5) * dens_n(5));
        } else if (dvr_style(5) == 3) {
          vr6_dcorr = 1 / (1 + exp(dvr_alpha(5) * dens_n(5) + dvr_beta(5)));
        } else if (dvr_style(5) == 4) {
          vr6_dcorr = 1 - (dens_n(5) / dvr_alpha(5));
          if (dvr_beta(5) != 0 && vr6_dcorr < -1.0) vr6_dcorr = 0.0; 
        }
      }
      
      // Adult fecundity
      if (dvr_yn(6)) {
        if (dvr_style(6) == 1) {
          vr7_dcorr = dvr_alpha(6) * exp(-1 * dvr_beta(6) * dens_n(6));
        } else if (dvr_style(6) == 2) {
          vr7_dcorr = dvr_alpha(6) / (1 + dvr_beta(6) * dens_n(6));
        } else if (dvr_style(6) == 3) {
          vr7_dcorr = 1 / (1 + exp(dvr_alpha(6) * dens_n(6) + dvr_beta(6)));
        } else if (dvr_style(6) == 4) {
          vr7_dcorr = 1 - (dens_n(6) / dvr_alpha(6));
          if (dvr_beta(6) != 0 && vr7_dcorr < -1.0) vr7_dcorr = 0.0; 
        }
      }
      
      // Juvenile survival
      if (dvr_yn(7)) {
        if (dvr_style(7) == 1) {
          vr8_dcorr = dvr_alpha(7) * exp(-1 * dvr_beta(7) * dens_n(7));
        } else if (dvr_style(7) == 2) {
          vr8_dcorr = dvr_alpha(7) / (1 + dvr_beta(7) * dens_n(7));
        } else if (dvr_style(7) == 3) {
          vr8_dcorr = 1 / (1 + exp(dvr_alpha(7) * dens_n(7) + dvr_beta(7)));
        } else if (dvr_style(7) == 4) {
          vr8_dcorr = 1 - (dens_n(7) / dvr_alpha(7));
          if (dvr_beta(7) != 0 && vr8_dcorr < -1.0) vr8_dcorr = 0.0; 
        }
      }
      
      // Juvenile observation
      if (dvr_yn(8)) {
        if (dvr_style(8) == 1) {
          vr9_dcorr = dvr_alpha(8) * exp(-1 * dvr_beta(8) * dens_n(8));
        } else if (dvr_style(8) == 2) {
          vr9_dcorr = dvr_alpha(8) / (1 + dvr_beta(8) * dens_n(8));
        } else if (dvr_style(8) == 3) {
          vr9_dcorr = 1 / (1 + exp(dvr_alpha(8) * dens_n(8) + dvr_beta(8)));
        } else if (dvr_style(8) == 4) {
          vr9_dcorr = 1 - (dens_n(8) / dvr_alpha(8));
          if (dvr_beta(8) != 0 && vr9_dcorr < -1.0) vr9_dcorr = 0.0; 
        }
      }
      
      // Juvenile sizea
      if (dvr_yn(9)) {
        if (dvr_style(9) == 1) {
          vr10_dcorr = dvr_alpha(9) * exp(-1 * dvr_beta(9) * dens_n(9));
        } else if (dvr_style(9) == 2) {
          vr10_dcorr = dvr_alpha(9) / (1 + dvr_beta(9) * dens_n(9));
        } else if (dvr_style(9) == 3) {
          vr10_dcorr = 1 / (1 + exp(dvr_alpha(9) * dens_n(9) + dvr_beta(9)));
        } else if (dvr_style(9) == 4) {
          vr10_dcorr = 1 - (dens_n(9) / dvr_alpha(9));
          if (dvr_beta(9) != 0 && vr10_dcorr < -1.0) vr10_dcorr = 0.0; 
        }
      }
      
      // Juvenile sizeb
      if (dvr_yn(10)) {
        if (dvr_style(10) == 1) {
          vr11_dcorr = dvr_alpha(10) * exp(-1 * dvr_beta(10) * dens_n(10));
        } else if (dvr_style(10) == 2) {
          vr11_dcorr = dvr_alpha(10) / (1 + dvr_beta(10) * dens_n(10));
        } else if (dvr_style(10) == 3) {
          vr11_dcorr = 1 / (1 + exp(dvr_alpha(10) * dens_n(10) + dvr_beta(10)));
        } else if (dvr_style(10) == 4) {
          vr11_dcorr = 1 - (dens_n(10) / dvr_alpha(10));
          if (dvr_beta(10) != 0 && vr11_dcorr < -1.0) vr11_dcorr = 0.0; 
        }
      }
      
      // Juvenile sizec
      if (dvr_yn(11)) {
        if (dvr_style(11) == 1) {
          vr12_dcorr = dvr_alpha(11) * exp(-1 * dvr_beta(11) * dens_n(11));
        } else if (dvr_style(11) == 2) {
          vr12_dcorr = dvr_alpha(11) / (1 + dvr_beta(11) * dens_n(11));
        } else if (dvr_style(11) == 3) {
          vr12_dcorr = 1 / (1 + exp(dvr_alpha(11) * dens_n(11) + dvr_beta(11)));
        } else if (dvr_style(11) == 4) {
          vr12_dcorr = 1 - (dens_n(11) / dvr_alpha(11));
          if (dvr_beta(11) != 0 && vr12_dcorr < -1.0) vr12_dcorr = 0.0; 
        }
      }
      
      // Juvenile reproduction
      if (dvr_yn(12)) {
        if (dvr_style(12) == 1) {
          vr13_dcorr = dvr_alpha(12) * exp(-1 * dvr_beta(12) * dens_n(12));
        } else if (dvr_style(12) == 2) {
          vr13_dcorr = dvr_alpha(12) / (1 + dvr_beta(12) * dens_n(12));
        } else if (dvr_style(12) == 3) {
          vr13_dcorr = 1 / (1 + exp(dvr_alpha(12) * dens_n(12) + dvr_beta(12)));
        } else if (dvr_style(12) == 4) {
          vr13_dcorr = 1 - (dens_n(12) / dvr_alpha(12));
          if (dvr_beta(12) != 0 && vr13_dcorr < -1.0) vr13_dcorr = 0.0; 
        }
      }
      
      // Juvenile maturity
      if (dvr_yn(13)) {
        if (dvr_style(13) == 1) {
          vr14_dcorr = dvr_alpha(13) * exp(-1 * dvr_beta(13) * dens_n(13));
        } else if (dvr_style(13) == 2) {
          vr14_dcorr = dvr_alpha(13) / (1 + dvr_beta(13) * dens_n(13));
        } else if (dvr_style(13) == 3) {
          vr14_dcorr = 1 / (1 + exp(dvr_alpha(13) * dens_n(13) + dvr_beta(13)));
        } else if (dvr_style(13) == 4) {
          vr14_dcorr = 1 - (dens_n(13) / dvr_alpha(13));
          if (dvr_beta(13) != 0 && vr14_dcorr < -1.0) vr14_dcorr = 0.0; 
        }
      }
      
      //Rcout << "mazurekd p ";
      
    }
    
    //Rcout << "mazurekd q ";
    
    // Matrix out collects conditional probabilities
    // It is a zero matrix with n rows and 7 columns: 0 surv, 1 obs, 2 repst,
    // 3 size, 4 size_b, 5 size_c, 6 matst, >6 are test variables
    if (err_check) {
      NumericMatrix zeroform(n, 7);
      std::fill(zeroform.begin(), zeroform.end(), 0.0); // Added for Linux issues
      out = zeroform;
      CharacterVector out_names = {"surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst"};
      colnames(out) = out_names;
    }
    NumericVector out_vec(7, 0.0);
    
    arma::mat survtransmat;
    arma::mat fectransmat;
    arma::sp_mat survtransmat_sp;
    arma::sp_mat fectransmat_sp;
    
    if (!sparse) {
      arma::mat survtransmat_pre(matrixdim, matrixdim, fill::zeros);
      arma::mat fectransmat_pre(matrixdim, matrixdim, fill::zeros);
      
      survtransmat = survtransmat_pre;
      fectransmat = fectransmat_pre;
    } else {
      arma::sp_mat survtransmat_sp_pre(matrixdim, matrixdim);
      arma::sp_mat fectransmat_sp_pre(matrixdim, matrixdim);
      
      survtransmat_sp = survtransmat_sp_pre;
      fectransmat_sp = fectransmat_sp_pre;
    }
    
    double fec_addedcoefs = sum(feccoefs);
    double jsurv_coefsadded = sum(jsurvcoefs);
    double mat_predicted {0.0};
    unsigned int k {0};
    
    //Rcout << "mazurekd r ";
    
    // Loop through each line of AllStages, calculating estimable elements
    for(int i = 0; i < n; i++) {
      out_vec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
      k = aliveandequal(i);
      mat_predicted = 0.0;
      
      if (err_check) out(i, 6) = 1.0; // Initialize maturity status prob for typical case
      
      //Rcout << "mazurekd s i: " << i << " k: " << k << " ";
      
      double anna1 {0.}; /////
      double anna2 {0.};
      double annb1 {0.};
      double annb2 {0.};
      double annc1 {0.};
      double annc2 {0.};
      
      Rcpp::NumericVector statusterms = {fl1(i), fl2n(i), sz1(i), sz2o(i),
        szb1(i), szb2o(i), szc1(i), szc2o(i), actualage2(i), inda1, inda2, indb1,
        indb2, indc1, indc2, dens, sz3(i), szb3(i), szc3(i), binwidth3(i),
        binbwidth3(i), bincwidth3(i), anna1, anna2, annb1, annb2, annc1, annc2};
      
      if (ovgivent(i) == -1 && indata(i) == 1 && stage2n(i) == stage2o(i)) {
        if ((mat2n(i) == 1 && mat3(i) == 1) || (mat2o(i) == 1 && mat3(i) == 1)) {
          
          // Adult survival transitions
          if (survdist < 5) {
            //Rcout << "mazurekd s1 ";
            out_vec(0) = preouterator(survproxy, survcoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, chosen_f2inda_cat,
              chosen_f1inda_cat, chosen_f2indb_cat, chosen_f1indb_cat,
              chosen_f2indc_cat, chosen_f1indc_cat, statusterms, survgroups2,
              survgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
              survind, survind_rownames, sizeindzi, sizeind_rownames_zi, false,
              survsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 1, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0);
            //Rcout << "mazurekd s2 ";
            out_vec(0) = out_vec(0) * vr1_dcorr;
          } else {
            out_vec(0) = survcoefs(0);
            out_vec(0) = out_vec(0) * vr1_dcorr;
          }
          if (err_check) out(i, 0) = out_vec(0);
          
          if (obsdist < 5) {
            //Rcout << "mazurekd s3 ";
            out_vec(1) = preouterator(obsproxy, obscoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, chosen_f2inda_cat,
              chosen_f1inda_cat, chosen_f2indb_cat, chosen_f1indb_cat,
              chosen_f2indc_cat, chosen_f1indc_cat, statusterms, obsgroups2,
              obsgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
              obsind, obsind_rownames, sizeindzi, sizeind_rownames_zi, false,
              obssigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 2, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0);
            //Rcout << "mazurekd s4 ";
            out_vec(1) = out_vec(1) * vr2_dcorr;
            
          } else {
            out_vec(1) = obscoefs(0);
            out_vec(1) = out_vec(1) * vr2_dcorr;
          }
          if (err_check) out(i, 1) = out_vec(1);
          
          if (ob3(i) == 1 || obsdist == 5) {
            
            if (sizedist < 5) {
              bool used_sizezero = false;
              if (sizezero && sz3(i) == 0) used_sizezero = sizezero;
              //Rcout << "mazurekd s5 ";
              out_vec(3) = preouterator(sizeproxy, sizecoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, sizegroups2, sizegroups1, sizegroups2zi, sizegroups1zi,
                sizeyearzi, sizepatchzi, sizeind, sizeind_rownames, sizeindzi,
                sizeind_rownames_zi, used_sizezero, sizesigma, grp2o(i), grp1(i),
                patchnumber, yearnumber, sizedist, 3, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
                sizetrunc);
              //Rcout << "mazurekd s6 ";
              out_vec(3) = out_vec(3) * vr3_dcorr;
              
            } else {
              out_vec(3) = 1.0;
              out_vec(3) = out_vec(3) * vr3_dcorr;
            }
            if (err_check) out(i, 3) = out_vec(3);
            
            if (sizebdist < 5) {
              bool used_sizebzero = false;
              if (sizebzero && szb3(i) == 0) used_sizebzero = sizebzero;
              //Rcout << "mazurekd s7 ";
              out_vec(4) = preouterator(sizebproxy, sizebcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, sizebgroups2, sizebgroups1, sizebgroups2zi,
                sizebgroups1zi, sizebyearzi, sizebpatchzi, sizebind,
                sizebind_rownames, sizebindzi, sizebind_rownames_zi, used_sizebzero,
                sizebsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist,
                4, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
                negfec, stage2n(i), nostages, sizebtrunc);
              //Rcout << "mazurekd s8 ";
              out_vec(4) = out_vec(4) * vr4_dcorr;
            } else {
              out_vec(4) = 1.0;
              out_vec(4) = out_vec(4) * vr4_dcorr;
            }
            if (err_check) out(i, 4) = out_vec(4);
            
            if (sizecdist < 5) {
              bool used_sizeczero = false;
              if (sizeczero && szc3(i) == 0) used_sizeczero = sizeczero;
              //Rcout << "mazurekd s9 ";
              out_vec(5) = preouterator(sizecproxy, sizeccoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, sizecgroups2, sizecgroups1, sizecgroups2zi,
                sizecgroups1zi, sizecyearzi, sizecpatchzi, sizecind,
                sizecind_rownames, sizecindzi, sizecind_rownames_zi, used_sizeczero,
                sizecsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist,
                5, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
                negfec, stage2n(i), nostages, sizectrunc);
              //Rcout << "mazurekd s10 ";
              out_vec(5) = out_vec(5) * vr5_dcorr;
            } else {
              out_vec(5) = 1.0;
              out_vec(5) = out_vec(5) * vr5_dcorr;
            }
            if (err_check) out(i, 5) = out_vec(5);
            
            if (repstdist < 5) {
              //Rcout << "mazurekd s11 ";
              out_vec(2) = preouterator(repstproxy, repstcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, repstgroups2, repstgroups1, dud_groups2zi, dud_groups1zi,
                dud_yearzi, dud_patchzi, repstind, repstind_rownames, sizeindzi,
                sizeind_rownames_zi, false, repstsigma, grp2o(i), grp1(i),
                patchnumber, yearnumber, 4, 6, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages, 0);
              //Rcout << "mazurekd s12 ";
              out_vec(2) = out_vec(2) * vr6_dcorr;
                
              if (fl3(i) == 0) {
                out_vec(2) = 1.0 - out_vec(2);
              }
            } else {
              if (fl3(i) == 0) {
                out_vec(2) = repstcoefs(0);
                out_vec(2) = out_vec(2) * vr6_dcorr;
                out_vec(2) = 1.0 - out_vec(2);
              } else if (fl3(i) == 1) {
                out_vec(2) = repstcoefs(0);
                out_vec(2) = out_vec(2) * vr6_dcorr;
              } else {
                out_vec(2) = 0.0;
              }
            }
            if (err_check) out(i, 2) = out_vec(2);
            
          } else {
            out_vec(1) = 1.0 - out_vec(1);
            out_vec(2) = 1.0;
            out_vec(3) = 1.0;
            out_vec(4) = 1.0;
            out_vec(5) = 1.0;
            out_vec(6) = 1.0;
            
            if (err_check) {
              out(i, 1) = out_vec(1);
              out(i, 2) = out_vec(2);
              out(i, 3) = out_vec(3);
              out(i, 4) = out_vec(4);
              out(i, 5) = out_vec(5);
              out(i, 6) = out_vec(6);
            }
          }
          if (!sparse) {
            survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
            out_vec(4) * out_vec(5) * out_vec(6);
          } else {
            survtransmat_sp(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
            out_vec(4) * out_vec(5) * out_vec(6);
          }
          
        } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurv_coefsadded != 0.0) {
          //Rcout << "mazurekd s13 ";
          // Juvenile to adult transitions
          if (jmatstdist < 5) {
            mat_predicted = preouterator(jmatstproxy, jmatstcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, jmatstgroups2, jmatstgroups1, dud_groups2zi,
              dud_groups1zi, dud_yearzi, dud_patchzi, jmatstind,
              jmatstind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              jmatstsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 21,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, 0);
            mat_predicted = mat_predicted * vr14_dcorr;
            
            if (mat3(i) > 0.5) {
              out_vec(6) = mat_predicted;
            } else {
              out_vec(6) = 1 - mat_predicted;
            }
          } else {
            if (mat3(i) > 0.5) {
              out_vec(6) = vr14_dcorr;
            } else {
              out_vec(6) = 1 - vr14_dcorr;
            }
          }
          if (err_check) out(i, 6) = out_vec(6);
          
          if (jsurvdist < 5) {
            out_vec(0) = preouterator(jsurvproxy, jsurvcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, jsurvgroups2, jsurvgroups1, dud_groups2zi,
              dud_groups1zi, dud_yearzi, dud_patchzi, jsurvind,
              jsurvind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              jsurvsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 8,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, 0);
            out_vec(0) = out_vec(0) * vr8_dcorr;
          } else {
            out_vec(0) = jsurvcoefs(0);
            out_vec(0) = out_vec(0) * vr8_dcorr;
          }
          if (err_check) out(i, 0) = out_vec(0);
          
          if (jobsdist < 5) {
            out_vec(1) = preouterator(jobsproxy, jobscoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, jobsgroups2, jobsgroups1, dud_groups2zi,
              dud_groups1zi, dud_yearzi, dud_patchzi, jobsind, jobsind_rownames,
              jsizeindzi, jsizeind_rownames_zi, false, jobssigma, grp2o(i),
              grp1(i), patchnumber, yearnumber, 4, 9, exp_tol, theta_tol,
              ipm_cdf, matrixformat, fecmod, repentry(i), negfec, stage2n(i),
              nostages, 0);
            out_vec(1) = out_vec(1) * vr9_dcorr;
          } else {
            out_vec(1) = jobscoefs(0);
            out_vec(1) = out_vec(1) * vr9_dcorr;
          }
          if (err_check) out(i, 1) = out_vec(1);
          
          if (ob3(i) == 1 || jobsdist == 5) {
            if (jsizedist < 5) {
              out_vec(3) = preouterator(jsizeproxy, jsizecoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, jsizegroups2, jsizegroups1, jsizegroups2zi,
                jsizegroups1zi, jsizeyearzi, jsizepatchzi, jsizeind,
                jsizeind_rownames, jsizeindzi, jsizeind_rownames_zi, jsizezero,
                jsizesigma, grp2o(i), grp1(i), patchnumber, yearnumber,
                sizedist, 10, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
                repentry(i), negfec, stage2n(i), nostages, jsizetrunc);
              out_vec(3) = out_vec(3) * vr10_dcorr;
            } else {
              out_vec(3) = 1.0;
              out_vec(3) = out_vec(3) * vr10_dcorr;
            }
            if (err_check) out(i, 3) = out_vec(3);
            
            if (jsizebdist < 5) {
              out_vec(4) = preouterator(jsizebproxy, jsizebcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, jsizebgroups2, jsizebgroups1, jsizebgroups2zi,
                jsizebgroups1zi, jsizebyearzi, jsizebpatchzi, jsizebind,
                jsizebind_rownames, jsizebindzi, jsizebind_rownames_zi,
                jsizebzero, jsizebsigma, grp2o(i), grp1(i), patchnumber,
                yearnumber, sizebdist, 11, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
                jsizebtrunc);
              out_vec(4) = out_vec(4) * vr11_dcorr;
            } else {
              out_vec(4) = 1.0;
              out_vec(4) = out_vec(4) * vr11_dcorr;
            }
            if (err_check) out(i, 4) = out_vec(4);
            
            if (jsizecdist < 5) {
              out_vec(5) = preouterator(jsizecproxy, jsizeccoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda,chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, jsizecgroups2, jsizecgroups1, jsizecgroups2zi,
                jsizecgroups1zi, jsizecyearzi, jsizecpatchzi, jsizecind,
                jsizecind_rownames, jsizecindzi, jsizecind_rownames_zi,
                jsizeczero, jsizecsigma, grp2o(i), grp1(i), patchnumber,
                yearnumber, sizecdist, 12, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
                jsizectrunc);
              out_vec(5) = out_vec(5) * vr12_dcorr;
            } else {
              out_vec(5) = 1.0;
              out_vec(5) = out_vec(5) * vr12_dcorr;
            }
            if (err_check) out(i, 5) = out_vec(5);
            
            if (jrepstdist < 5) {
              out_vec(2) = preouterator(jrepstproxy, jrepstcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda,
                chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
                chosen_r1indc, chosen_f2inda_cat, chosen_f1inda_cat,
                chosen_f2indb_cat, chosen_f1indb_cat, chosen_f2indc_cat,
                chosen_f1indc_cat, statusterms, jrepstgroups2, jrepstgroups1,
                dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
                jrepstind, jrepstind_rownames, jsizeindzi, jsizeind_rownames_zi,
                false, jrepstsigma, grp2o(i), grp1(i), patchnumber, yearnumber,
                4, 13, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
                repentry(i), negfec, stage2n(i), nostages, 0);
              out_vec(2) = out_vec(2) * vr13_dcorr;
                
              if (fl3(i) == 0) {
                out_vec(2) = 1.0 - out_vec(2);
              }
            } else {
              if (fl3(i) == 0) {
                out_vec(2) = jrepstcoefs(0);
                out_vec(2) = out_vec(2) * vr13_dcorr;
                out_vec(2) = 1.0 - out_vec(2);
              } else if (fl3(i) == 1) {
                out_vec(2) = jrepstcoefs(0);
                out_vec(2) = out_vec(2) * vr13_dcorr;
              } else {
                out_vec(2) = 0.0;
              }
            }
            if (err_check) out(i, 2) = out_vec(2);
            
          } else {
            out_vec(1) = 1.0 - out_vec(1);
            out_vec(2) = 1.0;
            out_vec(3) = 1.0;
            out_vec(4) = 1.0;
            out_vec(5) = 1.0;
            out_vec(6) = 1.0;
            
            if (err_check) {
              out(i, 1) = out_vec(1);
              out(i, 2) = out_vec(2);
              out(i, 3) = out_vec(3);
              out(i, 4) = out_vec(4);
              out(i, 5) = out_vec(5);
              out(i, 6) = out_vec(6);
            }
          }
          
          if (!sparse) {
            survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
              out_vec(4) * out_vec(5) * out_vec(6);
          } else {
            survtransmat_sp(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
              out_vec(4) * out_vec(5) * out_vec(6);
          }
          //Rcout << "mazurekd s14 ";
        }
      } else if (ovgivent(i) != -1.0) {
        // All other transitions
        if (!sparse) {
          survtransmat(k) = ovgivent(i);
        } else {
          survtransmat_sp(k) = ovgivent(i);
        }
      }
      
      // Fecundity calculation
      if (indata2n(i) == 1 && fec_addedcoefs != 0.0 && repentry(i) > 0) {
        if (fl2o(i) > 0.0 && ovgivenf(i) == -1.0) {
          
          if (!sparse) {
            fectransmat(k) = preouterator(fecproxy, feccoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, fecgroups2, fecgroups1, fecgroups2zi, fecgroups1zi,
              fecyearzi, fecpatchzi, fecind, fecind_rownames, fecindzi,
              fecind_rownames_zi, feczero, fecsigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, fecdist, 7, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
              fectrunc);
            fectransmat(k) = fectransmat(k) * vr7_dcorr;
          } else {
            fectransmat_sp(k) = preouterator(fecproxy, feccoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, fecgroups2, fecgroups1, fecgroups2zi, fecgroups1zi,
              fecyearzi, fecpatchzi, fecind, fecind_rownames, fecindzi,
              fecind_rownames_zi, feczero, fecsigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, fecdist, 7, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
              fectrunc);
            fectransmat_sp(k) = fectransmat_sp(k) * vr7_dcorr;
          }
          
        } else if (ovgivenf(i) != -1.0) {
          if (!sparse) {
            fectransmat(k) = ovgivenf(i);
            fectransmat(k) = fectransmat(k) * vr7_dcorr;
          } else {
            fectransmat_sp(k) = ovgivenf(i);
            fectransmat_sp(k) = fectransmat_sp(k) * vr7_dcorr;
          }
        }
      } else if (ovgivenf(i) != -1.0) {
        if (!sparse) {
          fectransmat(k) = ovgivenf(i);
          fectransmat(k) = fectransmat(k) * vr7_dcorr;
        } else {
          fectransmat_sp(k) = ovgivenf(i);
          fectransmat_sp(k) = fectransmat_sp(k) * vr7_dcorr;
        }
      }
    }
    
    //Rcout << "mazurekd t ";
    
    double ov_mult {0.0};
    if (replacementst > 0) {
      for (int i = 0; i < replacementst; i++) {
        
        repindex = replacetvec(i); // AllStages index
        properindex = aliveandequal(repindex);
        arma::uvec rightindex = find(index321 == ovestt(repindex));
        
        if (rightindex.n_elem > 0) {
          proxyindex = aliveandequal(rightindex(0));
          
          ov_mult = ovsurvmult(repindex);
          if (ov_mult < 0.0) ov_mult = 1.0;
          
          if (!sparse) {
            survtransmat(properindex) = survtransmat(proxyindex) * ov_mult;
          } else {
            survtransmat_sp(properindex) = survtransmat_sp(proxyindex) * ov_mult;
          }
        }
      }
    }
    
    //Rcout << "mazurekd u ";
    
    if (replacementsf > 0) {
      for (int i = 0; i < replacementsf; i++) {
        
        repindex = replacefvec(i); // AllStages index
        properindex = aliveandequal(repindex);
        arma::uvec rightindex = find(index321 == ovestf(repindex));
        
        if (rightindex.n_elem > 0) {
          proxyindex = aliveandequal(rightindex(0));
          
          ov_mult = ovfecmult(repindex);
          if (ov_mult < 0.0) ov_mult = 1.0;
          if (!sparse) {
            fectransmat(properindex) = fectransmat(proxyindex) * ov_mult;
          } else {
            fectransmat_sp(properindex) = fectransmat_sp(proxyindex) * ov_mult;
          }
        }
      }
    }
    
    //Rcout << "mazurekd v ";
    
    if (tmults_only_st > 0) {
      for (int i = 0; i < tmults_only_st; i++) {
        repindex = tmults_only(i);
        properindex = aliveandequal(repindex);
        ov_mult = ovsurvmult(repindex);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        if (!sparse) {
          survtransmat(properindex) = survtransmat(properindex) * ov_mult;
        } else {
          survtransmat_sp(properindex) = survtransmat_sp(properindex) * ov_mult;
        }
      }
    }
    
    //Rcout << "mazurekd w ";
    
    if (fmults_only_st > 0) {
      for (int i = 0; i < fmults_only_st; i++) {
        repindex = fmults_only(i);
        properindex = aliveandequal(repindex);
        ov_mult = ovfecmult(repindex);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        if (!sparse) {
          fectransmat(properindex) = fectransmat(properindex) * ov_mult;
        } else {
          fectransmat_sp(properindex) = fectransmat_sp(properindex) * ov_mult;
        }
      }
    }
    
    //Rcout << "mazurekd x ";
    
    // Final output
    List output;
    
    if (A_only) {
      List output_pre (2);
      
      if (!sparse) {
        arma::mat amatrix = survtransmat + fectransmat;
        output_pre(0) = amatrix;
        
        if (err_check) {
          output_pre(1) = out;
        } else {
          output_pre(1) = R_NilValue;
        }
      } else {
        arma::sp_mat amatrix = survtransmat_sp + fectransmat_sp;
        output_pre(0) = amatrix;
        
        if (err_check) {
          output_pre(1) = out;
        } else {
          output_pre(1) = R_NilValue;
        }
      }
      
      CharacterVector output_names = {"A", "out"};
      output_pre.attr("names") = output_names;
      
      output = output_pre;
    } else {
      List output_pre (4);
      
      if (!sparse) {
        arma::mat amatrix = survtransmat + fectransmat;
        output_pre(0) = amatrix;
        output_pre(1) = survtransmat;
        output_pre(2) = fectransmat;
        
        if (err_check) {
          output_pre(3) = out;
        } else {
          output_pre(3) = R_NilValue;
        }
      } else {
        arma::sp_mat amatrix_sp = survtransmat_sp + fectransmat_sp;
        output_pre(0) = amatrix_sp;
        output_pre(1) = survtransmat_sp;
        output_pre(2) = fectransmat_sp;
        
        if (err_check) {
          output_pre(3) = out;
        } else {
          output_pre(3) = R_NilValue;
        }
      }
      
      CharacterVector output_names = {"A", "U", "F", "out"};
      output_pre.attr("names") = output_names;
      
      output = output_pre;
    }
    
    //Rcout << "mazurekd y";
    return output;
  }
  
  //' Estimate All Elements of Function-based Leslie Population Projection Matrix
  //' 
  //' Function \code{mdabrowskiego()} swiftly calculates matrix elements in
  //' function-based Leslie population projection matrices.
  //' 
  //' @name mdabrowskiego
  //' 
  //' @param actualages An integer vector of all ages to be included in the
  //' matrices, in order.
  //' @param ageframe The modified stageframe used in matrix calculations.
  //' @param survproxy List of coefficients estimated in model of survival.
  //' @param fecproxy List of coefficients estimated in model of fecundity.
  //' @param f2_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param r2_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{a} at each time \emph{t} to be used in analysis.
  //' @param r1_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{a} at each time \emph{t}-1 to be used in analysis.
  //' @param r2_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{b} at each time \emph{t} to be used in analysis.
  //' @param r1_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{b} at each time \emph{t}-1 to be used in analysis.
  //' @param r2_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{c} at each time \emph{t} to be used in analysis.
  //' @param r1_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual random covariate
  //' \code{c} at each time \emph{t}-1 to be used in analysis.
  //' @param surv_dev A numeric value indicating the deviation to the linear
  //' model of survival input by the user.
  //' @param fec_dev A numeric value indicating the deviation to the linear
  //' model of fecundity input by the user.
  //' @param dens A numeric value equal to the density to be used in calculations.
  //' @param fecmod A scalar multiplier for fecundity.
  //' @param finalage The final age to be included in Leslie MPM estimation.
  //' @param negfec A logical value denoting whether to change negative estimated
  //' fecundity to \code{0}.
  //' @param yearnumber An integer specifying which time at time \emph{t} to
  //' develop matrices for. Must be in reference to the \code{listofyears} object
  //' developed in the \code{R} matrix estimator function.
  //' @param patchnumber An integer specifying which patch to develop matrices
  //' for. Must be in reference to the \code{listofyears} object developed in the
  //' \code{R} matrix estimator function.
  //' @param dens_vr A logical value indicating whether any vital rates are
  //' density dependent.
  //' @param dvr_yn A logical vector indicating whether each vital rate is density
  //' dependent.
  //' @param dvr_style An integer vector indicating the style of density
  //' dependence for each vital rate.
  //' @param dvr_alpha A numeric vector indicating the value of alpha to use in
  //' density dependence for each vital rate.
  //' @param dvr_beta A numeric vector indicating the value of beta to use in
  //' density dependence for each vital rate.
  //' @param dens_n A numeric vector corresponding to the population size to use
  //' in vital rate density dependence calculations.
  //' @param exp_tol A numeric value indicating the maximum limit for the
  //' \code{exp()} function to be used in vital rate calculations. Defaults to
  //' \code{700.0}.
  //' @param theta_tol A numeric value indicating a maximum value for theta in
  //' negative binomial probability density estimation. Defaults to
  //' \code{100000000.0}.
  //' @param sparse If \code{TRUE}, then only outputs matrices in sparse format.
  //' Defaults to \code{FALSE}.
  //' @param supplement An optional data frame edited and age-expanded showing
  //' supplemental transition information.
  //' 
  //' @return A list of 3 matrices, including the main MPM (A), the survival-
  //' transition matrix (U), and a fecundity matrix (F).
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List mdabrowskiego(IntegerVector actualages,
    const DataFrame& ageframe, const List& survproxy, const List& fecproxy,
    NumericVector f2_inda, NumericVector f1_inda, NumericVector f2_indb,
    NumericVector f1_indb, NumericVector f2_indc, NumericVector f1_indc,
    StringVector r2_inda, StringVector r1_inda, StringVector r2_indb,
    StringVector r1_indb, StringVector r2_indc, StringVector r1_indc,
    double surv_dev, double fec_dev, double dens, double fecmod,
    unsigned int finalage, bool negfec, int yearnumber, int patchnumber,
    bool dens_vr, LogicalVector dvr_yn, IntegerVector dvr_style,
    NumericVector dvr_alpha, NumericVector dvr_beta, NumericVector dens_n,
    double exp_tol = 700.0, double theta_tol = 100000000.0,
    bool sparse = false, Nullable<DataFrame> supplement = R_NilValue) {
    
    //Rcout << "mdabrowskiego a        ";
    
    // Determines the size of the matrix
    StringVector sf_agenames = as<StringVector>(ageframe["stage"]);
    IntegerVector sf_minage = as<IntegerVector>(ageframe["min_age"]);
    IntegerVector sf_maxage = as<IntegerVector>(ageframe["max_age"]);
    IntegerVector sf_repstatus = as<IntegerVector>(ageframe["repstatus"]);
    int noages = actualages.length();
    int start_age = min(actualages);
    
    bool cont = false;
    if (sf_maxage(noages - 1) == NA_INTEGER) {
      cont = true;
    }
    
    //Rcout << "mdabrowskiego b        ";
    
    // Supplement processing
    IntegerVector ov_age2;
    IntegerVector ov_estage2;
    NumericVector ov_givenrate;
    NumericVector ov_multiplier;
    IntegerVector ov_convtype;
    int supp_length {0};
    
    DataFrame supplement_;
    if (supplement.isNotNull()){
      supplement_ = as<DataFrame>(supplement);
      
      if (supplement_.containsElementNamed("age2")) {
        ov_age2 = clone(as<IntegerVector>(supplement_["age2"]));
        ov_estage2 = clone(as<IntegerVector>(supplement_["estage2"]));
        ov_givenrate = as<NumericVector>(supplement_["givenrate"]);
        ov_multiplier = as<NumericVector>(supplement_["multiplier"]);
        ov_convtype = as<IntegerVector>(supplement_["convtype"]);
        
        supp_length = static_cast<int>(ov_givenrate.length());
        
        for (int i = 0; i < supp_length; i++) {
          ov_age2(i) = ov_age2(i) - start_age;
          
          if (!IntegerVector::is_na(ov_estage2(i))) {
            ov_estage2(i) = ov_estage2(i) - start_age;
          }
        }
      }
    }
    
    //Rcout << "mdabrowskiego c        ";
    
    // Proxy model imports and settings
    NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
    NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
    
    bool feczero = as<bool>(fecproxy["zero_inflated"]);
    int survdist = as<int>(survproxy["dist"]);
    int fecdist = as<int>(fecproxy["dist"]);
    double fecsigma = as<double>(fecproxy["sigma"]);
    
    if (NumericVector::is_na(fecsigma)) {
      if (fecdist == 1) {
        fecsigma = 1.0;
      } else {
        fecsigma = 0.0;
      }
    }
    
    NumericMatrix vital_year = revelations_leslie(survproxy, fecproxy, 1);
    NumericMatrix vital_patch = revelations_leslie(survproxy, fecproxy, 2);
    
    NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
    NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
    NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
    NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
    NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
    NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
    NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
    NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
    
    NumericVector survind = flightoficarus(survproxy);
    NumericVector fecind = flightoficarus(fecproxy);
    NumericVector fecindzi = zero_flightoficarus(fecproxy);
    
    arma::imat rand_index = foi_index_leslie(survproxy, fecproxy);
    
    StringVector survind_rownames = bootson(survproxy);
    StringVector fecind_rownames = bootson(fecproxy);
    StringVector fecind_rownames_zi = zero_bootson(fecproxy);
    
    //Rcout << "mdabrowskiego d        ";
    
    // Determination of choices of fixed and random individual covariates
    double inda1 = f1_inda(0);
    double indb1 = f1_indb(0);
    double indc1 = f1_indc(0);
    double inda2 = f2_inda(0);
    double indb2 = f2_indb(0);
    double indc2 = f2_indc(0);
    
    String chosen_f2inda_cat("none");
    String chosen_f1inda_cat("none");
    String chosen_f2indb_cat("none");
    String chosen_f1indb_cat("none");
    String chosen_f2indc_cat("none");
    String chosen_f1indc_cat("none");
    
    String chosen_r2inda = r2_inda(0);
    String chosen_r1inda = r1_inda(0);
    String chosen_r2indb = r2_indb(0);
    String chosen_r1indb = r1_indb(0);
    String chosen_r2indc = r2_indc(0);
    String chosen_r1indc = r1_indc(0);
    
    //Rcout << "mdabrowskiego e        ";
    
    // The output matrices
    arma::mat survtransmat;
    arma::mat fectransmat;
    arma::sp_mat survtransmat_sp;
    arma::sp_mat fectransmat_sp;
    
    if (!sparse) {
      arma::mat survtransmat_chuck (noages, noages, fill::zeros);
      arma::mat fectransmat_chuck (noages, noages, fill::zeros);
      
      survtransmat = survtransmat_chuck;
      fectransmat = fectransmat_chuck;
    } else { 
      arma::sp_mat survtransmat_chuck (noages, noages);
      arma::sp_mat fectransmat_chuck (noages, noages);
      
      survtransmat_sp = survtransmat_chuck;
      fectransmat_sp = fectransmat_chuck;
    }
    
    //Rcout << "mdabrowskiego f        ";
    
    // Following loop runs through each age, and so runs through
    // each estimable element in the matrix
    double fec_addedcoefs = sum(feccoefs);
    for(int i = 0; i < noages; i++) {
      // Adult survival transitions
      
      double preout {0.0};
      
      if (survdist < 5) {
        //Rcout << "mdabrowskiego f2        ";
        
        // Random covariate processing
        double chosen_randcova2 {0.0};
        if (chosen_r2inda != "none") {
          for (int indcount = 0; indcount < rand_index(0, 0); indcount++) {
            if (chosen_r2inda == survind_rownames(indcount)) {
              chosen_randcova2 = survind(indcount);
            }
          }
        }
        double chosen_randcova1 {0.0};
        if (chosen_r1inda != "none") {
          int delectable_sum = rand_index(0, 0);
          for (int indcount = 0; indcount < rand_index(1, 0); indcount++) {
            if (chosen_r1inda == survind_rownames(indcount + delectable_sum)) {
              chosen_randcova1 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb2 {0.0};
        if (chosen_r2indb != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0);
          for (int indcount = 0; indcount < rand_index(2, 0); indcount++) {
            if (chosen_r2indb == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovb2 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb1 {0.0};
        if (chosen_r1indb != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0);
          for (int indcount = 0; indcount < rand_index(3, 0); indcount++) {
            if (chosen_r1indb == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovb1 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc2 {0.0};
        if (chosen_r2indc != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
            rand_index(3, 0);
          for (int indcount = 0; indcount < rand_index(4, 0); indcount++) {
            if (chosen_r2indc == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovc2 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc1 {0.0};
        if (chosen_r1indc != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
            rand_index(3, 0) + rand_index(4, 0);
          for (int indcount = 0; indcount < rand_index(5, 0); indcount++) {
            if (chosen_r1indc == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovc1 = survind(indcount + delectable_sum);
            }
          }
        }
        
        //Rcout << "mdabrowskiego f3        ";
        
        // Fixed factor covariate processing
        double chosen_fixcova2 {0.0};
        double chosen_fixcova1 {0.0};
        double chosen_fixcovb2 {0.0};
        double chosen_fixcovb1 {0.0};
        double chosen_fixcovc2 {0.0};
        double chosen_fixcovc1 {0.0};
        
        double mainsum = rimeotam(survcoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, static_cast<double>(actualages(i)), inda1, inda2, indb1, indb2,
          indc1, indc2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dens, false); // Final 6 zeroes are anncov
        
        //Rcout << "mdabrowskiego f4        ";
        
        preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
          chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
          chosen_randcovc1 + chosen_fixcova2 + chosen_fixcova1 +
          chosen_fixcovb2 + chosen_fixcovb1 + chosen_fixcovc2 +
          chosen_fixcovc1 + survgroups2(0) + survgroups1(0) + 
          vital_patch(patchnumber, 0) + vital_year(yearnumber, 0) + surv_dev);
        
        //Rcout << "mdabrowskiego f5        ";
        
        if (preout > exp_tol) preout = exp_tol; // Catches numbers too high
        if (i < (noages - 1)) {
          if (!sparse) {
            survtransmat(i+1, i) = exp(preout) / (1.0 + exp(preout));
          } else {
            survtransmat_sp(i+1, i) = exp(preout) / (1.0 + exp(preout));
          }
        } else {
          if (cont) {
            if (!sparse) {
              survtransmat(i, i) = exp(preout) / (1.0 + exp(preout));
            } else {
              survtransmat_sp(i, i) = exp(preout) / (1.0 + exp(preout));
            }
          }
        }
        //Rcout << "mdabrowskiego f6        ";
        
      } else {
        //Rcout << "mdabrowskiego f7        ";
        
        if (i < (noages - 1)) {
          if (!sparse) {
            survtransmat(i+1, i) = survcoefs(0);
          } else {
            survtransmat_sp(i+1, i) = survcoefs(0);
          }
        } else {
          if (!sparse) {
            survtransmat(i, i) = survcoefs(0);
          } else {
            survtransmat_sp(i, i) = survcoefs(0);
          }
        }
      }
      
        //Rcout << "mdabrowskiego f8        ";
        
      // This next block calculates fecundity
      if (fec_addedcoefs != 0.0) {
        if (sf_repstatus(i) == 1) {
          
        //Rcout << "mdabrowskiego f9        ";
        
          // Random covariate processing
          double chosen_randcova2 {0.0};
          if (chosen_r2inda != "none") {
            for (int indcount = 0; indcount < rand_index(0, 6); indcount++) {
              if (chosen_r2inda == fecind_rownames(indcount)) {
                chosen_randcova2 = fecind(indcount);
              }
            }
          }
          double chosen_randcova1 {0.0};
          if (chosen_r1inda != "none") {
            int delectable_sum = rand_index(0, 6);
            for (int indcount = 0; indcount < rand_index(1, 6); indcount++) {
              if (chosen_r1inda == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcova1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb2 {0.0};
          if (chosen_r2indb != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6);
            for (int indcount = 0; indcount < rand_index(2, 6); indcount++) {
              if (chosen_r2indb == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovb2 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb1 {0.0};
          if (chosen_r1indb != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6);
            for (int indcount = 0; indcount < rand_index(3, 6); indcount++) {
              if (chosen_r1indb == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovb1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc2 {0.0};
          if (chosen_r2indc != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
              rand_index(3, 6);
            for (int indcount = 0; indcount < rand_index(4, 6); indcount++) {
              if (chosen_r2indc == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovc2 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc1 {0.0};
          if (chosen_r1indc != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
              rand_index(3, 6) + rand_index(4, 6);
            for (int indcount = 0; indcount < rand_index(5, 6); indcount++) {
              if (chosen_r1indc == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovc1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcova2zi {0.0};
          if (chosen_r2inda != "none") {
            for (int indcount = 0; indcount < rand_index(0, 16); indcount++) {
              if (chosen_r2inda == fecind_rownames_zi(indcount)) {
                chosen_randcova2zi = fecindzi(indcount);
              }
            }
          }
          double chosen_randcova1zi {0.0};
          if (chosen_r1inda != "none") {
            int delectable_sum = rand_index(0, 16);
            for (int indcount = 0; indcount < rand_index(1, 16); indcount++) {
              if (chosen_r1inda == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcova1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb2zi {0.0};
          if (chosen_r2indb != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16);
            for (int indcount = 0; indcount < rand_index(2, 16); indcount++) {
              if (chosen_r2indb == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovb2zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb1zi {0.0};
          if (chosen_r1indb != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16);
            for (int indcount = 0; indcount < rand_index(3, 16); indcount++) {
              if (chosen_r1indb == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovb1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc2zi {0.0};
          if (chosen_r2indc != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
              rand_index(3, 16);
            for (int indcount = 0; indcount < rand_index(4, 16); indcount++) {
              if (chosen_r2indc == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovc2zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc1zi {0.0};
          if (chosen_r1indc != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
              rand_index(3, 16) + rand_index(4, 16);
            for (int indcount = 0; indcount < rand_index(5, 16); indcount++) {
              if (chosen_r1indc == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovc1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          //Rcout << "mdabrowskiego f10        ";
          
          
          // Fixed factor covariate processing
          double chosen_fixcova2 {0.0};
          double chosen_fixcova1 {0.0};
          double chosen_fixcovb2 {0.0};
          double chosen_fixcovb1 {0.0};
          double chosen_fixcovc2 {0.0};
          double chosen_fixcovc1 {0.0};
          
          double chosen_fixcova2zi {0.0};
          double chosen_fixcova1zi {0.0};
          double chosen_fixcovb2zi {0.0};
          double chosen_fixcovb1zi {0.0};
          double chosen_fixcovc2zi {0.0};
          double chosen_fixcovc1zi {0.0};
          
          double preoutx {0.0};
          
          //Rcout << "mdabrowskiego f11        ";
          
          if (fecdist < 4) {
            if (feczero) {
              
              double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, static_cast<double>(actualages(i)), inda1, inda2,
                indb1, indb2, indc1, indc2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dens,
                true); // Final 6 zeroes are anncov
              
              preoutx = (mainsum + chosen_randcova2zi + chosen_randcova1zi +
                chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
                chosen_randcovc1zi + chosen_fixcova2zi + chosen_fixcova1zi +
                chosen_fixcovb2zi + chosen_fixcovb1zi + chosen_fixcovc2zi +
                chosen_fixcovc1zi + fecgroups2zi(0) + fecgroups1zi(0) + 
                fecpatchzi(patchnumber) + fecyearzi(yearnumber) + fec_dev);
              
          //Rcout << "mdabrowskiego f11c        ";
          
            } else {
              
          //Rcout << "mdabrowskiego f11d        ";
          
              double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, static_cast<double>(actualages(i)), inda1, inda2,
                indb1, indb2, indc1, indc2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dens,
                false); // Final 6 zeroes are anncov
              
              preoutx = (mainsum + chosen_randcova2 + chosen_randcova1 +
                chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
                chosen_randcovc1 + chosen_fixcova2 + chosen_fixcova1 +
                chosen_fixcovb2 + chosen_fixcovb1 + chosen_fixcovc2 +
                chosen_fixcovc1 + fecgroups2(0) + fecgroups1(0) + 
                vital_patch(patchnumber, 1) + vital_year(yearnumber, 1) +
                fec_dev);
          //Rcout << "mdabrowskiego f11e        ";
          
            }
            
          //Rcout << "mdabrowskiego f11f        ";
          
            if (fecdist == 0 || fecdist == 1) {
              // Poisson and negative binomial fecundity
              
              if (feczero) {
                if (preoutx > exp_tol) preoutx = exp_tol;
                if (!sparse) {
                  fectransmat(0, i) = (exp(preoutx) / (1.0 + exp(preoutx))) * fecmod;
                } else {
                  fectransmat_sp(0, i) = (exp(preoutx) / (1.0 + exp(preoutx))) * fecmod;
                }
              } else {
                if (preoutx > exp_tol) preoutx = exp_tol;
                
                if (!sparse) {
                  fectransmat(0, i) = exp(preoutx) * fecmod;
                } else {
                  fectransmat_sp(0, i) = exp(preoutx) * fecmod;
                }
              }
          //Rcout << "mdabrowskiego f11g        ";
          
            } else if (fecdist == 2) {
          //Rcout << "mdabrowskiego f11h        ";
          
              // Gaussian fecundity
              
              if (!sparse) {
                fectransmat(0, i) = preoutx * fecmod;
                
                if (negfec && fectransmat(0, i) < 0.0) {
                  fectransmat(0, i) = 0.0;
                }
              } else { 
                fectransmat_sp(0, i) = preoutx * fecmod;
                
                if (negfec && fectransmat_sp(0, i) < 0.0) {
                  fectransmat_sp(0, i) = 0.0;
                }
              }
          //Rcout << "mdabrowskiego f11i        ";
          
            } else if (fecdist == 3) {
          //Rcout << "mdabrowskiego f11j        ";
          
              // Gamma fecundity
              if (!sparse) {
                fectransmat(0, i) = (1.0 / preoutx) * fecmod;
              } else {
                fectransmat_sp(0, i) = (1.0 / preoutx) * fecmod;
              }
            }
            
          //Rcout << "mdabrowskiego f11k        ";
          
          } else {
            if (!sparse) {
              fectransmat(0, i) = feccoefs(0);
            } else {
              fectransmat_sp(0, i) = feccoefs(0);
            }
          }
        }
      }
    }
    
    //Rcout << "mdabrowskiego g        ";
    
    // Supplement replacement portion
    int target_col {0};
    int target_row {0};
    int proxy_col {0};
    int proxy_row {0};
    
    for (int l = 0; l < supp_length; l++) {
      target_col = ov_age2(l);
      if (target_col >= (noages - 1) && cont) {
        target_col = noages - 1;
      } else if (target_col >= (noages - 1)) {
        throw Rcpp::exception("Some age2 values are too high.", false);
      }
      
      if (ov_convtype(l) == 1) {
        if (target_col >= (noages - 1) && cont) {
          target_row = target_col;
        } else {
          target_row = target_col + 1;
        }
        
        if (!NumericVector::is_na(ov_givenrate(l))) {
          if (!sparse) {
            survtransmat(target_row, target_col) = ov_givenrate(l);
          } else {
            survtransmat_sp(target_row, target_col) = ov_givenrate(l);
          }
        }
        if (!IntegerVector::is_na(ov_estage2(l))) {
          proxy_col = ov_estage2(l);
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_col = noages - 1;
          } else if (proxy_col >= (noages - 1)) {
            throw Rcpp::exception("Some estage2 values are too high.", false);
          }
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_row = proxy_col; // Was target_col
          } else {
            proxy_row = proxy_col + 1;
          }
          
          if (!sparse) {
            survtransmat(target_row, target_col) = survtransmat(proxy_row, proxy_col);
          } else {
            survtransmat_sp(target_row, target_col) = survtransmat_sp(proxy_row, proxy_col);
          }
        }
        if (!NumericVector::is_na(ov_multiplier(l))) {
          if (!sparse) {
            survtransmat(target_row, target_col) *= ov_multiplier(l);
          } else {
            survtransmat_sp(target_row, target_col) *= ov_multiplier(l);
          }
        }
      } else if (ov_convtype(l) == 2) {
        target_row = 0;
        
        if (!NumericVector::is_na(ov_givenrate(l))) {
          if (!sparse) {
            fectransmat(target_row, target_col) = ov_givenrate(l);
          } else {
            fectransmat_sp(target_row, target_col) = ov_givenrate(l);
          }
        }
        if (!IntegerVector::is_na(ov_estage2(l))) {
          proxy_col = ov_estage2(l);
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_col = noages - 1;
          } else if (proxy_col >= (noages - 1)) {
            throw Rcpp::exception("Some estage2 values are too high.", false);
          }
          
          proxy_row = 0;
          
          if (!sparse) {
            fectransmat(target_row, target_col) = fectransmat(proxy_row, proxy_col);
          } else {
            fectransmat_sp(target_row, target_col) = fectransmat_sp(proxy_row, proxy_col);
          }
        }
        if (!NumericVector::is_na(ov_multiplier(l))) {
          if (!sparse) {
            fectransmat(target_row, target_col) *= ov_multiplier(l);
          } else {
            fectransmat_sp(target_row, target_col) *= ov_multiplier(l);
          }
        }
      } else {
        target_row = 0;
        
        if (!NumericVector::is_na(ov_multiplier(l))) {
          if (!sparse) {
            fectransmat(target_row, target_col) *= ov_multiplier(l);
          } else {
            fectransmat_sp(target_row, target_col) *= ov_multiplier(l);
          }
        }
      }
    }
    
    //Rcout << "mdabrowskiego h        ";
    
    List output;
    if (!sparse) {
      arma::mat amatrix = survtransmat + fectransmat;
      output = List::create(Named("A") = amatrix);
    } else {
      arma::sp_mat amatrix_sp = survtransmat_sp + fectransmat_sp;
      output = List::create(Named("A") = amatrix_sp);
    }
    
    return output;
  }
  
  //' Create Element Index for Matrix Estimation with Trait Variants
  //' 
  //' Function \code{thenewpizzle()} creates a data frame object used by 
  //' function \code{invade3()} to alter elements in existing lefkoMat objects.
  //' 
  //' @name thenewpizzle
  //'
  //' @param StageFrame The stageframe object identifying the life history model
  //' being operationalized.
  //' @param trait_axis The trait axis data frame that will inform matrix
  //' alteration.
  //' @param firstage The first age to be used in the analysis. Should typically
  //' be \code{0} for pre-breeding and \code{1} for post-breeding life history
  //' models. If not building age-by-stage MPMs, then should be set to \code{0}.
  //' @param finalage The final age to be used in analysis. If not building
  //' age-by-stage MPMs, then should be set to \code{0}.
  //' @param format Indicates whether historical matrices should be in (\code{1})
  //' Ehrlen or (\code{2}) deVries format.
  //' @param style The style of analysis, where \code{0} is historical, \code{1}
  //' is ahistorical, and \code{2} is age-by-stage.
  //' @param filter An integer denoting whether to filter the output data frame to
  //' eliminate unusable rows, and if so, how to do so. Possible values: \code{0}:
  //' no filtering, \code{1}: filter out rows with \code{index321 == -1}, and
  //' \code{2}: filter out rows with \code{aliveandequal == -1}.
  //' @param mpm_only A Boolean value stipulating whether to allow only rows
  //' altering existing MPMs into the final stage expansion table.
  //' 
  //' @return The output is a large data frame describing every element to be
  //' altered in matrices.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List thenewpizzle(const DataFrame& StageFrame,
    const DataFrame& trait_axis, int firstage, int finalage, int format,
    int style, int filter) {
    
    int stageexpansion_size {0};
    
    IntegerVector tavariant = as<IntegerVector>(trait_axis["variant"]);
    StringVector tastage3 = as<StringVector>(trait_axis["stage3"]);
    StringVector tastage2 = as<StringVector>(trait_axis["stage2"]);
    StringVector tastage1 = as<StringVector>(trait_axis["stage1"]);
    IntegerVector taage3 = as<IntegerVector>(trait_axis["age3"]);
    IntegerVector taage2 = as<IntegerVector>(trait_axis["age2"]);
    StringVector taeststage3 = as<StringVector>(trait_axis["eststage3"]);
    StringVector taeststage2 = as<StringVector>(trait_axis["eststage2"]);
    StringVector taeststage1 = as<StringVector>(trait_axis["eststage1"]);
    IntegerVector taestage3 = as<IntegerVector>(trait_axis["estage3"]);
    IntegerVector taestage2 = as<IntegerVector>(trait_axis["estage2"]);
    arma::vec tagivenrate = as<arma::vec>(trait_axis["givenrate"]);
    arma::vec taoffset = as<arma::vec>(trait_axis["offset"]);
    arma::vec tamultiplier = as<arma::vec>(trait_axis["multiplier"]);
    arma::vec taconvtype = as<arma::vec>(trait_axis["convtype"]);
    arma::vec taconvt12 = as<arma::vec>(trait_axis["convtype_t12"]);
    arma::uvec tampmaltered = as<arma::uvec>(trait_axis["mpm_altered"]);
    int tarows = static_cast<int>(tagivenrate.n_elem);
    
    int totalages {(finalage - firstage) + 1};
    
    arma::vec taindex3 (tarows);
    arma::vec taindex2 (tarows);
    arma::vec taindex1 (tarows);
    arma::vec tanew3 (tarows);
    arma::vec tanew2 (tarows);
    arma::vec tanew1 (tarows);
    arma::vec taindexold321 (tarows);
    arma::vec taindexnew321 (tarows);
    arma::vec tanewgivenrate (tarows);
    arma::vec tanewoffset (tarows, fill::zeros);
    arma::vec tanewmultiplier (tarows, fill::zeros);
    arma::vec taconvtypeage (tarows);
    arma::uvec tanewmpmaltered (tarows, fill::zeros);
    taindex3.fill(-1.0);
    taindex2.fill(-1.0);
    taindex1.fill(-1.0);
    tanew3.fill(-1.0);
    tanew2.fill(-1.0);
    tanew1.fill(-1.0);
    taindexold321.fill(-1.0);
    taindexnew321.fill(-1.0);
    tanewgivenrate.fill(-1.0);
    taconvtypeage.fill(-1.0);
    
    arma::ivec newstageid = as<arma::ivec>(StageFrame["stage_id"]);
    StringVector origstageid = as<StringVector>(StageFrame["stage"]);
    arma::vec binsizectr = as<arma::vec>(StageFrame["sizebin_center"]);
    arma::vec repstatus = as<arma::vec>(StageFrame["repstatus"]);
    arma::vec obsstatus = as<arma::vec>(StageFrame["obsstatus"]);
    arma::vec immstatus = as<arma::vec>(StageFrame["immstatus"]);
    arma::vec matstatus = as<arma::vec>(StageFrame["matstatus"]);
    arma::vec indata = as<arma::vec>(StageFrame["indataset"]);
    arma::vec binsizewidth = as<arma::vec>(StageFrame["sizebin_width"]);
    arma::vec alive = as<arma::vec>(StageFrame["alive"]);
    arma::vec minage = as<arma::vec>(StageFrame["min_age"]);
    arma::vec maxage = as<arma::vec>(StageFrame["max_age"]);
    arma::vec group = as<arma::vec>(StageFrame["group"]);
    arma::vec almostborn = as<arma::vec>(StageFrame["almostborn"]);
    
    arma::vec binsizebctr = as<arma::vec>(StageFrame["sizebinb_center"]);
    arma::vec binsizecctr = as<arma::vec>(StageFrame["sizebinc_center"]);
    arma::vec binsizebwidth = as<arma::vec>(StageFrame["sizebinb_width"]);
    arma::vec binsizecwidth = as<arma::vec>(StageFrame["sizebinc_width"]);
    
    int sf_length = static_cast<int>(origstageid.length());
    
    //Rcout << "thenewpizzle a        ";
    
    // Create stage order vectors
    IntegerVector tastageorder3 (tarows); // Replace calls to stageorder below
    IntegerVector tastageorder2 (tarows);
    IntegerVector tastageorder1 (tarows);
    IntegerVector tanewstageid3 (tarows); // Replace calls to stage below
    IntegerVector tanewstageid2 (tarows);
    IntegerVector tanewstageid1 (tarows);
    
    IntegerVector tastageorder3est (tarows); // Identify proxy elements
    IntegerVector tastageorder2est (tarows);
    IntegerVector tastageorder1est (tarows);
    IntegerVector tanewstageid3est (tarows); 
    IntegerVector tanewstageid2est (tarows);
    IntegerVector tanewstageid1est (tarows);
    
    for (int i = 0; i < tarows; i++) {
      for (int j = 0; j < sf_length; j++) {
        if (tastage3(i) == origstageid(j)) {
          tastageorder3(i) = (j + 1);
          tanewstageid3(i) = newstageid(j);
        }
        if (tastage2(i) == origstageid(j)) {
          tastageorder2(i) = (j + 1);
          tanewstageid2(i) = newstageid(j);
        }
        if (tastage1(i) == origstageid(j)) {
          tastageorder1(i) = (j + 1);
          tanewstageid1(i) = newstageid(j);
        }
        
        if (taeststage3(i) == origstageid(j)) {
          tastageorder3est(i) = (j + 1);
          tanewstageid3est(i) = newstageid(j);
        }
        if (taeststage2(i) == origstageid(j)) {
          tastageorder2est(i) = (j + 1);
          tanewstageid2est(i) = newstageid(j);
        }
        if (taeststage1(i) == origstageid(j)) {
          tastageorder1est(i) = (j + 1);
          tanewstageid1est(i) = newstageid(j);
        }
      }
    }
    
    //Rcout << "thenewpizzle b        ";
    
    // Determine length of matrix map data frame
    int nostages = static_cast<int>(newstageid.n_elem);
    int nostages_nodead = nostages - 1;
    int nostages_nounborn = nostages;
    //int nostages_nodead_nounborn = nostages_nodead;
    int prior_stage = -1;
    arma::vec tarepentry_prior(nostages, fill::zeros);
    int totallength {0};
    
    //Rcout << "thenewpizzle c        ";
    //Rcout << "nostages: " << nostages << "         ";
    
    stageexpansion_size = tarows; // Might need to alter for age-by-stage or deVries historical
    
    if (style == 2) { // Age-by-stage
      totallength = (nostages * nostages * totalages * totalages);
    } else if (style == 1) { // Ahistorical & age-based
      totallength = (nostages * nostages);
    } else { // Historical
      if (format == 2) { // deVries format historical
        //nostages_nodead_nounborn = nostages - 2;
        prior_stage = nostages_nounborn;
        nostages_nounborn = nostages - 1;
        totallength = (2 * nostages_nounborn * nostages_nounborn * nostages_nounborn);
      } else { // Ehrlen format historical
        totallength = (nostages * (nostages_nodead * nostages_nodead));
      }
    }
    
    totallength = stageexpansion_size; // Need to check this - might not need former totallength calculations
    
    //Rcout << "thenewpizzle d        ";
    
    // Set up vectors that will be put together into matrix map data frame
    IntegerVector tanewage3 (totallength);
    IntegerVector tanewage2 (totallength);
    
    //arma::vec actualage (totallength, fill::zeros);
    arma::vec index321 (totallength); // No death transitions
    arma::vec index321d (totallength); // Death transitions included
    arma::vec index21 (totallength);
    arma::vec aliveequal (totallength);
    arma::vec aliveequal_proxy (totallength);
    arma::vec included (totallength, fill::zeros);
    index321.fill(-1.0);
    index321d.fill(-1.0);
    index21.fill(-1.0);
    aliveequal.fill(-1.0);
    aliveequal_proxy.fill(-1.0);
    
    arma::mat asadditions (totallength, 5, fill::zeros);
    
    arma::vec tanewconvtype (totallength, fill::ones);
    arma::vec tagivent (totallength);
    arma::vec tagivenf (totallength);
    arma::vec taestt (totallength);
    arma::vec taestf (totallength);
    arma::vec taoffsett (totallength, fill::zeros);
    arma::vec taoffsetf (totallength, fill::zeros);
    arma::vec tasurvmult (totallength);
    arma::vec tafecmult (totallength);
    arma::uvec tampmaltered_final (totallength, fill::zeros);
    tagivent.fill(-1.0);
    taestt.fill(-1.0);
    tagivenf.fill(-1.0);
    taestf.fill(-1.0);
    tasurvmult.fill(-1.0);
    tafecmult.fill(-1.0);
    
    //Rcout << "thenewpizzle e        ";
    
    // Change stage names to stage numbers per input stageframe for styles 0 and 1
    if (style < 3) {
      if (tarows > 0) {
        if (tarows > 1 || taconvtype(0) != -1.0) {
          for (int i = 0; i < tarows; i++) { // Loop across trait_axis rows
            //Rcout << "Entered thenewpizzle interior loop a (prep)        ";
            for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
              
              if (tastage3(i) == origstageid(j)) {
                taindex3(i) = newstageid(j);
              }
              
              if (tastage2(i) == origstageid(j)) {
                taindex2(i) = newstageid(j);
              }
              
              if (tastage1(i) == origstageid(j)) {
                taindex1(i) = newstageid(j);
              }
              
              if (taeststage3(i) == origstageid(j)) {
                tanew3(i) = newstageid(j);
              }
              
              if (taeststage2(i) == origstageid(j)) {
                tanew2(i) = newstageid(j);
              }
              
              if (taeststage1(i) == origstageid(j)) {
                tanew1(i) = newstageid(j);
              }
            } // j for loop
            
            if (style == 0 && format == 2) {
              //Rcout << "Entered thenewpizzle interior loop a1 (historical deVries)        ";
              
              if (taconvtype(i) > 1.0) { // Catches all changes to fecundity and reproductive multipliers
                tanewconvtype(i) = 2;
                taindexold321(i) = (taindex3(i) - 1) + (prior_stage * nostages) + 
                  ((taindex2(i) - 1) * nostages * nostages) + 
                  ((taindex1(i) - 1) * nostages * nostages * nostages);
                  
                taindexnew321(i) = (tanew3(i) - 1) + (prior_stage * nostages) + 
                  ((tanew2(i) - 1) * nostages * nostages) + 
                  ((tanew1(i) - 1) * nostages * nostages * nostages);
              } else if (taconvt12(i) == 2.0) { // Catches all survival terms with historical reproduction events
                taindexold321(i) = (taindex3(i) - 1) + ((taindex2(i) - 1) * nostages) + 
                  ((taindex2(i) - 1) * nostages * nostages) + 
                  (prior_stage * nostages * nostages * nostages);
                  
                taindexnew321(i) = (tanew3(i) - 1) + ((tanew2(i) - 1) * nostages) + 
                  ((tanew2(i) - 1) * nostages * nostages) + 
                  (prior_stage * nostages * nostages * nostages);
              } else { // Full survival transitions
                taindexold321(i) = (taindex3(i) - 1) + ((taindex2(i) - 1) * nostages) + 
                  ((taindex2(i) - 1) * nostages * nostages) + 
                  ((taindex1(i) - 1) * nostages * nostages * nostages);
                  
                taindexnew321(i) = (tanew3(i) - 1) + ((tanew2(i) - 1) * nostages) + 
                  ((tanew2(i) - 1) * nostages * nostages) + 
                  ((tanew1(i) - 1) * nostages * nostages * nostages);
              }
              if (taindexold321(i) < 0.0) taindexold321(i) = -1.0;
              if (taindexnew321(i) < 0.0) taindexnew321(i) = -1.0;
              
              if (!NumericVector::is_na(tagivenrate(i))) {
                if (taconvtype(i) == 1) {
                  tagivent(i) = tagivenrate(i);
                } else {
                  tagivenf(i) = tagivenrate(i);
                }
              }
              if (!NumericVector::is_na(taoffset(i))) {
                if (taoffset(i) != 0.) {
                  if (taconvtype(i) == 1) {
                    taoffsett(i) = taoffset(i);
                  } else {
                    taoffsetf(i) = taoffset(i);
                  }
                }
              }
              if (NumericVector::is_na(tamultiplier(i))) {
                if (taconvtype(i) == 1) {
                  tasurvmult(i) = tamultiplier(i);
                } else {
                  tafecmult(i) = tamultiplier(i);
                }
              }
              
              tanewmpmaltered(i) = tampmaltered(i);
              
              // This is from a different loop cycle, and might cause problems
              index321d(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages) + 
                ((tanewstageid2(i) - 1) * nostages * nostages) + 
                ((tanewstageid1(i) - 1) * nostages * nostages * nostages);
                  
              // Next index variable gives element in the final matrix
              aliveequal(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn) + 
                ((tanewstageid2(i) - 1) * nostages * nostages_nodead) + 
                ((tanewstageid1(i) - 1) * nostages_nounborn * nostages_nodead * nostages);
              
              // The next is an index for proxy elements
              aliveequal_proxy(i) = (tanewstageid3est(i) - 1) + 
                ((tanewstageid2est(i) - 1) * nostages_nounborn) + 
                ((tanewstageid2est(i) - 1) * nostages * nostages_nodead) + 
                ((tanewstageid1est(i) - 1) * nostages_nounborn * nostages_nodead * nostages);
              
              // Next two index variables used by ovreplace
              index321(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages) + 
                ((tanewstageid2(i) - 1) * nostages * nostages) + 
                ((tanewstageid1(i) - 1) * nostages * nostages * nostages);
                
              index21(i) = (tanewstageid2(i) - 1) + 
                ((tanewstageid1(i) - 1) * nostages);
            } else if (style == 0 && format == 1) {
              //Rcout << "Entered thenewpizzle interior loop a2 (historical Ehrlen)        ";
              
              if (taconvtype(i) > 1.0) tanewconvtype(i) = 2;
              
              taindexold321(i) = (taindex3(i) - 1) + ((taindex2(i) - 1) * nostages_nounborn) + 
                ((taindex2(i) - 1) * nostages_nounborn * nostages_nounborn) + 
                ((taindex1(i) - 1) * nostages_nounborn * nostages_nounborn * 
                  nostages_nounborn);
                
              taindexnew321(i) = (tanew3(i) - 1) + ((tanew2(i) - 1) * nostages_nodead) + 
                ((tanew2(i) - 1) * nostages_nounborn * nostages_nounborn) + 
                ((tanew1(i) - 1) * nostages_nounborn * nostages_nounborn * 
                  nostages_nounborn);
              
              if (taindexold321(i) < 0) taindexold321(i) = -1.0;
              if (taindexnew321(i) < 0) taindexnew321(i) = -1.0;
              
              if (!NumericVector::is_na(tagivenrate(i))) {
                if (taconvtype(i) == 1) {
                  tagivent(i) = tagivenrate(i);
                } else {
                  tagivenf(i) = tagivenrate(i);
                }
              }
              if (!NumericVector::is_na(taoffset(i))) {
                if (taoffset(i) != 0.) {
                  if (taconvtype(i) == 1) {
                    taoffsett(i) = taoffset(i);
                  } else {
                    taoffsetf(i) = taoffset(i);
                  }
                }
              }
              if (NumericVector::is_na(tamultiplier(i))) {
                if (taconvtype(i) == 1) {
                  tasurvmult(i) = tamultiplier(i);
                } else {
                  tafecmult(i) = tamultiplier(i);
                }
              }
              
              tanewmpmaltered(i) = tampmaltered(i);
              
              //Rcout << "Entered thenewpizzle interior loop a3        ";
              
              // The below is from the next loop cycle, so may cause problems
              index321d(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn * nostages_nounborn) + 
                ((tanewstageid1(i) - 1) * nostages_nounborn * nostages_nounborn * 
                  nostages_nounborn);
              
              aliveequal(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn * nostages_nounborn) + 
                ((tanewstageid1(i) - 1) * nostages_nounborn * nostages_nounborn * 
                  nostages_nounborn);
              
              aliveequal_proxy(i) = (tanewstageid3est(i) - 1) + 
                ((tanewstageid2est(i) - 1) * nostages_nounborn) + 
                ((tanewstageid2est(i) - 1) * nostages_nounborn * nostages_nounborn) + 
                ((tanewstageid1est(i) - 1) * nostages_nounborn * nostages_nounborn * 
                  nostages_nounborn);
              
              index321(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn) + 
                ((tanewstageid2(i) - 1) * nostages_nounborn * nostages_nounborn) + 
                ((tanewstageid1(i) - 1) * nostages_nounborn * nostages_nounborn * 
                  nostages_nounborn);
              index21(i) = (tanewstageid2(i) - 1) + ((tanewstageid1(i) - 1) * nostages);
            } else if (style == 1) { // Ahistorical
              //Rcout << "Entered thenewpizzle interior loop a4 (ahistorical)        ";
              
              if (taconvtype(i) > 1.0) tanewconvtype(i) = 2;
              
              taindexold321(i) = (taindex3(i) - 1) + ((taindex2(i) - 1) * nostages);
              taindexnew321(i) = (tanew3(i) - 1) + ((tanew2(i) - 1) * nostages);
              
              if (taindexold321(i) < 0) taindexold321(i) = -1.0;
              if (taindexnew321(i) < 0) taindexnew321(i) = -1.0;
              
              if (!NumericVector::is_na(tagivenrate(i))) {
                if (taconvtype(i) == 1) {
                  tagivent(i) = tagivenrate(i);
                } else {
                  tagivenf(i) = tagivenrate(i);
                }
              }
              if (!NumericVector::is_na(taoffset(i))) {
                if (taoffset(i) != 0.) {
                  if (taconvtype(i) == 1) {
                    taoffsett(i) = taoffset(i);
                  } else {
                    taoffsetf(i) = taoffset(i);
                  }
                }
              }
              if (!NumericVector::is_na(tamultiplier(i))) {
                if (taconvtype(i) == 1) {
                  tasurvmult(i) = tamultiplier(i);
                } else {
                  tafecmult(i) = tamultiplier(i);
                }
              }
              
              tanewmpmaltered(i) = tampmaltered(i);
              
              //Rcout << "Entered thenewpizzle interior loop a5        ";
              // The following is from the next loop cycle, and may cause problems
              aliveequal(i) = (tastageorder3(i) - 1) + 
                ((tastageorder2(i) - 1) * nostages);
              
              aliveequal_proxy(i) = (tastageorder3est(i) - 1) + 
                ((tastageorder2est(i) - 1) * nostages);
              
              index321(i) = (tanewstageid3(i) - 1) + 
                ((tanewstageid2(i) - 1) * nostages);
              index21(i) = (tanewstageid2(i) - 1);
            } else if (style == 2) { // age-by-stage
              //Rcout << "Entered thenewpizzle interior loop a6 (age-by-stage)        ";
              
              bool found_mpm_altered {false};
              
              if (taconvtype(i) > 1.0) tanewconvtype(i) = 2;
              
              int age3 = taage3(i);
              int age2 = taage2(i);
              int estage3 = taestage3(i);
              int estage2 = taestage2(i);
              
              tanewage3(i) = taage3(i);
              tanewage2(i) = taage2(i);
              
              if (tampmaltered(i) > 0) { // Check here
                if (tanewage3(i) < firstage || tanewage3(i) > finalage ||
                    tanewage2(i) < firstage || tanewage2(i) > finalage) {
                  throw Rcpp::exception("Some ages in trait_axis are outside of the possible range.", false);
                }
              }
              
              taindexold321(i) = taindex3(i) + ((age3 - firstage) * nostages) +
                (taindex2(i) * nostages * totalages) + 
                ((age2 - firstage) * nostages * nostages * totalages);
              
              if (!IntegerVector::is_na(taestage2(i)) && taestage2(i) != -1) {
                found_mpm_altered = true;
                int newage2 = taestage2(i);
                int newage3 = newage2 + 1;
                    
                taindexnew321(i) = tanew3(i) + ((newage3 - firstage) * nostages) +
                  (tanew2(i) * nostages * totalages) +
                  ((newage2 - firstage) * nostages * nostages * totalages);
              } else {
                taindexnew321(i) = tanew3(i) + ((age3 - firstage) * nostages) +
                  (tanew2(i) * nostages * totalages) +
                  ((age2 - firstage) * nostages * nostages * totalages);
              }
              
              if (!NumericVector::is_na(tagivenrate(i))) {
                found_mpm_altered = true;
                if (taconvtype(i) == 1) {
                  tagivent(i) = tagivenrate(i);
                } else {
                  tagivenf(i) = tagivenrate(i);
                }
              }
              if (!NumericVector::is_na(taoffset(i))) {
                found_mpm_altered = true;
                if (taoffset(i) != 0.) {
                  if (taconvtype(i) == 1) {
                    taoffsett(i) = taoffset(i);
                  } else {
                    taoffsetf(i) = taoffset(i);
                  }
                }
              }
              if (!NumericVector::is_na(tamultiplier(i))) {
                found_mpm_altered = true;
                if (taconvtype(i) == 1) {
                  tasurvmult(i) = tamultiplier(i);
                } else {
                  tafecmult(i) = tamultiplier(i);
                }
              }
              
              if (IntegerVector::is_na(taage3(i)) || IntegerVector::is_na(taage2(i))) {
                if (found_mpm_altered) {
                  throw Rcpp::exception("Age-by-stage MPMs require explicit ages in trait variants.",
                    false);
                }
              }
              
              tanewmpmaltered(i) = tampmaltered(i);
              
              //Rcout << "Entered thenewpizzle interior loop a7        ";
              // The following is from the next loop cycle, and may cause problems
              int time3 = tastageorder3(i) - 1;
              int time2n = tastageorder2(i) - 1;
              
              int time3est = tastageorder3est(i) - 1;
              int time2nest = tastageorder2est(i) - 1;
              
              int currentindex = time3 + ((age3 - firstage) * nostages) + 
                (time2n * nostages * totalages) +
                ((age2 - firstage) * nostages * nostages * totalages);
              int currentindex_proxy = time3est + ((estage3 - firstage) * nostages) + 
                (time2nest * nostages * totalages) +
                ((estage2 - firstage) * nostages * nostages * totalages);
              
              // Indexer order: (1st # age blocks) + (1st # stage cols) +
              // (1st # age rows) + stage in time 3
              index321(i) = currentindex;
              index21(i) = time2n + ((age2 - firstage) * nostages);
              
              // Identify elements with non-zero entries by element number in final matrix
              aliveequal(i) = currentindex;
              aliveequal_proxy(i) = currentindex_proxy;
            }
          } // i for loop
        } // tarows if statement
      }
    } // style if statement
    
    //Rcout << "thenewpizzle f         ";
    
    // Output formatting
    Rcpp::List output_longlist(23);
    arma::uvec used_indices;
    
    if (filter == 1) {
      used_indices = find(index321 != -1.0);
    } else if (filter == 2) {
      used_indices = find(aliveequal != -1.0);
    }
    
    //Rcout << "thenewpizzle g         ";
    
    int aliveequal_length = static_cast<int>(aliveequal.n_elem);
    
    output_longlist(0) = Rcpp::IntegerVector(tanewstageid3.begin(), tanewstageid3.end());
    output_longlist(1) = Rcpp::IntegerVector(tanewstageid2.begin(), tanewstageid2.end());
    output_longlist(2) = Rcpp::IntegerVector(tanewstageid1.begin(), tanewstageid1.end());
    
    output_longlist(3) = Rcpp::IntegerVector(tanewstageid3est.begin(), tanewstageid3est.end());
    output_longlist(4) = Rcpp::IntegerVector(tanewstageid2est.begin(), tanewstageid2est.end());
    output_longlist(5) = Rcpp::IntegerVector(tanewstageid1est.begin(), tanewstageid1est.end());
    
    output_longlist(6) = tanewage3;
    output_longlist(7) = tanewage2;
    
    output_longlist(8) = Rcpp::NumericVector(tagivent.begin(), tagivent.end());
    output_longlist(9) = Rcpp::NumericVector(taestt.begin(), taestt.end());
    output_longlist(10) = Rcpp::NumericVector(tagivenf.begin(), tagivenf.end());
    output_longlist(11) = Rcpp::NumericVector(taestf.begin(), taestf.end());
    output_longlist(12) = Rcpp::NumericVector(tasurvmult.begin(), tasurvmult.end());
    output_longlist(13) = Rcpp::NumericVector(tafecmult.begin(), tafecmult.end());
    
    output_longlist(14) = Rcpp::NumericVector(aliveequal.begin(), aliveequal.end());
    output_longlist(15) = Rcpp::NumericVector(aliveequal_proxy.begin(), aliveequal_proxy.end());
    output_longlist(16) = Rcpp::NumericVector(index321.begin(), index321.end());
    output_longlist(17) = Rcpp::NumericVector(index321d.begin(), index321d.end());
    output_longlist(18) = Rcpp::NumericVector(index21.begin(), index21.end());
    
    output_longlist(19) = Rcpp::NumericVector(taoffsett.begin(), taoffsett.end());
    output_longlist(20) = Rcpp::NumericVector(taoffsetf.begin(), taoffsetf.end());
    output_longlist(21) = Rcpp::NumericVector(tanewconvtype.begin(), tanewconvtype.end());
    output_longlist(22) = Rcpp::NumericVector(tanewmpmaltered.begin(), tanewmpmaltered.end());
    
    //Rcout << "thenewpizzle h         ";
    
    CharacterVector namevec = {"stage3", "stage2", "stage1", "eststage3",
      "eststage2", "eststage1", "age3", "age2", "tagiven_t", "taest_t",
      "tagiven_f", "taest_f", "tasurvmult", "tafecmult", "aliveandequal",
      "aliveandequal_proxy", "index321", "index321d", "index21", "taoffset_t",
      "taoffset_f", "taconvtype", "mpm_altered"};
    output_longlist.attr("names") = namevec;
    output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, aliveequal_length);
    output_longlist.attr("class") = "data.frame";
    
    return output_longlist;
  }

}

#endif
