#' @title Adaptive Dynamics and Community Matrix Model Projections
#' 
#' @description Runs projections of groups of matrix projection models (MPMs),
#' allowing density dependence mechanisms to work across MPMs. This package was
#' developed to run both adaptive dynamics simulations such as pairwise and
#' multiple invasion analysis, and community projections in which species are
#' represented by MPMs. All forms of MPMs are allowed, including integral
#' projection models (IPMs).
#' 
#' @details The lefko3 package provides three categories of functions:
#' 
#' 1. Core projection
#' 
#' 2. Function characterizing relationships among MPMs
#' 
#' 3. Functions describing, summarizing, or visualizing results
#' 
#' @details adapt3 also includes example datasets complete with sample code.
#' 
#' @aliases adapt3-package
"_PACKAGE"
#' @author Richard P. Shefferson <cdorm@g.ecc.u-tokyo.ac.jp>
#' @references Pending
#' @import Rcpp graphics
#' @importFrom graphics contour .filled.contour
#' @importFrom methods as is
#' @importFrom Rcpp evalCpp
#' @useDynLib adapt3
#' @name adapt3-package
NULL
