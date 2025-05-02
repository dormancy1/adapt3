#include <RcppArmadillo.h>
#define BOOST_DISABLE_ASSERTS
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(lefko3)]]

#include <LefkoUtils.h>
#include <AdaptUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;
using namespace AdaptUtils;





// Index of functions
// 1. DataFrame ta_skeleton  Create Skeleton Data Frame for Trait Variation for Invasion Analysis
// 2. List trait_axis  Create a Data Frame of Trait Data for Invasion Analysis



//' Create Skeleton Data Frame for Trait Variation for Invasion Analysis
//' 
//' Function \code{ta_skeleton()} creates a core data frame that can be modified
//' by users to provide the core variation in transition elements and vital
//' rates to use in invasion analysis. The resulting data frame should be used
//' as input in function \code{\link{invade3}()}.
//' 
//' @name ta_skeleton
//' 
//' @param rows The number of rows needed in the data frame. Defaults to 10.
//' 
//' @return A data frame of class \code{adaptAxis}, with the following columns:
//' \item{variant}{Denotes each variant in order, with each row corresponding to
//' a novel variant.}
//' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be
//' replaced.}
//' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
//' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be
//' replaced.}
//' \item{age3}{Age at occasion \emph{t}+1 in the transition to be replaced.}
//' \item{age2}{Age at occasion \emph{t} in the transition to be replaced.}
//' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
//' the transition designated by \code{stage3}, \code{stage2}, and 
//' \code{stage1}.}
//' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace the
//' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
//' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
//' the transition designated by \code{stage3}, \code{stage2}, and 
//' \code{stage1}.}
//' \item{estage3}{Age at occasion \emph{t}+1 in the transition to replace the
//' transition designated by \code{age2}.}
//' \item{estage2}{Age at occasion \emph{t} in the transition to replace the
//' transition designated by \code{age2}.}
//' \item{givenrate}{A constant to be used as the value of the transition.}
//' \item{offset}{A constant value to be added to the transition or proxy
//' transition.}
//' \item{multiplier}{A multiplier for proxy transitions or for fecundity.}
//' \item{convtype}{Designates whether the transition from occasion \emph{t} to
//' occasion \emph{t}+1 is a survival transition probability (1), a fecundity
//' rate (2), or a fecundity multiplier (3).}
//' \item{convtype_t12}{Designates whether the transition from occasion
//' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), a
//' fecundity rate (2).}
//' \item{surv_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for survival probability.}
//' \item{obs_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for observation probability.}
//' \item{size_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for primary size transition.}
//' \item{sizeb_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for secondary size transition.}
//' \item{sizec_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for tertiary size transition.}
//' \item{repst_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for reproduction probability.}
//' \item{fec_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for fecundity.}
//' \item{jsurv_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for juvenile survival probability.}
//' \item{jobs_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for juvenile observation probability.}
//' \item{jsize_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for juvenile primary size transition.}
//' \item{jsizeb_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for juvenile secondary size transition.}
//' \item{jsizec_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for juvenile tertiary size transition.}
//' \item{jrepst_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for juvenile reproduction probability.}
//' \item{jmatst_dev}{A numeric vector giving the deviations to the y-intercept of
//' the vital rate model for maturity status.}
//' 
//' @examples
//' 
//' current_traits <- ta_skeleton(4)
//' current_traits$stage3 <- c("Dorm", "Dorm", "Sdl1", NA)
//' current_traits$stage2 <- c("cut", "V0r", "rep", NA)
//' current_traits$convtype <- c(1, 1, 2, NA)
//' current_traits$offset <- c(0.1, 0.2, 0.3, NA)
//' current_traits$surv_dev <- c(NA, NA, NA, 0.1)
//' 
//' @export ta_skeleton
// [[Rcpp::export(ta_skeleton)]]
DataFrame ta_skeleton (int rows = 10) {
  IntegerVector row_id = seq(1, rows);
  CharacterVector stage3 (rows);
  CharacterVector stage2 (rows);
  CharacterVector stage1 (rows);
  IntegerVector age2 (rows, NA_INTEGER);
  NumericVector entered1 (rows, NA_REAL);
  
  for (int i = 0; i < rows; i++) {
    stage3(i) = NA_STRING;
    stage2(i) = NA_STRING;
    stage1(i) = NA_STRING;
  }
  
  DataFrame output = DataFrame::create(_["variant"] = row_id,
    _["stage3"] = stage3, _["stage2"] = stage2, _["stage1"] = stage1,
    _["age3"] = age2, _["age2"] = clone(age2), _["eststage3"] = clone(stage3),
    _["eststage2"] = clone(stage2), _["eststage1"] = clone(stage1),
    _["estage3"] = clone(age2),  _["estage2"] = clone(age2),
    _["givenrate"] = clone(entered1), _["offset"] = clone(entered1),
    _["multiplier"] = clone(entered1), _["convtype"] = clone(age2),
    _["convtype_t12"] = clone(age2), _["surv_dev"] = clone(entered1),
    _["obs_dev"] = clone(entered1), _["size_dev"] = clone(entered1),
    _["sizeb_dev"] = clone(entered1), _["sizec_dev"] = clone(entered1),
    _["repst_dev"] = clone(entered1), _["fec_dev"] = clone(entered1),
    _["jsurv_dev"] = clone(entered1), _["jobs_dev"] = clone(entered1),
    _["jsize_dev"] = clone(entered1), _["jsizeb_dev"] = clone(entered1),
    _["jsizec_dev"] = clone(entered1), _["jrepst_dev"] = clone(entered1),
    _["jmatst_dev"] = clone(entered1));
  
  CharacterVector newclasses = {"data.frame", "adaptAxis"};
  output.attr("class") = newclasses;
  
  return(output);
}

//' Create a Data Frame of Trait Data for Invasion Analysis
//' 
//' Function \code{trait_axis()} provides all necessary data for invasion
//' analysis. It lists the specific variations to MPMs for each variant run.
//' Variants can be given via overwritten matrix elements, proxy matrix
//' elements, additive offsets on matrix elements, matrix element multipliers,
//' and additive offsets to y-intercepts in vital rate models.
//' 
//' @name trait_axis
//' 
//' @param historical A single logical value indicating whether the MPMs
//' intended will be historical or ahistorical. Defaults to \code{TRUE}.
//' @param stagebased A single logical value indicating whether the MPM will be
//' stage-based or age-by-stage. Defaults to \code{TRUE}.
//' @param agebased A single logical value indicating whether the MPM will be
//' age-based or age-by-stage. Defaults to \code{FALSE}.
//' @param stageframe The stageframe used to produce the MPM. Required if
//' producing any stage-based or age-by-stage MPM. Must be omitted for purely
//' age-based MPMs.
//' @param stage3 String vector of stage names in occasion \emph{t}+1 in the
//' transition to be affected. Abbreviations for groups of stages are also
//' usable (see \code{Notes}). Required in all stage-based and age-by-stage
//' MPMs.
//' @param stage2 String vector of stage names in occasion \emph{t} in the
//' transition to be affected. Abbreviations for groups of stages are also
//' usable (see \code{Notes}). Required in all stage-based and age-by-stage
//' MPMs.
//' @param stage1 String vector of stage names in occasion \emph{t}-1 in the
//' transition to be affected. Only needed if a historical matrix is to be
//' produced. Abbreviations for groups of stages are also usable (see
//' \code{Notes}). Required for historical stage-based MPMs.
//' @param age3 An integer vector of the ages in occasion \emph{t}+1 to use in
//' transitions to be affected. Required for all age- and age-by-stage MPMs.
//' @param age2 An integer vector of the ages in occasion \emph{t} to use in
//' transitions to be affected. Required for all age- and age-by-stage MPMs.
//' @param eststage3 String vector of stage names to replace \code{stage3} in a
//' proxy transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage2 String vector of stage names to replace \code{stage2} in a
//' proxy transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage1 String vector of stage names to replace \code{stage1} in a
//' proxy historical transition. Only needed if a transition will be replaced by
//' another estimated transition, and the matrix to be estimated is historical
//' and stage-based. Stage \code{NotAlive} is also possible for raw hMPMs as a
//' means of handling the prior stage for individuals entering the population in
//' occasion \emph{t}.
//' @param estage3 Integer vector of age at time \emph{t}+1 to replace
//' \code{age3} in a proxy transition. Only needed if a transition will be
//' replaced by another estimated transition, and only in age-based and
//' age-by-stage MPMs.
//' @param estage2 Integer vector of age at time \emph{t} to replace \code{age2}
//' in a proxy transition. Only needed if a transition will be replaced by
//' another estimated transition, and only in age-based and age-by-stage MPMs.
//' @param givenrate A numeric vector of fixed rates or probabilities to replace
//' for the transition described by \code{stage3}, \code{stage2}, \code{stage1},
//' and/or \code{age2}.
//' @param offset A numeric vector of fixed numeric values to add to the
//' transitions described by \code{stage3}, \code{stage2}, \code{stage1}, and/or
//' \code{age2}.
//' @param multiplier A numeric vector of multipliers for the transition
//' described by \code{stage3}, \code{stage2}, \code{stage1}, and/or
//' \code{age2}, or for the proxy transitions described by \code{eststage3},
//' \code{eststage2}, \code{eststage1}, and/or \code{estage2}. Defaults to
//' \code{1}.
//' @param type Integer vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition;
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition; or \code{3}, \code{R}, or \code{r} for a fecundity set value /
//' general multiplier. If empty or not provided, then defaults to \code{1} for
//' survival transition.
//' @param type_t12 An optional integer vector denoting the kind of transition
//' between occasions \emph{t}-1 and \emph{t}. Only necessary if a historical
//' MPM in deVries format is desired. This should be entered as \code{1},
//' \code{S}, or \code{s} for a survival transition; or \code{2}, \code{F}, or
//' \code{f} for a fecundity transitions. Defaults to \code{1} for survival
//' transition, with impacts only on the construction of deVries-format hMPMs.
//' @param surv_dev An optional vector of numeric deviations to the y-intercept
//' of the survival model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param obs_dev An optional vector of numeric deviations to the y-intercept
//' of the observation model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param size_dev An optional vector of numeric deviations to the y-intercept
//' of the primary size model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param sizeb_dev An optional vector of numeric deviations to the y-intercept
//' of the secondary size model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param sizec_dev An optional vector of numeric deviations to the y-intercept
//' of the tertiary size model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param repst_dev An optional vector of numeric deviations to the y-intercept
//' of the reproduction model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param fec_dev An optional vector of numeric deviations to the y-intercept
//' of the fecundity model used in function-based MPM creation. Defaults to
//' \code{NA} for all values.
//' @param jsurv_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile survival model used in function-based MPM creation. Defaults
//' to \code{NA} for all values.
//' @param jobs_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile observation model used in function-based MPM creation.
//' Defaults to \code{NA} for all values.
//' @param jsize_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile primary size model used in function-based MPM creation.
//' Defaults to \code{NA} for all values.
//' @param jsizeb_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile secondary size model used in function-based MPM creation.
//' Defaults to \code{NA} for all values.
//' @param jsizec_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile tertiary size model used in function-based MPM creation.
//' Defaults to \code{NA} for all values.
//' @param jrepst_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile reproduction model used in function-based MPM creation.
//' Defaults to \code{NA} for all values.
//' @param jmatst_dev An optional vector of numeric deviations to the y-intercept
//' of the juvenile maturity model used in function-based MPM creation.
//' Defaults to \code{NA} for all values.
//' 
//' @return A data frame of class \code{adaptAxis}. This object can be used as
//' input in function \code{invade3()}.
//' 
//' Variables in this object include the following:
//' \item{variant}{Denotes each variant in order, with each row corresponding to
//' a novel variant.}
//' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be
//' replaced.}
//' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
//' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be
//' replaced.}
//' \item{age3}{Age at occasion \emph{t}+1 in the transition to be replaced.}
//' \item{age2}{Age at occasion \emph{t} in the transition to be replaced.}
//' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
//' the transition designated by \code{stage3}, \code{stage2}, and 
//' \code{stage1}.}
//' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace the
//' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
//' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
//' the transition designated by \code{stage3}, \code{stage2}, and 
//' \code{stage1}.}
//' \item{estage3}{Age at occasion \emph{t}+1 in the transition to replace the
//' transition designated by \code{age3}.}
//' \item{estage2}{Age at occasion \emph{t} in the transition to replace the
//' transition designated by \code{age2}.}
//' \item{givenrate}{A constant to be used as the value of the transition.}
//' \item{offset}{A constant value to be added to the transition or proxy
//' transition.}
//' \item{multiplier}{A multiplier for proxy transitions or for fecundity.}
//' \item{convtype}{Designates whether the transition from occasion \emph{t} to
//' occasion \emph{t}+1 is a survival transition probability (1), a fecundity
//' rate (2), or a fecundity multiplier (3).}
//' \item{convtype_t12}{Designates whether the transition from occasion
//' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), a
//' fecundity rate (2).}
//' \item{surv_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of survival.}
//' \item{obs_dev}{Numeric deviations to the y-intercept of the vital rate model
//' of observation.}
//' \item{size_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of primary size.}
//' \item{sizeb_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of secondary size.}
//' \item{sizec_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of tertiary size.}
//' \item{repst_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of reproduction.}
//' \item{fec_dev}{Numeric deviations to the y-intercept of the vital rate model
//' of fecundity.}
//' \item{jsurv_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile survival.}
//' \item{jobs_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile observation.}
//' \item{jsize_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile primary size.}
//' \item{jsizeb_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile secondary size.}
//' \item{jsizec_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile tertiary size.}
//' \item{jrepst_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile reproduction.}
//' \item{jmatst_dev}{Numeric deviations to the y-intercept of the vital rate
//' model of juvenile maturity.}
//' 
//' @section Notes:
//' Negative values are not allowed in \code{givenrate} and \code{multiplier}
//' input, but are allowed in \code{offset}, if values are to be subtracted from
//' specific estimated transitions. Stage entries should not be used for purely
//' age-based MPMs, and age entries should not be used for purely stage-based
//' MPMs.
//' 
//' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
//' abbreviations for groups of stages. Use \code{rep} if all reproductive
//' stages are to be used, \code{nrep} if all mature but non-reproductive stages
//' are to be used, \code{mat} if all mature stages are to be used, \code{immat}
//' if all immature stages are to be used, \code{prop} if all propagule stages
//' are to be used, \code{npr} if all non-propagule stages are to be used,
//' \code{obs} if all observable stages are to be used, \code{nobs} if all
//' unobservable stages are to be used, and leave empty or use \code{all} if all
//' stages in stageframe are to be used. Also use \code{groupX} to denote all
//' stages in group X (e.g. \code{group1} will use all stages in the respective
//' stageframe's group 1).
//' 
//' Type 3 conversions are referred to as fecundity set values, or general
//' fecundity multipliers. These set the transitions to be used as fecundity
//' transitions. Transitions set here will be interpreted as being generally
//' reproductive, meaning that the from and to stages will be used to determine
//' the general fecundity transitions to incorporate into stage-based MPMs,
//' while the age portion of the input will be used to incorporate the actual
//' multiplier(s) specified. If only stage transitions at certain ages are
//' expected to be the sole contributors to fecundity, then type 2 conversions
//' should also be included in the supplement (Type 1 and 2 conversions can be
//' purely age-specific, and do not set reproductive transitions in MPM
//' creation). For example, if all stage 2 to stage 3 transitions above age 2
//' yield fecundity, then stage 2 to stage 3 can be set to
//' \code{multiplier = 1.0} with \code{convtype = 3}, and the same transition
//' for \code{age2 = c(1, 2)} can be set to \code{multiplier = c(0, 0)}.
//' 
//' Several operations may be included per transition. Operations on the same
//' row of the resulting data frame are generally handled with given rate
//' substitutions first, then with proxy transitions, then by additive offsets,
//' and finally by multipliers. This order can be manipulated by ordering
//' operations across rows, with higher numbered rows in the data frame being
//' performed later.
//' 
//' @seealso \code{\link{ta_skeleton}()}
//' 
//' @examples
//' data(cypa_data)
//' 
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypa_data, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypa_ta <- trait_axis(stageframe = cypframe_raw,
//'   stage3 = c("P1", "P1", "P1", NA, NA, NA),
//'   stage2 = c("rep", "rep", "rep", NA, NA, NA),
//'   multiplier = c(0.5, 2.0, 10., NA, NA, NA), type = c(2, 2, 2, NA, NA, NA),
//'   obs_dev = c(NA, NA, NA, 0.5, 2.0, 50), fec_dev = c(NA, NA, NA, -1000, 0, 1000))
//' 
//' @export trait_axis
// [[Rcpp::export(trait_axis)]]
Rcpp::List trait_axis (Nullable<RObject> historical = R_NilValue,
  Nullable<RObject> stagebased = R_NilValue, Nullable<RObject> agebased = R_NilValue,
  Nullable<RObject> stageframe = R_NilValue,
  Nullable<RObject> stage3 = R_NilValue, Nullable<RObject> stage2 = R_NilValue,
  Nullable<RObject> stage1 = R_NilValue, Nullable<RObject> age3 = R_NilValue,
  Nullable<RObject> age2 = R_NilValue, Nullable<RObject> eststage3 = R_NilValue,
  Nullable<RObject> eststage2 = R_NilValue, Nullable<RObject> eststage1 = R_NilValue,
  Nullable<RObject> estage3 = R_NilValue, Nullable<RObject> estage2 = R_NilValue,
  Nullable<RObject> givenrate = R_NilValue, Nullable<RObject> offset = R_NilValue,
  Nullable<RObject> multiplier = R_NilValue, Nullable<RObject> type = R_NilValue,
  Nullable<RObject> type_t12 = R_NilValue, Nullable<RObject> surv_dev = R_NilValue,
  Nullable<RObject> obs_dev = R_NilValue, Nullable<RObject> size_dev = R_NilValue,
  Nullable<RObject> sizeb_dev = R_NilValue, Nullable<RObject> sizec_dev = R_NilValue,
  Nullable<RObject> repst_dev = R_NilValue, Nullable<RObject> fec_dev = R_NilValue,
  Nullable<RObject> jsurv_dev = R_NilValue, Nullable<RObject> jobs_dev = R_NilValue,
  Nullable<RObject> jsize_dev = R_NilValue, Nullable<RObject> jsizeb_dev = R_NilValue,
  Nullable<RObject> jsizec_dev = R_NilValue, Nullable<RObject> jrepst_dev = R_NilValue,
  Nullable<RObject> jmatst_dev = R_NilValue) {
  
  int wtf {-1};
  
  bool historical_bool {false};
  bool stagebased_bool {true};
  bool agebased_bool {false};
  
  bool trash_out {false};
  LefkoInputs::RObj_TF_input_check ("historical", "", historical_bool, trash_out,
    false, true, historical);
  LefkoInputs::RObj_TF_input_check ("stagebased", "", stagebased_bool, trash_out,
    false, true, stagebased);
  LefkoInputs::RObj_TF_input_check ("agebased", "", agebased_bool, trash_out,
    false, true, agebased);
  
  if (historical_bool && stagebased_bool && !agebased_bool) {
    wtf = 0;
  } else if (!historical_bool && stagebased_bool && !agebased_bool) {
    wtf = 1;
  } else if (!historical_bool && stagebased_bool && agebased_bool) {
    wtf = 2;
  } else if (!historical_bool && !stagebased_bool && agebased_bool) {
    wtf = 3;
  } else {
    throw Rcpp::exception("Unsupported MPM type.", false);
  }
  
  DataFrame stageframe_;
  int sf_yes {0};
  
  if (wtf < 3) {
    if (!stageframe.isNotNull()) {
      throw Rcpp::exception("Stageframe required for stage-based MPMs.", false);
    }
    
    if (is<DataFrame>(stageframe)) stageframe_ = as<DataFrame>(stageframe);
    StringVector sf_class = stageframe_.attr("class");
    
    String sf_error = "Please enter an object of class stageframe as input.";
    if (stageframe_.containsElementNamed("stage")) sf_yes++;
    if (stageframe_.containsElementNamed("min_age")) sf_yes++;
    if (stageframe_.containsElementNamed("max_age")) sf_yes++;
    if (stageframe_.containsElementNamed("group")) sf_yes++;
    
    for (int i = 0; i < static_cast<int>(sf_class.length()); i++) {
      if (sf_class(i) == "stageframe")  sf_yes++;
    }
    if (sf_yes < 5) throw Rcpp::exception(sf_error.get_cstring(), false);
  }
  
  StringVector stage3_;
  StringVector stage2_;
  StringVector stage1_;
  IntegerVector age3_;
  IntegerVector age2_;
  StringVector eststage3_;
  StringVector eststage2_;
  StringVector eststage1_;
  IntegerVector estage3_;
  IntegerVector estage2_;
  NumericVector givenrate_;
  NumericVector offset_;
  NumericVector multiplier_;
  IntegerVector type_;
  IntegerVector type_t12_;
  
  NumericVector surv_dev_;
  NumericVector obs_dev_;
  NumericVector size_dev_;
  NumericVector sizeb_dev_;
  NumericVector sizec_dev_;
  NumericVector repst_dev_;
  NumericVector fec_dev_;
  NumericVector jsurv_dev_;
  NumericVector jobs_dev_;
  NumericVector jsize_dev_;
  NumericVector jsizeb_dev_;
  NumericVector jsizec_dev_;
  NumericVector jrepst_dev_;
  NumericVector jmatst_dev_;
  
  int stage3_length {0};
  int stage2_length {0};
  int stage1_length {0};
  int age3_length {0};
  int age2_length {0};
  int eststage3_length {0};
  int eststage2_length {0};
  int eststage1_length {0};
  int estage3_length {0};
  int estage2_length {0};
  int type_length {0};
  
  StringVector all_stages;
  StringVector wildcard_names = {"all", "rep", "nrep", "mat", "immat", "prop",
    "npr", "notalive", "obs", "nobs"};;
  StringVector all_groups;
  
  if (wtf < 3) {
    all_stages = as<StringVector>(stageframe_["stage"]);
    
    StringVector all_groups_ = as<StringVector>(stageframe_["group"]);
    StringVector all_groups_temp (all_groups_.length());
    for (int i = 0; i < all_groups_.length(); i++) {
      String group_slot = "group";
      String group_term_added = all_groups_(i);
      group_slot += group_term_added;
      all_groups_temp(i) = group_slot;
    }
    all_groups = all_groups_temp;
  }
  
  if (stage3.isNotNull() && wtf != 3) {
    if (is<StringVector>(stage3)) {
      stage3_ = as<StringVector>(stage3);
      stage3_length = static_cast<int>(stage3_.length());
      
      for (int i = 0; i < stage3_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage3_(i)) || stage3_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage3_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage3_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage3_(i) == all_groups(j)) found_stage = true;
          }
          if (stage3_(i) == "AlmostBorn") {
            throw Rcpp::exception("Stage AlmostBorn can only be used in arguments stage1 and eststage1.", false);
          }
          if (!found_stage) {
            String ems = "Some entered stage names (stage3) do not ";
            String ems1 = "match expectation.";
            ems += ems1;
            
            throw Rcpp::exception(ems.get_cstring(), false);
          }
        }
      }
      
    } else {
      throw Rcpp::exception("Please enter argument stage3 as text.", 
        false);
    }
  } else if (wtf != 3) {
    throw Rcpp::exception("Stage information (stage3) for transitions is required.",
      false);
  }
  
  if (stage2.isNotNull() && wtf != 3) {
    if (is<StringVector>(stage2)) {
      stage2_ = as<StringVector>(stage2);
      stage2_length = static_cast<int>(stage2_.length());
      
      for (int i = 0; i < stage2_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage2_(i)) || stage2_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage2_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage2_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage2_(i) == all_groups(j)) found_stage = true;
          }
          if (stage2_(i) == "AlmostBorn") {
            throw Rcpp::exception("Stage AlmostBorn can only be used in arguments stage1 and eststage1.",
              false);
          }
          if (!found_stage) {
            String ems = "Some entered stage names (stage2) do not ";
            String ems1 = "match expectation.";
            ems += ems1;
            
            throw Rcpp::exception(ems.get_cstring(), false);
          }
        }
      }
      
    } else {
      throw Rcpp::exception("Please enter argument stage2 as text.", false);
    }
  } else if (wtf != 3) {
    throw Rcpp::exception("Stage information (stage2) for transitions is required.",
      false);
  }
  
  if (stage1.isNotNull() && wtf == 0) {
    if (is<StringVector>(stage1)) {
      stage1_ = as<StringVector>(stage1);
      stage1_length = static_cast<int>(stage1_.length());
      
      for (int i = 0; i < stage1_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage1_(i)) || stage1_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage1_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage1_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage1_(i) == all_groups(j)) found_stage = true;
          }
          if (stage1_(i) == "AlmostBorn") found_stage = true;
          if (!found_stage) {
            String ems = "Some entered stage names (stage1) do not ";
            String ems1 = "match expectation.";
            ems += ems1;
            
            throw Rcpp::exception(ems.get_cstring(), false);
          }
        }
      }
      
    } else if (is<LogicalVector>(stage1)) {
      if (wtf == 0) {
        throw Rcpp::exception("Stage information (stage1) for transitions is required.",
          false);
      }
      
      StringVector stage1_temp (stage2_length, NA_STRING);
      stage1_ = stage1_temp;
      stage1_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter argument stage1 as text.", 
        false);
    }
  } else {
    if (wtf == 0) {
      throw Rcpp::exception("Stage information (stage1) for transitions is required.",
        false);
    } else if (wtf != 3) {
      StringVector stage1_temp (stage2_length, NA_STRING);
      stage1_ = stage1_temp;
      stage1_length = stage2_length;
    }
  }
  
  if (age3.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(age3) || is<NumericVector>(age3)) {
      age3_ = as<IntegerVector>(age3);
      age3_length = static_cast<int>(age3_.length());
      
      arma::ivec a3_arma = as<arma::ivec>(age3_);
      arma::uvec neg_tester = find(a3_arma < 0 && a3_arma > -2147483648);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered age3 values are negative.");
      }
    } else if (is<LogicalVector>(age3)) { 
      IntegerVector age3_temp (stage2_length, NA_INTEGER);
      age3_ = age3_temp;
      age3_length = stage2_length;
    } else {
      throw Rcpp::exception("Please enter argument age3 in integer format.",
        false);
    }
  } else {
    IntegerVector age3_temp (stage2_length, NA_INTEGER);
    age3_ = age3_temp;
    age3_length = stage3_length;
  }
  
  if (age2.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(age2) || is<NumericVector>(age2)) {
      age2_ = as<IntegerVector>(age2);
      age2_length = static_cast<int>(age2_.length());
      
      arma::ivec a2_arma = as<arma::ivec>(age2_);
      arma::uvec neg_tester = find(a2_arma < 0 && a2_arma > -2147483648);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered age2 values are negative.");
      }
    } else if (is<LogicalVector>(age2)) { 
      IntegerVector age2_temp (stage2_length, NA_INTEGER);
      age2_ = age2_temp;
      age2_length = stage2_length;
    } else {
      throw Rcpp::exception("Please enter argument age2 in integer format.",
        false);
    }
  } else {
    IntegerVector age2_temp (stage2_length, NA_INTEGER);
    age2_ = age2_temp;
    age2_length = stage2_length;
  }
  
  if (eststage3.isNotNull() && wtf != 3) {
    if (is<StringVector>(eststage3)) {
      eststage3_ = as<StringVector>(eststage3);
      eststage3_length = static_cast<int>(eststage3_.length());
      
      for (int i = 0; i < eststage3_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage3_(i)) || eststage3_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage3_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage3_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage3_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage3_(i) == "AlmostBorn") {
            throw Rcpp::exception("Stage AlmostBorn can only be used in arguments stage1 and eststage1.", false);
          }
          if (!found_stage) {
            String ems = "Some entered stage names (eststage3) do not ";
            String ems1 = "match expectation.";
            ems += ems1;
            
            throw Rcpp::exception(ems.get_cstring(), false);
          }
        } else found_stage = true;
        
      }
      
    } else if (is<LogicalVector>(eststage3)) {
      StringVector eststage3_temp (stage3_length, NA_STRING);
      eststage3_ = eststage3_temp;
      eststage3_length = stage3_length;
      
    } else {
      throw Rcpp::exception("Please enter argument eststage3 as text.", false);
    }
  } else if (stage3_length != 0) {
    StringVector eststage3_temp (stage3_length, NA_STRING);
    eststage3_ = eststage3_temp;
    eststage3_length = stage3_length;
  }
  
  if (eststage2.isNotNull() && wtf != 3) {
    if (is<StringVector>(eststage2)) {
      eststage2_ = as<StringVector>(eststage2);
      eststage2_length = static_cast<int>(eststage2_.length());
      
      for (int i = 0; i < eststage2_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage2_(i)) || eststage2_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage2_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage2_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage2_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage2_(i) == "AlmostBorn") {
            throw Rcpp::exception("Stage AlmostBorn can only be used in arguments stage1 and eststage1.", false);
          }
          if (!found_stage) {
            String ems = "Some entered stage names (eststage2) do not ";
            String ems1 = "match expectation.";
            ems += ems1;
            
            throw Rcpp::exception(ems.get_cstring(), false);
          }
        } else found_stage = true;
        
      }
      
    } else if (is<LogicalVector>(eststage2)) {
      StringVector eststage2_temp (stage2_length, NA_STRING);
      eststage2_ = eststage2_temp;
      eststage2_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter argument eststage2 as text.", false);
    }
  } else if (stage2_length != 0) {
    StringVector eststage2_temp (stage2_length, NA_STRING);
    eststage2_ = eststage2_temp;
    eststage2_length = stage2_length;
  }
  
  if (eststage1.isNotNull() && wtf == 0) {
    if (is<StringVector>(eststage1)) {
      eststage1_ = as<StringVector>(eststage1);
      eststage1_length = static_cast<int>(eststage1_.length());
      
      for (int i = 0; i < eststage1_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage1_(i)) || eststage1_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage1_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage1_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage1_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage1_(i) == "AlmostBorn") found_stage = true;
          if (!found_stage) {
            String ems = "Some entered stage names (eststage1) do not ";
            String ems1 = "match expectation.";
            ems += ems1;
            
            throw Rcpp::exception(ems.get_cstring(), false);
          }
        } else found_stage = true;
        
      }
      
    } else if (is<LogicalVector>(eststage1)) {
      StringVector eststage1_temp (stage2_length, NA_STRING);
      eststage1_ = eststage1_temp;
      eststage1_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter argument eststage1 as text.", false);
    }
  } else if (stage2_length != 0) {
    StringVector eststage1_temp (stage2_length, NA_STRING);
    eststage1_ = eststage1_temp;
    eststage1_length = stage2_length;
  }
  
  if (estage3.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(estage3) || is<NumericVector>(estage3)) {
      estage3_ = as<IntegerVector>(estage3);
      estage3_length = static_cast<int>(estage3_.length());
      
      arma::ivec ea3_arma = as<arma::ivec>(estage3_);
      arma::uvec neg_tester = find(ea3_arma < 0 && ea3_arma > -2147483648);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered estage3 values are negative.");
      }
    } else if (is<LogicalVector>(estage3)) {
      LogicalVector estage3_temp = as<LogicalVector>(estage3);
      estage3_length = static_cast<int>(estage3_temp.length());
      
      IntegerVector estage3_new (estage3_length, NA_INTEGER);
      estage3_ = estage3_new;
      
    } else {
      throw Rcpp::exception("Please enter argument estage3 in integer format.",
        false);
    }
  } else {
    if (stage2_length > 0) {
      IntegerVector estage3_temp (stage2_length, NA_INTEGER);
      estage3_ = estage3_temp;
      estage3_length = stage2_length;
    } else if (age2_length > 0) {
      IntegerVector estage3_temp (age2_length, NA_INTEGER);
      estage3_ = estage3_temp;
      estage3_length = age2_length;
    }
  }
  
  if (estage2.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(estage2) || is<NumericVector>(estage2)) {
      estage2_ = as<IntegerVector>(estage2);
      estage2_length = static_cast<int>(estage2_.length());
      
      arma::ivec ea2_arma = as<arma::ivec>(estage2_);
      arma::uvec neg_tester = find(ea2_arma < 0 && ea2_arma > -2147483648);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered estage2 values are negative.");
      }
    } else if (is<LogicalVector>(estage2)) {
      LogicalVector estage2_temp = as<LogicalVector>(estage2);
      estage2_length = static_cast<int>(estage2_temp.length());
      
      IntegerVector estage2_new (estage2_length, NA_INTEGER);
      estage2_ = estage2_new;
      
    } else {
      throw Rcpp::exception("Please enter argument estage2 in integer format.",
        false);
    }
  } else {
    if (stage2_length > 0) {
      IntegerVector estage2_temp (stage2_length, NA_INTEGER);
      estage2_ = estage2_temp;
      estage2_length = stage2_length;
    } else if (age2_length > 0) {
      IntegerVector estage2_temp (age2_length, NA_INTEGER);
      estage2_ = estage2_temp;
      estage2_length = age2_length;
    }
  }
  
  LefkoInputs::numeric_vectorizer(givenrate_, givenrate, "givenrate", stage2_length, age2_length, false, 0);
  LefkoInputs::numeric_vectorizer(offset_, offset, "offset", stage2_length, age2_length, false, 0);
  LefkoInputs::numeric_vectorizer(multiplier_, multiplier, "multiplier", stage2_length, age2_length, false, 0);
  
  LefkoInputs::integer_vectorizer (type_, type, "type", stage2_length, age2_length, 1, 3,
    true, false, 0);
  LefkoInputs::integer_vectorizer (type_t12_, type_t12, "type_t12", stage2_length,
    age2_length, 1, 2, true, false, 0);
  
  if (wtf < 3 && type_length != 0  && type_length != stage2_length) { 
    throw Rcpp::exception("All input vectors must be of the same length.", false);
  }
  
  if (wtf == 0) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    
    if (stage2_length != stage3_length || stage2_length != stage1_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage1_length != 0 && stage1_length != eststage1_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage2_length != 0 && stage2_length != eststage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage3_length != 0 && stage3_length != eststage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
  } else if (wtf == 1) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    
    if (stage2_length != stage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage2_length != 0 && stage2_length != eststage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage3_length != 0 && stage3_length != eststage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
  } else if (wtf == 2) {
    if (stage2_length == 0) throw Rcpp::exception("Argument stage2 is required.", false);
    if (age3_length == 0) throw Rcpp::exception("Argument age3 is required.", false);
    if (age2_length == 0) throw Rcpp::exception("Argument age2 is required.", false);
    
    if (stage2_length != stage3_length || stage2_length != age2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (estage3_length != 0 && age2_length != estage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    if (estage2_length != 0 && age2_length != estage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  
  } else if (wtf == 3) {
    if (age2_length == 0) throw Rcpp::exception("Argument age2 is required.", false);
    if (age3_length == 0) throw Rcpp::exception("Argument age2 is required.", false);
    
    if (estage3_length != 0 && age2_length != estage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    if (estage2_length != 0 && age2_length != estage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  }
    
  if (stage3_length == 0) { 
    StringVector s3_temp (age2_length, NA_STRING);
    stage3_ = s3_temp;
    stage3_length = age2_length;
  }
  if (stage2_length == 0) { 
    StringVector s2_temp (age2_length, NA_STRING);
    stage2_ = s2_temp;
    stage2_length = age2_length;
  }
  if (stage1_length == 0) { 
    StringVector s1_temp (stage2_length, NA_STRING);
    stage1_ = s1_temp;
    stage1_length = stage2_length;
  }
  if (eststage3_length == 0) { 
    StringVector s3_temp (stage2_length, NA_STRING);
    eststage3_ = s3_temp;
    eststage3_length = stage2_length;
  }
  if (eststage2_length == 0) { 
    StringVector s2_temp (stage2_length, NA_STRING);
    eststage2_ = s2_temp;
    eststage2_length = stage2_length;
  }
  if (eststage1_length == 0) { 
    StringVector s1_temp (stage2_length, NA_STRING);
    eststage1_ = s1_temp;
    eststage1_length = stage2_length;
  }
  
  LefkoInputs::numeric_vectorizer(surv_dev_, surv_dev, "surv_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(obs_dev_, obs_dev, "obs_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(size_dev_, size_dev, "size_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(sizeb_dev_, sizeb_dev, "sizeb_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(sizec_dev_, sizec_dev, "sizec_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(repst_dev_, repst_dev, "repst_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(fec_dev_, fec_dev, "fec_dev", stage2_length, age2_length, false);
  
  LefkoInputs::numeric_vectorizer(jsurv_dev_, jsurv_dev, "jsurv_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(jobs_dev_, jobs_dev, "jobs_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(jsize_dev_, jsize_dev, "jsize_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(jsizeb_dev_, jsizeb_dev, "jsizeb_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(jsizec_dev_, jsizec_dev, "jsizec_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(jrepst_dev_, jrepst_dev, "jrepst_dev", stage2_length, age2_length, false);
  LefkoInputs::numeric_vectorizer(jmatst_dev_, jmatst_dev, "jmatst_dev", stage2_length, age2_length, false);
  
  IntegerVector new_variants = seq(1, stage2_length);
  
  List newtraitaxis (30);
  
  newtraitaxis(0) = new_variants;
  newtraitaxis(1) = stage3_;
  newtraitaxis(2) = stage2_;
  newtraitaxis(3) = stage1_;
  newtraitaxis(4) = age3_;
  newtraitaxis(5) = age2_;
  newtraitaxis(6) = eststage3_;
  newtraitaxis(7) = eststage2_;
  newtraitaxis(8) = eststage1_;
  newtraitaxis(9) = estage3_;
  newtraitaxis(10) = estage2_;
  newtraitaxis(11) = givenrate_;
  newtraitaxis(12) = offset_;
  newtraitaxis(13) = multiplier_;
  newtraitaxis(14) = type_;
  newtraitaxis(15) = type_t12_;
  newtraitaxis(16) = surv_dev_;
  newtraitaxis(17) = obs_dev_;
  newtraitaxis(18) = size_dev_;
  newtraitaxis(19) = sizeb_dev_;
  newtraitaxis(20) = sizec_dev_;
  newtraitaxis(21) = repst_dev_;
  newtraitaxis(22) = fec_dev_;
  newtraitaxis(23) = jsurv_dev_;
  newtraitaxis(24) = jobs_dev_;
  newtraitaxis(25) = jsize_dev_;
  newtraitaxis(26) = jsizeb_dev_;
  newtraitaxis(27) = jsizec_dev_;
  newtraitaxis(28) = jrepst_dev_;
  newtraitaxis(29) = jmatst_dev_;
  
  StringVector ta_names = {"variant", "stage3", "stage2", "stage1", "age3",
    "age2", "eststage3", "eststage2", "eststage1", "estage3", "estage2",
    "givenrate", "offset", "multiplier", "convtype", "convtype_t12", "surv_dev",
    "obs_dev", "size_dev", "sizeb_dev", "sizec_dev", "repst_dev", "fec_dev",
    "jsurv_dev", "jobs_dev", "jsize_dev", "jsizeb_dev", "jsizec_dev",
    "jrepst_dev", "jmatst_dev"};
  StringVector ta_class = {"data.frame", "adaptAxis"};
  
  newtraitaxis.attr("class") = ta_class;
  newtraitaxis.attr("names") = ta_names;
  newtraitaxis.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stage2_length);
  
  return newtraitaxis;
}

