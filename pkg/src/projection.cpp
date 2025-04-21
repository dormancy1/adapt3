#include <RcppArmadillo.h>
#define BOOST_DISABLE_ASSERTS
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(lefko3)]]

#include <RcppArmadilloExtensions/sample.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <LefkoUtils.h>
#include <AdaptUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;
using namespace AdaptUtils;





// Index of functions
// 1. List cleanup3  Clean Up RObject Inputs to Projection Functions
// 2. void project3_pre_core  Engine Projecting Multiple Existing MPMs With or Without Density Dependence
// 3. void project3_fb_core  Engine Projecting Multiple Function-based MPMs With or Without Density Dependence
// 4. List project3  Project Multiple MPMs With or Without Density Dependence
// 5. List cleanup3_inv  Clean Up RObject Inputs for Invasion Analysis
// 6. void invade3_pre_core  Engine Running Invasion Analyses of Existing MPMs
// 7. void invade3_fb_core  Engine Running Invasion Analyses of Function-based MPMs
// 8. List invade3  Run Pairwise and Multiple Invasion Analysis



//' Clean Up RObject Inputs to Projection Functions
//' 
//' This function takes RObject inputs in the core projection functions, and
//' uses them to create the strict inputs for projection.
//' 
//' @name cleanup3
//' 
//' @param mpms An optional list of MPMs. Each MPM must be of class
//' \code{lefkoMat}.
//' @param vrms An optional list of \code{vrm_input} objects, each corresponding
//' to a distinct MPM that will be created during projection. Each
//' \code{vrm_input} object requires its own stageframe, entered in the same
//' order via argument \code{stageframes}.
//' @param stageframes An optional list of stageframes, corresponding in number
//' and order to the MPMs in argument \code{vrms}. Each stageframe must be of
//' class \code{stageframe}.
//' @param supplements An optional list of data frames of class \code{lefkoSD}
//' that provide supplemental data that should be incorporated into
//' function-based MPMs. If used, then should be the same number of data frames
//' as the number of MPMs provided in the list for argument \code{vrms}. MPMs
//' that do not need supplemental data should be entered as \code{NULL} in this
//' list. See \code{\link[lefko3]{supplemental}()} for details.
//' @param format An optional integer vector indicating the kind of
//' function-based MPM to create for each \code{vrm_input} object entered in
//' argument \code{vrms}. Possible choices include: \code{1}, Ehrlen-format
//' historical MPM; \code{2}, deVries-format historical MPM; \code{3},
//' ahistorical MPM (default); \code{4}, age-by-stage MPM; and \code{5}, Leslie
//' (age-based) MPM.
//' @param firstage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the starting ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Typically,
//' the starting age in such MPMs should be set to \code{0} if post-breeding and
//' \code{1} if pre-breeding. All other MPMs should be set to \code{0}. Do not
//' use if no MPM has age structure. 
//' @param finalage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the final ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Do not use if no
//' MPM has age structure.
//' @param fecage_min An optional integer vector used for function-based Leslie
//' MPMs giving the first age at which organisms can be reproductive in such
//' MPMs. Use only if at least one MPM is a function-based Leslie MPM. Defaults
//' to the values given in \code{firstage}.
//' @param fecage_max An optional integer vector used for function-based Leslie
//' MPMs giving the final age at which organisms can be reproductive in such
//' MPMs. Use only if at least one MPM is a function-based Leslie MPM. Defaults
//' to the values given in \code{finalage}.
//' @param cont An optional vector used for function-based Leslie and
//' age-by-stage MPMs stating whether the MPM should should include a stasis
//' transition within the final age. This should be used only when an organism
//' can maintain the demographic characteristics of the final described age
//' after reaching that age. Can be entered as a logical vector or an integer
//' vector. MPMs without age structure should be entered as \code{0} or
//' \code{FALSE}. Do not use if no MPM has age structure.
//' @param fecmod An optional vector used for function-based MPMs giving scalar
//' multipliers for fecundity terms, when two fecundity variables are used for a
//' collective fecundity per individual. Each entry refers to each 
//' \code{vrm_input} object in argument \code{vrms}, in the same order.
//' @param starts An optional list of \code{lefkoSV} objects, which are data
//' frames providing the starting numbers of individuals of each stage. If
//' provided, then one is needed per MPM. If not provided, then all projections
//' start with a single individual of each stage per MPM.
//' @param patches An optional string vector with length equal to the number of
//' MPMs, detailing the name of each patch to project for each MPM, in order.
//' Only a single pop-patch may be projected for each MPM given. A value of
//' \code{NA} can be supplied to indicate that the population-level matrices
//' should be projected (if argument \code{mpms} is used and a population-level
//' set of matrices exist), or that the first patch noted should be used.
//' Defaults to the population-level set or the first patch, depending on
//' whether the former exists.
//' @param years An optional term corresponding either to a single integer
//' vector of time \code{t} values, if all MPMs will use the same time \code{t}
//' or set of time \code{t}'s, or a list of such vectors with each vector
//' corresponding to each MPM in order. In the latter case, a vector composed of
//' a single \code{NA} value is interpreted to mean that all time \code{t}
//' values in the MPM should be utilized. If a vector shorter than \code{times}
//' is supplied, then this vector will be cycled.
//' @param tweights An optional list composed of numeric vectors or matrices
//' denoting the probabilities of choosing each matrix in each MPM in a
//' stochastic projection. If an element of the list is a matrix, then a
//' first-order Markovian environment is assumed, in which the probability of
//' choosing a specific annual matrix depends on which annual matrix is
//' currently chosen. If an element of the list is a vector, then the choice of
//' annual matrix is assumed to be independent of the current matrix. Defaults
//' to equal weighting among matrices. If used, then one element per MPM is
//' required, with equal weighting assumed for any element set to \code{NULL}.
//' @param density An optional list of data frames of class \code{lefkoDens},
//' which provide details for density dependence in MPM elements and have been
//' created with function \code{\link[lefko3]{density_input}()}. If used, then
//' one such data frame per MPM is required. MPMs to be run without density
//' dependence should be set to \code{NULL}.
//' @param entry_time An optional integer vector giving the entry time for each
//' MPM into the projection. Defaults to a zero vector with the length of the
//' number of MPMs, as given either by argument \code{mpms} or \code{vrms}.
//' @param density_vr An optional list of data frames of class
//' \code{lefkoDensVR}, which provide details for density dependence in vital
//' rate models and have been created with function
//' \code{link[lefko3]{density_vr}()}. If used, then one such data frame per MPM
//' is required. MPMs to be run without vital describing density dependence
//' relationships in vital rates should be set to \code{NULL}. Can only be used
//' with function-based projections.
//' @param sp_density An optional argument for use with \code{vrm_input} objects
//' that specifies the spatial density to be used in each time step. If used,
//' may either be a numeric vector giving a single spatial density for each
//' \code{vrm_input} object entered in argument \code{vrms} (in this case, the
//' value of spatial density given for each \code{vrm_input} object will be held
//' constant through the projection), or a list of as many numeric vectors as
//' \code{vrm_input} objects, with the length of each vector giving the spatial
//' density at each time step. If vectors are shorter than specified in 
//' \code{times}, then these values will be cycled.
//' @param ind_terms An optional argument providing values of individual or
//' environmental covariate values for \code{vrm_input} objects used in
//' function-based projection. Can be set either to a single data frame with 3
//' columns giving values for up to 3 covariates across time (rows give the time
//' order of these values), or a list of as many such data frames as
//' \code{vrm_input} objects. In the latter case, \code{vrm_input} objects that
//' do not use such covariates should have the associated element set to
//' \code{NULL}. Unused terms within each data frame must be set to \code{0}
//' (use of \code{NA} will produce errors.) If the number of rows is less than
//' \code{times}, then these values will be cycled.
//' @param dev_terms An optional list of data frames, one for each
//' \code{vrm_input} object. Each should include 14 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival,
//' 9: juvenile observation, 10: juvenile primary size, 11: juvenile secondary
//' size, 12: juvenile tertiary size, 13: juvenile reproduction, and
//' 14: juvenile maturity transition. Unused terms must be set to \code{0} (use
//' of \code{NA} will produce errors). Single or small numbers of values per
//' vital rate model are also allowed, and if the number of rows is less than
//' \code{times}, then the terms will be cycled.
//' @param fb_sparse A logical vector indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE} for
//' each MPM.
//' @param equivalence An optional numeric vector or list of numeric vectors. If
//' a numeric Vector, then must have the same number of elements as the number
//' of MPMs, with each element giving the effect of an individual of each
//' MPM relative to a reference individual. If a list, then the list should be
//' composed of as many numeric vectors as MPMs, with each vector giving the
//' effect of each individual in each stage relative to a reference individual.
//' Numeric entries used in these vectors can be thought of as Lotka-Volterra
//' interaction terms, such as are used in multiple species competition models.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param theta_tol A numeric value used to indicate a maximum value to theta as
//' used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' @param prep_mats An integer value for use when creating function-based MPM
//' projections. If using \code{vrms} input instead of \code{mpms} input, then
//' this argument determines how many matrices should be used as a limit to
//' develop matrices prior to running the projection. See \code{Notes} for
//' further details.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent and density
//' independent simulations. Defaults to \code{0}, which does not enforce
//' substochasticity. Alternatively, \code{1} forces all survival-transition
//' elements to range from 0.0 to 1.0, and forces fecundity to be non-negative;
//' and \code{2} forces all column rows in the survival-transition matrices to
//' total no more than 1.0, in addition to the actions outlined for option
//' \code{1}. Both settings \code{1} and \code{2} change negative fecundity
//' elements to \code{0.0}.
//' @param force_fb A logical value indicating whether to force function-based
//' MPMs to be developed at each time step even if fewer than \code{prep_mats}.
//' Defaults to \code{FALSE}.
//' @param stochastic A Boolean value indicating to perform a temporally
//' stochastic projection.
//' @param err_check A logical value indicating whether to add extra terms to
//' the output.
//' 
//' @return A list of R-defined objects, including vectors, lists, integers, and
//' data frames, for use in later stages of analysis.
//' 
//' @keywords internal
//' @noRd
Rcpp::List cleanup3(Nullable<RObject> mpms = R_NilValue,
  Nullable<RObject> vrms = R_NilValue, Nullable<RObject> stageframes = R_NilValue,
  Nullable<RObject> supplements = R_NilValue, Nullable<RObject> format = R_NilValue,
  Nullable<RObject> firstage = R_NilValue, Nullable<RObject> finalage = R_NilValue,
  Nullable<RObject> fecage_min = R_NilValue, Nullable<RObject> fecage_max = R_NilValue,
  Nullable<RObject> cont = R_NilValue, Nullable<RObject> fecmod = R_NilValue,
  Nullable<RObject> starts = R_NilValue, Nullable<RObject> patches = R_NilValue,
  Nullable<RObject> years = R_NilValue, Nullable<RObject> tweights = R_NilValue,
  Nullable<RObject> density = R_NilValue, Nullable<RObject> entry_time = R_NilValue,
  Nullable<RObject> density_vr = R_NilValue, Nullable<RObject> sp_density = R_NilValue,
  Nullable<RObject> ind_terms = R_NilValue, Nullable<RObject> dev_terms = R_NilValue,
  Nullable<RObject> fb_sparse = R_NilValue, Nullable<RObject> equivalence = R_NilValue,
  double exp_tol = 700.0, double theta_tol = 100000000.0, int prep_mats = 20,
  int substoch = 0, bool force_fb = false, bool stochastic = false,
  bool err_check = false) {
  
  List mpm_list;
  List vrm_list;
  List A_list;
  List stageframe_list;
  List stageframe_list_fb;
  List supplement_list;
  List supplement_list_fb;
  List repmatrix_list;
  List hstages_list;
  List agestages_list;
  List start_list;
  List year_list;
  //List patch_list;
  List labels_list;
  List tweights_list;
  List density_list;
  List dens_index_list;
  List density_vr_list;
  List ind_terms_num_list;
  List ind_terms_cat_list;
  List dev_terms_list;
  List sp_density_list;
  List equivalence_list;
  List allstages_all;
  List allmodels_all;
  
  IntegerVector stagecounts;
  LogicalVector sparse_vec;
  IntegerVector format_vec;
  IntegerVector firstage_vec;
  IntegerVector finalage_vec;
  IntegerVector cont_vec;
  NumericVector fecmod_vec;
  IntegerVector fecage_min_vec;
  IntegerVector fecage_max_vec;
  IntegerVector matrowcounts;
  CharacterVector patch_vec;
  IntegerVector total_years_vec;
  IntegerVector tweights_type_vec;
  IntegerVector dens_yn_vec;
  IntegerVector entry_time_vec;
  IntegerVector dens_vr_yn_vec;
  IntegerVector sp_density_num_vec;
  IntegerVector dev_terms_num_vec;
  IntegerVector inda_terms_num_vec;
  IntegerVector indb_terms_num_vec;
  IntegerVector indc_terms_num_vec;
  IntegerVector inda_terms_cat_vec;
  IntegerVector indb_terms_cat_vec;
  IntegerVector indc_terms_cat_vec; 
  NumericVector equivalence_vec;
  
  DataFrame labels;
  
  int mpm_count {0};
  int vrm_count {0};
  int total_mpms {0};
  int stageframe_count {0};
  int stageframe_notNull_count {0};
  int supplement_count {0};
  int sparse_vec_count {0};
  int start_count {0};
  int found_fleslie {0};
  int tweights_count {0};
  int density_count {0};
  int entry_time_count {0};
  int density_vr_count {0};
  int equivalence_count {0};

  bool preexisting {false};
  bool funcbased {false};
  bool pure_fleslie {false};
  bool entry_time_vec_use {false};
  bool stages_not_equal {false};
  bool fb_override {false};
  
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Argument substoch must equal 0, 1, or 2.", false);
  }
  
  //Rcout << "cleanup3 A    ";
  
  if (mpms.isNotNull()) {
    if (vrms.isNotNull() || stageframes.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "stageframes", "projecting existing MPMms", 24);
    }
    
    if (is<List>(mpms)) { 
      mpm_list = as<List>(mpms);
      mpm_count = static_cast<int>(mpm_list.length());
      
    } else {
      AdaptUtils::pop_error2("mpms", "a list of lefkoMat objects", "", 1);
    }
    
    List found_stageframes (mpm_count);
    List A_list_pre (mpm_count);
    LogicalVector sparse_vec_temp (mpm_count);
    
    for (int i = 0; i < mpm_count; i++) {
      if (is<List>(mpm_list(i))) {
        List chosen_mpm = as<List>(mpm_list(i));
        
        if (!chosen_mpm.hasAttribute("class")) {
          AdaptUtils::pop_error2("mpms", "a list of lefkoMat objects", "", 1);
        }
        CharacterVector chosen_mpm_class = chosen_mpm.attr("class");
        
        bool found_lMat {false};
        for (int j = 0; j < static_cast<int>(chosen_mpm_class.length()); j++) {
          if (chosen_mpm_class(j) == "lefkoMat") found_lMat = true;
        }
        
        if (!found_lMat) {
          AdaptUtils::pop_error2("mpms", "a list of lefkoMat objects", "", 1);
        }
        
        DataFrame found_stageframe = as<DataFrame>(chosen_mpm["ahstages"]);
        found_stageframes(i) = found_stageframe;
        
        List current_A_list = as<List>(chosen_mpm["A"]);
        A_list_pre(i) = current_A_list;
        
        if (is<S4>(current_A_list(0))) {
          sparse_vec_temp(i) = true;
        }
        
      } else {
        AdaptUtils::pop_error2("mpms", "a list of lefkoMat objects", "", 1);
      }
    }
    
    A_list = A_list_pre;
    stageframe_list = found_stageframes;
    stageframe_count = mpm_count;
    sparse_vec = sparse_vec_temp;
    sparse_vec_count = mpm_count;
    preexisting = true;
  }
  
  //Rcout << "cleanup3 B    ";
  
  if (vrms.isNotNull()) {
    if (is<List>(vrms)) {
      vrm_list = as<List>(vrms);
      vrm_count = static_cast<int>(vrm_list.length());
      
      if (format.isNotNull()) {
        if (is<NumericVector>(format) || is<IntegerVector>(format)) {
          format_vec = as<IntegerVector>(format);
          int format_count = static_cast<int>(format_vec.length());
          
          if (format_count != vrm_count) {
            AdaptUtils::pop_error2("format", "vrm_input objects", "vrms", 2);
          }
          
          for (int i = 0; i < vrm_count; i++) {
            if (IntegerVector::is_na(format_vec(i))) {
              AdaptUtils::pop_error2("NA values", "format", "", 25);
            }
            if (format_vec(i) != 5 && !stageframes.isNotNull()) {
              AdaptUtils::pop_error2("stageframes", "run function-based projections", "", 26);
            } else if (format_vec(i) == 5 && !stageframes.isNotNull()) {
              found_fleslie++;
            }
          }
          if (found_fleslie == vrm_count) pure_fleslie = true;
        } else AdaptUtils::pop_error2("format", "an integer vector", "", 1);
      } else if (funcbased) {
        IntegerVector format_vec_pre (vrm_count, 3);
        format_vec = format_vec_pre;
        
        if (!stageframes.isNotNull()) {
          AdaptUtils::pop_error2("stageframes", "run function-based projections", "", 26);
        }
      }
    } else AdaptUtils::pop_error2("vrms", "a list of vrm_input objects", "", 1);
    
    if (is<List>(stageframes)) {
      stageframe_list_fb = as<List>(stageframes);
      stageframe_count = static_cast<int>(stageframe_list_fb.length());
      
      if (stageframe_count != vrm_count) {
        throw Rcpp::exception("A stageframe is required for each vrm_input object.", false);
      }
    } else {
      bool throw_error {false};
      for (int i = 0; i < vrm_count; i++) {
        if (format_vec(i) != 5) throw_error = true;
      }
      if (throw_error) AdaptUtils::pop_error2("stageframes", "a list of stageframe objects", "", 1);
      pure_fleslie = true;
    }
    
    for (int i = 0; i < vrm_count; i++) {
      if (is<List>(vrm_list(i))) {
        List chosen_vrm = as<List>(vrm_list(i));
        
        if (!chosen_vrm.hasAttribute("class")) {
          AdaptUtils::pop_error2("vrms", "a list of vrm_input objects", "", 1);
        }
        CharacterVector chosen_vrm_class = chosen_vrm.attr("class");
        
        bool found_vrmi {false};
        for (int j = 0; j < static_cast<int>(chosen_vrm_class.length()); j++) {
          if (chosen_vrm_class(j) == "vrm_input") found_vrmi = true;
        }
        
        if (!found_vrmi) AdaptUtils::pop_error2("vrms", "a list of vrm_input objects", "", 1);
      } else AdaptUtils::pop_error2("vrms", "a list of vrm_input objects", "", 1);
    }
    
    for (int i = 0; i < stageframe_count; i++) {
      if (is<DataFrame>(stageframe_list_fb(i))) {
        DataFrame chosen_stageframe = as<DataFrame>(stageframe_list_fb(i));
        
        if (!chosen_stageframe.hasAttribute("class")) {
          AdaptUtils::pop_error2("stageframes", "a list of stageframe objects", "", 1);
        }
        CharacterVector chosen_stageframe_class = chosen_stageframe.attr("class");
        
        bool found_stageframe {false};
        for (int j = 0; j < static_cast<int>(chosen_stageframe_class.length()); j++) {
          if (chosen_stageframe_class(j) == "stageframe") found_stageframe = true;
        }
        
        if (!found_stageframe) {
          AdaptUtils::pop_error2("stageframes", "a list of stageframe objects", "", 1);
        }
        stageframe_notNull_count++;
      } else if (stageframe_list_fb(i) == R_NilValue) {
        if (format_vec(i) != 5) {
          throw Rcpp::exception("All non-Leslie MPMs need stageframes.", false);
        }
      } else {
        AdaptUtils::pop_error2("stageframes", "a list of stageframe objects", "", 1);
      }
    }
    
    if ((vrm_count - found_fleslie) != stageframe_notNull_count && !pure_fleslie) {
      throw Rcpp::exception("Each vrm_input object must have its own stageframe.",
        false);
    }
    
    funcbased = true;
  }
  
  //Rcout << "cleanup3 C    ";
  
  if (!preexisting && !funcbased) {
    throw Rcpp::exception("Cannot proceed without either argument mpms, or arguments vrms and stageframes set.",
      false);
  } else if (preexisting && funcbased) {
    throw Rcpp::exception("Cannot proceed with argument mpms, vrms, and stageframes set.",
      false);
  }
  
  total_mpms = mpm_count + vrm_count;
  
  if (supplements.isNotNull()) {
    if (!funcbased) {
      throw Rcpp::exception("Argument supplements can only be used with argument vrms.",
        false);
    }
    
    if (is<List>(supplements)) {
      supplement_list_fb = as<List>(supplements);
      supplement_count = static_cast<int>(supplement_list_fb.length());
      if (supplement_count != vrm_count) {
        AdaptUtils::pop_error2("vrms", "supplements", "lists of the same length", 27);
      }
      
      for (int i = 0; i < supplement_count; i++) {
        if (is<DataFrame>(supplement_list_fb(i))) {
          DataFrame chosen_supplement = as<DataFrame>(supplement_list_fb(i));
          
          if (chosen_supplement.hasAttribute("class")) {
            CharacterVector chosen_supplement_class = chosen_supplement.attr("class");
            bool found_lefkoSD {false};
            
            for (int j = 0; j < static_cast<int>(chosen_supplement_class.length()); j++) {
              if (chosen_supplement_class(j) == "lefkoSD") found_lefkoSD = true;
            }
            if (!found_lefkoSD) {
              AdaptUtils::pop_error2("supplements", "a list of lefkoSD objects and NULL values only", "", 1);
            }
          } else {
            AdaptUtils::pop_error2("supplements", "a list of lefkoSD objects and NULL values only", "", 1);
          }
        } else if (supplement_list_fb(i) != R_NilValue) { 
          AdaptUtils::pop_error2("supplements", "a list of lefkoSD objects and NULL values only", "", 1);
        }
      }
      
    } else {
      AdaptUtils::pop_error2("supplements", "a list of lefkoSD objects and NULL values only", "", 1);
    }
  } else {
    int used_mpm_count = mpm_count;
    if (vrm_count > mpm_count) used_mpm_count = vrm_count;
    
    if (funcbased) {
      List supplement_list_fb_pre (used_mpm_count);
      supplement_list_fb = supplement_list_fb_pre;
    }
  }
  
  //Rcout << "cleanup3 D    ";
  
  if (format.isNotNull()) {
    if (!funcbased) {
      AdaptUtils::pop_error2("vrms", "use argument format", "", 26);
    }
    
  } else if (preexisting) {
    IntegerVector format_vec_pre (mpm_count, 3);
    
    for (int i = 0; i < mpm_count; i++) {
      List chosen_mpm = as<List>(mpm_list(i));
      RObject hstages_element = as<RObject>(chosen_mpm["hstages"]);
      RObject agestages_element = as<RObject>(chosen_mpm["agestages"]);
      
      if (!is<LogicalVector>(hstages_element)) {
        if (is<DataFrame>(hstages_element)) {
          DataFrame hst_input = as<DataFrame>(hstages_element);
          int hst_cols = hst_input.length();
          
          if (hst_cols > 1) format_vec_pre(i) = 1;
        }
      }
      
      if (!is<LogicalVector>(agestages_element) && format_vec_pre(i) == 3) {
        if (is<DataFrame>(agestages_element)) {
          DataFrame ast_input = as<DataFrame>(agestages_element);
          int ast_cols = ast_input.length();
          
          if (ast_cols > 1) format_vec_pre(i) = 4;
        }
      }
    }
    format_vec = format_vec_pre;
  } else if (funcbased) {
    IntegerVector format_vec_pre (vrm_count, 3);
    format_vec = format_vec_pre;
  }
  
  //Rcout << "cleanup3 E    ";
  
  // firstage, finalage, and cont processing for age-by-stage and Leslie MPMs
  if (firstage.isNotNull() || finalage.isNotNull() || cont.isNotNull()) {
    if (!funcbased) {
      AdaptUtils::pop_error2("vrms", "use arguments firstage, finalage, and cont", "", 26);
    }
    
    bool found_age_MPM {false};
    for (int i = 0; i < static_cast<int>(format_vec.length()); i++) {
      if (format_vec(i) > 3) found_age_MPM = true;
    }
    
    if (!found_age_MPM) {
      AdaptUtils::pop_error2("firstage, finalage, and cont", "age-based, function-based MPMs", "", 28);
    }
  }
  
  if (firstage.isNotNull()) {
    if (is<IntegerVector>(firstage) || is<NumericVector>(firstage)) { 
      firstage_vec = as<IntegerVector>(firstage);
      
      int firstage_vec_length = static_cast<int>(firstage_vec.length());
      if (firstage_vec_length != vrm_count && firstage_vec_length != mpm_count) {
        AdaptUtils::pop_error2("firstage", "number of MPMs to project", "", 29);
      }
      
      for (int i = 0; i < firstage_vec_length; i++) {
        if (firstage_vec(i) < 0) {
          AdaptUtils::pop_error2("firstage", "", "", 30);
        }
        
        if (firstage_vec(i) > 0 && format_vec(i) < 4) {
          throw Rcpp::exception("Entries in argument firstage must equal 0 for MPMs without age structure.", false);
        }
      }
    } else AdaptUtils::pop_error2("firstage", "an integer vector", "", 1);
  } else {
    if (preexisting) {
      IntegerVector firstage_vec_pre (mpm_count);
      firstage_vec = firstage_vec_pre;
    } else {
      IntegerVector firstage_vec_pre (vrm_count);
      firstage_vec = firstage_vec_pre;
    }
  }
  
  if (finalage.isNotNull()) {
    if (is<IntegerVector>(finalage) || is<NumericVector>(finalage)) { 
      finalage_vec = as<IntegerVector>(finalage);
      
      int finalage_vec_length = static_cast<int>(finalage_vec.length());
      if (finalage_vec_length != vrm_count && finalage_vec_length != mpm_count) {
        AdaptUtils::pop_error2("finalage", "number of MPMs to project", "", 29);
      }
      
      for (int i = 0; i < finalage_vec_length; i++) {
        if (finalage_vec(i) < 0) {
          AdaptUtils::pop_error2("finalage", "", "", 30);
        }
        
        if (finalage_vec(i) > 0 && format_vec(i) < 4) {
          throw Rcpp::exception("Entries in argument finalage must equal 0 for MPMs without age structure.", false);
        }
        
        if (finalage_vec(i) < firstage_vec(i)) {
          throw Rcpp::exception("Entries in argument finalage may not be less than respective entries in argument firstage.", false);
        }
      }
    } else AdaptUtils::pop_error2("finalage", "an integer vector", "", 1);
  } else {
    if (preexisting) {
      IntegerVector finalage_vec_pre (mpm_count);
      finalage_vec = finalage_vec_pre;
    } else {
      IntegerVector finalage_vec_pre (vrm_count);
      finalage_vec = finalage_vec_pre;
    }
  }
  
  //Rcout << "cleanup3 F    ";
  
  if (fecage_min.isNotNull()) {
    if (is<IntegerVector>(fecage_min) || is<NumericVector>(fecage_min)) {
      IntegerVector fecage_min_prevec = as<IntegerVector>(fecage_min);
      
      int fecage_min_length = static_cast<int>(fecage_min_prevec.length());
      
      if (fecage_min_length == vrm_count) {
        fecage_min_vec = fecage_min_prevec;
      } else if (fecage_min_length == 1) {
        fecage_min_vec = rep(fecage_min_prevec(0), vrm_count);
      } else AdaptUtils::pop_error2("fecage_min", "vrm_input objects", "vrms", 2);
    } else AdaptUtils::pop_error2("fecage_min", "an integer vector", "", 1);
  } else {
    fecage_min_vec = firstage_vec;
  }
  
  if (fecage_max.isNotNull()) {
    if (is<IntegerVector>(fecage_max) || is<NumericVector>(fecage_max)) {
      IntegerVector fecage_max_prevec = as<IntegerVector>(fecage_max);
      
      int fecage_max_length = static_cast<int>(fecage_max_prevec.length());
      
      if (fecage_max_length == vrm_count) {
        fecage_max_vec = fecage_max_prevec;
      } else if (fecage_max_length == 1) {
        fecage_max_vec = rep(fecage_max_prevec(0), vrm_count);
      } else AdaptUtils::pop_error2("fecage_max", "vrm_input objects", "vrms", 2);
    } else AdaptUtils::pop_error2("fecage_max", "an integer vector", "", 1);
  } else {
    fecage_max_vec = finalage_vec;
  }
  
  if (fecmod.isNotNull()) {
    if (is<NumericVector>(fecmod) && funcbased) {
      NumericVector fecmod_pre = as<NumericVector>(fecmod);
      
      if (static_cast<int>(fecmod_pre.length()) != vrm_count) {
        AdaptUtils::pop_error2("fecmod", "number of vrm inputs", "", 29);
      }
      
      fecmod_vec = fecmod_pre;
      
    } else if (!funcbased) {
      AdaptUtils::pop_error2("vrms", "use argument fecmod", "", 26);
    }
  } else if (funcbased) {
    NumericVector fecmod_vec_pre (vrm_count, 1.0);
    
    fecmod_vec = fecmod_vec_pre;
  }
  
  //Rcout << "cleanup3 G    ";
  
  if (cont.isNotNull()) {
    if (is<IntegerVector>(cont) || is<NumericVector>(cont)) { 
      IntegerVector cont_vec_temp = as<IntegerVector>(cont);
      
      if (static_cast<int>(cont_vec_temp.length() == 1)) {
        IntegerVector cont_maxed_out (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          cont_maxed_out(i) = cont_vec_temp(0);
        }
        cont_vec = cont_maxed_out;
      } else cont_vec = cont_vec_temp;
      
      int cont_vec_length = static_cast<int>(cont_vec.length());
      if (cont_vec_length != vrm_count) {
        AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
      }
      
      for (int i = 0; i < cont_vec_length; i++) {
        if (cont_vec(i) < 0 || cont_vec(i) > 1) {
          throw Rcpp::exception("Entries in argument cont must equal 0 or 1.",
            false);
        }
        
        if (cont_vec(i) > 0 && format_vec(i) < 4) {
          throw Rcpp::exception("Entries in argument cont must equal 0 for MPMs without age structure.", 
            false);
        }
        
        if (IntegerVector::is_na(cont_vec(i))) {
          cont_vec(i) = 0;
        }
      }
    } if (is<LogicalVector>(cont)) {
      LogicalVector cont_vec_temp = as<LogicalVector>(cont);
      
      LogicalVector cont_vec_log;
      
      if (static_cast<int>(cont_vec_temp.length() == 1)) {
        LogicalVector cont_maxed_out (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          cont_maxed_out(i) = cont_vec_temp(0);
        }
        cont_vec_log = cont_maxed_out;
      } else cont_vec_log = cont_vec_temp;
      
      int cont_vec_log_length = static_cast<int>(cont_vec_log.length());
      if (cont_vec_log_length != vrm_count) {
        AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
      }
      
      IntegerVector cont_vec_pre (cont_vec_log_length);
      for (int i = 0; i < cont_vec_log_length; i++) {
        if (LogicalVector::is_na(cont_vec_log(i))) {
          cont_vec_pre(i) = 0;
        } else if (cont_vec_log(i)) {
          cont_vec_pre(i) = 1;
        } else if (!cont_vec_log(i)) {
          cont_vec_pre(i) = 0;
        }
      }
      
      cont_vec = cont_vec_pre;
      
    } else {
      AdaptUtils::pop_error2("cont", "an integer vector", "", 1);
    }
  } else {
    if (preexisting) {
      IntegerVector cont_vec_pre (mpm_count);
      cont_vec = cont_vec_pre;
    } else {
      IntegerVector cont_vec_pre (vrm_count);
      cont_vec = cont_vec_pre;
    }
  }
  
  //Rcout << "cleanup3 H    ";
  
  // Altered stageframe processing
  if (funcbased) {
    List stageframe_list_pre (vrm_count);
    List supplement_list_pre (vrm_count);
    List repmatrix_list_pre (vrm_count);
    
    List hstages_list_pre (vrm_count);
    List agestages_list_pre (vrm_count);
    IntegerVector matrows_pre (vrm_count);
    
    for (int i = 0; i < vrm_count; i++) {
      DataFrame chosen_stageframe;
      
      RObject supp_trial = RObject(supplement_list_fb(i));
      DataFrame chosen_supplement;
      bool trial_supp_null {false};
      
      if (is<DataFrame>(supp_trial)) {
        chosen_supplement = as<DataFrame>(supplement_list_fb(i));
      } else {
        trial_supp_null = true;
      }
      
      if (format_vec(i) < 5) {
        chosen_stageframe = as<DataFrame>(stageframe_list_fb(i));
        
        bool agemat = false;
        bool historical = false;
        int ehrlen {1};
        //int style {0};
        //int filter {1};
        
        if (format_vec(i) == 2) ehrlen = 2;
        //if (format_vec(i) == 3) style = 1;
        if (format_vec(i) == 4) {
          agemat = true;
          //style = 2;
          //filter = 2;
        }
        if (format_vec(i) < 3) historical = true;
        
        List melchett;
        
        if (!trial_supp_null) {
          melchett = LefkoMats::sf_reassess_internal(chosen_stageframe,
            chosen_supplement, R_NilValue, R_NilValue, agemat, historical,
            ehrlen);
        } else {
          melchett = LefkoMats::sf_reassess_internal(chosen_stageframe,
            R_NilValue, R_NilValue, R_NilValue, agemat, historical, ehrlen);
        }
        DataFrame new_stageframe = as<DataFrame>(melchett["stageframe"]);
        
        stageframe_list_pre(i) = new_stageframe;
        
        arma::mat new_repmatrix = as<arma::mat>(melchett["repmatrix"]);
        repmatrix_list_pre(i) = new_repmatrix;
        
        if (format_vec(i) < 4) {
          DataFrame new_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          if (new_ovtable_temp.containsElementNamed("stage3")) {
            supplement_list_pre(i) = new_ovtable_temp;
          } else {
            StringVector nsst3 = {};
            IntegerVector nsa2 = {};
            NumericVector nsgr = {};
            
            DataFrame intro_ovtable = DataFrame::create(_["stage3"] = nsst3,
              _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
              _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
              _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
              _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
              _["multiplier"] = clone(nsgr), _["convtype"] = clone(nsa2),
              _["convtype_t12"] = clone(nsa2), _["pop"] = clone(nsst3),
              _["patch"] = clone(nsst3), _["year2"] = clone(nsst3));
            supplement_list_pre(i) = intro_ovtable;
          }
        } else {
          DataFrame new_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          if (new_ovtable_temp.containsElementNamed("stage3")) {
            supplement_list_pre(i) = LefkoMats::age_expanded(new_ovtable_temp,
              firstage_vec(i), finalage_vec(i));
          } else {
            StringVector nsst3 = {};
            IntegerVector nsa2 = {};
            NumericVector nsgr = {};
            
            DataFrame intro_ovtable = DataFrame::create(_["stage3"] = nsst3,
              _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
              _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
              _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
              _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
              _["multiplier"] = clone(nsgr), _["convtype"] = clone(nsa2),
              _["convtype_t12"] = clone(nsa2), _["pop"] = clone(nsst3),
              _["patch"] = clone(nsst3), _["year2"] = clone(nsst3));
            supplement_list_pre(i) = intro_ovtable;
          }
        }
        
        DataFrame chosen_stageframe_pre = clone(as<DataFrame>(stageframe_list_pre(i)));
        
        IntegerVector removal_row = {static_cast<int>(chosen_stageframe_pre.nrows())};
        StringVector removal_var = {"stage_id"};
        DataFrame chosen_stageframe = LefkoUtils::df_remove(chosen_stageframe_pre,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        if (format_vec(i) < 3) {
          DataFrame hstages_temp;
          hst_maker(hstages_temp, chosen_stageframe, format_vec(i));
          hstages_list_pre(i) = hstages_temp;
          
          int found_numrows = static_cast<int>(hstages_temp.nrows());
          matrows_pre(i) = found_numrows;
        } else if (format_vec(i) == 4) {
          DataFrame agestages_temp = age_maker(chosen_stageframe,
            firstage_vec(i), finalage_vec(i));
          agestages_list_pre(i) = agestages_temp;
          
          int found_numrows = static_cast<int>(agestages_temp.nrows());
          matrows_pre(i) = found_numrows;
          
        } else {
          int found_numrows = static_cast<int>(chosen_stageframe.nrows());
          matrows_pre(i) = found_numrows;
        }
        
      } else {
        bool cont_used {false};
        if (cont_vec(i) > 0) cont_used = true;
        
        DataFrame melchett = LefkoMats::sf_leslie(firstage_vec(i),
          finalage_vec(i), fecage_min_vec(i), fecage_max_vec(i), cont_used);
        DataFrame new_stageframe = melchett;
        
        DataFrame new_ovtable;
        if (!trial_supp_null) {
          new_ovtable = LefkoMats::age_expanded(chosen_supplement,
            firstage_vec(i), finalage_vec(i));
        }
        stageframe_list_pre(i) = new_stageframe;
        supplement_list_pre(i) = new_ovtable;
        
        stageframe_count++;
      }
    }
    stageframe_list = stageframe_list_pre;
    supplement_list = supplement_list_pre;
    repmatrix_list = repmatrix_list_pre;
    
    hstages_list = hstages_list_pre;
    agestages_list = agestages_list_pre;
    matrowcounts = matrows_pre;
    
  } else {
    List hstages_list_pre (mpm_count);
    List agestages_list_pre (mpm_count);
    IntegerVector matrows_pre (mpm_count);
    
    for (int i = 0; i < mpm_count; i++) {
      List chosen_mpm = as<List>(mpm_list(i));
      if (format_vec(i) < 3) {
        DataFrame chosen_hstages = as<DataFrame>(chosen_mpm["hstages"]);
        hstages_list_pre(i) = chosen_hstages;
        
        int found_numrows = static_cast<int>(chosen_hstages.nrows());
        matrows_pre(i) = found_numrows;
      } else if (format_vec(i) == 4) {
        DataFrame chosen_agestages = as<DataFrame>(chosen_mpm["agestages"]);
        agestages_list_pre(i) = chosen_agestages;
        
        int found_numrows = static_cast<int>(chosen_agestages.nrows());
        matrows_pre(i) = found_numrows;
        
        IntegerVector all_ages_agestages = chosen_agestages["age"];
        int min_age_agestages = min(all_ages_agestages);
        int max_age_agestages = max(all_ages_agestages);
        
        firstage_vec(i) = min_age_agestages;
        finalage_vec(i) = max_age_agestages;
        
      } else {
        DataFrame chosen_ahstages = as<DataFrame>(chosen_mpm["ahstages"]);
        
        int found_numrows = static_cast<int>(chosen_ahstages.nrows());
        matrows_pre(i) = found_numrows;
      }
    }
    hstages_list = hstages_list_pre;
    agestages_list = agestages_list_pre;
    matrowcounts = matrows_pre;
  }
  
  //Rcout << "cleanup3 I    ";
  
  // start vector
  if (starts.isNotNull()) {
    if (is<List>(starts)) {
      List start_list_pre = as<List>(starts);
      start_count = static_cast<int>(start_list_pre.length());
      IntegerVector stagecounts_pre (stageframe_count);
      
      List start_list_new (stageframe_count);
      
      for (int i = 0; i < start_count; i++) {
        if (is<DataFrame>(start_list_pre(i))) {
          DataFrame chosen_start = as<DataFrame>(start_list_pre(i));
        
          if (!chosen_start.hasAttribute("class")) {
            AdaptUtils::pop_error2("starts", "a list of lefkoSV objects", "", 1);
          }
          CharacterVector chosen_start_class = chosen_start.attr("class");
          
          bool found_lSt {false};
          for (int j = 0; j < static_cast<int>(chosen_start_class.length()); j++) {
            if (chosen_start_class(j) == "lefkoSV") found_lSt = true;
          }
          
          if (!found_lSt) {
            AdaptUtils::pop_error2("starts", "a list of lefkoSV objects", "", 1);
          }
          
          DataFrame chosen_stageframe;
          if (format_vec(i) < 3) {
            if (preexisting) {
              List chosen_mpm = as<List>(mpm_list(i));
              chosen_stageframe = as<DataFrame>(chosen_mpm["hstages"]);
            } else {
              chosen_stageframe = as<DataFrame>(hstages_list(i));
            }
          } else if (format_vec(i) == 4) {
            if (preexisting) {
              List chosen_mpm = as<List>(mpm_list(i));
              chosen_stageframe = as<DataFrame>(chosen_mpm["agestages"]);
            } else {
              chosen_stageframe = as<DataFrame>(agestages_list(i));
            }
          } else {
            chosen_stageframe = as<DataFrame>(stageframe_list(i));
          }
          
          int scp = static_cast<int>(chosen_stageframe.nrows());
          if (format_vec(i) == 3 && funcbased) scp--;
          stagecounts_pre(i) = scp;
          
          arma::vec start_vec (scp, fill::zeros);
          arma::uvec start_elems = as<arma::uvec>(chosen_start["row_num"]);
          start_elems = start_elems - 1;
          arma::vec start_values = as<arma::vec>(chosen_start["value"]);
          
          if (static_cast<int>(start_elems.max()) > (scp - 1)) {
            throw Rcpp::exception("lefkoStart object includes element indices too high for associated MPM.",
              false);
          }
          
          for (int j = 0; j < static_cast<int>(start_elems.n_elem); j++) {
            start_vec(start_elems(j)) = start_values(j);
          }
          
          start_list_new(i) = start_vec;
        } else {
          AdaptUtils::pop_error2("starts", "a list of lefkoSV objects", "", 1);
        }
      }
      
      stagecounts = stagecounts_pre;
      start_list = start_list_new;
      
      if (start_count != mpm_count && start_count != vrm_count) {
        throw Rcpp::exception("In argument starts, list length entered must equal number of MPMs.",
          false);
      }
    } else {
      AdaptUtils::pop_error2("starts", "a list of lefkoSV objects", "", 1);
    }
  } else {
    // Construct default list of start vectors (1 indiv / stage)
    IntegerVector stagecounts_pre (stageframe_count);
    List start_list_pre (stageframe_count);
    
    for (int i = 0; i < stageframe_count; i++) {
      DataFrame chosen_stageframe;
      
      if (format_vec(i) == 3 || format_vec(i) == 5) {
        chosen_stageframe = as<DataFrame>(stageframe_list(i));
      } else if (format_vec(i) < 3) {
        chosen_stageframe = as<DataFrame>(hstages_list(i));
      } else if (format_vec(i) == 4) {
        chosen_stageframe = as<DataFrame>(agestages_list(i));
      }
      
      int scp = static_cast<int>(chosen_stageframe.nrows());
      if (format_vec(i) == 3 && funcbased) scp--;
      stagecounts_pre(i) = scp;
      
      arma::vec start_vec (scp, fill::ones);
      start_list_pre(i) = start_vec;
    }
    
    stagecounts = stagecounts_pre;
    start_list = start_list_pre;
    start_count = stageframe_count;
  }
  
  //Rcout << "cleanup3 J    ";
  
  // patches vector
  if (patches.isNotNull()) {
    //List patch_list_pre (total_mpms);
    List labels_list_pre (total_mpms);
    
    if (is<CharacterVector>(patches) || is<LogicalVector>(patches)) {
      CharacterVector patch_vec_pre = as<CharacterVector>(patches);
      
      if (preexisting) {
        CharacterVector mpm_patch_vec (mpm_count);
        for (int i = 0; i < mpm_count; i++) {
          CharacterVector label_portion (2);
          
          List chosen_mpm = as<List>(mpm_list(i));
          DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
          CharacterVector mpm_labels_pop = as<CharacterVector>(mpm_labels["pop"]);
          CharacterVector mpm_labels_patch = as<CharacterVector>(mpm_labels["patch"]);
          
          if (!CharacterVector::is_na(patch_vec_pre(i))) {
            IntegerVector found_indices = index_l3(mpm_labels_patch, patch_vec_pre(i));
            if (found_indices.length() == 0) {
              throw Rcpp::exception("Some values in vector patches do not exist in entered MPMs.", 
                false);
            }
            mpm_patch_vec(i) = patch_vec_pre(i);
            
            int key_labels_index = found_indices(0);
            label_portion(0) = mpm_labels_pop(key_labels_index);
            label_portion(1) = mpm_labels_patch(key_labels_index);
          } else {
            int chosen_patch_index {0};
            for (int j = 0; j < static_cast<int>(mpm_labels_patch.length()); j++) {
              if (mpm_labels_patch(j) == "0" || CharacterVector::is_na(mpm_labels_patch(j))) {
                chosen_patch_index = j;
                break;
              }
            }
            mpm_patch_vec(i) = mpm_labels_patch(chosen_patch_index);
            
            label_portion(0) = mpm_labels_pop(chosen_patch_index);
            label_portion(1) = mpm_labels_patch(chosen_patch_index);
            
          }
          labels_list_pre(i) = label_portion;
          //patch_list_pre(i) = patch_vec;
        }
        
        patch_vec = mpm_patch_vec;
        
      } else if (funcbased) {
        CharacterVector vrm_patch_vec (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          CharacterVector label_portion (2);
          label_portion(0) = "pop1";
          
          List chosen_vrm = as<List>(vrm_list(i));
          DataFrame vrm_patchframe = as<DataFrame>(chosen_vrm["patch_frame"]);
          CharacterVector vrm_patchframe_patches = as<CharacterVector>(vrm_patchframe["patches"]);
          
          if (!CharacterVector::is_na(patch_vec_pre(i))) {
            IntegerVector found_indices = index_l3(vrm_patchframe_patches, patch_vec_pre(i));
            if (found_indices.length() == 0) {
              throw Rcpp::exception("Some values in vector patches do not exist in entered MPMs.", 
                false);
            }
            vrm_patch_vec(i) = patch_vec_pre(i);
            
            label_portion(1) = patch_vec_pre(i);
          } else { 
            vrm_patch_vec(i) = vrm_patchframe_patches(0);
            
            label_portion(1) = vrm_patchframe_patches(0);
          }
          labels_list_pre(i) = label_portion;
          //patch_list_pre(i) = patch_vec;
        }
        patch_vec = vrm_patch_vec;
      }
      
      labels_list = labels_list_pre;
    }
    //patch_list = patch_list_pre;
  } else {
    List labels_list_pre (total_mpms);
    
    if (preexisting) {
      CharacterVector patch_vec_pre (mpm_count);
      
      for (int i = 0; i < mpm_count; i++) {
        CharacterVector label_portion (2);
        
        List chosen_mpm = as<List>(mpm_list(i));
        DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
        CharacterVector mpm_labels_pop = as<CharacterVector>(mpm_labels["pop"]);
        CharacterVector mpm_labels_patch = as<CharacterVector>(mpm_labels["patch"]);
        
        int chosen_patch_index {0};
        for (int j = 0; j < static_cast<int>(mpm_labels_patch.length()); j++) {
          if (mpm_labels_patch(j) == "0" || CharacterVector::is_na(mpm_labels_patch(j))) {
            chosen_patch_index = j;
            break;
          }
        }
        patch_vec_pre(i) = mpm_labels_patch(chosen_patch_index);
        
        label_portion(0) = mpm_labels_pop(chosen_patch_index);
        label_portion(1) = mpm_labels_patch(chosen_patch_index);
        labels_list_pre(i) = label_portion;
      }
      patch_vec = patch_vec_pre;
    } else if (funcbased) {
      CharacterVector patch_vec_pre (vrm_count);
      
      for (int i = 0; i < vrm_count; i++) {
        CharacterVector label_portion (2);
        label_portion(0) = "pop1";
        
        List chosen_vrm = as<List>(vrm_list(i));
        DataFrame vrm_patchframe = as<DataFrame>(chosen_vrm["patch_frame"]);
        CharacterVector vrm_patchframe_patches = as<CharacterVector>(vrm_patchframe["patches"]);
        
        patch_vec_pre(i) = vrm_patchframe_patches(0);
        
        label_portion(1) = vrm_patchframe_patches(0);
        labels_list_pre(i) = label_portion;
      }
      patch_vec = patch_vec_pre;
    }
    
    labels_list = labels_list_pre;
  }
  
  //Rcout << "cleanup3 K    ";
  
  // label construction
  {
    CharacterVector labels_pops (total_mpms);
    CharacterVector labels_patches (total_mpms);
    IntegerVector labels_mpms = seq(1, total_mpms);
    
    for (int i = 0; i < total_mpms; i++) {
      CharacterVector current_label = as<CharacterVector>(labels_list(i));
      labels_pops(i) = current_label(0);
      labels_patches(i) = current_label(1);
    }
    
    labels = DataFrame::create(_["mpm"] = labels_mpms, _["pop"] = labels_pops,
      _["patch"] = labels_patches);
  }
  
  //Rcout << "cleanup3 L    ";
  
  // years vector
  if (years.isNotNull()) {
    if (is<NumericVector>(years) || is<CharacterVector>(years)) {
      CharacterVector year_vec = as<CharacterVector>(years);
      
      if (preexisting) {
        List year_list_pre (mpm_count);
        IntegerVector total_years_vec_pre (mpm_count);
        
        for (int i = 0; i < mpm_count; i++) {
          List chosen_mpm = as<List>(mpm_list(i));
          DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
          CharacterVector mpm_labels_vars = mpm_labels.attr("names");
          IntegerVector mpm_labels_y2_var = index_l3(mpm_labels_vars, "year2");
          
          if (mpm_labels_y2_var.length() == 0) {
            throw Rcpp::exception("Some MPMs appear to be mean MPMs.", false);
          } else {
            CharacterVector mpm_labels_year2 = as<CharacterVector>(mpm_labels["year2"]);
            CharacterVector mly2_diffs = setdiff(year_vec, mpm_labels_year2);
            
            int found_total_years = static_cast<int>(mpm_labels_year2.length());
            total_years_vec_pre(i) = found_total_years;
            
            if (mly2_diffs.length() > 0) { 
              throw Rcpp::exception("Some entered values in years do not exist in some MPMs.",
                false);
            }
          }
          year_list_pre(i) = year_vec;
        }
        
        year_list = year_list_pre;
        total_years_vec = total_years_vec_pre;
        
      } else if (funcbased) {
        List year_list_pre (vrm_count);
        IntegerVector total_years_vec_pre (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          List chosen_vrm = as<List>(vrm_list(i));
          DataFrame vrm_yearframe = as<DataFrame>(chosen_vrm["year_frame"]);
          
          CharacterVector vrm_yearframe_years = as<CharacterVector>(vrm_yearframe["years"]);
          CharacterVector vyy2_diffs = setdiff(year_vec, vrm_yearframe_years);
          
          int found_total_years = static_cast<int>(vrm_yearframe_years.length());
          total_years_vec_pre(i) = found_total_years;
            
          if (vyy2_diffs.length() > 0) { 
            throw Rcpp::exception("Some entered values in years do not exist in some MPMs.",
              false);
          }
          
          year_list_pre(i) = year_vec;
        }
        
        year_list = year_list_pre;
        total_years_vec = total_years_vec_pre;
      }
    } else if (is<List>(years)) {
      List year_list_pre = as<List>(years);
      
      if (preexisting) {
        List mpm_year_list (mpm_count);
        IntegerVector total_years_vec_pre (mpm_count);
        
        for (int i = 0; i < mpm_count; i++) {
          if (!is<NumericVector>(year_list_pre(i)) && !is<CharacterVector>(year_list_pre(i))) {
            throw Rcpp::exception("Elements in list years must be numeric vectors.", false);
          }
          
          CharacterVector year_vec = as<CharacterVector>(year_list_pre(i));
          
          List chosen_mpm = as<List>(mpm_list(i));
          DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
          CharacterVector mpm_labels_vars = mpm_labels.attr("names");
          IntegerVector mpm_labels_y2_var = index_l3(mpm_labels_vars, "year2");
          
          if (mpm_labels_y2_var.length() == 0) {
            throw Rcpp::exception("Some MPMs appear to be mean MPMs.", false);
          } else {
            CharacterVector mpm_labels_year2 = as<CharacterVector>(mpm_labels["year2"]);
            
            if (!CharacterVector::is_na(year_vec(0))) {
              CharacterVector mly2_diffs = setdiff(year_vec, mpm_labels_year2);
              CharacterVector years_unique = sort_unique(mpm_labels_year2);
              
              mpm_year_list(i) = year_vec;
              if (mly2_diffs.length() > 0) { 
                throw Rcpp::exception("Some entered values in years do not exist in some MPMs.",
                  false);
              }
              
              int found_total_years = static_cast<int>(years_unique.length());
              total_years_vec_pre(i) = found_total_years;
            } else {
              CharacterVector years_unique = sort_unique(mpm_labels_year2);
              mpm_year_list(i) = years_unique;
              
              int found_total_years = static_cast<int>(years_unique.length());
              total_years_vec_pre(i) = found_total_years;
            }
          }
        }
        year_list = mpm_year_list;
        total_years_vec = total_years_vec_pre;
        
      } else if (funcbased) {
        List vrm_year_list (vrm_count);
        IntegerVector total_years_vec_pre (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          if (!is<NumericVector>(year_list_pre(i)) && !is<CharacterVector>(year_list_pre(i))) {
            throw Rcpp::exception("Elements in list years must be numeric vectors.",
              false);
          }
          
          CharacterVector year_vec = as<CharacterVector>(year_list_pre(i));
          
          List chosen_vrm = as<List>(vrm_list(i));
          DataFrame vrm_yearframe = as<DataFrame>(chosen_vrm["year_frame"]);
          
          CharacterVector vrm_yearframe_years = as<CharacterVector>(vrm_yearframe["years"]);
          
          int found_total_years = static_cast<int>(vrm_yearframe_years.length());
          total_years_vec_pre(i) = found_total_years;
            
          if (!CharacterVector::is_na(year_vec(0))) {
            CharacterVector vyy2_diffs = setdiff(year_vec, vrm_yearframe_years);
            
            vrm_year_list(i) = year_vec;
            if (vyy2_diffs.length() > 0) { 
              throw Rcpp::exception("Some entered values in years do not exist in some MPMs.",
                false);
            }
          } else {
            vrm_year_list(i) = vrm_yearframe_years;
          }
        }
        
        year_list = vrm_year_list;
        total_years_vec = total_years_vec_pre;
      }
    }
  } else {
    CharacterVector year_vec (1, NA_STRING);
    
    if (preexisting) {
      IntegerVector total_years_vec_pre (mpm_count);
      
      List year_list_pre (mpm_count);
      
      for (int i = 0; i < mpm_count; i++) {
        List chosen_mpm = as<List>(mpm_list(i));
        DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
        CharacterVector mpm_labels_vars = mpm_labels.attr("names");
        IntegerVector mpm_labels_y2_var = index_l3(mpm_labels_vars, "year2");
        
        if (mpm_labels_y2_var.length() == 0) {
          total_years_vec_pre(i) = 1;
          year_list_pre(i) = year_vec;
        } else {
          CharacterVector mpm_labels_year2 = as<CharacterVector>(mpm_labels["year2"]);
          
          CharacterVector mly2_unique = sort_unique(mpm_labels_year2);
          year_list_pre(i) = mly2_unique;
          
          int found_total_years = static_cast<int>(mly2_unique.length());
          total_years_vec_pre(i) = found_total_years;
        }
      }
      year_list = year_list_pre;
      total_years_vec = total_years_vec_pre;
      
    } else if (funcbased) {
      IntegerVector total_years_vec_pre (vrm_count);
      
      List year_list_pre (vrm_count);
      
      for (int i = 0; i < vrm_count; i++) {
        List chosen_vrm = as<List>(vrm_list(i));
        DataFrame vrm_yearframe = as<DataFrame>(chosen_vrm["year_frame"]);
        
        CharacterVector vrm_yearframe_years = as<CharacterVector>(vrm_yearframe["years"]);
        
        int found_total_years = static_cast<int>(vrm_yearframe_years.length());
        total_years_vec_pre(i) = found_total_years;
        
        year_list_pre(i) = vrm_yearframe_years;
      }
      
      year_list = year_list_pre;
      total_years_vec = total_years_vec_pre;
    }
  }
  
  //Rcout << "cleanup3 M    ";
  
  // tweights list
  if (tweights.isNotNull()) {
    int assumed_mpms = mpm_count;
    if (vrm_count > mpm_count) assumed_mpms = vrm_count;
    
    IntegerVector tweights_type_vec_temp (assumed_mpms);
    
    if (is<List>(tweights)) {
      tweights_list = as<List>(tweights);
      tweights_count = static_cast<int>(tweights_list.length());
      
      if (tweights_count != mpm_count && tweights_count != vrm_count) {
        throw Rcpp::exception("Argument tweights must have as many elements as MPMs.",
          false);
      }
      
      for (int i = 0; i < tweights_count; i++) {
        if (Rf_isMatrix(tweights_list(i))) {
          NumericMatrix chosen_matrix = as<NumericMatrix>(tweights_list(i));
          int mat_rows = chosen_matrix.nrow();
          int mat_cols = chosen_matrix.ncol();
          
          if (mat_rows != mat_cols) {
            throw Rcpp::exception("Matrices in argument tweights must be square.", false);
          }
          
          if (mat_rows != total_years_vec(i)) {
            throw Rcpp::exception("Matrices in argument tweights must account for all years.",
              false);
          }
          tweights_type_vec_temp(i) = 2;
          
        } else if (is<NumericVector>(tweights_list(i))) {
          NumericVector chosen_vector = as<NumericVector>(tweights_list(i));
          
          if (static_cast<int>(chosen_vector.length()) != total_years_vec(i)) {
            throw Rcpp::exception("Vectors in argument tweights must account for all years.",
              false);
          }
          tweights_type_vec_temp(i) = 1;
          
        } else {
          AdaptUtils::pop_error2("tweights", "a list of numeric vectors or matrices", "", 1);
        }
        tweights_type_vec = tweights_type_vec_temp;
      }
    } else {
      AdaptUtils::pop_error2("tweights", "a list of numeric vectors or matrices", "", 1);
    }
  } else {
    if (preexisting) {
      List tweights_list_pre (mpm_count);
      
      for (int i = 0; i < mpm_count; i++) {
        tweights_list_pre(i) = R_NilValue;
      }
      tweights_list = tweights_list_pre;
      tweights_count = mpm_count;
      
      IntegerVector tweights_type_vec_temp (mpm_count);
      tweights_type_vec = tweights_type_vec_temp;
    } else if (funcbased) {
      List tweights_list_pre (vrm_count);
      
      for (int i = 0; i < vrm_count; i++) {
        tweights_list_pre(i) = R_NilValue;
      }
      tweights_list = tweights_list_pre;
      tweights_count = vrm_count;
      
      IntegerVector tweights_type_vec_temp (vrm_count);
      tweights_type_vec = tweights_type_vec_temp;
    }
  }
  
  //Rcout << "cleanup3 N    ";
  
  // density list
  if (density.isNotNull()) {
    if (is<List>(density)) {
      density_list = as<List>(density);
      density_count = static_cast<int>(density_list.length());
      if (density_count != vrm_count && density_count != mpm_count) {
        throw Rcpp::exception("Argument density must be list of same length as number of MPMs.",
          false);
      }
      IntegerVector dens_yn_vec_temp (density_count);
      List dens_index_list_pre (density_count);
      //List hstages_list_fb_pre (density_count);
      
      for (int i = 0; i < density_count; i++) {
        if (is<DataFrame>(density_list(i))) {
          DataFrame chosen_density = as<DataFrame>(density_list(i));
          
          if (chosen_density.hasAttribute("class")) {
            CharacterVector chosen_density_class = chosen_density.attr("class");
            bool found_lefkoDens {false};
            
            for (int j = 0; j < static_cast<int>(chosen_density_class.length()); j++) {
              if (chosen_density_class(j) == "lefkoDens") found_lefkoDens = true;
            }
            if (!found_lefkoDens) {
              AdaptUtils::pop_error2("density", "a list of lefkoDens objects and NULL values", "", 1);
            }
            
            CharacterVector dl_stage1 = as<CharacterVector>(chosen_density["stage1"]);
            IntegerVector dl_age2 = as<IntegerVector>(chosen_density["age2"]);
            
            if (format_vec(i) < 3) {
              if (is<LogicalVector>(chosen_density["stage1"])) {
                throw Rcpp::exception("Argument density requires real stage1 entries other than NA if MPMs are historical.", false);
              }
              for (int j = 0; j < static_cast<int>(dl_stage1.length()); j++) {
                if (CharacterVector::is_na(dl_stage1(j))) {
                  throw Rcpp::exception("Argument density requires real stage1 entries other than NA if MPMs are historical.", false);
                }
              }
            } else if (format_vec(i) > 3) {
              if (is<LogicalVector>(chosen_density["age2"])) {
                throw Rcpp::exception("Argument density requires real stage1 entries other than NA if MPMs are historical.", false);
              }
              for (int j = 0; j < static_cast<int>(dl_age2.length()); j++) {
                if (IntegerVector::is_na(dl_age2(j)) || LogicalVector::is_na(dl_age2(j))) {
                  throw Rcpp::exception("Argument density requires real age2 entries other than NA if MPMs are age-by-stage.", false);
                }
              }
            }
            
            dens_yn_vec_temp(i) = 1;
          } else {
            AdaptUtils::pop_error2("density", "a list of lefkoDens objects and NULL values", "", 1);
          }
          
          Rcpp::StringVector di_stage3 = as<StringVector>(chosen_density["stage3"]);
          Rcpp::StringVector di_stage2 = as<StringVector>(chosen_density["stage2"]);
          Rcpp::StringVector di_stage1 = as<StringVector>(chosen_density["stage1"]);
          int di_size = di_stage3.length();
          
          if (format_vec(i) < 3) {
            DataFrame hstages = as<DataFrame>(hstages_list(i));
            
            StringVector stage3 = as<StringVector>(hstages["stage_2"]);
            StringVector stage2r = as<StringVector>(hstages["stage_1"]);
            StringVector stage2c = as<StringVector>(hstages["stage_2"]);
            StringVector stage1 = as<StringVector>(hstages["stage_1"]);
            int hst_size = stage3.length();
            
            arma::uvec hst_3(hst_size, fill::zeros);
            arma::uvec hst_2r(hst_size, fill::zeros);
            arma::uvec hst_2c(hst_size, fill::zeros);
            arma::uvec hst_1(hst_size, fill::zeros);
            
            arma::uvec di_stage32_id(di_size, fill::zeros);
            arma::uvec di_stage21_id(di_size, fill::zeros);
            arma::uvec di_index(di_size, fill::zeros);
            
            for (int j = 0; j < di_size; j++) { // Loop through each density_input line
              for (int k = 0; k < hst_size; k++) {
                if (di_stage3(j) == stage3(k)) {
                  hst_3(k) = 1;
                } else {
                  hst_3(k) = 0;
                }
              }
              
              for (int k = 0; k < hst_size; k++) {
                if (di_stage2(j) == stage2r(k)) {
                  hst_2r(k) = 1;
                } else {
                  hst_2r(k) = 0;
                }
              }
              
              for (int k = 0; k < hst_size; k++) {
                if (di_stage2(j) == stage2c(k)) {
                  hst_2c(k) = 1;
                } else {
                  hst_2c(k) = 0;
                }
              }
              
              for (int k = 0; k < hst_size; k++) {
                if (di_stage1(j) == stage1(k)) {
                  hst_1(k) = 1;
                } else {
                  hst_1(k) = 0;
                }
              }
              
              arma::uvec find_hst3 = find(hst_3);
              arma::uvec find_hst2r = find(hst_2r);
              arma::uvec find_hst2c = find(hst_2c);
              arma::uvec find_hst1 = find(hst_1);
              
              arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
              arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
              
              if (static_cast<int>(pop_32.n_elem) == 0 || static_cast<int>(pop_21.n_elem) == 0) {
                throw Rcpp::exception("Some stages in argument density could not be found.", 
                  false);
              }
              di_stage32_id(j) = pop_32(0);
              di_stage21_id(j) = pop_21(0);
              di_index(j) = pop_32(0) + (pop_21(0) * hst_size);
              
              hst_3.zeros();
              hst_2r.zeros();
              hst_2c.zeros();
              hst_1.zeros();
            }
            
            List dens_index_list_mpm = Rcpp::List::create(_["index32"] = di_stage32_id,
              _["index21"] = di_stage21_id, _["index321"] = di_index);
            dens_index_list_pre(i) = dens_index_list_mpm;
            
          } else if (format_vec(i) == 4 ) { 
            IntegerVector di_age2 = as<IntegerVector>(chosen_density["age2"]);
            DataFrame agestages = as<DataFrame>(agestages_list(i));
            
            StringVector stage3 = as<StringVector>(agestages["stage"]);
            StringVector stage2 = as<StringVector>(agestages["stage"]);
            IntegerVector age2 = as<IntegerVector>(agestages["age"]);
            int agst_size = stage3.length();
            
            arma::uvec agst_s3(agst_size, fill::zeros);
            arma::uvec agst_a3(agst_size, fill::zeros);
            arma::uvec agst_s2(agst_size, fill::zeros);
            arma::uvec agst_a2(agst_size, fill::zeros);
            
            arma::uvec di_s3a3_id(di_size, fill::zeros);
            arma::uvec di_s2a2_id(di_size, fill::zeros);
            arma::uvec di_index(di_size, fill::zeros);
            
            for (int j = 0; j < di_size; j++) { // Loop through each density_input line
              for (int k = 0; k < agst_size; k++) {
                if (di_stage3(j) == stage3(k)) {
                  agst_s3(k) = 1;
                } else {
                  agst_s3(k) = 0;
                }
              }
              
              for (int k = 0; k < agst_size; k++) {
                if (di_stage2(j) == stage2(k)) {
                  agst_s2(k) = 1;
                } else {
                  agst_s2(k) = 0;
                }
              }
              
              for (int k = 0; k < agst_size; k++) {
                if (di_age2(j) < finalage_vec(i)) {
                  if (di_age2(j) == age2(k)) {
                    agst_a2(k) = 1;
                    
                    for (int l = 0; l < agst_size; l++) {
                      if ((di_age2(j) + 1) == age2(l)) {
                        agst_a3(l) = 1;
                      } else {
                        agst_a3(l) = 0;
                      }
                    }
                  } else {
                    agst_a2(k) = 0;
                  }
                } else {
                  if (di_age2(j) == age2(k)) {
                    agst_a2(k) = 1;
                    agst_a3(k) = 1;
                  } else {
                    agst_a2(k) = 0;
                    agst_a3(k) = 0;
                  }
                }
              }
              
              arma::uvec find_agst_s3 = find(agst_s3);
              arma::uvec find_agst_s2 = find(agst_s2);
              arma::uvec find_agst_a3 = find(agst_a3);
              arma::uvec find_agst_a2 = find(agst_a2);
              
              arma::uvec pop_32 = intersect(find_agst_s3, find_agst_a3);
              arma::uvec pop_21 = intersect(find_agst_s2, find_agst_a2);
              
              if (static_cast<int>(pop_32.n_elem) == 0 || static_cast<int>(pop_21.n_elem) == 0) {
                throw Rcpp::exception("Some age-stages in argument density could not be found.", 
                  false);
              }
              di_s3a3_id(j) = pop_32(0);
              di_s2a2_id(j) = pop_21(0);
              di_index(j) = pop_32(0) + (pop_21(0) * agst_size);
              
              agst_s3.zeros();
              agst_s2.zeros();
              agst_a3.zeros();
              agst_a2.zeros();
            }
            
            List dens_index_list_mpm = Rcpp::List::create(_["index32"] = di_s3a3_id,
              _["index21"] = di_s2a2_id, _["index321"] = di_index);
            dens_index_list_pre(i) = dens_index_list_mpm;
            
          } else {
            DataFrame stageframe;
            
            if (preexisting) {
              List chosen_mpm = as<List>(mpm_list(i));
              stageframe = as<DataFrame>(chosen_mpm["ahstages"]);
            } else {
              stageframe = as<DataFrame>(stageframe_list(i));
            }
            
            StringVector stage3 = as<StringVector>(stageframe["stage"]);
            StringVector stage2 = as<StringVector>(stageframe["stage"]);
            int ahst_size = stage3.length();
            if (funcbased) ahst_size--;
            
            arma::uvec ahst_3(ahst_size, fill::zeros);
            arma::uvec ahst_2(ahst_size, fill::zeros);
            
            arma::uvec di_stage32_id(di_size, fill::zeros);
            arma::uvec di_stage21_id(di_size, fill::zeros);
            arma::uvec di_index(di_size, fill::zeros);
            
            for (int j = 0; j < di_size; j++) { // Loop through each density_input
              for (int k = 0; k < ahst_size; k++) {
                if (di_stage3(j) == stage3(k)) {
                  ahst_3(k) = 1;
                } else {
                  ahst_3(k) = 0;
                }
              }
              
              for (int k = 0; k < ahst_size; k++) {
                if (di_stage2(j) == stage2(k)) {
                  ahst_2(k) = 1;
                } else {
                  ahst_2(k) = 0;
                }
              }
              
              arma::uvec find_ahst3 = find(ahst_3);
              arma::uvec find_ahst2 = find(ahst_2);
              di_stage32_id(j) = find_ahst3(0);
              di_stage21_id(j) = find_ahst2(0);
              di_index(j) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
              
              ahst_3.zeros();
              ahst_2.zeros();
            }
            
            List dens_index_list_mpm = Rcpp::List::create(_["index3"] = di_stage32_id,
              _["index2"] = di_stage21_id, _["index321"] = di_index);
            dens_index_list_pre(i) = dens_index_list_mpm;
          }
          
          
          arma::uvec dyn_style = as<arma::uvec>(chosen_density["style"]);
          arma::vec dyn_alpha = as<arma::vec>(chosen_density["alpha"]);
          arma::vec dyn_beta = as<arma::vec>(chosen_density["beta"]);
          
          for (int j = 0; j < static_cast<int>(dyn_style.n_elem); j++) {
            if (dyn_style(j) < 1 || dyn_style(j) > 4) {
              String eat_my_shorts = "Some density inputs are stated as yielding density ";
              eat_my_shorts += "dependence but not in an accepted style.";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
            
            if (dyn_style(j) == 1) {
              if (dyn_beta(j) > exp_tol) {
                Rf_warningcall(R_NilValue,
                  "Beta used in Ricker function may be too high. Results may be unpredictable.");
                
              } else if (dyn_beta(j) < (-1.0 * exp_tol)) {
                Rf_warningcall(R_NilValue,
                  "Beta used in Ricker function may be too high. Results may be unpredictable.");
                
              }
              
            } else if (dyn_style(j) == 3) {
              double summed_stuff = dyn_alpha(j) + dyn_beta(j);
              
              if (summed_stuff > exp_tol) {
                Rf_warningcall(R_NilValue,
                  "Alpha and beta used in Usher function may be too high.");
                
              } else if (summed_stuff < (-1.0 * exp_tol)) {
                Rf_warningcall(R_NilValue,
                  "Alpha and beta used in Usher function may be too high.");
              }
            }
          }
          
        } else if (density_list(i) != R_NilValue) { 
          AdaptUtils::pop_error2("density", "a list of lefkoDens objects and NULL values", "", 1);
        }
      }
      dens_index_list = dens_index_list_pre;
      dens_yn_vec = dens_yn_vec_temp;
      
    } else {
      AdaptUtils::pop_error2("density", "a list of lefkoDens objects and NULL values", "", 1);
    }
  } else {
    if (preexisting) {
      List density_list_pre (mpm_count);
      
      for (int i = 0; i < mpm_count; i++) {
        density_list_pre(i) = R_NilValue;
      }
      density_list = density_list_pre;
      density_count = mpm_count;
      
      IntegerVector dens_yn_vec_temp (mpm_count);
      dens_yn_vec = dens_yn_vec_temp;
      
    } else if (funcbased) {
      List density_list_pre (vrm_count);
      
      for (int i = 0; i < vrm_count; i++) {
        density_list_pre(i) = R_NilValue;
      }
      density_list = density_list_pre;
      density_count = vrm_count;
      
      IntegerVector dens_yn_vec_temp (vrm_count);
      dens_yn_vec = dens_yn_vec_temp;
    }
  }
  
  //Rcout << "cleanup3 O    ";
  
  // entry time vector
  if (entry_time.isNotNull()) {
    if (is<NumericVector>(entry_time) || is<IntegerVector>(entry_time)) {
      entry_time_vec = as<IntegerVector>(entry_time);
      entry_time_count = static_cast<int>(entry_time_vec.length());
      
      if (entry_time_count != mpm_count && entry_time_count != vrm_count) {
        throw Rcpp::exception("Argument entry_time must be integer vector with same length as number of MPMs.", false);
      }
      
      int etv_sum = sum(entry_time_vec);
      if (etv_sum > 0) entry_time_vec_use = true;
      
    } else {
      throw Rcpp::exception("Argument entry_time must be integer vector with same length as number of MPMs.", false);
    }
  } else {
    if (preexisting) {
      IntegerVector entry_time_vec_temp (mpm_count);
      entry_time_vec = entry_time_vec_temp;
    } else {
      IntegerVector entry_time_vec_temp (vrm_count);
      entry_time_vec = entry_time_vec_temp;
    }
  }
  
  //Rcout << "cleanup3 P    ";
  
  // vrms-only arguments
  if (!funcbased) {
    if (density_vr.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "use argument density_vr", "", 26);
    }
    if (sp_density.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "use argument sp_density", "", 26);
    }
    if (ind_terms.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "use argument ind_terms", "", 26);
    }
    if (dev_terms.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "use argument ind_terms", "", 26);
    }
    if (fb_sparse.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "use argument fb_sparse", "", 26);
    }
    
  } else {
    // density_vr list
    IntegerVector dvr_yn_count (vrm_count);
    
    if (density_vr.isNotNull()) {
      if (is<List>(density_vr)) {
        density_vr_list = as<List>(density_vr);
        density_vr_count = static_cast<int>(density_vr_list.length());
        
        if (density_vr_count != vrm_count) {
          throw Rcpp::exception("Arguments vrms and density_vr must be lists of same length.",
            false);
        }
        
        for (int i = 0; i < density_vr_count; i++) {
          if (is<DataFrame>(density_vr_list(i))) {
            DataFrame chosen_density_vr = as<DataFrame>(density_vr_list(i));
            dvr_yn_count(i) = 1;
            
            if (chosen_density_vr.hasAttribute("class")) {
              CharacterVector chosen_density_vr_class = chosen_density_vr.attr("class");
              bool found_lefkoDensVR {false};
              
              for (int j = 0; j < static_cast<int>(chosen_density_vr_class.length()); j++) {
                if (chosen_density_vr_class(j) == "lefkoDensVR") found_lefkoDensVR = true;
              }
              if (!found_lefkoDensVR) {
                AdaptUtils::pop_error2("density_vr", "a list of lefkoDensVR objects and NULL values", "", 1);
              }
            } else {
              AdaptUtils::pop_error2("density_vr", "a list of lefkoDensVR objects and NULL values", "", 1);
            }
          } else if (density_vr_list(i) != R_NilValue) { 
            AdaptUtils::pop_error2("density_vr", "a list of lefkoDensVR objects and NULL values", "", 1);
          }
        }
      } else {
        AdaptUtils::pop_error2("density_vr", "a list of lefkoDensVR objects and NULL values", "", 1);
      }
    } else {
      if (preexisting) {
        List density_vr_list_pre (mpm_count);
        
        for (int i = 0; i < mpm_count; i++) {
          density_vr_list_pre(i) = R_NilValue;
        }
        density_vr_list = density_vr_list_pre;
        density_vr_count = mpm_count;
        
      } else if (funcbased) {
        List density_vr_list_pre (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          density_vr_list_pre(i) = R_NilValue;
        }
        density_vr_list = density_vr_list_pre;
        density_vr_count = vrm_count;
      }
    }
    dens_vr_yn_vec = dvr_yn_count;
    
    // sp_density list
    IntegerVector spd_num_count (vrm_count);
    
    if (sp_density.isNotNull()) {
      if (is<NumericVector>(sp_density)) {
        NumericVector sp_density_temp = as<NumericVector>(sp_density);
        int sp_density_count = static_cast<int>(sp_density_temp.length());
        
        if (sp_density_count != vrm_count) {
          throw Rcpp::exception("Length of vector sp_density must equal number of vrm_input objects.", false);
        }
        
        List sp_density_initial_list (vrm_count);
        for (int i = 0; i < vrm_count; i++) {
          NumericVector single_value {static_cast<double>(sp_density_temp(i))};
          sp_density_initial_list(i) = single_value;
          spd_num_count(i) = 1;
        }
      } else if (is<List>(sp_density)) {
        sp_density_list = as<List>(sp_density);
        int sp_density_count = static_cast<int>(sp_density_list.length());
        
        if (sp_density_count != vrm_count) {
          throw Rcpp::exception("Length of list sp_density must equal number of vrm_input objects.", false);
        }
        
        // Check list elements
        for (int i = 0; i < vrm_count; i++) {
          if (is<NumericVector>(sp_density_list(i))) {
            NumericVector sp_density_current_vec = as<NumericVector>(sp_density_list(i));
            int sdcv_length = static_cast<int>(sp_density_current_vec.length());
            
            spd_num_count(i) = sdcv_length;
          } else if (is<LogicalVector>(sp_density_list(i))) {
            LogicalVector sp_density_current_vec = as<LogicalVector>(sp_density_list(i));
            int sdcv_length = static_cast<int>(sp_density_current_vec.length());
            
            if (sdcv_length != 1 || !LogicalVector::is_na(sp_density_current_vec(0))) {
              throw Rcpp::exception("Some list elements in argument sp_density are not valid format.", false);
            }
          } else if (sp_density_list(i) != R_NilValue) { 
            throw Rcpp::exception("Some list elements in argument sp_density are not valid format.", false);
          }
        }
        
      } else {
        throw Rcpp::exception("Input in argument sp_density is not valid.",
          false);
      }
    } else {
      List sp_density_list_temp (vrm_count);
      
      for (int i = 0; i < vrm_count; i++) {
        NumericVector sp_temp {0.0};
        sp_density_list_temp(i) = sp_temp;
      }
    }
    sp_density_num_vec = spd_num_count;
    
    // ind_terms list
    IntegerVector inda_num_count (vrm_count);
    IntegerVector indb_num_count (vrm_count);
    IntegerVector indc_num_count (vrm_count);
    IntegerVector inda_cat_count (vrm_count);
    IntegerVector indb_cat_count (vrm_count);
    IntegerVector indc_cat_count (vrm_count);
    
    if (ind_terms.isNotNull()) {
      if (is<DataFrame>(ind_terms)) {
        DataFrame ind_terms_df = as<DataFrame>(ind_terms);
        int idt_df_size = static_cast<int>(ind_terms_df.size());
        if (idt_df_size != 3) {
          throw Rcpp::exception("Data frame ind_terms should have 3 columns only.",
            false);
        }
        
        if ((!is<NumericVector>(ind_terms_df(0)) && !is<CharacterVector>(ind_terms_df(0))) ||
          (!is<NumericVector>(ind_terms_df(1)) && !is<CharacterVector>(ind_terms_df(1))) ||
          (!is<NumericVector>(ind_terms_df(2)) && !is<CharacterVector>(ind_terms_df(2)))) {
            throw Rcpp::exception("Data frame in argument ind_terms must be either numeric or character values.", false);
        }
        
        List idt_num_pre (vrm_count);
        List idt_cat_pre (vrm_count);
        
        int idt_df_nrows = static_cast<int>(ind_terms_df.nrows());
        for (int i = 0; i < vrm_count; i++) {
          List current_idt_cat (3);
          List current_idt_num (3);
          
          for (int j = 0; j < 3; j++) {
            if (is<CharacterVector>(ind_terms_df(j))) {
              CharacterVector current_idt_cat_col = as<CharacterVector>(ind_terms_df(j));
              NumericVector current_idt_num_col (idt_df_nrows);
              
              current_idt_cat(j) = current_idt_cat_col;
              current_idt_num(j) = current_idt_num_col;
              
              if (j == 0) {
                inda_num_count(i) = 0;
                inda_cat_count(i) = idt_df_nrows;
              } else if (j == 1) {
                indb_num_count(i) = 0;
                indb_cat_count(i) = idt_df_nrows;
              } else {
                indc_num_count(i) = 0;
                indc_cat_count(i) = idt_df_nrows;
              }
              
            } else if (is<NumericVector>(ind_terms_df(j))) {
              CharacterVector single_none {"none"};
              CharacterVector current_idt_cat_col = rep(single_none, idt_df_nrows);
              NumericVector current_idt_num_col = as<NumericVector>(ind_terms_df(j));
              
              current_idt_cat(j) = current_idt_cat_col;
              current_idt_num(j) = current_idt_num_col;
              
              if (j == 0) {
                inda_num_count(i) = idt_df_nrows;
                inda_cat_count(i) = 0;
              } else if (j == 1) {
                indb_num_count(i) = idt_df_nrows;
                indb_cat_count(i) = 0;
              } else {
                indc_num_count(i) = idt_df_nrows;
                indc_cat_count(i) = 0;
              }
            }
          }
          
          idt_num_pre(i) = current_idt_num;
          idt_cat_pre(i) = current_idt_cat;
        }
        
        ind_terms_num_list = idt_num_pre;
        ind_terms_cat_list = idt_cat_pre;
        
      } else if (is<List>(ind_terms)) {
        List ind_terms_list = as<List>(ind_terms);
        int idt_list_length = static_cast<int>(ind_terms_list.length());
        
        if (idt_list_length != vrm_count) {
          throw Rcpp::exception("List in argument ind_terms must have as many elements as vrm_input objects.", false);
        }
        
        List idt_num_pre (vrm_count);
        List idt_cat_pre (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          DataFrame ind_terms_df = as<DataFrame>(ind_terms_list(i));
          int idt_df_size = static_cast<int>(ind_terms_df.size());
          if (idt_df_size != 3) {
            throw Rcpp::exception("All data frames in argument ind_terms must have 3 columns.",
              false);
          }
          
          if ((!is<NumericVector>(ind_terms_df(0)) && !is<CharacterVector>(ind_terms_df(0))) ||
          (!is<NumericVector>(ind_terms_df(1)) && !is<CharacterVector>(ind_terms_df(1))) ||
          (!is<NumericVector>(ind_terms_df(2)) && !is<CharacterVector>(ind_terms_df(2)))) {
              throw Rcpp::exception("All data frames in argument ind_terms must be numeric or character values.", false);
          }
          
          int idt_df_nrows = static_cast<int>(ind_terms_df.nrows());
          
          List current_idt_cat (3);
          List current_idt_num (3);
          
          for (int j = 0; j < 3; j++) {
            if (is<CharacterVector>(ind_terms_df(j))) {
              CharacterVector current_idt_cat_col = as<CharacterVector>(ind_terms_df(j));
              NumericVector current_idt_num_col (idt_df_nrows);
              
              current_idt_cat(j) = current_idt_cat_col;
              current_idt_num(j) = current_idt_num_col;
              
              if (j == 0) {
                inda_num_count(i) = 0;
                inda_cat_count(i) = idt_df_nrows;
              } else if (j == 1) {
                indb_num_count(i) = 0;
                indb_cat_count(i) = idt_df_nrows;
              } else {
                indc_num_count(i) = 0;
                indc_cat_count(i) = idt_df_nrows;
              }
              
            } else if (is<NumericVector>(ind_terms_df(j))) {
              CharacterVector single_none {"none"};
              CharacterVector current_idt_cat_col = rep(single_none, idt_df_nrows);
              NumericVector current_idt_num_col = as<NumericVector>(ind_terms_df(j));
              
              current_idt_cat(j) = current_idt_cat_col;
              current_idt_num(j) = current_idt_num_col;
              
              if (j == 0) {
                inda_num_count(i) = idt_df_nrows;
                inda_cat_count(i) = 0;
              } else if (j == 1) {
                indb_num_count(i) = idt_df_nrows;
                indb_cat_count(i) = 0;
              } else {
                indc_num_count(i) = idt_df_nrows;
                indc_cat_count(i) = 0;
              }
            }
          }
          
          idt_num_pre(i) = current_idt_num;
          idt_cat_pre(i) = current_idt_cat;
        }
        
        ind_terms_num_list = idt_num_pre;
        ind_terms_cat_list = idt_cat_pre;
        
      } else {
        throw Rcpp::exception("Input in argument ind_terms is not valid.",
          false);
      }
    } else {
      List ind_terms_num_list_pre (vrm_count);
      List ind_terms_cat_list_pre (vrm_count);
      
      NumericVector region_A = NumericVector::create(0.);
      CharacterVector region_B {"none"};
      
      DataFrame region_A_df = DataFrame::create(_["A"] = region_A,
        _["B"] = clone(region_A), _["C"] = clone(region_A));
      DataFrame region_B_df = DataFrame::create(_["A"] = region_B,
        _["B"] = region_B, _["C"] = region_B);
      
      for (int i = 0; i < vrm_count; i++) {
        ind_terms_num_list_pre(i) = region_A_df;
        ind_terms_cat_list_pre(i) = region_B_df;
      }
      
      ind_terms_num_list = ind_terms_num_list_pre;
      ind_terms_cat_list = ind_terms_cat_list_pre;
    }
    inda_terms_num_vec = inda_num_count;
    indb_terms_num_vec = indb_num_count;
    indc_terms_num_vec = indc_num_count;
    inda_terms_cat_vec = inda_cat_count;
    indb_terms_cat_vec = indb_cat_count;
    indc_terms_cat_vec = indc_cat_count;
    
    // dev_terms list
    IntegerVector dev_num_count (vrm_count);
    
    if (dev_terms.isNotNull()) {
      if (is<List>(dev_terms)) {
        dev_terms_list = as<List>(dev_terms);
        int dvt_list_length = static_cast<int>(dev_terms_list.length());
        
        if (dvt_list_length != vrm_count) {
          throw Rcpp::exception("List in argument dev_terms must have as many elements as vrm_input objects.", false);
        }
        
        for (int i = 0; i < vrm_count; i++) {
          if (!is<DataFrame>(dev_terms_list(i))) {
            throw Rcpp::exception("List in argument dev_terms must be composed of data frames.",
              false);
          }
          
          DataFrame dev_terms_current_df = as<DataFrame>(dev_terms_list(i));
          int dvtc_df_size = static_cast<int>(dev_terms_current_df.size());
          int dvtc_df_nrows = static_cast<int>(dev_terms_current_df.nrows());
          
          if (dvtc_df_size != 14) {
            throw Rcpp::exception("Data frames in argument dev_terms must have 14 columns.",
              false);
          }
          
          for (int j = 0; j < 14; j++) {
            if (!is<NumericVector>(dev_terms_current_df(j))) {
              throw Rcpp::exception("Data frames in argument dev_terms must be composed of numeric variables.", false);
            }
          }
          
          dev_num_count(i) = dvtc_df_nrows;
        }
      } else {
        AdaptUtils::pop_error2("dev_terms", "a list of data frames", "", 1);
      }
    }
    dev_terms_num_vec = dev_num_count;
    
    // fb_sparse
    if (fb_sparse.isNotNull()) {
      if (is<LogicalVector>(fb_sparse)) {
        sparse_vec = as<LogicalVector>(fb_sparse);
        sparse_vec_count = static_cast<int>(sparse_vec.length());
        
        if (sparse_vec_count != vrm_count) {
          throw Rcpp::exception("Argument fb_sparse must be a logical vector of same length as list vrms.", false);
        }
        
        for (int i = 0; i < sparse_vec_count; i++) {
          if (LogicalVector::is_na(sparse_vec(i))) {
            throw Rcpp::exception("No NA values are allowed in argument fb_sparse.",
              false);
          }
        }
      }
    } else {
      LogicalVector fb_sparse_temp (vrm_count);
      sparse_vec = fb_sparse_temp;
    }
  } // End of vrm-only section
  
  //Rcout << "cleanup3 Q    ";
  
  // equivalence interpretation
  if (equivalence.isNotNull()) {
    if (is<NumericVector>(equivalence)) {
      equivalence_vec = as<NumericVector>(equivalence);
      
      int trial_count = mpm_count;
      if (vrm_count > mpm_count) trial_count = vrm_count;
      
      equivalence_count = static_cast<int>(equivalence_vec.length());
      
      for (int i = 0; i < equivalence_count; i++) {
        if (equivalence_vec(i) < 0.0) {
          AdaptUtils::pop_error2("equivalence", "", "", 30);
        } else if (NumericVector::is_na(equivalence_vec(i))) {
          throw Rcpp::exception("No NA values are allowed in argument equivalence.",
            false);
        }
      }
      
    } else if (is<List>(equivalence)) {
      List equivalence_list_temp = as<List>(equivalence);
      stages_not_equal = true;
      
      int trial_count = mpm_count;
      if (vrm_count > mpm_count) trial_count = vrm_count;
      
      equivalence_count = static_cast<int>(equivalence_list_temp.length());
      if (equivalence_count != trial_count) {
        throw Rcpp::exception("There must be as many elements in argument equivalence as MPMs.",
          false);
      }
      
      List equivalence_list_pre (equivalence_count);
      for (int i = 0; i < equivalence_count; i++) {
        if (is<NumericVector>(equivalence_list_temp(i))) {
          NumericVector trial_equivalence = as<NumericVector>(equivalence_list(i));
          int trial_eq_length = static_cast<int>(trial_equivalence.length());
          
          if (trial_eq_length != matrowcounts(i)) {
            throw Rcpp::exception("Numeric vectors in argument equivalence must account for all MPM rows.", false);
          }
          
          NumericVector equivalence_list_vec = as<NumericVector>(equivalence_list_temp(i));
          
          for (int j = 0; j < static_cast<int>(equivalence_list_vec.length()); j++) {
            if (equivalence_list_vec(j) < 0.0) {
              AdaptUtils::pop_error2("equivalence", "", "", 30);
            } else if (NumericVector::is_na(equivalence_list_vec(j))) {
              throw Rcpp::exception("No NA values are allowed in argument equivalence.",
                false);
            }
          }
          equivalence_list_pre(i) = equivalence_list_vec;
          
        } else if (is<DataFrame>(equivalence_list_temp(i))) {
          DataFrame eq_list_df = as<DataFrame>(equivalence_list_temp(i));
          if (!eq_list_df.hasAttribute("class")) {
            throw Rcpp::exception("Argument equivalence should include data frames of class adaptEq, or numeric vectors.", false);
          }
          CharacterVector eq_list_df_class = eq_list_df.attr("class");
          bool found_adaptEq {false};
          for (int j = 0; j < static_cast<int>(eq_list_df_class.length()); j++) {
            if (eq_list_df_class(j) == "adaptEq") found_adaptEq = true;
          }
          if (!found_adaptEq) {
            throw Rcpp::exception("Argument equivalence should include data frames of class adaptEq, or numeric vectors.", false);
          }
          
          IntegerVector eq_s2 = as<IntegerVector>(eq_list_df["stage_id_2"]);
          IntegerVector eq_s1 = as<IntegerVector>(eq_list_df["stage_id_1"]);
          IntegerVector eq_a2 = as<IntegerVector>(eq_list_df["age2"]);
          IntegerVector eq_rn = clone(as<IntegerVector>(eq_list_df["row_num"]));
          NumericVector eq_val = as<NumericVector>(eq_list_df["value"]);
          
          eq_rn = eq_rn - 1;
          
          if (format_vec(i) < 3) {
            if (IntegerVector::is_na(eq_s1(0))) {
              throw Rcpp::exception("Enter stage pairs in adaptEq objects used for historical MPMs.", 
                false);
            }
            if (IntegerVector::is_na(eq_s2(0))) {
              throw Rcpp::exception("Entries in column stage2 of adaptEq objects cannot be empty except in Leslie MPMs.", false);
            }
          } else if (format_vec(i) > 3) {
            if (IntegerVector::is_na(eq_a2(0))) {
              throw Rcpp::exception("Enter ages in adaptEq objects used for age-by-stage MPMs.",
                false);
            }
            if (format_vec(i) == 4) {
              if (IntegerVector::is_na(eq_s2(0))) {
                throw Rcpp::exception("Entries in column stage2 of adaptEq objects cannot be empty except in Leslie MPMs.", false);
              }
            }
          } else {
            if (IntegerVector::is_na(eq_s2(0))) {
              throw Rcpp::exception("Entries in column stage2 of adaptEq objects cannot be empty except in Leslie MPMs.", false);
            }
          }
          
          if (max(eq_rn) > matrowcounts(i)) {
            throw Rcpp::exception("Highest row numbers in an entered adaptEq object are too high.", 
              false);
          }
          
          if (min(eq_val) < 0.0) {
            throw Rcpp::exception("Entered equivalence values cannot be negative.",
              false);
          }
          
          NumericVector current_eq (matrowcounts(i), 1.0);
          for (int j = 0; j < static_cast<int>(eq_rn.length()); j++) {
            current_eq(eq_rn(j)) = eq_val(j);
          }
          
          equivalence_list_pre(i) = current_eq;
          
        }
      }
      equivalence_list = equivalence_list_pre;
    } else {
      throw Rcpp::exception("Argument equivalence should be either a numeric vector or a list of such vectors.", false);
    }
  } else {
    if (preexisting) {
      equivalence_count = mpm_count;
    } else {
      equivalence_count = vrm_count;
    }
    NumericVector equivalance_vec_pre (equivalence_count, 1.0);
    equivalence_vec = equivalance_vec_pre;
  }
  
  //Rcout << "cleanup3 R    ";
  
  // process stageframe, supplement, repmatrix, and allstages list for fbMPMs
  if (funcbased) {
    // Create function-based MPMs and assign them to mpm_list
    List allstages_all_pre (vrm_count);
    List allmodels_all_pre (vrm_count);
    
    for (int i = 0; i < vrm_count; i++) {
      List current_vrm = as<List>(vrm_list(i));
      DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
      DataFrame current_supplement = as<DataFrame>(supplement_list(i));
      
      arma::mat current_repmatrix;
      if (format_vec(i) < 5) current_repmatrix = as<arma::mat>(repmatrix_list(i));
      
      int ehrlen_format {1}; // This will need to be dealt with differently later
      
      int mpm_style {1};
      int filter_style {1};
      if (format_vec(i) < 3) {
        mpm_style = 0;
      } else if (format_vec(i) == 4) {
        mpm_style = 2;
        filter_style = 2;
      }
      
      DataFrame current_mpm_allstages;
      if (format_vec(i) < 5) {
        current_mpm_allstages = theoldpizzle(current_stageframe,
          current_supplement, current_repmatrix, firstage_vec(i), finalage_vec(i),
          ehrlen_format, mpm_style, cont_vec(i), filter_style); // Last term removes unused rows & cols
      } else {
        DataFrame leslie_allstages = clone(as<DataFrame>(stageframe_list(i)));
        current_mpm_allstages = leslie_allstages;
      }
      allstages_all_pre(i) = current_mpm_allstages;
      
      if (format_vec(i) < 5) {
        DataFrame chosen_stageframe_pre = clone(as<DataFrame>(stageframe_list(i)));
        
        IntegerVector removal_row = {static_cast<int>(chosen_stageframe_pre.nrows())};
        StringVector removal_var = {"stage_id"};
        DataFrame chosen_stageframe = LefkoUtils::df_remove(chosen_stageframe_pre,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        stageframe_list(i) = chosen_stageframe;
      }
      
      // vrm_input processing
      // Move model summaries to appropriate RObjects
      RObject current_surv_model;
      RObject current_obs_model;
      RObject current_size_model;
      RObject current_sizeb_model;
      RObject current_sizec_model;
      RObject current_repst_model;
      RObject current_fec_model;
      RObject current_jsurv_model;
      RObject current_jobs_model;
      RObject current_jsize_model;
      RObject current_jsizeb_model;
      RObject current_jsizec_model;
      RObject current_jrepst_model;
      RObject current_jmatst_model;
      DataFrame current_paramnames;
      
      DataFrame vrm_frame = as<DataFrame>(current_vrm["vrm_frame"]);
      DataFrame year_frame = as<DataFrame>(current_vrm["year_frame"]);
      DataFrame patch_frame = as<DataFrame>(current_vrm["patch_frame"]);
      DataFrame group2_frame = as<DataFrame>(current_vrm["group2_frame"]);
      DataFrame group1_frame = as<DataFrame>(current_vrm["group1_frame"]);
      DataFrame dist_frame = as<DataFrame>(current_vrm["dist_frame"]);
      NumericVector st_frame = as<NumericVector>(current_vrm["st_frame"]);
      
      CharacterVector main_effect_1 = as<CharacterVector>(vrm_frame["main_effect_1"]);
      CharacterVector effects_names = clone(main_effect_1);
      
      CharacterVector main_effect_2;
      if (main_effect_1.length() > 20) {
        main_effect_2 = as<CharacterVector>(vrm_frame["main_effect_2"]);
        
        for (int i = 0; i < main_effect_1.length(); i++) {
          if (i > 16) {
            effects_names(i) += ":";
            effects_names(i) += main_effect_2(i);
          }
        }
      }
        
      CharacterVector year_names = as<CharacterVector>(year_frame["years"]);
      CharacterVector patch_names = as<CharacterVector>(patch_frame["patches"]);
      CharacterVector group_names = as<CharacterVector>(group2_frame["groups"]);
      
      bool zi_yn = false;
      int vrm_length = vrm_frame.length();
      
      NumericVector surv_num = as<NumericVector>(vrm_frame["surv"]);
      NumericVector obs_num = as<NumericVector>(vrm_frame["obs"]);
      NumericVector sizea_num = as<NumericVector>(vrm_frame["sizea"]);
      NumericVector sizeb_num = as<NumericVector>(vrm_frame["sizeb"]);
      NumericVector sizec_num = as<NumericVector>(vrm_frame["sizec"]);
      NumericVector repst_num = as<NumericVector>(vrm_frame["repst"]);
      NumericVector fec_num = as<NumericVector>(vrm_frame["fec"]);
      NumericVector jsurv_num = as<NumericVector>(vrm_frame["jsurv"]);
      NumericVector jobs_num = as<NumericVector>(vrm_frame["jobs"]);
      NumericVector jsizea_num = as<NumericVector>(vrm_frame["jsizea"]);
      NumericVector jsizeb_num = as<NumericVector>(vrm_frame["jsizeb"]);
      NumericVector jsizec_num = as<NumericVector>(vrm_frame["jsizec"]);
      NumericVector jrepst_num = as<NumericVector>(vrm_frame["jrepst"]);
      NumericVector jmatst_num = as<NumericVector>(vrm_frame["jmatst"]);
      
      NumericVector surv_year = as<NumericVector>(year_frame["surv"]);
      NumericVector obs_year = as<NumericVector>(year_frame["obs"]);
      NumericVector sizea_year = as<NumericVector>(year_frame["sizea"]);
      NumericVector sizeb_year = as<NumericVector>(year_frame["sizeb"]);
      NumericVector sizec_year = as<NumericVector>(year_frame["sizec"]);
      NumericVector repst_year = as<NumericVector>(year_frame["repst"]);
      NumericVector fec_year = as<NumericVector>(year_frame["fec"]);
      NumericVector jsurv_year = as<NumericVector>(year_frame["jsurv"]);
      NumericVector jobs_year = as<NumericVector>(year_frame["jobs"]);
      NumericVector jsizea_year = as<NumericVector>(year_frame["jsizea"]);
      NumericVector jsizeb_year = as<NumericVector>(year_frame["jsizeb"]);
      NumericVector jsizec_year = as<NumericVector>(year_frame["jsizec"]);
      NumericVector jrepst_year = as<NumericVector>(year_frame["jrepst"]);
      NumericVector jmatst_year = as<NumericVector>(year_frame["jmatst"]);
      
      NumericVector surv_patch = as<NumericVector>(patch_frame["surv"]);
      NumericVector obs_patch = as<NumericVector>(patch_frame["obs"]);
      NumericVector sizea_patch = as<NumericVector>(patch_frame["sizea"]);
      NumericVector sizeb_patch = as<NumericVector>(patch_frame["sizeb"]);
      NumericVector sizec_patch = as<NumericVector>(patch_frame["sizec"]);
      NumericVector repst_patch = as<NumericVector>(patch_frame["repst"]);
      NumericVector fec_patch = as<NumericVector>(patch_frame["fec"]);
      NumericVector jsurv_patch = as<NumericVector>(patch_frame["jsurv"]);
      NumericVector jobs_patch = as<NumericVector>(patch_frame["jobs"]);
      NumericVector jsizea_patch = as<NumericVector>(patch_frame["jsizea"]);
      NumericVector jsizeb_patch = as<NumericVector>(patch_frame["jsizeb"]);
      NumericVector jsizec_patch = as<NumericVector>(patch_frame["jsizec"]);
      NumericVector jrepst_patch = as<NumericVector>(patch_frame["jrepst"]);
      NumericVector jmatst_patch = as<NumericVector>(patch_frame["jmatst"]);
      
      NumericVector surv_group2 = as<NumericVector>(group2_frame["surv"]);
      NumericVector obs_group2 = as<NumericVector>(group2_frame["obs"]);
      NumericVector sizea_group2 = as<NumericVector>(group2_frame["sizea"]);
      NumericVector sizeb_group2 = as<NumericVector>(group2_frame["sizeb"]);
      NumericVector sizec_group2 = as<NumericVector>(group2_frame["sizec"]);
      NumericVector repst_group2 = as<NumericVector>(group2_frame["repst"]);
      NumericVector fec_group2 = as<NumericVector>(group2_frame["fec"]);
      NumericVector jsurv_group2 = as<NumericVector>(group2_frame["jsurv"]);
      NumericVector jobs_group2 = as<NumericVector>(group2_frame["jobs"]);
      NumericVector jsizea_group2 = as<NumericVector>(group2_frame["jsizea"]);
      NumericVector jsizeb_group2 = as<NumericVector>(group2_frame["jsizeb"]);
      NumericVector jsizec_group2 = as<NumericVector>(group2_frame["jsizec"]);
      NumericVector jrepst_group2 = as<NumericVector>(group2_frame["jrepst"]);
      NumericVector jmatst_group2 = as<NumericVector>(group2_frame["jmatst"]);
      
      NumericVector surv_group1 = as<NumericVector>(group1_frame["surv"]);
      NumericVector obs_group1 = as<NumericVector>(group1_frame["obs"]);
      NumericVector sizea_group1 = as<NumericVector>(group1_frame["sizea"]);
      NumericVector sizeb_group1 = as<NumericVector>(group1_frame["sizeb"]);
      NumericVector sizec_group1 = as<NumericVector>(group1_frame["sizec"]);
      NumericVector repst_group1 = as<NumericVector>(group1_frame["repst"]);
      NumericVector fec_group1 = as<NumericVector>(group1_frame["fec"]);
      NumericVector jsurv_group1 = as<NumericVector>(group1_frame["jsurv"]);
      NumericVector jobs_group1 = as<NumericVector>(group1_frame["jobs"]);
      NumericVector jsizea_group1 = as<NumericVector>(group1_frame["jsizea"]);
      NumericVector jsizeb_group1 = as<NumericVector>(group1_frame["jsizeb"]);
      NumericVector jsizec_group1 = as<NumericVector>(group1_frame["jsizec"]);
      NumericVector jrepst_group1 = as<NumericVector>(group1_frame["jrepst"]);
      NumericVector jmatst_group1 = as<NumericVector>(group1_frame["jmatst"]);
        
      StringVector distribs = as<StringVector>(dist_frame["dist"]);
      String surv_dist = distribs(0);
      String obs_dist = distribs(1);
      String sizea_dist = distribs(2);
      String sizeb_dist = distribs(3);
      String sizec_dist = distribs(4);
      String repst_dist = distribs(5);
      String fec_dist = distribs(6);
      String jsurv_dist = distribs(7);
      String jobs_dist = distribs(8);
      String jsizea_dist = distribs(9);
      String jsizeb_dist = distribs(10);
      String jsizec_dist = distribs(11);
      String jrepst_dist = distribs(12);
      String jmatst_dist = distribs(13);
      
      double sizea_st = st_frame(2);
      double sizeb_st = st_frame(3);
      double sizec_st = st_frame(4);
      double fec_st = st_frame(6);
      double jsizea_st = st_frame(9);
      double jsizeb_st = st_frame(10);
      double jsizec_st = st_frame(11);
      
      NumericVector sizea_zi;
      NumericVector sizeb_zi;
      NumericVector sizec_zi;
      NumericVector fec_zi;
      NumericVector jsizea_zi;
      NumericVector jsizeb_zi;
      NumericVector jsizec_zi;
      
      NumericVector year_sizea_zi;
      NumericVector year_sizeb_zi;
      NumericVector year_sizec_zi;
      NumericVector year_fec_zi;
      NumericVector year_jsizea_zi;
      NumericVector year_jsizeb_zi;
      NumericVector year_jsizec_zi;
      
      NumericVector patch_sizea_zi;
      NumericVector patch_sizeb_zi;
      NumericVector patch_sizec_zi;
      NumericVector patch_fec_zi;
      NumericVector patch_jsizea_zi;
      NumericVector patch_jsizeb_zi;
      NumericVector patch_jsizec_zi;
      
      NumericVector group2_sizea_zi;
      NumericVector group2_sizeb_zi;
      NumericVector group2_sizec_zi;
      NumericVector group2_fec_zi;
      NumericVector group2_jsizea_zi;
      NumericVector group2_jsizeb_zi;
      NumericVector group2_jsizec_zi;
      
      NumericVector group1_sizea_zi;
      NumericVector group1_sizeb_zi;
      NumericVector group1_sizec_zi;
      NumericVector group1_fec_zi;
      NumericVector group1_jsizea_zi;
      NumericVector group1_jsizeb_zi;
      NumericVector group1_jsizec_zi;
      
      NumericVector dud_zi;
      
      if (vrm_length > 16) {
        zi_yn = true;
        
        sizea_zi = as<NumericVector>(vrm_frame["sizea_zi"]);
        sizeb_zi = as<NumericVector>(vrm_frame["sizeb_zi"]);
        sizec_zi = as<NumericVector>(vrm_frame["sizec_zi"]);
        fec_zi = as<NumericVector>(vrm_frame["fec_zi"]);
        jsizea_zi = as<NumericVector>(vrm_frame["jsizea_zi"]);
        jsizeb_zi = as<NumericVector>(vrm_frame["jsizeb_zi"]);
        jsizec_zi = as<NumericVector>(vrm_frame["jsizec_zi"]);
        
        year_sizea_zi = as<NumericVector>(year_frame["sizea_zi"]);
        year_sizeb_zi = as<NumericVector>(year_frame["sizeb_zi"]);
        year_sizec_zi = as<NumericVector>(year_frame["sizec_zi"]);
        year_fec_zi = as<NumericVector>(year_frame["fec_zi"]);
        year_jsizea_zi = as<NumericVector>(year_frame["jsizea_zi"]);
        year_jsizeb_zi = as<NumericVector>(year_frame["jsizeb_zi"]);
        year_jsizec_zi = as<NumericVector>(year_frame["jsizec_zi"]);
        
        patch_sizea_zi = as<NumericVector>(patch_frame["sizea_zi"]);
        patch_sizeb_zi = as<NumericVector>(patch_frame["sizeb_zi"]);
        patch_sizec_zi = as<NumericVector>(patch_frame["sizec_zi"]);
        patch_fec_zi = as<NumericVector>(patch_frame["fec_zi"]);
        patch_jsizea_zi = as<NumericVector>(patch_frame["jsizea_zi"]);
        patch_jsizeb_zi = as<NumericVector>(patch_frame["jsizeb_zi"]);
        patch_jsizec_zi = as<NumericVector>(patch_frame["jsizec_zi"]);
        
        group2_sizea_zi = as<NumericVector>(group2_frame["sizea_zi"]);
        group2_sizeb_zi = as<NumericVector>(group2_frame["sizeb_zi"]);
        group2_sizec_zi = as<NumericVector>(group2_frame["sizec_zi"]);
        group2_fec_zi = as<NumericVector>(group2_frame["fec_zi"]);
        group2_jsizea_zi = as<NumericVector>(group2_frame["jsizea_zi"]);
        group2_jsizeb_zi = as<NumericVector>(group2_frame["jsizeb_zi"]);
        group2_jsizec_zi = as<NumericVector>(group2_frame["jsizec_zi"]);
        
        group1_sizea_zi = as<NumericVector>(group1_frame["sizea_zi"]);
        group1_sizeb_zi = as<NumericVector>(group1_frame["sizeb_zi"]);
        group1_sizec_zi = as<NumericVector>(group1_frame["sizec_zi"]);
        group1_fec_zi = as<NumericVector>(group1_frame["fec_zi"]);
        group1_jsizea_zi = as<NumericVector>(group1_frame["jsizea_zi"]);
        group1_jsizeb_zi = as<NumericVector>(group1_frame["jsizeb_zi"]);
        group1_jsizec_zi = as<NumericVector>(group1_frame["jsizec_zi"]);
      }
      
      CharacterVector indcova_names;
      CharacterVector indcovb_names;
      CharacterVector indcovc_names;
      
      NumericVector surv_indcova2;
      NumericVector surv_indcovb2;
      NumericVector surv_indcovc2;
      NumericVector obs_indcova2;
      NumericVector obs_indcovb2;
      NumericVector obs_indcovc2;
      NumericVector sizea_indcova2;
      NumericVector sizea_indcovb2;
      NumericVector sizea_indcovc2;
      NumericVector sizeb_indcova2;
      NumericVector sizeb_indcovb2;
      NumericVector sizeb_indcovc2;
      NumericVector sizec_indcova2;
      NumericVector sizec_indcovb2;
      NumericVector sizec_indcovc2;
      NumericVector repst_indcova2;
      NumericVector repst_indcovb2;
      NumericVector repst_indcovc2;
      NumericVector fec_indcova2;
      NumericVector fec_indcovb2;
      NumericVector fec_indcovc2;
      NumericVector jsurv_indcova2;
      NumericVector jsurv_indcovb2;
      NumericVector jsurv_indcovc2;
      NumericVector jobs_indcova2;
      NumericVector jobs_indcovb2;
      NumericVector jobs_indcovc2;
      NumericVector jsizea_indcova2;
      NumericVector jsizea_indcovb2;
      NumericVector jsizea_indcovc2;
      NumericVector jsizeb_indcova2;
      NumericVector jsizeb_indcovb2;
      NumericVector jsizeb_indcovc2;
      NumericVector jsizec_indcova2;
      NumericVector jsizec_indcovb2;
      NumericVector jsizec_indcovc2;
      NumericVector jrepst_indcova2;
      NumericVector jrepst_indcovb2;
      NumericVector jrepst_indcovc2;
      NumericVector jmatst_indcova2;
      NumericVector jmatst_indcovb2;
      NumericVector jmatst_indcovc2;
      
      NumericVector sizea_indcova2_zi;
      NumericVector sizea_indcovb2_zi;
      NumericVector sizea_indcovc2_zi;
      NumericVector sizeb_indcova2_zi;
      NumericVector sizeb_indcovb2_zi;
      NumericVector sizeb_indcovc2_zi;
      NumericVector sizec_indcova2_zi;
      NumericVector sizec_indcovb2_zi;
      NumericVector sizec_indcovc2_zi;
      NumericVector fec_indcova2_zi;
      NumericVector fec_indcovb2_zi;
      NumericVector fec_indcovc2_zi;
      NumericVector jsizea_indcova2_zi;
      NumericVector jsizea_indcovb2_zi;
      NumericVector jsizea_indcovc2_zi;
      NumericVector jsizeb_indcova2_zi;
      NumericVector jsizeb_indcovb2_zi;
      NumericVector jsizeb_indcovc2_zi;
      NumericVector jsizec_indcova2_zi;
      NumericVector jsizec_indcovb2_zi;
      NumericVector jsizec_indcovc2_zi;
      
      NumericVector surv_indcova1;
      NumericVector surv_indcovb1;
      NumericVector surv_indcovc1;
      NumericVector obs_indcova1;
      NumericVector obs_indcovb1;
      NumericVector obs_indcovc1;
      NumericVector sizea_indcova1;
      NumericVector sizea_indcovb1;
      NumericVector sizea_indcovc1;
      NumericVector sizeb_indcova1;
      NumericVector sizeb_indcovb1;
      NumericVector sizeb_indcovc1;
      NumericVector sizec_indcova1;
      NumericVector sizec_indcovb1;
      NumericVector sizec_indcovc1;
      NumericVector repst_indcova1;
      NumericVector repst_indcovb1;
      NumericVector repst_indcovc1;
      NumericVector fec_indcova1;
      NumericVector fec_indcovb1;
      NumericVector fec_indcovc1;
      NumericVector jsurv_indcova1;
      NumericVector jsurv_indcovb1;
      NumericVector jsurv_indcovc1;
      NumericVector jobs_indcova1;
      NumericVector jobs_indcovb1;
      NumericVector jobs_indcovc1;
      NumericVector jsizea_indcova1;
      NumericVector jsizea_indcovb1;
      NumericVector jsizea_indcovc1;
      NumericVector jsizeb_indcova1;
      NumericVector jsizeb_indcovb1;
      NumericVector jsizeb_indcovc1;
      NumericVector jsizec_indcova1;
      NumericVector jsizec_indcovb1;
      NumericVector jsizec_indcovc1;
      NumericVector jrepst_indcova1;
      NumericVector jrepst_indcovb1;
      NumericVector jrepst_indcovc1;
      NumericVector jmatst_indcova1;
      NumericVector jmatst_indcovb1;
      NumericVector jmatst_indcovc1;
      
      NumericVector sizea_indcova1_zi;
      NumericVector sizea_indcovb1_zi;
      NumericVector sizea_indcovc1_zi;
      NumericVector sizeb_indcova1_zi;
      NumericVector sizeb_indcovb1_zi;
      NumericVector sizeb_indcovc1_zi;
      NumericVector sizec_indcova1_zi;
      NumericVector sizec_indcovb1_zi;
      NumericVector sizec_indcovc1_zi;
      NumericVector fec_indcova1_zi;
      NumericVector fec_indcovb1_zi;
      NumericVector fec_indcovc1_zi;
      NumericVector jsizea_indcova1_zi;
      NumericVector jsizea_indcovb1_zi;
      NumericVector jsizea_indcovc1_zi;
      NumericVector jsizeb_indcova1_zi;
      NumericVector jsizeb_indcovb1_zi;
      NumericVector jsizeb_indcovc1_zi;
      NumericVector jsizec_indcova1_zi;
      NumericVector jsizec_indcovb1_zi;
      NumericVector jsizec_indcovc1_zi;
      
      int modelsuite_length = current_vrm.length();
      CharacterVector modelsuite_names = current_vrm.attr("names");
      
      for (int i = 0; i < modelsuite_length; i++) {
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova2_frame")) {
          DataFrame indcova2_frame = as<DataFrame>(current_vrm["indcova2_frame"]);
          
          indcova_names = indcova2_frame["indcova"];
          
          surv_indcova2 = indcova2_frame["surv"];
          obs_indcova2 = indcova2_frame["obs"];
          sizea_indcova2 = indcova2_frame["sizea"];
          sizeb_indcova2 = indcova2_frame["sizeb"];
          sizec_indcova2 = indcova2_frame["sizec"];
          repst_indcova2 = indcova2_frame["repst"];
          fec_indcova2 = indcova2_frame["fec"];
          
          jsurv_indcova2 = indcova2_frame["jsurv"];
          jobs_indcova2 = indcova2_frame["jobs"];
          jsizea_indcova2 = indcova2_frame["jsizea"];
          jsizeb_indcova2 = indcova2_frame["jsizeb"];
          jsizec_indcova2 = indcova2_frame["jsizec"];
          jrepst_indcova2 = indcova2_frame["jrepst"];
          jmatst_indcova2 = indcova2_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcova2_zi = indcova2_frame["sizea_zi"];
            sizeb_indcova2_zi = indcova2_frame["sizeb_zi"];
            sizec_indcova2_zi = indcova2_frame["sizec_zi"];
            fec_indcova2_zi = indcova2_frame["fec_zi"];
            jsizea_indcova2_zi = indcova2_frame["jsizea_zi"];
            jsizeb_indcova2_zi = indcova2_frame["jsizeb_zi"];
            jsizec_indcova2_zi = indcova2_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova1_frame")) {
          DataFrame indcova1_frame = as<DataFrame>(current_vrm["indcova1_frame"]);
          
          indcova_names = indcova1_frame["indcova"];
          
          surv_indcova1 = indcova1_frame["surv"];
          obs_indcova1 = indcova1_frame["obs"];
          sizea_indcova1 = indcova1_frame["sizea"];
          sizeb_indcova1 = indcova1_frame["sizeb"];
          sizec_indcova1 = indcova1_frame["sizec"];
          repst_indcova1 = indcova1_frame["repst"];
          fec_indcova1 = indcova1_frame["fec"];
          
          jsurv_indcova1 = indcova1_frame["jsurv"];
          jobs_indcova1 = indcova1_frame["jobs"];
          jsizea_indcova1 = indcova1_frame["jsizea"];
          jsizeb_indcova1 = indcova1_frame["jsizeb"];
          jsizec_indcova1 = indcova1_frame["jsizec"];
          jrepst_indcova1 = indcova1_frame["jrepst"];
          jmatst_indcova1 = indcova1_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcova1_zi = indcova1_frame["sizea_zi"];
            sizeb_indcova1_zi = indcova1_frame["sizeb_zi"];
            sizec_indcova1_zi = indcova1_frame["sizec_zi"];
            fec_indcova1_zi = indcova1_frame["fec_zi"];
            jsizea_indcova1_zi = indcova1_frame["jsizea_zi"];
            jsizeb_indcova1_zi = indcova1_frame["jsizeb_zi"];
            jsizec_indcova1_zi = indcova1_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb2_frame")) {
          DataFrame indcovb2_frame = as<DataFrame>(current_vrm["indcovb2_frame"]);
          
          indcovb_names = indcovb2_frame["indcovb"];
          
          surv_indcovb2 = indcovb2_frame["surv"];
          obs_indcovb2 = indcovb2_frame["obs"];
          sizea_indcovb2 = indcovb2_frame["sizea"];
          sizeb_indcovb2 = indcovb2_frame["sizeb"];
          sizec_indcovb2 = indcovb2_frame["sizec"];
          repst_indcovb2 = indcovb2_frame["repst"];
          fec_indcovb2 = indcovb2_frame["fec"];
          
          jsurv_indcovb2 = indcovb2_frame["jsurv"];
          jobs_indcovb2 = indcovb2_frame["jobs"];
          jsizea_indcovb2 = indcovb2_frame["jsizea"];
          jsizeb_indcovb2 = indcovb2_frame["jsizeb"];
          jsizec_indcovb2 = indcovb2_frame["jsizec"];
          jrepst_indcovb2 = indcovb2_frame["jrepst"];
          jmatst_indcovb2 = indcovb2_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovb2_zi = indcovb2_frame["sizea_zi"];
            sizeb_indcovb2_zi = indcovb2_frame["sizeb_zi"];
            sizec_indcovb2_zi = indcovb2_frame["sizec_zi"];
            fec_indcovb2_zi = indcovb2_frame["fec_zi"];
            jsizea_indcovb2_zi = indcovb2_frame["jsizea_zi"];
            jsizeb_indcovb2_zi = indcovb2_frame["jsizeb_zi"];
            jsizec_indcovb2_zi = indcovb2_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb1_frame")) {
          DataFrame indcovb1_frame = as<DataFrame>(current_vrm["indcovb1_frame"]);
          
          indcovb_names = indcovb1_frame["indcovb"];
          
          surv_indcovb1 = indcovb1_frame["surv"];
          obs_indcovb1 = indcovb1_frame["obs"];
          sizea_indcovb1 = indcovb1_frame["sizea"];
          sizeb_indcovb1 = indcovb1_frame["sizeb"];
          sizec_indcovb1 = indcovb1_frame["sizec"];
          repst_indcovb1 = indcovb1_frame["repst"];
          fec_indcovb1 = indcovb1_frame["fec"];
          
          jsurv_indcovb1 = indcovb1_frame["jsurv"];
          jobs_indcovb1 = indcovb1_frame["jobs"];
          jsizea_indcovb1 = indcovb1_frame["jsizea"];
          jsizeb_indcovb1 = indcovb1_frame["jsizeb"];
          jsizec_indcovb1 = indcovb1_frame["jsizec"];
          jrepst_indcovb1 = indcovb1_frame["jrepst"];
          jmatst_indcovb1 = indcovb1_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovb1_zi = indcovb1_frame["sizea_zi"];
            sizeb_indcovb1_zi = indcovb1_frame["sizeb_zi"];
            sizec_indcovb1_zi = indcovb1_frame["sizec_zi"];
            fec_indcovb1_zi = indcovb1_frame["fec_zi"];
            jsizea_indcovb1_zi = indcovb1_frame["jsizea_zi"];
            jsizeb_indcovb1_zi = indcovb1_frame["jsizeb_zi"];
            jsizec_indcovb1_zi = indcovb1_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc2_frame")) {
          DataFrame indcovc2_frame = as<DataFrame>(current_vrm["indcovc2_frame"]);
          
          indcovc_names = indcovc2_frame["indcovc"];
          
          surv_indcovc2 = indcovc2_frame["surv"];
          obs_indcovc2 = indcovc2_frame["obs"];
          sizea_indcovc2 = indcovc2_frame["sizea"];
          sizeb_indcovc2 = indcovc2_frame["sizeb"];
          sizec_indcovc2 = indcovc2_frame["sizec"];
          repst_indcovc2 = indcovc2_frame["repst"];
          fec_indcovc2 = indcovc2_frame["fec"];
          
          jsurv_indcovc2 = indcovc2_frame["jsurv"];
          jobs_indcovc2 = indcovc2_frame["jobs"];
          jsizea_indcovc2 = indcovc2_frame["jsizea"];
          jsizeb_indcovc2 = indcovc2_frame["jsizeb"];
          jsizec_indcovc2 = indcovc2_frame["jsizec"];
          jrepst_indcovc2 = indcovc2_frame["jrepst"];
          jmatst_indcovc2 = indcovc2_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovc2_zi = indcovc2_frame["sizea_zi"];
            sizeb_indcovc2_zi = indcovc2_frame["sizeb_zi"];
            sizec_indcovc2_zi = indcovc2_frame["sizec_zi"];
            fec_indcovc2_zi = indcovc2_frame["fec_zi"];
            jsizea_indcovc2_zi = indcovc2_frame["jsizea_zi"];
            jsizeb_indcovc2_zi = indcovc2_frame["jsizeb_zi"];
            jsizec_indcovc2_zi = indcovc2_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc1_frame")) {
          DataFrame indcovc1_frame = as<DataFrame>(current_vrm["indcovc1_frame"]);
          
          indcovc_names = indcovc1_frame["indcovc"];
          
          surv_indcovc1 = indcovc1_frame["surv"];
          obs_indcovc1 = indcovc1_frame["obs"];
          sizea_indcovc1 = indcovc1_frame["sizea"];
          sizeb_indcovc1 = indcovc1_frame["sizeb"];
          sizec_indcovc1 = indcovc1_frame["sizec"];
          repst_indcovc1 = indcovc1_frame["repst"];
          fec_indcovc1 = indcovc1_frame["fec"];
          
          jsurv_indcovc1 = indcovc1_frame["jsurv"];
          jobs_indcovc1 = indcovc1_frame["jobs"];
          jsizea_indcovc1 = indcovc1_frame["jsizea"];
          jsizeb_indcovc1 = indcovc1_frame["jsizeb"];
          jsizec_indcovc1 = indcovc1_frame["jsizec"];
          jrepst_indcovc1 = indcovc1_frame["jrepst"];
          jmatst_indcovc1 = indcovc1_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovc1_zi = indcovc1_frame["sizea_zi"];
            sizeb_indcovc1_zi = indcovc1_frame["sizeb_zi"];
            sizec_indcovc1_zi = indcovc1_frame["sizec_zi"];
            fec_indcovc1_zi = indcovc1_frame["fec_zi"];
            jsizea_indcovc1_zi = indcovc1_frame["jsizea_zi"];
            jsizeb_indcovc1_zi = indcovc1_frame["jsizeb_zi"];
            jsizec_indcovc1_zi = indcovc1_frame["jsizec_zi"];
          }
        }
      }
      
      CharacterVector list_names = {"fixed_slopes", "year_slopes", "patch_slopes",
        "group2_slopes", "dist", "zi", "fixed_zi", "year_zi", "patch_zi",
        "group2_zi", "indcova_names", "indcova2_slopes", "indcova2_zi",
        "indcovb_names", "indcovb2_slopes", "indcovb2_zi", "indcovc_names",
        "indcovc2_slopes", "indcovc2_zi", "year_names", "patch_names",
        "group_names", "main_effect_1", "main_effect_2", "sigma_theta",
        "effects_names", "group1_slopes", "group1_zi", "indcova1_slopes",
        "indcovb1_slopes", "indcovc1_slopes", "indcova1_zi", "indcovb1_zi",
        "indcovc1_zi"};
      
      List surv_list(34);
      surv_list(0) = surv_num;
      surv_list(1) = surv_year;
      surv_list(2) = surv_patch;
      surv_list(3) = surv_group2;
      surv_list(4) = surv_dist;
      surv_list(5) = false;
      surv_list(6) = dud_zi;
      surv_list(7) = dud_zi;
      surv_list(8) = dud_zi;
      surv_list(9) = dud_zi;
      surv_list(10) = indcova_names;
      surv_list(11) = surv_indcova2;
      surv_list(12) = dud_zi;
      surv_list(13) = indcovb_names;
      surv_list(14) = surv_indcovb2;
      surv_list(15) = dud_zi;
      surv_list(16) = indcovc_names;
      surv_list(17) = surv_indcovc2;
      surv_list(18) = dud_zi;
      surv_list(19) = year_names;
      surv_list(20) = patch_names;
      surv_list(21) = group_names;
      surv_list(22) = main_effect_1;
      surv_list(23) = main_effect_2;
      surv_list(24) = 1.0;
      surv_list(25) = effects_names;
      surv_list(26) = surv_group1;
      surv_list(27) = dud_zi;
      surv_list(28) = surv_indcova1;
      surv_list(29) = surv_indcovb1;
      surv_list(30) = surv_indcovc1;
      surv_list(31) = dud_zi;
      surv_list(32) = dud_zi;
      surv_list(33) = dud_zi;
      
      List obs_list(34);
      obs_list(0) = obs_num;
      obs_list(1) = obs_year;
      obs_list(2) = obs_patch;
      obs_list(3) = obs_group2;
      obs_list(4) = obs_dist;
      obs_list(5) = false;
      obs_list(6) = dud_zi;
      obs_list(7) = dud_zi;
      obs_list(8) = dud_zi;
      obs_list(9) = dud_zi;
      obs_list(10) = indcova_names;
      obs_list(11) = obs_indcova2;
      obs_list(12) = dud_zi;
      obs_list(13) = indcovb_names;
      obs_list(14) = obs_indcovb2;
      obs_list(15) = dud_zi;
      obs_list(16) = indcovc_names;
      obs_list(17) = obs_indcovc2;
      obs_list(18) = dud_zi;
      obs_list(19) = year_names;
      obs_list(20) = patch_names;
      obs_list(21) = group_names;
      obs_list(22) = main_effect_1;
      obs_list(23) = main_effect_2;
      obs_list(24) = 1.0;
      obs_list(25) = effects_names;
      obs_list(26) = obs_group1;
      obs_list(27) = dud_zi;
      obs_list(28) = obs_indcova1;
      obs_list(29) = obs_indcovb1;
      obs_list(30) = obs_indcovc1;
      obs_list(31) = dud_zi;
      obs_list(32) = dud_zi;
      obs_list(33) = dud_zi;
      
      List sizea_list(34);
      sizea_list(0) = sizea_num;
      sizea_list(1) = sizea_year;
      sizea_list(2) = sizea_patch;
      sizea_list(3) = sizea_group2;
      sizea_list(4) = sizea_dist;
      sizea_list(5) = zi_yn;
      sizea_list(6) = sizea_zi;
      sizea_list(7) = year_sizea_zi;
      sizea_list(8) = patch_sizea_zi;
      sizea_list(9) = group2_sizea_zi;
      sizea_list(10) = indcova_names;
      sizea_list(11) = sizea_indcova2;
      sizea_list(12) = sizea_indcova2_zi;
      sizea_list(13) = indcovb_names;
      sizea_list(14) = sizea_indcovb2;
      sizea_list(15) = sizea_indcovb2_zi;
      sizea_list(16) = indcovc_names;
      sizea_list(17) = sizea_indcovc2;
      sizea_list(18) = sizea_indcovc2_zi;
      sizea_list(19) = year_names;
      sizea_list(20) = patch_names;
      sizea_list(21) = group_names;
      sizea_list(22) = main_effect_1;
      sizea_list(23) = main_effect_2;
      sizea_list(24) = sizea_st;
      sizea_list(25) = effects_names;
      sizea_list(26) = sizea_group1;
      sizea_list(27) = group1_sizea_zi;
      sizea_list(28) = sizea_indcova1;
      sizea_list(29) = sizea_indcovb1;
      sizea_list(30) = sizea_indcovc1;
      sizea_list(31) = sizea_indcova1_zi;
      sizea_list(32) = sizea_indcovb1_zi;
      sizea_list(33) = sizea_indcovc1_zi;
      
      List sizeb_list(34);
      sizeb_list(0) = sizeb_num;
      sizeb_list(1) = sizeb_year;
      sizeb_list(2) = sizeb_patch;
      sizeb_list(3) = sizeb_group2;
      sizeb_list(4) = sizeb_dist;
      sizeb_list(5) = zi_yn;
      sizeb_list(6) = sizeb_zi;
      sizeb_list(7) = year_sizeb_zi;
      sizeb_list(8) = patch_sizeb_zi;
      sizeb_list(9) = group2_sizeb_zi;
      sizeb_list(10) = indcova_names;
      sizeb_list(11) = sizeb_indcova2;
      sizeb_list(12) = sizeb_indcova2_zi;
      sizeb_list(13) = indcovb_names;
      sizeb_list(14) = sizeb_indcovb2;
      sizeb_list(15) = sizeb_indcovb2_zi;
      sizeb_list(16) = indcovc_names;
      sizeb_list(17) = sizeb_indcovc2;
      sizeb_list(18) = sizeb_indcovc2_zi;
      sizeb_list(19) = year_names;
      sizeb_list(20) = patch_names;
      sizeb_list(21) = group_names;
      sizeb_list(22) = main_effect_1;
      sizeb_list(23) = main_effect_2;
      sizeb_list(24) = sizeb_st;
      sizeb_list(25) = effects_names;
      sizeb_list(26) = sizeb_group1;
      sizeb_list(27) = group1_sizeb_zi;
      sizeb_list(28) = sizeb_indcova1;
      sizeb_list(29) = sizeb_indcovb1;
      sizeb_list(30) = sizeb_indcovc1;
      sizeb_list(31) = sizeb_indcova1_zi;
      sizeb_list(32) = sizeb_indcovb1_zi;
      sizeb_list(33) = sizeb_indcovc1_zi;
      
      List sizec_list(34);
      sizec_list(0) = sizec_num;
      sizec_list(1) = sizec_year;
      sizec_list(2) = sizec_patch;
      sizec_list(3) = sizec_group2;
      sizec_list(4) = sizec_dist;
      sizec_list(5) = zi_yn;
      sizec_list(6) = sizec_zi;
      sizec_list(7) = year_sizec_zi;
      sizec_list(8) = patch_sizec_zi;
      sizec_list(9) = group2_sizec_zi;
      sizec_list(10) = indcova_names;
      sizec_list(11) = sizec_indcova2;
      sizec_list(12) = sizec_indcova2_zi;
      sizec_list(13) = indcovb_names;
      sizec_list(14) = sizec_indcovb2;
      sizec_list(15) = sizec_indcovb2_zi;
      sizec_list(16) = indcovc_names;
      sizec_list(17) = sizec_indcovc2;
      sizec_list(18) = sizec_indcovc2_zi;
      sizec_list(19) = year_names;
      sizec_list(20) = patch_names;
      sizec_list(21) = group_names;
      sizec_list(22) = main_effect_1;
      sizec_list(23) = main_effect_2;
      sizec_list(24) = sizec_st;
      sizec_list(25) = effects_names;
      sizec_list(26) = sizec_group1;
      sizec_list(27) = group1_sizec_zi;
      sizec_list(28) = sizec_indcova1;
      sizec_list(29) = sizec_indcovb1;
      sizec_list(30) = sizec_indcovc1;
      sizec_list(31) = sizec_indcova1_zi;
      sizec_list(32) = sizec_indcovb1_zi;
      sizec_list(33) = sizec_indcovc1_zi;
      
      List repst_list(34);
      repst_list(0) = repst_num;
      repst_list(1) = repst_year;
      repst_list(2) = repst_patch;
      repst_list(3) = repst_group2;
      repst_list(4) = repst_dist;
      repst_list(5) = false;
      repst_list(6) = dud_zi;
      repst_list(7) = dud_zi;
      repst_list(8) = dud_zi;
      repst_list(9) = dud_zi;
      repst_list(10) = indcova_names;
      repst_list(11) = repst_indcova2;
      repst_list(12) = dud_zi;
      repst_list(13) = indcovb_names;
      repst_list(14) = repst_indcovb2;
      repst_list(15) = dud_zi;
      repst_list(16) = indcovc_names;
      repst_list(17) = repst_indcovc2;
      repst_list(18) = dud_zi;
      repst_list(19) = year_names;
      repst_list(20) = patch_names;
      repst_list(21) = group_names;
      repst_list(22) = main_effect_1;
      repst_list(23) = main_effect_2;
      repst_list(24) = 1.0;
      repst_list(25) = effects_names;
      repst_list(26) = repst_group1;
      repst_list(27) = dud_zi;
      repst_list(28) = repst_indcova1;
      repst_list(29) = repst_indcovb1;
      repst_list(30) = repst_indcovc1;
      repst_list(31) = dud_zi;
      repst_list(32) = dud_zi;
      repst_list(33) = dud_zi;
      
      List fec_list(34);
      fec_list(0) = fec_num;
      fec_list(1) = fec_year;
      fec_list(2) = fec_patch;
      fec_list(3) = fec_group2;
      fec_list(4) = fec_dist;
      fec_list(5) = zi_yn;
      fec_list(6) = fec_zi;
      fec_list(7) = year_fec_zi;
      fec_list(8) = patch_fec_zi;
      fec_list(9) = group2_fec_zi;
      fec_list(10) = indcova_names;
      fec_list(11) = fec_indcova2;
      fec_list(12) = fec_indcova2_zi;
      fec_list(13) = indcovb_names;
      fec_list(14) = fec_indcovb2;
      fec_list(15) = fec_indcovb2_zi;
      fec_list(16) = indcovc_names;
      fec_list(17) = fec_indcovc2;
      fec_list(18) = fec_indcovc2_zi;
      fec_list(19) = year_names;
      fec_list(20) = patch_names;
      fec_list(21) = group_names;
      fec_list(22) = main_effect_1;
      fec_list(23) = main_effect_2;
      fec_list(24) = fec_st;
      fec_list(25) = effects_names;
      fec_list(26) = fec_group1;
      fec_list(27) = group1_fec_zi;
      fec_list(28) = fec_indcova1;
      fec_list(29) = fec_indcovb1;
      fec_list(30) = fec_indcovc1;
      fec_list(31) = fec_indcova1_zi;
      fec_list(32) = fec_indcovb1_zi;
      fec_list(33) = fec_indcovc1_zi;
      
      List jsurv_list(34);
      jsurv_list(0) = jsurv_num;
      jsurv_list(1) = jsurv_year;
      jsurv_list(2) = jsurv_patch;
      jsurv_list(3) = jsurv_group2;
      jsurv_list(4) = jsurv_dist;
      jsurv_list(5) = false;
      jsurv_list(6) = dud_zi;
      jsurv_list(7) = dud_zi;
      jsurv_list(8) = dud_zi;
      jsurv_list(9) = dud_zi;
      jsurv_list(10) = indcova_names;
      jsurv_list(11) = jsurv_indcova2;
      jsurv_list(12) = dud_zi;
      jsurv_list(13) = indcovb_names;
      jsurv_list(14) = jsurv_indcovb2;
      jsurv_list(15) = dud_zi;
      jsurv_list(16) = indcovc_names;
      jsurv_list(17) = jsurv_indcovc2;
      jsurv_list(18) = dud_zi;
      jsurv_list(19) = year_names;
      jsurv_list(20) = patch_names;
      jsurv_list(21) = group_names;
      jsurv_list(22) = main_effect_1;
      jsurv_list(23) = main_effect_2;
      jsurv_list(24) = 1.0;
      jsurv_list(25) = effects_names;
      jsurv_list(26) = jsurv_group1;
      jsurv_list(27) = dud_zi;
      jsurv_list(28) = jsurv_indcova1;
      jsurv_list(29) = jsurv_indcovb1;
      jsurv_list(30) = jsurv_indcovc1;
      jsurv_list(31) = dud_zi;
      jsurv_list(32) = dud_zi;
      jsurv_list(33) = dud_zi;
      
      List jobs_list(34);
      jobs_list(0) = jobs_num;
      jobs_list(1) = jobs_year;
      jobs_list(2) = jobs_patch;
      jobs_list(3) = jobs_group2;
      jobs_list(4) = jobs_dist;
      jobs_list(5) = false;
      jobs_list(6) = dud_zi;
      jobs_list(7) = dud_zi;
      jobs_list(8) = dud_zi;
      jobs_list(9) = dud_zi;
      jobs_list(10) = indcova_names;
      jobs_list(11) = jobs_indcova2;
      jobs_list(12) = dud_zi;
      jobs_list(13) = indcovb_names;
      jobs_list(14) = jobs_indcovb2;
      jobs_list(15) = dud_zi;
      jobs_list(16) = indcovc_names;
      jobs_list(17) = jobs_indcovc2;
      jobs_list(18) = dud_zi;
      jobs_list(19) = year_names;
      jobs_list(20) = patch_names;
      jobs_list(21) = group_names;
      jobs_list(22) = main_effect_1;
      jobs_list(23) = main_effect_2;
      jobs_list(24) = 1.0;
      jobs_list(25) = effects_names;
      jobs_list(26) = jobs_group1;
      jobs_list(27) = dud_zi;
      jobs_list(28) = jobs_indcova1;
      jobs_list(29) = jobs_indcovb1;
      jobs_list(30) = jobs_indcovc1;
      jobs_list(31) = dud_zi;
      jobs_list(32) = dud_zi;
      jobs_list(33) = dud_zi;
      
      List jsizea_list(34);
      jsizea_list(0) = jsizea_num;
      jsizea_list(1) = jsizea_year;
      jsizea_list(2) = jsizea_patch;
      jsizea_list(3) = jsizea_group2;
      jsizea_list(4) = jsizea_dist;
      jsizea_list(5) = zi_yn;
      jsizea_list(6) = jsizea_zi;
      jsizea_list(7) = year_jsizea_zi;
      jsizea_list(8) = patch_jsizea_zi;
      jsizea_list(9) = group2_jsizea_zi;
      jsizea_list(10) = indcova_names;
      jsizea_list(11) = jsizea_indcova2;
      jsizea_list(12) = jsizea_indcova2_zi;
      jsizea_list(13) = indcovb_names;
      jsizea_list(14) = jsizea_indcovb2;
      jsizea_list(15) = jsizea_indcovb2_zi;
      jsizea_list(16) = indcovc_names;
      jsizea_list(17) = jsizea_indcovc2;
      jsizea_list(18) = jsizea_indcovc2_zi;
      jsizea_list(19) = year_names;
      jsizea_list(20) = patch_names;
      jsizea_list(21) = group_names;
      jsizea_list(22) = main_effect_1;
      jsizea_list(23) = main_effect_2;
      jsizea_list(24) = jsizea_st;
      jsizea_list(25) = effects_names;
      jsizea_list(26) = jsizea_group1;
      jsizea_list(27) = group1_jsizea_zi;
      jsizea_list(28) = jsizea_indcova1;
      jsizea_list(29) = jsizea_indcovb1;
      jsizea_list(30) = jsizea_indcovc1;
      jsizea_list(31) = jsizea_indcova1_zi;
      jsizea_list(32) = jsizea_indcovb1_zi;
      jsizea_list(33) = jsizea_indcovc1_zi;
      
      List jsizeb_list(34);
      jsizeb_list(0) = jsizeb_num;
      jsizeb_list(1) = jsizeb_year;
      jsizeb_list(2) = jsizeb_patch;
      jsizeb_list(3) = jsizeb_group2;
      jsizeb_list(4) = jsizeb_dist;
      jsizeb_list(5) = zi_yn;
      jsizeb_list(6) = jsizeb_zi;
      jsizeb_list(7) = year_jsizeb_zi;
      jsizeb_list(8) = patch_jsizeb_zi;
      jsizeb_list(9) = group2_jsizeb_zi;
      jsizeb_list(10) = indcova_names;
      jsizeb_list(11) = jsizeb_indcova2;
      jsizeb_list(12) = jsizeb_indcova2_zi;
      jsizeb_list(13) = indcovb_names;
      jsizeb_list(14) = jsizeb_indcovb2;
      jsizeb_list(15) = jsizeb_indcovb2_zi;
      jsizeb_list(16) = indcovc_names;
      jsizeb_list(17) = jsizeb_indcovc2;
      jsizeb_list(18) = jsizeb_indcovc2_zi;
      jsizeb_list(19) = year_names;
      jsizeb_list(20) = patch_names;
      jsizeb_list(21) = group_names;
      jsizeb_list(22) = main_effect_1;
      jsizeb_list(23) = main_effect_2;
      jsizeb_list(24) = jsizeb_st;
      jsizeb_list(25) = effects_names;
      jsizeb_list(26) = jsizeb_group1;
      jsizeb_list(27) = group1_jsizeb_zi;
      jsizeb_list(28) = jsizeb_indcova1;
      jsizeb_list(29) = jsizeb_indcovb1;
      jsizeb_list(30) = jsizeb_indcovc1;
      jsizeb_list(31) = jsizeb_indcova1_zi;
      jsizeb_list(32) = jsizeb_indcovb1_zi;
      jsizeb_list(33) = jsizeb_indcovc1_zi;
      
      List jsizec_list(34);
      jsizec_list(0) = jsizec_num;
      jsizec_list(1) = jsizec_year;
      jsizec_list(2) = jsizec_patch;
      jsizec_list(3) = jsizec_group2;
      jsizec_list(4) = jsizec_dist;
      jsizec_list(5) = zi_yn;
      jsizec_list(6) = jsizec_zi;
      jsizec_list(7) = year_jsizec_zi;
      jsizec_list(8) = patch_jsizec_zi;
      jsizec_list(9) = group2_jsizec_zi;
      jsizec_list(10) = indcova_names;
      jsizec_list(11) = jsizec_indcova2;
      jsizec_list(12) = jsizec_indcova2_zi;
      jsizec_list(13) = indcovb_names;
      jsizec_list(14) = jsizec_indcovb2;
      jsizec_list(15) = jsizec_indcovb2_zi;
      jsizec_list(16) = indcovc_names;
      jsizec_list(17) = jsizec_indcovc2;
      jsizec_list(18) = jsizec_indcovc2_zi;
      jsizec_list(19) = year_names;
      jsizec_list(20) = patch_names;
      jsizec_list(21) = group_names;
      jsizec_list(22) = main_effect_1;
      jsizec_list(23) = main_effect_2;
      jsizec_list(24) = jsizec_st;
      jsizec_list(25) = effects_names;
      jsizec_list(26) = jsizec_group1;
      jsizec_list(27) = group1_jsizec_zi;
      jsizec_list(28) = jsizec_indcova1;
      jsizec_list(29) = jsizec_indcovb1;
      jsizec_list(30) = jsizec_indcovc1;
      jsizec_list(31) = jsizec_indcova1_zi;
      jsizec_list(32) = jsizec_indcovb1_zi;
      jsizec_list(33) = jsizec_indcovc1_zi;
      
      List jrepst_list(34);
      jrepst_list(0) = jrepst_num;
      jrepst_list(1) = jrepst_year;
      jrepst_list(2) = jrepst_patch;
      jrepst_list(3) = jrepst_group2;
      jrepst_list(4) = jrepst_dist;
      jrepst_list(5) = false;
      jrepst_list(6) = dud_zi;
      jrepst_list(7) = dud_zi;
      jrepst_list(8) = dud_zi;
      jrepst_list(9) = dud_zi;
      jrepst_list(10) = indcova_names;
      jrepst_list(11) = jrepst_indcova2;
      jrepst_list(12) = dud_zi;
      jrepst_list(13) = indcovb_names;
      jrepst_list(14) = jrepst_indcovb2;
      jrepst_list(15) = dud_zi;
      jrepst_list(16) = indcovc_names;
      jrepst_list(17) = jrepst_indcovc2;
      jrepst_list(18) = dud_zi;
      jrepst_list(19) = year_names;
      jrepst_list(20) = patch_names;
      jrepst_list(21) = group_names;
      jrepst_list(22) = main_effect_1;
      jrepst_list(23) = main_effect_2;
      jrepst_list(24) = 1.0;
      jrepst_list(25) = effects_names;
      jrepst_list(26) = jrepst_group1;
      jrepst_list(27) = dud_zi;
      jrepst_list(28) = jrepst_indcova1;
      jrepst_list(29) = jrepst_indcovb1;
      jrepst_list(30) = jrepst_indcovc1;
      jrepst_list(31) = dud_zi;
      jrepst_list(32) = dud_zi;
      jrepst_list(33) = dud_zi;
      
      List jmatst_list(34);
      jmatst_list(0) = jmatst_num;
      jmatst_list(1) = jmatst_year;
      jmatst_list(2) = jmatst_patch;
      jmatst_list(3) = jmatst_group2;
      jmatst_list(4) = jmatst_dist;
      jmatst_list(5) = false;
      jmatst_list(6) = dud_zi;
      jmatst_list(7) = dud_zi;
      jmatst_list(8) = dud_zi;
      jmatst_list(9) = dud_zi;
      jmatst_list(10) = indcova_names;
      jmatst_list(11) = jmatst_indcova2;
      jmatst_list(12) = dud_zi;
      jmatst_list(13) = indcovb_names;
      jmatst_list(14) = jmatst_indcovb2;
      jmatst_list(15) = dud_zi;
      jmatst_list(16) = indcovc_names;
      jmatst_list(17) = jmatst_indcovc2;
      jmatst_list(18) = dud_zi;
      jmatst_list(19) = year_names;
      jmatst_list(20) = patch_names;
      jmatst_list(21) = group_names;
      jmatst_list(22) = main_effect_1;
      jmatst_list(23) = main_effect_2;
      jmatst_list(24) = 1.0;
      jmatst_list(25) = effects_names;
      jmatst_list(26) = jmatst_group1;
      jmatst_list(27) = dud_zi;
      jmatst_list(28) = jmatst_indcova1;
      jmatst_list(29) = jmatst_indcovb1;
      jmatst_list(30) = jmatst_indcovc1;
      jmatst_list(31) = dud_zi;
      jmatst_list(32) = dud_zi;
      jmatst_list(33) = dud_zi;
      
      current_surv_model = surv_list;
      current_obs_model = obs_list;
      current_size_model = sizea_list;
      current_sizeb_model = sizeb_list;
      current_sizec_model = sizec_list;
      current_repst_model = repst_list;
      current_fec_model = fec_list;
      
      current_jsurv_model = jsurv_list;
      current_jobs_model = jobs_list;
      current_jsize_model = jsizea_list;
      current_jsizeb_model = jsizeb_list;
      current_jsizec_model = jsizec_list;
      current_jrepst_model = jrepst_list;
      current_jmatst_model = jmatst_list;
      
      current_surv_model.attr("names") = list_names;
      current_obs_model.attr("names") = list_names;
      current_size_model.attr("names") = list_names;
      current_sizeb_model.attr("names") = list_names;
      current_sizec_model.attr("names") = list_names;
      current_repst_model.attr("names") = list_names;
      current_fec_model.attr("names") = list_names;
      current_jsurv_model.attr("names") = list_names;
      current_jobs_model.attr("names") = list_names;
      current_jsize_model.attr("names") = list_names;
      current_jsizeb_model.attr("names") = list_names;
      current_jsizec_model.attr("names") = list_names;
      current_jrepst_model.attr("names") = list_names;
      current_jmatst_model.attr("names") = list_names;
      
      DataFrame c_paramnames = paramnames_skeleton(true);
      CharacterVector modelparams = as<CharacterVector>(c_paramnames["modelparams"]);
      CharacterVector mainparams = as<CharacterVector>(c_paramnames["mainparams"]);
      CharacterVector parameter_names = as<CharacterVector>(c_paramnames["parameter_names"]);
      
      bool current_check = false;
      for (int i = 0; i < modelparams.length(); i++) {
        for (int j = 0; j < 17; j++) {
          current_check = stringcompare_hard(as<std::string>(mainparams(i)), 
            as<std::string>(main_effect_1(j)));
          
          if (current_check) modelparams(i) = main_effect_1(j);
        }
      }
      
      current_paramnames = DataFrame::create(_["parameter_names"] = parameter_names,
        _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
      CharacterVector current_mainyears = as<CharacterVector>(year_list(i));
      unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
      
      CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
      CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
      
      DataFrame indcova2_frame = as<DataFrame>(current_vrm["indcova2_frame"]);
      DataFrame indcovb2_frame = as<DataFrame>(current_vrm["indcovb2_frame"]);
      DataFrame indcovc2_frame = as<DataFrame>(current_vrm["indcovc2_frame"]);
      CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
      CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
      CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
      
      List surv_proxy = LefkoUtils::modelextract(current_surv_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List obs_proxy = LefkoUtils::modelextract(current_obs_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List size_proxy = LefkoUtils::modelextract(current_size_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List sizeb_proxy = LefkoUtils::modelextract(current_sizeb_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List sizec_proxy = LefkoUtils::modelextract(current_sizec_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List repst_proxy = LefkoUtils::modelextract(current_repst_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List fec_proxy = LefkoUtils::modelextract(current_fec_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      
      List jsurv_proxy = LefkoUtils::modelextract(current_jsurv_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List jobs_proxy = LefkoUtils::modelextract(current_jobs_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List jsize_proxy = LefkoUtils::modelextract(current_jsize_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List jsizeb_proxy = LefkoUtils::modelextract(current_jsizeb_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List jsizec_proxy = LefkoUtils::modelextract(current_jsizec_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List jrepst_proxy = LefkoUtils::modelextract(current_jrepst_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      List jmatst_proxy = LefkoUtils::modelextract(current_jmatst_model,
        current_paramnames, current_mainyears, current_mainpatches,
        current_maingroups, current_mainindcova, current_mainindcovb,
        current_mainindcovc, true);
      
      List current_vrm_extract (15);
      current_vrm_extract(0) = surv_proxy;
      current_vrm_extract(1) = obs_proxy;
      current_vrm_extract(2) = size_proxy;
      current_vrm_extract(3) = sizeb_proxy;
      current_vrm_extract(4) = sizec_proxy;
      current_vrm_extract(5) = repst_proxy;
      current_vrm_extract(6) = fec_proxy;
      current_vrm_extract(7) = jsurv_proxy;
      current_vrm_extract(8) = jobs_proxy;
      current_vrm_extract(9) = jsize_proxy;
      current_vrm_extract(10) = jsizeb_proxy;
      current_vrm_extract(11) = jsizec_proxy;
      current_vrm_extract(12) = jrepst_proxy;
      current_vrm_extract(13) = jmatst_proxy;
      current_vrm_extract(14) = current_paramnames;
      
      allmodels_all_pre(i) = current_vrm_extract;
      
    }
    allstages_all = allstages_all_pre;
    allmodels_all = allmodels_all_pre;
  }
  
  //Rcout << "cleanup3 S    ";
  
  // matrix creation for fbMPMs under override option
  if (funcbased) { // Also set this up to check indcovs
    bool check_for_override {true};
    
    if (max(total_years_vec) > prep_mats) check_for_override = false;
    if (max(sp_density_num_vec) > prep_mats) check_for_override = false;
    if (max(dev_terms_num_vec) > prep_mats) check_for_override = false;
    if (max(inda_terms_num_vec) > prep_mats) check_for_override = false;
    if (max(indb_terms_num_vec) > prep_mats) check_for_override = false;
    if (max(indc_terms_num_vec) > prep_mats) check_for_override = false;
    if (max(inda_terms_cat_vec) > prep_mats) check_for_override = false;
    if (max(indb_terms_cat_vec) > prep_mats) check_for_override = false;
    if (max(indc_terms_cat_vec) > prep_mats) check_for_override = false;
    
    IntegerVector found_maxes = {max(total_years_vec), max(sp_density_num_vec),
      max(dev_terms_num_vec), max(inda_terms_num_vec), max(indb_terms_num_vec),
      max(indc_terms_num_vec), max(inda_terms_cat_vec), max(indb_terms_cat_vec),
      max(indc_terms_cat_vec)};
    
    IntegerVector unique_maxes = sort_unique(found_maxes);
    if (static_cast<int>(unique_maxes.length()) > 2) {
      check_for_override = false;
    } else if (static_cast<int>(unique_maxes.length()) == 2) {
      if (unique_maxes(0) != 0) check_for_override = false;
    }
    
    if (check_for_override && !force_fb && !stochastic) fb_override = true;
  }
  
  //Rcout << "cleanup3 T    ";
  //Rcout << "vrm_count: " << vrm_count << " ";
  
  List errcheck_mpmout_vrm (vrm_count);
  
  if (fb_override) {
    // Create function-based MPMs and assign them to mpm_list
    //Rcout << "Entered fb_override    ";
    
    List mpm_list_pre (vrm_count);
    List A_list_pre (vrm_count);
    
    IntegerVector year_counter (vrm_count);
    
    for (int i = 0; i < vrm_count; i++) {
      
      //Rcout << "cleanup3 T1 i: " << i << "    ";
      
      List current_vrm_extract = as<List>(allmodels_all(i));
      List current_vrm_unextract = as<List>(vrm_list(i));
      DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
      
      int ehrlen_format {1}; // This will need to be dealt with differently later
      
      int mpm_style {1};
      if (format_vec(i) < 3) {
        mpm_style = 0;
      } else if (format_vec(i) == 4) {
        mpm_style = 2;
      }
      
      DataFrame current_mpm_allstages = as<DataFrame>(allstages_all(i));
      
      List surv_proxy = as<List>(current_vrm_extract(0));
      List obs_proxy = as<List>(current_vrm_extract(1));
      List size_proxy = as<List>(current_vrm_extract(2));
      List sizeb_proxy = as<List>(current_vrm_extract(3));
      List sizec_proxy = as<List>(current_vrm_extract(4));
      List repst_proxy = as<List>(current_vrm_extract(5));
      List fec_proxy = as<List>(current_vrm_extract(6));
      List jsurv_proxy = as<List>(current_vrm_extract(7));
      List jobs_proxy = as<List>(current_vrm_extract(8));
      List jsize_proxy = as<List>(current_vrm_extract(9));
      List jsizeb_proxy = as<List>(current_vrm_extract(10));
      List jsizec_proxy = as<List>(current_vrm_extract(11));
      List jrepst_proxy = as<List>(current_vrm_extract(12));
      List jmatst_proxy = as<List>(current_vrm_extract(13));
      DataFrame current_paramnames = as<DataFrame>(current_vrm_extract(14));
      
      //Rcout << "cleanup3 T2    ";
      
      CharacterVector current_mainyears = as<CharacterVector>(year_list(i));
      unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
      
      //Rcout << "cleanup3 current_mainyears: " << current_mainyears << "          ";
      
      DataFrame group2_frame = as<DataFrame>(current_vrm_unextract["group2_frame"]);
      CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
      
      //CharacterVector current_chosenpatches = patch_vec(i);
      //unsigned int no_chosenpatches = static_cast<unsigned int>(current_chosenpatches.length());
      DataFrame patch_frame = as<DataFrame>(current_vrm_unextract["patch_frame"]);
      CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
      
      //Rcout << "cleanup3 current_mainpatches: " << current_mainpatches << "          ";
      
      DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
      DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
      DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
      CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
      CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
      CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
      
      //Rcout << "cleanup3 T3    ";
      
      //int year_counter {0};
      //int patch_counter {0};
      IntegerVector inda_num_terms_counter (vrm_count);
      IntegerVector indb_num_terms_counter (vrm_count);
      IntegerVector indc_num_terms_counter (vrm_count);
      IntegerVector inda_cat_terms_counter (vrm_count);
      IntegerVector indb_cat_terms_counter (vrm_count);
      IntegerVector indc_cat_terms_counter (vrm_count);
      IntegerVector inda_num_terms_previous (vrm_count);
      IntegerVector indb_num_terms_previous (vrm_count);
      IntegerVector indc_num_terms_previous (vrm_count);
      IntegerVector inda_cat_terms_previous (vrm_count);
      IntegerVector indb_cat_terms_previous (vrm_count);
      IntegerVector indc_cat_terms_previous (vrm_count);
      IntegerVector dev_num_counter (vrm_count);
      IntegerVector sp_density_counter (vrm_count);
      
      //Rcout << "cleanup3 T4    ";
      
      IntegerVector found_calls = {total_years_vec(i), sp_density_num_vec(i),
        dev_terms_num_vec(i), inda_terms_num_vec(i), indb_terms_num_vec(i),
        indc_terms_num_vec(i), inda_terms_cat_vec(i), indb_terms_cat_vec(i),
        indc_terms_cat_vec(i)};
      int needed_calls = max(found_calls);
      
      List current_building_mpm (needed_calls);
      List errcheck_mpmout_vrm_time (needed_calls);
      CharacterVector labels_year2_terms (needed_calls);
      CharacterVector labels_patch_terms (needed_calls);
      
      //Rcout << "cleanup3 T5    ";
      
      for (int j = 0; j < needed_calls; j++) { // time loop
        
        //Rcout << "cleanup3 T6 j:" << j << "    ";
        
        // Counter resets
        if (year_counter(i) == no_mainyears || j == 0) year_counter(i) = 0;
        //if (patch_counter == no_chosenpatches) patch_counter = 0;
        int yearnumber = year_counter(i);
        int patchnumber {0};
        
        CharacterVector current_year = as<CharacterVector>(current_mainyears(year_counter(i)));
        for (int z = 0; z < static_cast<int>(current_mainpatches.length()); z++) {
          if (LefkoUtils::stringcompare_simple(String(patch_vec(i)),
              String(current_mainpatches(z)), false)) patchnumber = z;
        }
        //patch_counter++;
        //Rcout << "cleanup3 yearnumber: " << yearnumber << "          ";
        //Rcout << "cleanup3 patchnumber: " << patchnumber << "          ";
        
        if (inda_num_terms_counter(i) >= inda_terms_num_vec(i)) inda_num_terms_counter(i) = 0;
        if (indb_num_terms_counter(i) >= indb_terms_num_vec(i)) indb_num_terms_counter(i) = 0;
        if (indc_num_terms_counter(i) >= indc_terms_num_vec(i)) indc_num_terms_counter(i) = 0;
        if (inda_cat_terms_counter(i) >= inda_terms_cat_vec(i)) inda_cat_terms_counter(i) = 0;
        if (indb_cat_terms_counter(i) >= indb_terms_cat_vec(i)) indb_cat_terms_counter(i) = 0;
        if (indc_cat_terms_counter(i) >= indc_terms_cat_vec(i)) indc_cat_terms_counter(i) = 0;
        
        List current_ind_terms_num = ind_terms_num_list(i);
        List current_ind_terms_cat = ind_terms_cat_list(i);
        
        NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
        NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
        NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
        CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
        CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
        CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
        
        NumericVector f2_inda = {f_inda_full(inda_num_terms_counter(i))};
        NumericVector f1_inda = {f_inda_full(inda_num_terms_previous(i))};
        NumericVector f2_indb = {f_indb_full(indb_num_terms_counter(i))};
        NumericVector f1_indb = {f_indb_full(indb_num_terms_previous(i))};
        NumericVector f2_indc = {f_indc_full(indc_num_terms_counter(i))};
        NumericVector f1_indc = {f_indc_full(indc_num_terms_previous(i))};
        CharacterVector r2_inda = as<CharacterVector>(r_inda_full(inda_cat_terms_counter(i)));
        CharacterVector r1_inda = as<CharacterVector>(r_inda_full(inda_cat_terms_previous(i)));
        CharacterVector r2_indb = as<CharacterVector>(r_indb_full(indb_cat_terms_counter(i)));
        CharacterVector r1_indb = as<CharacterVector>(r_indb_full(indb_cat_terms_previous(i)));
        CharacterVector r2_indc = as<CharacterVector>(r_indc_full(indc_cat_terms_counter(i)));
        CharacterVector r1_indc = as<CharacterVector>(r_indc_full(indc_cat_terms_previous(i)));
        
        NumericVector dv_terms (14);
        if (dev_terms_num_vec(i) > 0) {
          DataFrame used_dv_df = as<DataFrame>(dev_terms_list(i));
        
          if (dev_num_counter(i) >= dev_terms_num_vec(i)) dev_num_counter(i) = 0;
          
          for (int j = 0; j < 14; j++) {
            dv_terms(j) = used_dv_df(dev_num_counter(i), j);
          }
          dev_num_counter(i) = dev_num_counter(i) + 1;
        }
        
        bool dvr_bool {false};
        
        LogicalVector dvr_yn = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        IntegerVector dvr_style (14);
        IntegerVector dvr_time_delay (14);
        NumericVector dvr_alpha (14);
        NumericVector dvr_beta (14);
        
        if (dens_vr_yn_vec(i) > 0) {
          dvr_bool = true;
          
          DataFrame current_dvr = as<DataFrame>(density_vr_list(i));
          LogicalVector true_dvr_yn = as<LogicalVector>(current_dvr(1));
          IntegerVector true_dvr_style = as<IntegerVector>(current_dvr(2));
          IntegerVector true_dvr_time_delay = as<IntegerVector>(current_dvr(3));
          NumericVector true_dvr_alpha = as<NumericVector>(current_dvr(4));
          NumericVector true_dvr_beta = as<NumericVector>(current_dvr(5));
          
          dvr_yn = true_dvr_yn;
          dvr_style = true_dvr_style;
          dvr_time_delay = true_dvr_time_delay;
          dvr_alpha = true_dvr_alpha;
          dvr_beta = true_dvr_beta;
        }
        
        double maxsize {0.0};
        double maxsizeb {0.0};
        double maxsizec {0.0};
        
        if (format_vec(i) < 5) {
          DataFrame current_allstages = as<DataFrame>(allstages_all(i));
          
          NumericVector size3 = as<NumericVector>(current_allstages["size3"]);
          NumericVector size2n = as<NumericVector>(current_allstages["size2n"]);
          NumericVector size2o = as<NumericVector>(current_allstages["size2o"]);
          NumericVector sizeb3 = as<NumericVector>(current_allstages["sizeb3"]);
          NumericVector sizeb2n = as<NumericVector>(current_allstages["sizeb2n"]);
          NumericVector sizeb2o = as<NumericVector>(current_allstages["sizeb2o"]);
          NumericVector sizec3 = as<NumericVector>(current_allstages["sizec3"]);
          NumericVector sizec2n = as<NumericVector>(current_allstages["sizec2n"]);
          NumericVector sizec2o = as<NumericVector>(current_allstages["sizec2o"]);
          
          NumericVector maxveca = {max(size3), max(size2n), max(size2o)};
          NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
          NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
          
          maxsize = max(maxveca);
          maxsizeb = max(maxvecb);
          maxsizec = max(maxvecc);
        }
        
        double dens_sp {1.0};
        
        if (sp_density_num_vec(i) > 0) {
          if (sp_density_counter(i) >= sp_density_num_vec(i)) sp_density_counter(i) = 0;
          
          NumericVector current_sp_density = as<NumericVector>(sp_density_list(i));
          dens_sp = current_sp_density(sp_density_counter(i));
          
          sp_density_counter(i) = sp_density_counter(i) + 1;
        }
        
        NumericVector dens_n (14, 1.0); // This needs to be updated with the actual pop size
        
        List current_mpm;
        if (format_vec(i) < 5) {
          current_mpm = AdaptMats::mazurekd(current_mpm_allstages,
            current_stageframe, format_vec(i), surv_proxy, obs_proxy,
            size_proxy, sizeb_proxy, sizec_proxy, repst_proxy, fec_proxy,
            jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy, jsizec_proxy,
            jrepst_proxy, jmatst_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
            f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
            r1_indc, dv_terms, dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
            dens_n, dens_sp, fecmod_vec(i), maxsize, maxsizeb, maxsizec,
            static_cast<unsigned int>(firstage_vec(i)),
            static_cast<unsigned int>(finalage_vec(i)), true, yearnumber,
            patchnumber, exp_tol, theta_tol, true, err_check, sparse_vec(i),
            true);
            
          if (err_check) {
            NumericMatrix mpm_out = as<NumericMatrix>(current_mpm["out"]);
            errcheck_mpmout_vrm_time(j) = mpm_out; 
          }
        } else {
          IntegerVector all_ages = seq(firstage_vec(i), finalage_vec(i));
          DataFrame current_supplement;
          if (supplement_list(i) == R_NilValue) {
            current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
              surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
              f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
              r1_indc, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(i),
              finalage_vec(i), true, yearnumber, patchnumber, dvr_bool, dvr_yn,
              dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
              sparse_vec(i));
            
          } else {
            current_supplement = as<DataFrame>(supplement_list(i));
            
            current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
              surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
              f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
              r1_indc, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(i),
              finalage_vec(i), true, yearnumber, patchnumber, dvr_bool, dvr_yn,
              dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
              sparse_vec(i), current_supplement);
          }
        }
        arma::mat pulled_matrix = as<arma::mat>(current_mpm(0));
        current_building_mpm(j) = pulled_matrix;
        
        labels_year2_terms(j) = current_year(0);
        labels_patch_terms(j) = current_mainpatches(patchnumber);
        
        year_counter(i) = year_counter(i) + 1;
      } // time loop (j)
      if (err_check) errcheck_mpmout_vrm(i) = errcheck_mpmout_vrm_time;
      //Rcout << "cleanup3 T7    ";
      
      List A_mats = current_building_mpm;
      DataFrame current_labels = DataFrame::create(_["patch"] = labels_patch_terms,
        _["year2"] = labels_year2_terms);
      
      List trial_mpm = List::create(_["A"] = A_mats, _["labels"] = current_labels);
      
      A_list_pre(i) = A_mats;
      
      mpm_list_pre(i) = trial_mpm;
      
    } // vrm loop (i)
    A_list = A_list_pre;
    mpm_list = mpm_list_pre;
    mpm_count = vrm_count;
  }
  
  //Rcout << "cleanup3 U    ";
  
  // Output processing
  List out_list (72);
  out_list(0) = mpm_list;
  out_list(1) = mpm_count;
  out_list(2) = vrm_list;
  out_list(3) = vrm_count;
  out_list(4) = A_list;
  out_list(5) = stageframe_list;
  out_list(6) = stageframe_list_fb;
  out_list(7) = stageframe_count;
  out_list(8) = supplement_list;
  out_list(9) = supplement_list_fb;
  out_list(10) = supplement_count;
  out_list(11) = repmatrix_list;
  out_list(12) = sparse_vec;
  out_list(13) = sparse_vec_count;
  out_list(14) = format_vec;
  out_list(15) = found_fleslie;
  out_list(16) = stageframe_notNull_count;
  out_list(17) = preexisting;
  out_list(18) = funcbased;
  out_list(19) = pure_fleslie;
  out_list(20) = firstage_vec;
  out_list(21) = finalage_vec;
  out_list(22) = cont_vec;
  out_list(23) = fecmod_vec;
  out_list(24) = fecage_min_vec;
  out_list(25) = fecage_max_vec;
  out_list(26) = hstages_list;
  out_list(27) = agestages_list;
  out_list(28) = matrowcounts;
  out_list(29) = stagecounts;
  out_list(30) = start_list;
  out_list(31) = start_count;
  out_list(32) = labels_list;
  out_list(33) = labels;
  out_list(34) = patch_vec;
  out_list(35) = year_list;
  out_list(36) = total_years_vec;
  out_list(37) = tweights_list;
  out_list(38) = tweights_count;
  out_list(39) = tweights_type_vec;
  out_list(40) = density_list;
  out_list(41) = dens_index_list;
  out_list(42) = dens_yn_vec;
  out_list(43) = density_count;
  out_list(44) = entry_time_vec;
  out_list(45) = entry_time_count;
  out_list(46) = entry_time_vec_use;
  out_list(47) = density_vr_list;
  out_list(48) = ind_terms_num_list;
  out_list(49) = ind_terms_cat_list;
  out_list(50) = dev_terms_list;
  out_list(51) = dens_vr_yn_vec;
  out_list(52) = sp_density_num_vec;
  out_list(53) = dev_terms_num_vec;
  out_list(54) = inda_terms_num_vec;
  out_list(55) = indb_terms_num_vec;
  out_list(56) = indc_terms_num_vec;
  out_list(57) = inda_terms_cat_vec;
  out_list(58) = indb_terms_cat_vec;
  out_list(59) = indc_terms_cat_vec;
  out_list(60) = sparse_vec;
  out_list(61) = density_vr_count;
  out_list(62) = sparse_vec_count;
  out_list(63) = sp_density_list;
  out_list(64) = equivalence_list;
  out_list(65) = equivalence_vec;
  out_list(66) = equivalence_count;
  out_list(67) = stages_not_equal;
  out_list(68) = allstages_all;
  out_list(69) = allmodels_all;
  out_list(70) = fb_override;
  out_list(71) = errcheck_mpmout_vrm;
  
  CharacterVector out_list_names = {"mpm_list", "mpm_count", "vrm_list",
    "vrm_count", "A_list", "stageframe_list", "stageframe_list_fb",
    "stageframe_count", "supplement_list", "supplement_list_fb",
    "supplement_count", "repmatrix_list", "sparse_vec", "sparse_vec_count",
    "format_vec", "found_fleslie", "stageframe_notNull_count", "preexisting",
    "funcbased", "pure_fleslie", "firstage_vec", "finalage_vec", "cont_vec",
    "fecmod_vec", "fecage_min_vec", "fecage_max_vec", "hstages_list",
    "agestages_list", "matrowcounts", "stagecounts", "start_list",
    "start_count", "labels_list", "labels", "patch_vec", "year_list",
    "total_years_vec", "tweights_list", "tweights_count", "tweights_type_vec",
    "density_list", "dens_index_list", "dens_yn_vec", "density_count",
    "entry_time_vec", "entry_time_count", "entry_time_vec_use",
    "density_vr_list", "ind_terms_num_list", "ind_terms_cat_list",
    "dev_terms_list", "dens_vr_yn_vec", "sp_density_num_vec",
    "dev_terms_num_vec", "inda_terms_num_vec", "indb_terms_num_vec",
    "indc_terms_num_vec", "inda_terms_cat_vec", "indb_terms_cat_vec",
    "indc_terms_cat_vec", "sparse_vec", "density_vr_count", "sparse_vec_count",
    "sp_density_list", "equivalence_list", "equivalence_vec",
    "equivalence_count", "stages_not_equal", "allstages_all", "allmodels_all",
    "fb_override", "errcheck_mpmout_vrm"};
  out_list.attr("names") = out_list_names;
  
  return(out_list);
}

//' Engine Projecting Multiple Existing MPMs With or Without Density Dependence
//' 
//' Function \code{project3_pre_core} is the main function running processing of
//' projections with existing MPMs supplied by the user.
//' 
//' @name project3_pre_core
//' 
//' @param N_out The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param extreme_mpm_list A multi-level list output if
//' \code{err_check = "extreme"}.
//' @param mpm_list A list of MPMs in \code{lefkoMat} format.
//' @param A_list A list of lists of \code{A} matrices covering all entered
//' MPMs.
//' @param tweights_list A list of tweights vectors covering all MPMs.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param vrm_list A list of \code{vrm_input} objects.
//' @param stageframe_list A list of stageframe objects covering all MPMs.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param allstages_all A list of data frames giving the \code{allstages}
//' stage expansions giving indices for all matrix elements.
//' @param supplement_list A list of supplements in \code{lefkoSD} format.
//' @param year_list A list of vectors giving the main years used in each MPM.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates for each MPM.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates for each MPM.
//' @param dev_terms_list List of deviations for vital rate models in all MPMs.
//' @param density_vr_list List of \code{lefkoDensVR} objects holding density
//' relationships for all 14 vital rate models, for all MPMs.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param density_list A list of data frames of class \code{lefkoDens} for all
//' MPMs.
//' @param dens_index_list A list of data frames giving indices for density
//' dependent transitions for each MPM.
//' @param equivalence_list An optional numeric vector, list of numeric vectors,
//' data frame of class \code{adaptEq}, or list of data frames of class
//' \code{adaptEq}. If a numeric vector, then must have the same number of
//' elements as the number of MPMs, with each element giving the effect of an
//' individual of each MPM relative to a reference individual. If a list of
//' vectors, then the list should be composed of as many numeric vectors as
//' MPMs, with each vector giving the effect of each individual in each stage
//' relative to a reference individual. Data frames of class \code{adaptEq}, and
//' lists of such data frames, can be made with function
//' \code{\link{equiv_input}()}. Numeric entries used in these vectors can be
//' thought of as Lotka-Volterra interaction terms, such as are used in multiple
//' species competition models.
//' @param dev_terms_num_vec A vector giving the number of vital rate deviations
//' in each MPM. Used to create function-based matrices under override option.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms per MPM. Used to create function-based matrices under override option.
//' @param firstage_vec An integer vector containing the first age values for
//' all MPMs in order.
//' @param finalage_vec  An integer vector containing the final age values for
//' all MPMs in order.
//' @param stagecounts An integer vector containing the number of stages in each
//' MPM.
//' @param entry_time_vec An integer vector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param format_vec An integer vector giving the MPM format of each MPM.
//' @param inda_terms_num_vec An integer vector giving the number of numeric
//' terms given in individual covariate a.
//' @param indb_terms_num_vec An integer vector giving the number of numeric
//' terms given in individual covariate b.
//' @param indc_terms_num_vec An integer vector giving the number of numeric
//' terms given in individual covariate c.
//' @param inda_terms_cat_vec An integer vector giving the number of factor
//' terms given in individual covariate a.
//' @param indb_terms_cat_vec An integer vector giving the number of factor
//' terms given in individual covariate b.
//' @param indc_terms_cat_vec An integer vector giving the number of factor
//' terms given in individual covariate c.
//' @param dens_yn_vec A vector stating whether density dependence is used in
//' each MPM, given through \code{lefkoDens} objects.
//' @param dens_vr_yn_vec A vector stating whether density dependence is used in
//' each MPM, given through \code{lefkoDensVR} objects.
//' @param tweights_type_vec An integer vector giving the style of
//' \code{tweights} used in each MPM.
//' @param fecmod_vec A numeric vector giving the fecmod values for all MPMs.
//' @param sparse_vec A logical vector stating whether each MPM is in sparse
//' matrix format or not.
//' @param patch_vec A string vector giving the name of each patch used in each
//' MPM projection.
//' @param vrm_count An integer giving the number of vrm_inputs supplied.
//' @param mpm_count An integer giving the number of existing MPMs supplied.
//' @param nreps An integer giving the number of replicates to perform.
//' @param times An integer giving the amount of time steps to run the
//' projection for.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param stochastic A Boolean value indicating to perform a temporally
//' stochastic projection.
//' @param entry_time_vec_use A Boolean value indicating whether entry
//' times differ among MPMs.
//' @param err_check A Boolean value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A Boolean value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' 
//' @return The first three arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
void project3_pre_core (List& N_out, List& comm_out, List& extreme_mpm_out,
  const List mpm_list, const List A_list, const List tweights_list,
  const List start_list, const List vrm_list, const List stageframe_list,
  const List allmodels_all, const List allstages_all, const List supplement_list,
  const List year_list, const List ind_terms_num_list,
  const List ind_terms_cat_list, const List dev_terms_list,
  const List density_vr_list, const List sp_density_list,
  const List density_list, const List dens_index_list,
  const List equivalence_list, const IntegerVector dev_terms_num_vec,
  const IntegerVector sp_density_num_vec, const IntegerVector firstage_vec,
  const IntegerVector finalage_vec, const IntegerVector stagecounts,
  const IntegerVector entry_time_vec, const IntegerVector format_vec,
  const IntegerVector inda_terms_num_vec, const IntegerVector indb_terms_num_vec,
  const IntegerVector indc_terms_num_vec, const IntegerVector inda_terms_cat_vec,
  const IntegerVector indb_terms_cat_vec, const IntegerVector indc_terms_cat_vec,
  const IntegerVector dens_yn_vec, const IntegerVector dens_vr_yn_vec,
  const IntegerVector tweights_type_vec, const NumericVector fecmod_vec,
  const LogicalVector sparse_vec, const CharacterVector patch_vec,
  const int vrm_count, const int mpm_count, const int nreps, const int times,
  const int substoch, const double exp_tol, const double theta_tol,
  const bool integeronly, const bool stages_not_equal, const bool stochastic,
  const bool entry_time_vec_use, const bool err_check,
  const bool err_check_extreme) {
  
  // start_list
  
  // Matrix order set up and creation of zero stage vectors
  //Rcout << "Entered project3_pre_core          ";
  
  List matrix_choice_list (mpm_count);
  List zvl (mpm_count);
  for (int i = 0; i < mpm_count; i++) {
    
    //Rcout << "project3_pre_core a  i (mpm): " << i << "          ";
    
    IntegerVector chosen_matrix_vec;
    
    List chosen_mpm = as<List>(mpm_list(i));
    DataFrame chosen_labels = as<DataFrame>(chosen_mpm["labels"]);
    CharacterVector chosen_labels_names = as<CharacterVector>(chosen_labels.attr("names"));
    IntegerVector clm_y2 = index_l3(chosen_labels_names, "year2");
    if (static_cast<int>(clm_y2.length()) == 0) {
      throw Rcpp::exception("Value for argument year not found.", false);
    }
    
    CharacterVector mpm_labels_patch = as<CharacterVector>(chosen_labels["patch"]);
    IntegerVector pvis = index_l3(mpm_labels_patch, patch_vec(i));
    //Rcout << "mpm_labels_patch: " << mpm_labels_patch << endl;
    //Rcout << "pvis (mpm labels chosen for patch match): " << pvis << endl;
    
    if (static_cast<int>(pvis.length()) == 0) {
      throw Rcpp::exception("Value for argument patch not found.", false);
    }
    
    IntegerVector chosen_mats;
    
    if (clm_y2.length() > 0) {
      CharacterVector mpm_labels_year2 = as<CharacterVector>(chosen_labels["year2"]);
      
      CharacterVector chosen_years = as<CharacterVector>(year_list(i));
      int chosen_years_length = static_cast<int>(chosen_years.length());
      
      IntegerVector yvis = index_l3(mpm_labels_year2, chosen_years(0)); // Do we need to drop the 0?
      if (chosen_years_length > 1) {
        for (int j = 1; j < chosen_years_length; j++) {
          IntegerVector yvis_next = index_l3(mpm_labels_year2, chosen_years(j));
          IntegerVector yvis_append = concat_int(yvis, yvis_next);
          
          yvis = sort_unique(yvis_append);
        }
      }
      
      IntegerVector chosen_mats_pre = intersect(pvis, yvis);
      chosen_mats = chosen_mats_pre;
      
      //Rcout << "mpm_labels_year2: " << mpm_labels_year2 << endl;
      //Rcout << "chosen_years: " << chosen_years << endl;
      //Rcout << "yvis: " << yvis << endl;
      //Rcout << "chosen_mats_pre: " << chosen_mats_pre << endl;
    } else {
      IntegerVector chosen_mats_pre = sort_unique(pvis);
      chosen_mats = chosen_mats_pre;
      //Rcout << "No years found, so going with all associated with patch" << endl;
      //Rcout << "chosen_mats_pre: " << chosen_mats_pre << endl;
    }
    
    //Rcout << "project3_pre_core b    ";
    
    matrix_choice_list(i) = chosen_mats;
  }
  
  // Year order determination
  List comm_out_pre (mpm_count);
  List used_times (mpm_count);
  
  //Rcout << "project3_pre_core c    ";
  
  for (int i = 0; i < mpm_count; i++) {
    List pop_reps (nreps);
    
    arma::mat pops_out_pre (stagecounts(i), (times + 1), fill::zeros);
    
    IntegerVector chosen_mats = as<IntegerVector>(matrix_choice_list(i));
    int chosen_mats_length = static_cast<int>(chosen_mats.length());
    
    List used_times_mpm (nreps);
    
    for (int j = 0; j < nreps; j++) {
      IntegerVector years_topull;
      
      if (!stochastic) {
        IntegerVector years_topull_pre (times);
        
        int mat_tracker {0};
        for (int k = 0; k < times; k++) {
          if (k >= entry_time_vec(i)) {
            if (mat_tracker >= chosen_mats_length) mat_tracker = 0;
            
            years_topull_pre(k) = chosen_mats(mat_tracker);
            mat_tracker++;
          }
        }
        years_topull = years_topull_pre;
        
      } else {
        if (tweights_type_vec(i) == 0) {
          NumericVector twinput (chosen_mats_length,
            (1.0 / static_cast<double>(chosen_mats_length)));
          
          IntegerVector years_topull_pre = Rcpp::RcppArmadillo::sample(chosen_mats,
            (times - entry_time_vec(i)), true, twinput);
          
          IntegerVector years_topull_almost (times);
          for (int copy_elem = 0; copy_elem < (times - entry_time_vec(i)); copy_elem++) {
            years_topull_almost(copy_elem + entry_time_vec(i)) = years_topull_pre(copy_elem);
          }
          years_topull = years_topull_almost;
        } else if (tweights_type_vec(i) == 1) {
          NumericVector twinput = as<NumericVector>(tweights_list(i));
          NumericVector twinput_st = twinput / sum(twinput);
          
          IntegerVector years_topull_pre = Rcpp::RcppArmadillo::sample(chosen_mats,
            (times - entry_time_vec(i)), true, twinput_st);
          
          IntegerVector years_topull_almost (times);
          for (int copy_elem = 0; copy_elem < (times - entry_time_vec(i)); copy_elem++) {
            years_topull_almost(copy_elem + entry_time_vec(i)) = years_topull_pre(copy_elem);
          }
          years_topull = years_topull_almost;
        } else if (tweights_type_vec(i) == 2) {
          arma::ivec chosen_mats_arma = as<arma::ivec>(chosen_mats);
          arma::mat twinput_mat = as<arma::mat>(tweights_list(i));
          arma::vec twinput = twinput_mat.col(0);
          twinput = twinput / sum(twinput);
          
          IntegerVector years_topull_pre (times);
          NumericVector twinput_setup (chosen_mats_length, (1.0 / 
            static_cast<double>(chosen_mats_length)));
          arma::ivec first_choice = Rcpp::RcppArmadillo::sample(chosen_mats_arma,
            times, true, twinput_setup);
          years_topull_pre(entry_time_vec(i)) = chosen_mats(first_choice(0));
          
          for (int k = 1; k < (times - entry_time_vec(i)); k++) {
            arma::ivec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(chosen_mats_arma,
              1, true, twinput);
            years_topull_pre(k + entry_time_vec(i)) = theprophecy_piecemeal(0);
              
            arma::uvec tnotb_preassigned = find(chosen_mats_arma == theprophecy_piecemeal(0));
            twinput = twinput_mat.col(static_cast<int>(tnotb_preassigned(0)));
            twinput = twinput / sum(twinput);
          }
          years_topull = years_topull_pre;
          
        } else {
          throw Rcpp::exception("tweights_type_vec error.", false);
        }
      }
      used_times_mpm(j) = years_topull;
      pop_reps(j) = pops_out_pre;
    }
    used_times(i) = used_times_mpm;
    comm_out_pre(i) = pop_reps;
  }
  
  //Rcout << "project3_pre_core d    ";
  
  // Main projection
  List N_out_pre (nreps);
  List extreme_mpm_reps (nreps);
  
  for (int i = 0; i < nreps; i ++) {
    List running_popvecs = clone(start_list);
    NumericMatrix N_mpm (mpm_count, (times + 1));
    List extreme_mpm_reps_times (times);
    
    //Rcout << "project3_pre_core e i (rep): " << i << "      ";
    
    for (int j = 0; j < times; j++) {
      if (j % 10 == 0){
        Rcpp::checkUserInterrupt();
      }
      
      //Rcout << "project3_pre_core f j (time): " << j << "      ";
      
      List extreme_mpm_reps_times_mpms (mpm_count);
      
      for (int k = 0; k < mpm_count; k++) {
        //Rcout << "project3_pre_core g k (mpm): " << k << "      ";
        
        List reps_out = comm_out_pre(k);
        arma::mat pops_out = as<arma::mat>(reps_out(i));
        
        if (j > (entry_time_vec(k) - 1)) {
          List current_used_times_mpm = as<List>(used_times(k));
          IntegerVector current_times_vec = as<IntegerVector>(current_used_times_mpm(i));
          
          arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs(k));
          
          if (j == entry_time_vec(k)) {
            pops_out.col(j) = running_popvec_mpm;
            
            double N_current = accu(running_popvec_mpm);
            N_mpm(k, j) = N_current;
          }
          
          List chosen_As = as<List>(A_list(k));
          
          if (dens_yn_vec(k) == 0) {
            //Rcout << "project3_pre_core h no density       ";
            
            if (sparse_vec(k) == 0) {
              arma::mat current_A = as<arma::mat>(chosen_As(current_times_vec(j)));
              running_popvec_mpm = current_A * running_popvec_mpm;
              if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = current_A;
              
            } else {
              arma::sp_mat current_A = as<arma::sp_mat>(chosen_As(current_times_vec(j)));
              running_popvec_mpm = current_A * running_popvec_mpm;
              if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = current_A;
            }
          } else {
            //Rcout << "project3_pre_core i density       ";
            
            DataFrame used_density_input = as<DataFrame>(density_list(k));
            DataFrame used_density_index_input = as<DataFrame>(dens_index_list(k));
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            if (j > used_delay) { // Change to allow different delay Ns for different entries
              if (!stages_not_equal) {
                //Rcout << "project3_pre_core j density stages equal       ";
                
                NumericVector delay_issue = N_mpm(_, (j + 1 - used_delay));
                
                double delay_N_sum = sum(delay_issue);
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, chosen_As, delay_N_sum,
                  static_cast<int>(current_times_vec(j)), integeronly,
                  substoch, used_density_input, used_density_index_input,
                  false, static_cast<bool>(sparse_vec(k)),
                  static_cast<bool>(sparse_vec(k)), false, err_check_extreme);
                
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = new_projmat;
                
              } else {
                //Rcout << "project3_pre_core k density stages not equal       ";
                
                double delay_N_sum {0.0};
                
                if (j > 0) {
                  for (int l = 0; l < mpm_count; l++) {
                    List current_pop_list = as<List>(comm_out_pre(l));
                    arma::mat delay_pop = as<arma::mat>(current_pop_list(i));
                    arma::vec delay_pop_vec = delay_pop.col(j + 1 - used_delay);
                    arma::vec current_equiv_vec = as<arma::vec>(equivalence_list(l));
                    arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                    double delay_pop_N = accu(adjusted_delay_pop_vec);
                    
                    delay_N_sum += delay_pop_N;
                  }
                }
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, chosen_As, delay_N_sum,
                  static_cast<int>(current_times_vec(j)), integeronly,
                  substoch, used_density_input, used_density_index_input,
                  false, static_cast<bool>(sparse_vec(k)),
                  static_cast<bool>(sparse_vec(k)), false, err_check_extreme);
                
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = new_projmat;
              }
            } else {
              //Rcout << "project3_pre_core l density initial time       ";
              
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              proj3dens_ad(new_popvec, new_projmat, new_projsp,
                running_popvec_mpm, chosen_As, 0.0,
                static_cast<int>(current_times_vec(j)), integeronly, substoch,
                used_density_input, used_density_index_input, false,
                static_cast<bool>(sparse_vec(k)),
                static_cast<bool>(sparse_vec(k)), false, err_check_extreme);
              
              running_popvec_mpm = new_popvec;
              if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = new_projmat;
            }
          }
          
          if (integeronly) running_popvec_mpm = floor(running_popvec_mpm);
          double N_current = accu(running_popvec_mpm);
          N_mpm(k, (j + 1)) = N_current;
          
          running_popvecs(k) = running_popvec_mpm;
          pops_out.col(j + 1) = running_popvec_mpm;
        }
        
        reps_out(i) = pops_out;
        comm_out_pre(k) = reps_out;
      }
      if (err_check_extreme) extreme_mpm_reps_times(j) = extreme_mpm_reps_times_mpms;
    }
    comm_out = comm_out_pre;
    N_out_pre(i) = N_mpm;
    if (err_check_extreme) extreme_mpm_reps(i) = extreme_mpm_reps_times;
  }
  N_out = N_out_pre;
  if (err_check_extreme) extreme_mpm_out = extreme_mpm_reps;
}

//' Engine Projecting Multiple Function-based MPMs With or Without Density Dependence
//' 
//' Function \code{project3_fb_core} is the main function running processing of
//' projections in which matrices must be created at each time step.
//' 
//' @name project3_fb_core
//' 
//' @param N_out The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param extreme_mpm_out A multi-level list output if
//' \code{err_check = "extreme"}.
//' @param errcheck_fb_out A multi-level list of out matrices from matrix
//' calculation, output if \code{err_check = TRUE}.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param vrm_list A list of unextracted \code{vrm_input} objects holding model
//' information for all MPMs to be created.
//' @param tweights_list A list of tweights vectors covering all MPMs.
//' @param stageframe_list A list of stageframe objects covering all MPMs.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param allstages_all A list of data frames giving the \code{allstages}
//' stage expansions giving indices for all matrix elements.
//' @param supplement_list A list of supplements in \code{lefkoSD} format.
//' @param year_list A list of vectors giving the main years used in each MPM.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates for each MPM.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates for each MPM.
//' @param dev_terms_list List of deviations for vital rate models in all MPMs.
//' @param density_vr_list List of \code{lefkoDensVR} objects holding density
//' relationships for all 14 vital rate models, for all MPMs.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param density_list A list of data frames of class \code{lefkoDens} for all
//' MPMs.
//' @param dens_index_list A list of data frames giving indices for density
//' dependent transitions for each MPM.
//' @param equivalence_list An optional numeric vector, list of numeric vectors,
//' data frame of class \code{adaptEq}, or list of data frames of class
//' \code{adaptEq}. If a numeric vector, then must have the same number of
//' elements as the number of MPMs, with each element giving the effect of an
//' individual of each MPM relative to a reference individual. If a list of
//' vectors, then the list should be composed of as many numeric vectors as
//' MPMs, with each vector giving the effect of each individual in each stage
//' relative to a reference individual. Data frames of class \code{adaptEq}, and
//' lists of such data frames, can be made with function
//' \code{\link{equiv_input}()}. Numeric entries used in these vectors can be
//' thought of as Lotka-Volterra interaction terms, such as are used in multiple
//' species competition models.
//' @param dev_terms_num_vec A vector giving the number of vital rate deviations
//' in each MPM. Used to create function-based matrices under override option.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms per MPM. Used to create function-based matrices under override option.
//' @param firstage_vec An integer vector containing the first age values for
//' all MPMs in order.
//' @param finalage_vec  An integer vector containing the final age values for
//' all MPMs in order.
//' @param stagecounts An integer vector containing the number of stages in each
//' MPM.
//' @param entry_time_vec An integer vector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param format_vec An integer vector giving the MPM format of each MPM.
//' @param inda_terms_num_vec An integer vector giving the number of numeric
//' terms given in individual covariate a.
//' @param indb_terms_num_vec An integer vector giving the number of numeric
//' terms given in individual covariate b.
//' @param indc_terms_num_vec An integer vector giving the number of numeric
//' terms given in individual covariate c.
//' @param inda_terms_cat_vec An integer vector giving the number of factor
//' terms given in individual covariate a.
//' @param indb_terms_cat_vec An integer vector giving the number of factor
//' terms given in individual covariate b.
//' @param indc_terms_cat_vec An integer vector giving the number of factor
//' terms given in individual covariate c.
//' @param dens_yn_vec A vector stating whether density dependence is used in
//' each MPM, given through \code{lefkoDens} objects.
//' @param dens_vr_yn_vec A vector stating whether density dependence is used in
//' each MPM, given through \code{lefkoDensVR} objects.
//' @param tweights_type_vec An integer vector giving the style of
//' \code{tweights} used in each MPM.
//' @param fecmod_vec A numeric vector giving the fecmod values for all MPMs.
//' @param sparse_vec A logical vector stating whether each MPM is in sparse
//' matrix format or not.
//' @param patch_vec A string vector giving the name of each patch used in each
//' MPM projection.
//' @param vrm_count An integer giving the number of vrm_inputs supplied.
//' @param nreps An integer giving the number of replicates to perform.
//' @param times An integer giving the amount of time steps to run the
//' projection for.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param stochastic A Boolean value indicating to perform a temporally
//' stochastic projection.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' 
//' @return The first three arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
void project3_fb_core (List& N_out, List& comm_out, List& extreme_mpm_out,
  List& errcheck_fb_out, const List start_list, const List vrm_list,
  const List tweights_list, const List stageframe_list, const List allmodels_all,
  const List allstages_all, const List supplement_list, const List year_list,
  const List ind_terms_num_list, const List ind_terms_cat_list,
  const List dev_terms_list, const List density_vr_list,
  const List sp_density_list, const List density_list,
  const List dens_index_list, const List equivalence_list,
  const IntegerVector dev_terms_num_vec, const IntegerVector sp_density_num_vec,
  const IntegerVector firstage_vec, const IntegerVector finalage_vec,
  const IntegerVector stagecounts, const IntegerVector entry_time_vec,
  const IntegerVector format_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_yn_vec,
  const IntegerVector dens_vr_yn_vec, const IntegerVector tweights_type_vec,
  const NumericVector fecmod_vec, const LogicalVector sparse_vec,
  const CharacterVector patch_vec, const int vrm_count, const int nreps,
  const int times, const int substoch, const double exp_tol,
  const double theta_tol, const bool integeronly, const bool stages_not_equal,
  const bool stochastic, const bool err_check, const bool err_check_extreme) {
  
  // start_list
  // density_vr_list
  // dens_vr_yn_vec

  //Rcout << "Entered project3_fb_core          ";
  
  int year_counter {0};
  IntegerVector inda_num_terms_counter (vrm_count);
  IntegerVector indb_num_terms_counter (vrm_count);
  IntegerVector indc_num_terms_counter (vrm_count);
  IntegerVector inda_cat_terms_counter (vrm_count);
  IntegerVector indb_cat_terms_counter (vrm_count);
  IntegerVector indc_cat_terms_counter (vrm_count);
  IntegerVector inda_num_terms_previous (vrm_count);
  IntegerVector indb_num_terms_previous (vrm_count);
  IntegerVector indc_num_terms_previous (vrm_count);
  IntegerVector inda_cat_terms_previous (vrm_count);
  IntegerVector indb_cat_terms_previous (vrm_count);
  IntegerVector indc_cat_terms_previous (vrm_count);
  IntegerVector dev_num_counter (vrm_count);
  IntegerVector sp_density_counter (vrm_count);
  
  List N_out_pre (nreps);
  List comm_out_pre (vrm_count);
  List stochastic_year_vector (vrm_count);
  
  // Stochastic yearnumber vector generation
  if (stochastic) {
    for (int i = 0; i < vrm_count; i++) {
      CharacterVector current_mainyears = as<CharacterVector>(year_list(i));
      unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
      IntegerVector current_mainyears_int = seq(0, (no_mainyears - 1));
      
      List stochastic_year_vector_reps (nreps);
     
      for (int j = 0; j < nreps; j++) {
        IntegerVector years_topull;
        
        if (tweights_type_vec(i) == 0) {
          NumericVector twinput (no_mainyears, (1.0 / static_cast<double>(no_mainyears)));
          
          IntegerVector years_topull_pre = Rcpp::RcppArmadillo::sample(current_mainyears_int,
            (times - entry_time_vec(i)), true, twinput);
          
          IntegerVector years_topull_almost (times);
          for (int copy_elem = 0; copy_elem < (times - entry_time_vec(i)); copy_elem++) {
            years_topull_almost(copy_elem + entry_time_vec(i)) = years_topull_pre(copy_elem);
          }
          years_topull = years_topull_almost;
        } else if (tweights_type_vec(i) == 1) {
          NumericVector twinput = as<NumericVector>(tweights_list(i));
          NumericVector twinput_st = twinput / sum(twinput);
          
          IntegerVector years_topull_pre = Rcpp::RcppArmadillo::sample(current_mainyears_int,
            (times - entry_time_vec(i)), true, twinput_st);
          
          IntegerVector years_topull_almost (times);
          for (int copy_elem = 0; copy_elem < (times - entry_time_vec(i)); copy_elem++) {
            years_topull_almost(copy_elem + entry_time_vec(i)) = years_topull_pre(copy_elem);
          }
          years_topull = years_topull_almost;
        } else if (tweights_type_vec(i) == 2) {
          arma::ivec chosen_mats_arma = as<arma::ivec>(current_mainyears_int);
          arma::mat twinput_mat = as<arma::mat>(tweights_list(i));
          arma::vec twinput = twinput_mat.col(0);
          twinput = twinput / sum(twinput);
          
          IntegerVector years_topull_pre (times);
          NumericVector twinput_setup (no_mainyears, (1.0 / static_cast<double>(no_mainyears)));
          arma::ivec first_choice = Rcpp::RcppArmadillo::sample(chosen_mats_arma,
            times, true, twinput_setup);
          years_topull_pre(entry_time_vec(i)) = current_mainyears_int(first_choice(0));
          
          for (int k = 1; k < (times - entry_time_vec(i)); k++) {
            arma::ivec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(chosen_mats_arma,
              1, true, twinput);
            years_topull_pre(k + entry_time_vec(i)) = theprophecy_piecemeal(0);
            
            arma::uvec tnotb_preassigned = find(chosen_mats_arma == theprophecy_piecemeal(0));
            twinput = twinput_mat.col(static_cast<int>(tnotb_preassigned(0)));
            twinput = twinput / sum(twinput);
          }
          years_topull = years_topull_pre;
          
        } else {
          throw Rcpp::exception("tweights_type_vec error.", false);
        }
        stochastic_year_vector_reps(j) = years_topull;
      }
      stochastic_year_vector(i) = stochastic_year_vector_reps;
    }
  }
  
  for (int i = 0; i < vrm_count; i++) {
    
    //Rcout << "project3_fb_core a i: " << i << "          ";
    
    List pop_reps (nreps);
    
    arma::mat pops_out_pre (stagecounts(i), (times + 1), fill::zeros);
    for (int j = 0; j < nreps; j++) {
      pop_reps(j) = pops_out_pre;
    }
    
    comm_out_pre(i) = pop_reps;
  }
  
  List extreme_mpm_reps (nreps);
  List errcheck_mpmout_rep (nreps);
  
  // Main projection
  for (int current_rep = 0; current_rep < nreps; current_rep++) {
    //Rcout << "\nMain projection start" << endl;
    //Rcout << "project3_fb_core b current_rep: " << current_rep << "          ";
    
    List running_popvecs = clone(start_list);
    NumericMatrix N_vrm (vrm_count, (times + 1));
    List extreme_mpm_reps_times (times);
    List errcheck_mpmout_rep_time (times); 
    
    IntegerVector year_counter (vrm_count);
    
    for (int current_time = 0; current_time < times; current_time++) {
      //Rcout << "project3_fb_core c current_time: " << current_time << "          ";
      
      Rcpp::checkUserInterrupt();
      
      List extreme_mpm_reps_times_vrms (vrm_count);
      List errcheck_mpmout_rep_time_vrm (vrm_count); 
      
      for (int i = 0; i < vrm_count; i++) {
        //Rcout << "\nproject3_fb_core d i (current vrm): " << i << "          ";
        
        List reps_out = comm_out_pre(i);
        arma::mat pops_out = as<arma::mat>(reps_out(current_rep));
        
        if (current_time > (entry_time_vec(i) - 1)) {
          arma::vec running_popvec_vrm = as<arma::vec>(running_popvecs(i));
          
          if (current_time == entry_time_vec(i)) {
            pops_out.col(current_time) = running_popvec_vrm;
            
            double N_current = accu(running_popvec_vrm);
            N_vrm(i, current_time) = N_current;
          }
          
          List current_vrm_extract = as<List>(allmodels_all(i));
          List current_vrm_unextract = as<List>(vrm_list(i));
          DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
          int ehrlen_format {1}; // This will need to be dealt with differently later
          
          int mpm_style {1};
          if (format_vec(i) < 3) {
            mpm_style = 0;
          } else if (format_vec(i) == 4) {
            mpm_style = 2;
          }
          
          DataFrame current_mpm_allstages = as<DataFrame>(allstages_all(i));
          
          List surv_proxy = as<List>(current_vrm_extract(0));
          List obs_proxy = as<List>(current_vrm_extract(1));
          List size_proxy = as<List>(current_vrm_extract(2));
          List sizeb_proxy = as<List>(current_vrm_extract(3));
          List sizec_proxy = as<List>(current_vrm_extract(4));
          List repst_proxy = as<List>(current_vrm_extract(5));
          List fec_proxy = as<List>(current_vrm_extract(6));
          List jsurv_proxy = as<List>(current_vrm_extract(7));
          List jobs_proxy = as<List>(current_vrm_extract(8));
          List jsize_proxy = as<List>(current_vrm_extract(9));
          List jsizeb_proxy = as<List>(current_vrm_extract(10));
          List jsizec_proxy = as<List>(current_vrm_extract(11));
          List jrepst_proxy = as<List>(current_vrm_extract(12));
          List jmatst_proxy = as<List>(current_vrm_extract(13));
          DataFrame current_paramnames = as<DataFrame>(current_vrm_extract(14));
          
          CharacterVector current_mainyears = as<CharacterVector>(year_list(i));
          unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
          
          DataFrame group2_frame = as<DataFrame>(current_vrm_unextract["group2_frame"]);
          CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
          
          DataFrame patch_frame = as<DataFrame>(current_vrm_unextract["patch_frame"]);
          CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
          
          //Rcout << "project3_fb_core e current_mainyears: " << current_mainyears << "          ";
          //Rcout << "project3_fb_core f current_mainpatches: " << current_mainpatches << "          ";
          
          // Not sure if we need the next bit
          DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
          DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
          DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
          CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
          CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
          CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
          
          // Counter resets
          if (year_counter(i) >= no_mainyears || current_time == 0) year_counter(i) = 0;
          int yearnumber {0};
          
          if (!stochastic) {
            yearnumber = year_counter(i);
          } else {
            List stochastic_times_list = as<List>(stochastic_year_vector(i));
            IntegerVector stochastic_times = as<IntegerVector>(stochastic_times_list(current_rep));
            yearnumber = stochastic_times(current_time);
          }
          
          CharacterVector current_year = as<CharacterVector>(current_mainyears(yearnumber));
          //Rcout << "project3_fb_core g current_year: " << current_year << "          ";
          
          IntegerVector pvis = index_l3(current_mainpatches, patch_vec(i));
          if (static_cast<int>(pvis.length()) == 0) {
            throw Rcpp::exception("Value for argument patch not found.", false);
          }
          int patchnumber = pvis(0);
          //Rcout << "project3_fb_core h patch_vec: " << patch_vec << "          ";
          //Rcout << "project3_fb_core i i (the ith element in patch_vec, corresponding to vrm, is chosen): " << i << "          ";
          //Rcout << "project3_fb_core j patchnumber: " << patchnumber << "          ";
          
          if (inda_num_terms_counter(i) >= inda_terms_num_vec(i)) inda_num_terms_counter(i) = 0;
          if (indb_num_terms_counter(i) >= indb_terms_num_vec(i)) indb_num_terms_counter(i) = 0;
          if (indc_num_terms_counter(i) >= indc_terms_num_vec(i)) indc_num_terms_counter(i) = 0;
          if (inda_cat_terms_counter(i) >= inda_terms_cat_vec(i)) inda_cat_terms_counter(i) = 0;
          if (indb_cat_terms_counter(i) >= indb_terms_cat_vec(i)) indb_cat_terms_counter(i) = 0;
          if (indc_cat_terms_counter(i) >= indc_terms_cat_vec(i)) indc_cat_terms_counter(i) = 0;
          
          List current_ind_terms_num = ind_terms_num_list(i);
          List current_ind_terms_cat = ind_terms_cat_list(i);
          
          NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
          NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
          NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
          CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
          CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
          CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
          
          NumericVector f2_inda = {f_inda_full(inda_num_terms_counter(i))};
          NumericVector f1_inda = {f_inda_full(inda_num_terms_previous(i))};
          NumericVector f2_indb = {f_indb_full(indb_num_terms_counter(i))};
          NumericVector f1_indb = {f_indb_full(indb_num_terms_previous(i))};
          NumericVector f2_indc = {f_indc_full(indc_num_terms_counter(i))};
          NumericVector f1_indc = {f_indc_full(indc_num_terms_previous(i))};
          CharacterVector r2_inda = as<CharacterVector>(r_inda_full(inda_cat_terms_counter(i)));
          CharacterVector r1_inda = 
            as<CharacterVector>(r_inda_full(inda_cat_terms_previous(i)));
          CharacterVector r2_indb = as<CharacterVector>
            (r_indb_full(indb_cat_terms_counter(i)));
          CharacterVector r1_indb = as<CharacterVector>
            (r_indb_full(indb_cat_terms_previous(i)));
          CharacterVector r2_indc = as<CharacterVector>
            (r_indc_full(indc_cat_terms_counter(i)));
          CharacterVector r1_indc = 
            as<CharacterVector>(r_indc_full(indc_cat_terms_previous(i)));
          
          NumericVector dv_terms (14);
          if (dev_terms_num_vec(i) > 0) {
            DataFrame used_dv_df = as<DataFrame>(dev_terms_list(i));
          
            if (dev_num_counter(i) >= dev_terms_num_vec(i)) dev_num_counter(i) = 0;
            
            for (int j = 0; j < 14; j++) {
              dv_terms(j) = used_dv_df(dev_num_counter(i), j);
            }
            dev_num_counter(i) = dev_num_counter(i) + 1;
          }
          
          bool dvr_bool {false};
          
          LogicalVector dvr_yn = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          IntegerVector dvr_style (14);
          IntegerVector dvr_time_delay (14);
          NumericVector dvr_alpha (14);
          NumericVector dvr_beta (14);
          
          if (dens_vr_yn_vec(i) > 0) {
            dvr_bool = true;
            
            DataFrame current_dvr = as<DataFrame>(density_vr_list(i));
            LogicalVector true_dvr_yn = as<LogicalVector>(current_dvr(1));
            IntegerVector true_dvr_style = as<IntegerVector>(current_dvr(2));
            IntegerVector true_dvr_time_delay = as<IntegerVector>(current_dvr(3));
            NumericVector true_dvr_alpha = as<NumericVector>(current_dvr(4));
            NumericVector true_dvr_beta = as<NumericVector>(current_dvr(5));
            
            dvr_yn = true_dvr_yn;
            dvr_style = true_dvr_style;
            dvr_time_delay = true_dvr_time_delay;
            dvr_alpha = true_dvr_alpha;
            dvr_beta = true_dvr_beta;
          }
          
          double maxsize {0.0};
          double maxsizeb {0.0};
          double maxsizec {0.0};
          
          if (format_vec(i) < 5) {
            DataFrame current_allstages = as<DataFrame>(allstages_all(i));
            
            NumericVector size3 = as<NumericVector>(current_allstages["size3"]);
            NumericVector size2n = as<NumericVector>(current_allstages["size2n"]);
            NumericVector size2o = as<NumericVector>(current_allstages["size2o"]);
            NumericVector sizeb3 = as<NumericVector>(current_allstages["sizeb3"]);
            NumericVector sizeb2n = as<NumericVector>(current_allstages["sizeb2n"]);
            NumericVector sizeb2o = as<NumericVector>(current_allstages["sizeb2o"]);
            NumericVector sizec3 = as<NumericVector>(current_allstages["sizec3"]);
            NumericVector sizec2n = as<NumericVector>(current_allstages["sizec2n"]);
            NumericVector sizec2o = as<NumericVector>(current_allstages["sizec2o"]);
            
            NumericVector maxveca = {max(size3), max(size2n), max(size2o)};
            NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
            NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
            
            maxsize = max(maxveca);
            maxsizeb = max(maxvecb);
            maxsizec = max(maxvecc);
          }
          
          double dens_sp {1.0};
          
          if (sp_density_num_vec(i) > 0) {
            if (sp_density_counter(i) >= sp_density_num_vec(i)) sp_density_counter(i) = 0;
            
            NumericVector current_sp_density = as<NumericVector>(sp_density_list(i));
            dens_sp = current_sp_density(sp_density_counter(i));
            
            sp_density_counter(i) = sp_density_counter(i) + 1;
          }
          
          NumericVector dens_n (14, 1.0); // This needs to be updated with the actual pop size
          
          //Rcout << "project3_fb_core k          ";
          
          List current_mpm;
          if (format_vec(i) < 5) {
            current_mpm = AdaptMats::mazurekd(current_mpm_allstages, current_stageframe,
              format_vec(i), surv_proxy, obs_proxy, size_proxy, sizeb_proxy,
              sizec_proxy, repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy,
              jsize_proxy, jsizeb_proxy, jsizec_proxy, jrepst_proxy, jmatst_proxy,
              f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc, r2_inda,
              r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, dv_terms, dvr_bool,
              dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n, dens_sp,
              fecmod_vec(i), maxsize, maxsizeb, maxsizec,
              static_cast<unsigned int>(firstage_vec(i)),
              static_cast<unsigned int>(finalage_vec(i)), true, yearnumber,
              patchnumber, exp_tol, theta_tol, true, err_check, sparse_vec(i),
              true);
            
            if (err_check) {
              NumericMatrix mpm_out = as<NumericMatrix>(current_mpm["out"]);
              errcheck_mpmout_rep_time_vrm(i) = mpm_out; 
            }
              
          } else {
            IntegerVector all_ages = seq(firstage_vec(i), finalage_vec(i));
            DataFrame current_supplement;
            if (supplement_list(i) == R_NilValue) {
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe, surv_proxy,
                fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc,
                r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, dv_terms(0),
                dv_terms(6), dens_sp, fecmod_vec(i), finalage_vec(i), true,
                yearnumber, patchnumber, dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
                dens_n, exp_tol, theta_tol, sparse_vec(i));
              
            } else {
              current_supplement = as<DataFrame>(supplement_list(i));
              
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe, surv_proxy,
                fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc,
                r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, dv_terms(0),
                dv_terms(6), dens_sp, fecmod_vec(i), finalage_vec(i), true,
                yearnumber, patchnumber, dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
                dens_n, exp_tol, theta_tol, sparse_vec(i), current_supplement);
            }
          }
          if (dens_yn_vec(i) == 0) {
            if (!sparse_vec(i)) {
              arma::mat current_A = as<arma::mat>(current_mpm(0));
              
              running_popvec_vrm = current_A * running_popvec_vrm;
              if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = current_A;
            } else {
              arma::sp_mat current_A = as<arma::sp_mat>(current_mpm(0));
              
              running_popvec_vrm = current_A * running_popvec_vrm;
              if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = current_A;
            }
          } else {
            //Rcout << "project3_fb_core l          ";
            
            DataFrame used_density_input = as<DataFrame>(density_list(i));
            DataFrame used_density_index_input = as<DataFrame>(dens_index_list(i));
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            if (current_time > used_delay) { // Changed to allow different delay Ns
              if (!stages_not_equal) {
                NumericVector delay_issue = N_vrm(_, (current_time + 1 - used_delay));
                
                double delay_N_sum = sum(delay_issue);
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_vrm, current_mpm, delay_N_sum, 0, integeronly,
                  substoch, used_density_input, used_density_index_input, false,
                  static_cast<bool>(sparse_vec(i)),
                  static_cast<bool>(sparse_vec(i)), false, err_check_extreme);
                
                running_popvec_vrm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = new_projmat;
              } else {
                double delay_N_sum {0.0};
                
                if (current_time > 0) {
                  for (int l = 0; l < vrm_count; l++) {
                    List current_pop_list = as<List>(comm_out_pre(l));
                    arma::mat delay_pop = as<arma::mat>(current_pop_list(current_rep));
                    arma::vec delay_pop_vec = delay_pop.col(current_time + 1 - used_delay);
                    arma::vec current_equiv_vec = as<arma::vec>(equivalence_list(l));
                    arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                    double delay_pop_N = accu(adjusted_delay_pop_vec);
                    
                    delay_N_sum += delay_pop_N;
                  }
                }
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_vrm, current_mpm, delay_N_sum, 0, integeronly,
                  substoch, used_density_input, used_density_index_input, false,
                  static_cast<bool>(sparse_vec(i)),
                  static_cast<bool>(sparse_vec(i)), false, err_check_extreme);
                
                running_popvec_vrm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = new_projmat;
              }
            } else {
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              proj3dens_ad(new_popvec, new_projmat, new_projsp,
                running_popvec_vrm, current_mpm, 0.0, 0, integeronly, substoch,
                used_density_input, used_density_index_input, false,
                static_cast<bool>(sparse_vec(i)), static_cast<bool>(sparse_vec(i)),
                false, err_check_extreme);
              
              running_popvec_vrm = new_popvec;
              if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = new_projmat;
            }
          }
          
          //Rcout << "project3_fb_core m          ";
          
          if (integeronly) running_popvec_vrm = floor(running_popvec_vrm);
          double N_current = accu(running_popvec_vrm);
          N_vrm(i, (current_time + 1)) = N_current;
          
          inda_num_terms_previous(i) = static_cast<int>(inda_num_terms_counter(i));
          indb_num_terms_previous(i) = static_cast<int>(indb_num_terms_counter(i));
          indc_num_terms_previous(i) = static_cast<int>(indc_num_terms_counter(i));
          inda_cat_terms_previous(i) = static_cast<int>(inda_cat_terms_counter(i));
          indb_cat_terms_previous(i) = static_cast<int>(indb_cat_terms_counter(i));
          indc_cat_terms_previous(i) = static_cast<int>(indc_cat_terms_counter(i));
          
          inda_num_terms_counter(i) = inda_num_terms_counter(i) + 1;
          indb_num_terms_counter(i) = indb_num_terms_counter(i) + 1;
          indc_num_terms_counter(i) = indc_num_terms_counter(i) + 1;
          inda_cat_terms_counter(i) = inda_cat_terms_counter(i) + 1;
          indb_cat_terms_counter(i) = indb_cat_terms_counter(i) + 1;
          indc_cat_terms_counter(i) = indc_cat_terms_counter(i) + 1;
          
          running_popvecs(i) = running_popvec_vrm;
          pops_out.col(current_time + 1) = running_popvec_vrm;
          
        }
        year_counter(i) = year_counter(i) + 1;
        
        reps_out(current_rep) = pops_out;
        comm_out_pre(i) = reps_out;
      } // vrm loop
      
      if (err_check_extreme) extreme_mpm_reps_times(current_time) = extreme_mpm_reps_times_vrms;
      if (err_check) errcheck_mpmout_rep_time(current_time) = errcheck_mpmout_rep_time_vrm;
    } // current_time loop
    comm_out = comm_out_pre;
    N_out_pre(current_rep) = N_vrm;
    if (err_check_extreme) extreme_mpm_reps(current_rep) = extreme_mpm_reps_times;
    if (err_check) errcheck_mpmout_rep(current_rep) = errcheck_mpmout_rep_time;
    
    //Rcout << "project3_fb_core n          ";
    
  } // current_rep loop
  N_out = N_out_pre;
  if (err_check_extreme) extreme_mpm_out = extreme_mpm_reps;
  if (err_check) errcheck_fb_out = errcheck_mpmout_rep;
}

//' Project Multiple MPMs With or Without Density Dependence
//' 
//' Function \code{project3} uses pre-existing or function-based MPMs to run
//' community projection simulations, in which different populations are run as
//' separate MPMs. Density dependence can be used with individual equivalence
//' vectors specifying Lotka-Volterra coefficients to adjust overall population
//' sizes to make them comparable.
//' 
//' @name project3
//' 
//' @param mpms An optional list of MPMs. Each MPM must be of class
//' \code{lefkoMat}.
//' @param vrms An optional list of \code{vrm_input} objects, each corresponding
//' to a distinct MPM that will be created during projection. Each
//' \code{vrm_input} object requires its own stageframe, entered in the same
//' order via argument \code{stageframes}.
//' @param stageframes An optional list of stageframes, corresponding in number
//' and order to the MPMs in argument \code{vrms}. Each stageframe must be of
//' class \code{stageframe}.
//' @param supplements An optional list of data frames of class \code{lefkoSD}
//' that provide supplemental data that should be incorporated into
//' function-based MPMs. If used, then should be the same number of data frames
//' as the number of MPMs provided in the list for argument \code{vrms}. MPMs
//' that do not need supplemental data should be entered as \code{NULL} in this
//' list. See \code{\link[lefko3]{supplemental}()} for details.
//' @param equivalence An optional numeric vector, list of numeric vectors,
//' data frame of class \code{adaptEq}, or list of data frames of class
//' \code{adaptEq}. If a numeric vector, then must have the same number of
//' elements as the number of MPMs, with each element giving the effect of an
//' individual of each MPM relative to a reference individual. If a list of
//' vectors, then the list should be composed of as many numeric vectors as
//' MPMs, with each vector giving the effect of each individual in each stage
//' relative to a reference individual. Data frames of class \code{adaptEq}, and
//' lists of such data frames, can be made with function
//' \code{\link{equiv_input}()}. Numeric entries used in these vectors can be
//' thought of as Lotka-Volterra interaction terms, such as are used in multiple
//' species competition models.
//' @param starts An optional list of \code{lefkoSV} objects, which are data
//' frames providing the starting numbers of individuals of each stage. If
//' provided, then one is needed per MPM. If not provided, then all projections
//' start with a single individual of each stage per MPM.
//' @param years An optional term corresponding either to a single integer vector
//' of time \code{t} values, if all MPMs will use the same time \code{t} or set
//' of time \code{t}'s, or a list of such vectors with each vector corresponding
//' to each MPM in order. In the latter case, a vector composed of a single
//' \code{NA} value is interpreted to mean that all time \code{t} values in the
//' MPM should be utilized. If a vector shorter than \code{times} is supplied,
//' then this vector will be cycled.
//' @param patches An optional string vector with length equal to the number of
//' MPMs, detailing the name of each patch to project for each MPM, in order.
//' Only a single pop-patch may be projected for each MPM given. A value of
//' \code{NA} can be supplied to indicate that the population-level matrices
//' should be projected (if argument \code{mpms} is used and a population-level
//' set of matrices exist), or that the first patch noted should be used.
//' Defaults to the population-level set or the first patch, depending on
//' whether the former exists.
//' @param tweights An optional list composed of numeric vectors or matrices
//' denoting the probabilities of choosing each matrix in each MPM in a
//' stochastic projection. If an element of the list is a matrix, then a
//' first-order Markovian environment is assumed, in which the probability of
//' choosing a specific annual matrix depends on which annual matrix is
//' currently chosen. If an element of the list is a vector, then the choice of
//' annual matrix is assumed to be independent of the current matrix. Defaults
//' to equal weighting among matrices. If used, then one element per MPM is
//' required, with equal weighting assumed for any element set to \code{NULL}.
//' @param format An optional integer vector indicating the kind of
//' function-based MPM to create for each \code{vrm_input} object entered in
//' argument \code{vrms}. Possible choices include: \code{1}, Ehrlen-format
//' historical MPM; \code{2}, deVries-format historical MPM; \code{3},
//' ahistorical MPM (default); \code{4}, age-by-stage MPM; and \code{5}, Leslie
//' (age-based) MPM.
//' @param entry_time An optional integer vector giving the entry time for each
//' MPM into the projection. Defaults to a zero vector with the length of the
//' number of MPMs, as given either by argument \code{mpms} or \code{vrms}.
//' @param sp_density An optional argument for use with \code{vrm_input} objects
//' that specifies the spatial density to be used in each time step. If used,
//' may either be a numeric vector giving a single spatial density for each
//' \code{vrm_input} object entered in argument \code{vrms} (in this case, the
//' value of spatial density given for each \code{vrm_input} object will be held
//' constant through the projection), or a list of as many numeric vectors as
//' \code{vrm_input} objects, with the length of each vector giving the spatial
//' density at each time step. If vectors are shorter than specified in 
//' \code{times}, then these values will be cycled.
//' @param ind_terms An optional argument providing values of individual or
//' environmental covariate values for \code{vrm_input} objects used in
//' function-based projection. Can be set either to a single data frame with 3
//' columns giving values for up to 3 covariates across time (rows give the time
//' order of these values), or a list of as many such data frames as
//' \code{vrm_input} objects. In the latter case, \code{vrm_input} objects that
//' do not use such covariates should have the associated element set to
//' \code{NULL}. Unused terms within each data frame must be set to \code{0}
//' (use of \code{NA} will produce errors.) If the number of rows is less than
//' \code{times}, then these values will be cycled.
//' @param dev_terms An optional list of data frames, one for each
//' \code{vrm_input} object. Each should include 14 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival,
//' 9: juvenile observation, 10: juvenile primary size, 11: juvenile secondary
//' size, 12: juvenile tertiary size, 13: juvenile reproduction, and
//' 14: juvenile maturity transition. Unused terms must be set to \code{0} (use
//' of \code{NA} will produce errors). Single or small numbers of values per
//' vital rate model are also allowed, and if the number of rows is less than
//' \code{times}, then the terms will be cycled.
//' @param fb_sparse A logical vector indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE} for
//' each MPM.
//' @param firstage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the starting ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Typically,
//' the starting age in such MPMs should be set to \code{0} if post-breeding and
//' \code{1} if pre-breeding. All other MPMs should be set to \code{0}. Do not
//' use if no MPM has age structure. 
//' @param finalage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the final ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Do not use if no
//' MPM has age structure.
//' @param fecage_min An optional integer vector used for function-based Leslie
//' MPMs giving the first age at which organisms can be reproductive in such
//' MPMs. Use only if at least one MPM is a function-based Leslie MPM. Defaults
//' to the values given in \code{firstage}.
//' @param fecage_max An optional integer vector used for function-based Leslie
//' MPMs giving the final age at which organisms can be reproductive in such
//' MPMs. Use only if at least one MPM is a function-based Leslie MPM. Defaults
//' to the values given in \code{finalage}.
//' @param cont An optional vector used for function-based Leslie and
//' age-by-stage MPMs stating whether the MPM should should include a stasis
//' transition within the final age. This should be used only when an organism
//' can maintain the demographic characteristics of the final described age
//' after reaching that age. Can be entered as a logical vector or an integer
//' vector. MPMs without age structure should be entered as \code{0} or
//' \code{FALSE}. Do not use if no MPM has age structure.
//' @param fecmod An optional vector used for function-based MPMs giving scalar
//' multipliers for fecundity terms, when two fecundity variables are used for a
//' collective fecundity per individual. Each entry refers to each 
//' \code{vrm_input} object in argument \code{vrms}, in the same order.
//' @param density An optional list of data frames of class \code{lefkoDens},
//' which provide details for density dependence in MPM elements and have been
//' created with function \code{\link[lefko3]{density_input}()}. If used, then
//' one such data frame per MPM is required. MPMs to be run without density
//' dependence should be set to \code{NULL}.
//' @param density_vr An optional list of data frames of class
//' \code{lefkoDensVR}, which provide details for density dependence in vital
//' rate models and have been created with function
//' \code{link[lefko3]{density_vr}()}. If used, then one such data frame per MPM
//' is required. MPMs to be run without vital describing density dependence
//' relationships in vital rates should be set to \code{NULL}. Can only be used
//' with function-based projections.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking. Can also be set to the text value
//' \code{"extreme"}, in which case all \code{err_check} output plus a multiple
//' level list with each MPM used in each time step will be output.
//' @param stochastic A logical value indicating whether the projection will be
//' run as a temporally stochastic projection. Defaults to \code{FALSE}.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion in each MPM to the
//' nearest integer. Defaults to \code{FALSE}.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent and density
//' independent simulations. Defaults to \code{0}, which does not enforce
//' substochasticity. Alternatively, \code{1} forces all survival-transition
//' elements to range from 0.0 to 1.0, and forces fecundity to be non-negative;
//' and \code{2} forces all column rows in the survival-transition matrices to
//' total no more than 1.0, in addition to the actions outlined for option
//' \code{1}. Both settings \code{1} and \code{2} change negative fecundity
//' elements to \code{0.0}.
//' @param nreps The number of replicate projections. Defaults to \code{1}.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' \code{10000}.
//' @param prep_mats An integer value for use when creating function-based MPM
//' projections. If using \code{vrms} input instead of \code{mpms} input, then
//' this argument determines how many matrices should be used as a limit to
//' develop matrices prior to running the projection. See \code{Notes} for
//' further details.
//' @param force_fb A logical value indicating whether to force function-based
//' MPMs to be developed at each time step even if fewer than \code{prep_mats}.
//' Defaults to \code{FALSE}.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param theta_tol A numeric value used to indicate a maximum value to theta as
//' used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' 
//' @return A list of class \code{adaptProj}, with the following elements:
//' \item{comm_out}{A two-level list with the top level list having number of
//' elements equal to the number of MPMs used as input, and the lower level
//' corresponding to the number of replicates. Each element of the lower level
//' list is a data frame showing the number of individuals in each stage at each
//' time. Rows and columns in the data frames correspond to stages and time
//' steps, respectively.}
//' \item{N_out}{A list with the number of elements equal to the number of
//' replicates. Each element within this list is data frame showing the number
//' of individuals of each species or genotype alive at each time. The number of
//' rows are equal to the number of MPMs used, and the columns correspond to the
//' time steps.}
//' \item{stageframe_list}{A list in which each element is the stageframe for
//' each MPM used.}
//' \item{hstages_list}{A list giving the used \code{hstages} data frames, which
//' identify the correct stage pairing for each row / column in each
//' historical MPM utilized.}
//' \item{agestages_list}{A list giving the used \code{agestages} data frames,
//' which identify the correct age-stage pairing for each row / column in each
//' age-by-stage MPM utilized.}
//' \item{labels}{A small data frame giving the the population and patch
//' identities for each MPM entered.}
//' \item{err_check}{An optional list composed of an additional six lists, each
//' of which has the number of elements equal to the number of MPMs utilized.
//' List output include \code{allstages_all}, which gives the indices of
//' estimatedtransitions in MPMs constructed by function \code{project3()} from
//' input vital rate models; \code{allmodels_all}, which provides all vital rate
//' models as decomposed and interpreted by function \code{project3()};
//' \code{equivalence_list}, which gives the stage equivalence for density
//' calculations across MPMs; \code{density_list}, which gives the
//' \code{density} inputs utilized; \code{dens_index_list}, which provides
//' indices used to identify matrix elements for density dependence; and
//' \code{density_vr_list}, which gives the \code{density_vr} inputs utilized.}
//' 
//' @section Notes:
//' 
//' This function has been optimized in the function-based approach such that
//' if there are relatively few matrices required per MPM to run the projection
//' forward, then these matrices will be made prior to running the projection.
//' This approach saves time, but only if there are relatively few unique
//' matrices required for each MPM. If many or only unique MPMs are required at
//' each time step, then the matrices will be made on the fly during the
//' projection itself. Such a situation will most likely occur if each time
//' step requires a new matrix resulting from a unique individual covariate
//' value, or if the \code{density_vr} argument is used. The key argument
//' determining this behavior is \code{prep_mats}, which provides the maximum
//' limit for the number of matrices required per MPM in order to create
//' matrices prior to projection.
//' 
//' @examples
//' 
//' library(lefko3)
//' data(cypdata)
//' 
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
//' cycaraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//'   
//' cyparaw_v1 <- verticalize3(data = cypa_data, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' cyp_supp_list1 <- list(cypsupp2r, cypsupp2r)
//' 
//' cycamatrix2r <- rlefko2(data = cycaraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cypamatrix2r <- rlefko2(data = cyparaw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r_alt,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cyp_mpm_list <- list(cycamatrix2r, cypamatrix2r)
//' 
//' cyca2_start <- start_input(cycamatrix2r, stage2 = c("SD", "P1", "P2"),
//'   value = c(500, 100, 200))
//' cypa2_start <- start_input(cypamatrix2r, stage2 = c("SD", "P1", "P2"),
//'   value = c(5000, 1000, 2000))
//' cyp_start_list <- list(cyca2_start, cypa2_start)
//' 
//' cyp2_dv <- density_input(cypamatrix2r, stage3 = c("SD", "P1"),
//'   stage2 = c("rep", "rep"), style = c(1, 1), alpha = c(0.5, 1.2),
//'   beta = c(1.0, 2.0), type = c(2, 1))
//' cyp_dv_list <- list(cyp2_dv, cyp2_dv)
//' 
//' cyp_comm_proj <- project3(mpms = cyp_mpm_list, starts = cyp_start_list,
//'   density = cyp_dv_list, times = 10)
//' 
//' @export project3
// [[Rcpp::export(project3)]]
List project3 (Nullable<RObject> mpms  = R_NilValue,
  Nullable<RObject> vrms = R_NilValue, Nullable<RObject> stageframes  = R_NilValue,
  Nullable<RObject> supplements = R_NilValue, Nullable<RObject> equivalence = R_NilValue,
  Nullable<RObject> starts = R_NilValue, Nullable<RObject> years = R_NilValue,
  Nullable<RObject> patches = R_NilValue, Nullable<RObject> tweights = R_NilValue,
  
  Nullable<RObject> format = R_NilValue, Nullable<RObject> entry_time = R_NilValue,
  Nullable<RObject> sp_density = R_NilValue, Nullable<RObject> ind_terms = R_NilValue,
  Nullable<RObject> dev_terms = R_NilValue, Nullable<RObject> fb_sparse = R_NilValue,
  
  Nullable<RObject> firstage = R_NilValue, Nullable<RObject> finalage = R_NilValue,
  Nullable<RObject> fecage_min = R_NilValue, Nullable<RObject> fecage_max = R_NilValue,
  Nullable<RObject> cont = R_NilValue, Nullable<RObject> fecmod = R_NilValue,
  
  Nullable<RObject> density = R_NilValue, Nullable<RObject> density_vr = R_NilValue,
  Nullable<RObject> err_check = R_NilValue,
  
  bool stochastic = false, bool integeronly = false, int substoch = 0,
  int nreps = 1, int times = 10000, int prep_mats = 20, bool force_fb = false,
  double exp_tol = 700.0, double theta_tol = 100000000.0) {
  
  bool preexisting {false}; // Are preexisting MPMs being used?
  bool funcbased {false}; // Will function-based MPMs be created?
  bool fb_override {false}; // Will fbMPMs be created but run via code for preexisting MPMs?
  bool entry_time_vec_use {false}; // Are any elements in entry_time greater than 0?
  bool stages_not_equal {false}; // Are equivalence vectors supplied separating even stages?
  bool pure_fleslie {false}; // Are all function-based MPMs Leslie MPMs?
  
  int mpm_count {0};
  int vrm_count {0};
  int total_mpms {0};  // This includes all MPMs and VRMs
  int stageframe_count {0};
  int stageframe_notNull_count {0};
  int supplement_count {0};
  int equivalence_count {0};
  int start_count {0};
  int tweights_count {0};
  int density_count {0};
  int entry_time_count {0};
  int density_vr_count {0};
  int sparse_vec_count {0};
  int found_fleslie {0};
  
  bool inda_char {false};
  bool indb_char {false};
  bool indc_char {false};
  bool err_check_bool {false};
  bool err_check_extreme {false};
  
  // err_check processing
  LefkoInputs::RObj_TF_input_check("err_check", "extreme", err_check_bool,
    err_check_extreme, true, true, err_check);
  
  // Age-by-stage and Leslie MPM settings
  IntegerVector firstage_vec;
  IntegerVector finalage_vec;
  IntegerVector cont_vec;
  
  // Leslie MPM only
  IntegerVector fecage_min_vec;
  IntegerVector fecage_max_vec;
  
  // Main lists
  List mpm_list;
  List A_list;
  List vrm_list;
  List stageframe_list;
  List stageframe_list_fb; // List ending in _fb are only used in function-based cases
  List supplement_list;
  List supplement_list_fb;
  List repmatrix_list;
  List equivalence_list;
  List hstages_list;
  List agestages_list;
  List start_list;
  List year_list;
  List tweights_list;
  List density_list;
  List dens_index_list;  // Holds element index vectors for density_frames in density_list
  List density_vr_list;
  List sp_density_list;
  List ind_terms_num_list;
  List ind_terms_cat_list;
  List dev_terms_list;
  List allstages_all; // Used in fbMPM processing
  List allmodels_all; // Used in fbMPM processirg
  List extreme_mpm_out; // err_check processing only
  List err_check_fb_out; // err_check processing only
  
  LogicalVector sparse_vec; // Vector accounting for whether MPMs are / should be in sparse format
  IntegerVector stagecounts; // # stages in each MPM
  IntegerVector matrowcounts; // # rows in each MPM
  NumericVector equivalence_vec; // equivalence vector if !stages_not_equal
  CharacterVector patch_vec; // choice of patch in each MPM
  IntegerVector tweights_type_vec; // tweights input as vector (1) or matrix (2) or null (0)
  IntegerVector total_years_vec; // total # years in each MPM
  IntegerVector format_vec; // MPM format (1:Ehrlen; 2:deVries; 3:ahist; 4:age-stage; 5: Leslie)
  IntegerVector entry_time_vec; // times of entry for each MPM
  IntegerVector dens_yn_vec; // density input for each MPM (0 = no, 1 = yes)
  IntegerVector dens_vr_yn_vec; // density_vr input for each MPM (0 = no, 1 = yes)
  IntegerVector sp_density_num_vec; // # of spatial density terms per MPM
  IntegerVector inda_terms_num_vec; // # of indcova (double) times per MPM
  IntegerVector indb_terms_num_vec; // # of indcovb (double) times per MPM
  IntegerVector indc_terms_num_vec; // # of indcovc (double) times per MPM
  IntegerVector inda_terms_cat_vec; // # of indcova (cat) times per MPM
  IntegerVector indb_terms_cat_vec; // # of indcovb (cat) times per MPM
  IntegerVector indc_terms_cat_vec; // # of indcovc (cat) times per MPM
  IntegerVector dev_terms_num_vec; // dev_terms rows input for each MPM
  NumericVector fecmod_vec; // fecundity multipliers for multiple offspring stages
  List fb_mpmout;
  
  List comm_out; // List of matrices of population vectors (top lvl: mpm, lower lvl: reps, mats: stages x times)
  List N_out;  // List of pop size matrices (top level: reps, mats: mpm by times)
  DataFrame labels;  // Data frame to provide order of MPMs
  List labels_list;  // List to hold data for data frame labels
  
  //Rcout << "a ";
  
  List cleaned_input = cleanup3(mpms, vrms, stageframes, supplements, format,
    firstage, finalage, fecage_min, fecage_max, cont, fecmod, starts, patches,
    years, tweights, density, entry_time, density_vr, sp_density, ind_terms,
    dev_terms, fb_sparse, equivalence, exp_tol, theta_tol, prep_mats, substoch,
    force_fb, stochastic, err_check_bool);
  
  //Rcout << "a1 ";
  
  mpm_list = as<List>(cleaned_input(0));
  mpm_count = static_cast<int>(cleaned_input(1));
  vrm_list = as<List>(cleaned_input(2));
  vrm_count = static_cast<int>(cleaned_input(3));
  A_list = as<List>(cleaned_input(4));
  stageframe_list = as<List>(cleaned_input(5));
  stageframe_list_fb = as<List>(cleaned_input(6));
  stageframe_count = static_cast<int>(cleaned_input(7));
  supplement_list = as<List>(cleaned_input(8));
  supplement_list_fb = as<List>(cleaned_input(9));
  supplement_count = static_cast<int>(cleaned_input(10));
  repmatrix_list = as<List>(cleaned_input(11));
  sparse_vec = as<LogicalVector>(cleaned_input(12));
  sparse_vec_count = static_cast<int>(cleaned_input(13));
  format_vec = as<IntegerVector>(cleaned_input(14));
  found_fleslie = static_cast<int>(cleaned_input(15));
  stageframe_notNull_count = static_cast<int>(cleaned_input(16));
  preexisting = static_cast<bool>(cleaned_input(17));
  funcbased = static_cast<bool>(cleaned_input(18));
  pure_fleslie = static_cast<bool>(cleaned_input(19));
  firstage_vec = as<IntegerVector>(cleaned_input(20));
  finalage_vec = as<IntegerVector>(cleaned_input(21));
  cont_vec = as<IntegerVector>(cleaned_input(22));
  fecmod_vec = as<NumericVector>(cleaned_input(23));
  fecage_min_vec = as<IntegerVector>(cleaned_input(24));
  fecage_max_vec = as<IntegerVector>(cleaned_input(25));
  hstages_list = as<List>(cleaned_input(26));
  agestages_list = as<List>(cleaned_input(27));
  matrowcounts = as<IntegerVector>(cleaned_input(28));
  stagecounts = as<IntegerVector>(cleaned_input(29));
  start_list = as<List>(cleaned_input(30));
  start_count = static_cast<int>(cleaned_input(31));
  labels_list = as<List>(cleaned_input(32));
  labels = as<DataFrame>(cleaned_input(33));
  patch_vec = as<CharacterVector>(cleaned_input(34));
  year_list = as<List>(cleaned_input(35));
  total_years_vec = as<IntegerVector>(cleaned_input(36));
  tweights_list = as<List>(cleaned_input(37));
  tweights_count = static_cast<int>(cleaned_input(38));
  tweights_type_vec = as<IntegerVector>(cleaned_input(39));
  density_list = as<List>(cleaned_input(40));
  dens_index_list = as<List>(cleaned_input(41));
  dens_yn_vec = as<IntegerVector>(cleaned_input(42));
  density_count = static_cast<int>(cleaned_input(43));
  entry_time_vec = as<IntegerVector>(cleaned_input(44));
  entry_time_count = static_cast<int>(cleaned_input(45));
  entry_time_vec_use = static_cast<bool>(cleaned_input(46));
  density_vr_list = as<List>(cleaned_input(47));
  ind_terms_num_list = as<List>(cleaned_input(48));
  ind_terms_cat_list = as<List>(cleaned_input(49));
  dev_terms_list = as<List>(cleaned_input(50));
  dens_vr_yn_vec = as<IntegerVector>(cleaned_input(51));
  sp_density_num_vec = as<IntegerVector>(cleaned_input(52));
  dev_terms_num_vec = as<IntegerVector>(cleaned_input(53));
  inda_terms_num_vec = as<IntegerVector>(cleaned_input(54));
  indb_terms_num_vec = as<IntegerVector>(cleaned_input(55));
  indc_terms_num_vec = as<IntegerVector>(cleaned_input(56));
  inda_terms_cat_vec = as<IntegerVector>(cleaned_input(57));
  indb_terms_cat_vec = as<IntegerVector>(cleaned_input(58));
  indc_terms_cat_vec = as<IntegerVector>(cleaned_input(59)); 
  sparse_vec = as<LogicalVector>(cleaned_input(60));
  density_vr_count = static_cast<int>(cleaned_input(61));
  sparse_vec_count = static_cast<int>(cleaned_input(62));
  sp_density_list = as<List>(cleaned_input(63));
  equivalence_list = as<List>(cleaned_input(64));
  equivalence_vec = as<NumericVector>(cleaned_input(65));
  equivalence_count = static_cast<int>(cleaned_input(66));
  stages_not_equal = static_cast<bool>(cleaned_input(67));
  allstages_all = as<List>(cleaned_input(68));
  allmodels_all = as<List>(cleaned_input(69));
  fb_override = static_cast<bool>(cleaned_input(70));
  fb_mpmout = as<List>(cleaned_input(71));
  
  List final_out_matrices;
  total_mpms = mpm_count + vrm_count;
  
  //Rcout << "b ";
  
  // Projection runs
  
  // nreps
  if (preexisting || fb_override) {
    project3_pre_core (N_out, comm_out, extreme_mpm_out, mpm_list, A_list,
      tweights_list, start_list, vrm_list, stageframe_list, allmodels_all,
      allstages_all, supplement_list, year_list, ind_terms_num_list,
      ind_terms_cat_list, dev_terms_list, density_vr_list, sp_density_list,
      density_list, dens_index_list, equivalence_list, dev_terms_num_vec,
      sp_density_num_vec, firstage_vec, finalage_vec, stagecounts,
      entry_time_vec, format_vec, inda_terms_num_vec, indb_terms_num_vec,
      indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
      indc_terms_cat_vec, dens_yn_vec, dens_vr_yn_vec, tweights_type_vec,
      fecmod_vec, sparse_vec, patch_vec, vrm_count, mpm_count, nreps, times,
      substoch, exp_tol, theta_tol, integeronly, stages_not_equal, stochastic, 
      entry_time_vec_use, err_check_bool, err_check_extreme);
    final_out_matrices = fb_mpmout;
  } else if (funcbased && !fb_override) {
    project3_fb_core (N_out, comm_out, extreme_mpm_out, err_check_fb_out,
      start_list, vrm_list, tweights_list, stageframe_list, allmodels_all, allstages_all,
      supplement_list, year_list, ind_terms_num_list, ind_terms_cat_list,
      dev_terms_list, density_vr_list, sp_density_list, density_list,
      dens_index_list, equivalence_list, dev_terms_num_vec, sp_density_num_vec,
      firstage_vec, finalage_vec, stagecounts, entry_time_vec, format_vec,
      inda_terms_num_vec, indb_terms_num_vec, indc_terms_num_vec,
      inda_terms_cat_vec, indb_terms_cat_vec, indc_terms_cat_vec, dens_yn_vec,
      dens_vr_yn_vec, tweights_type_vec, fecmod_vec, sparse_vec, patch_vec, vrm_count, nreps,
      times, substoch, exp_tol, theta_tol, integeronly, stages_not_equal,
      stochastic, err_check_bool, err_check_extreme);
    final_out_matrices = err_check_fb_out;
  } // funcbased processing
  
  //Rcout << "m ";
  
  int out_dim = 6;
  if (err_check_bool) out_dim++;
  
  List output (out_dim);
  output(0) = comm_out; // Needed in final output
  output(1) = N_out; // Needed in final output
  output(2) = stageframe_list; // Needed in final output
  output(3) = hstages_list;
  output(4) = agestages_list;
  output(5) = labels;
  
  if (err_check_extreme) {
    
    List output_errcheck (8);
    output_errcheck(0) = allstages_all;
    output_errcheck(1) = allmodels_all;
    output_errcheck(2) = equivalence_list;
    output_errcheck(3) = density_list;
    output_errcheck(4) = dens_index_list;
    output_errcheck(5) = density_vr_list;
    output_errcheck(6) = final_out_matrices;
    output_errcheck(7) = extreme_mpm_out;
    
    CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
      "equivalence_list", "density_list", "dens_index_list", "density_vr_list",
      "fb_mpm_out_matrices", "modified_mpms"};
    output_errcheck.attr("names") = output_errcheck_names;
    
    output(6) = output_errcheck;
    
    CharacterVector output_main_names = {"comm_out", "N_out", "stageframe_list",
      "hstages_list", "agestages_list", "labels", "err_check"};
    output.attr("names") = output_main_names;
    
  } else if (err_check_bool) {
    
    List output_errcheck (7);
    output_errcheck(0) = allstages_all;
    output_errcheck(1) = allmodels_all;
    output_errcheck(2) = equivalence_list;
    output_errcheck(3) = density_list;
    output_errcheck(4) = dens_index_list;
    output_errcheck(5) = density_vr_list;
    output_errcheck(6) = final_out_matrices;
    
    CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
      "equivalence_list", "density_list", "dens_index_list", "density_vr_list",
      "fb_mpm_out_matrices"};
    output_errcheck.attr("names") = output_errcheck_names;
    
    output(6) = output_errcheck;
    
    CharacterVector output_main_names = {"comm_out", "N_out", "stageframe_list",
      "hstages_list", "agestages_list", "labels", "err_check"};
    output.attr("names") = output_main_names;
    
  } else {
    CharacterVector output_main_names = {"comm_out", "N_out", "stageframe_list",
      "hstages_list", "agestages_list", "labels"};
    output.attr("names") = output_main_names;
  }
  
  output.attr("class") = "adaptProj";
  
  return output;
}

//' Clean Up RObject Inputs for Invasion Analysis
//' 
//' This function takes RObject inputs in the core projection functions, and
//' uses them to create the strict inputs for projection.
//' 
//' @name cleanup3_inv
//' 
//' @param mpm An MPM of class \code{lefkoMat}, for use if using existing MPMs.
//' @param vrm A \code{vrm_input} object corresponding to a distinct MPM that
//' will be created during analysis. Requires a stageframe, entered in argument
//' \code{stageframe}.
//' @param stageframe A stageframe defining stages and the life cycle for the
//' entered object in argument \code{vrms}. Must be of class \code{stageframe}.
//' @param supplement An optional data frame of class \code{lefkoSD} providing
//' supplemental data that should be incorporated into function-based MPMs. See
//' \code{\link[lefko3]{supplemental}()} for details. Use only with argument
//' \code{vrm}.
//' @param format An optional integer indicating the kind of function-based MPM
//' to create, if argument \code{vrm} is provided. Possible choices include:
//' \code{1}, Ehrlen-format historical MPM; \code{2}, deVries-format historical
//' MPM; \code{3}, ahistorical MPM (default); \code{4}, age-by-stage MPM; and
//' \code{5}, Leslie (age-based) MPM. Defaults to \code{3}.
//' @param firstage An optional integer used for function-based Leslie and
//' age-by-stage MPMs giving the starting ages in such MPMs. Use only if MPM is
//' both function-based and has age structure. Typically, the starting age in
//' such MPMs should be set to \code{0} if post-breeding and \code{1} if
//' pre-breeding. All other MPMs should be set to \code{0}. Do not use if not
//' using age structure. 
//' @param finalage An optional integer used for function-based Leslie and
//' age-by-stage MPMs giving the final ages in such MPMs. Use only if MPM is
//' both function-based and has age structure. Do not use if not using age
//' structure.
//' @param fecage_min An optional integer used for function-based Leslie MPMs
//' giving the first age at which organisms can be reproductive. Use only for
//' function-based Leslie MPMs. Defaults to the values given in \code{firstage}.
//' @param fecage_max An integer value used for function-based Leslie MPMs
//' giving the final age at which organisms can be reproductive. Use only for
//' function-based Leslie MPMs. Defaults to the values given in \code{finalage}.
//' @param cont An optional vector used for function-based Leslie and
//' age-by-stage MPMs stating whether the MPM should should include a stasis
//' transition within the final age. This should be used only when an organism
//' can maintain the demographic characteristics of the final described age
//' after reaching that age. Can be entered as a logical value or an integer. 
//' MPMs without age structure should be entered as \code{0} or \code{FALSE}.
//' Do not use if not using age structure.
//' @param fecmod An optional value used for function-based MPMs giving scalar
//' multipliers for fecundity terms, when two fecundity variables are used for a
//' collective fecundity per individual.
//' @param start An optional \code{lefkoSV} object, which is a data frame
//' providing the starting numbers of individuals of each stage. If not
//' provided, then all projections start with a single individual of each stage.
//' @param patch An optional string giving the name of the patch to project.
//' Only a single pop-patch may be projected. A value of \code{NA} can be
//' supplied to indicate that the population-level matrices should be projected
//' (if argument \code{mpm} is used and a population-level set of matrices
//' exist), or that the first patch noted should be used. Defaults to the
//' population-level set or the first patch, depending on whether the former
//' exists.
//' @param years An optional term corresponding to a single integer vector of
//' time \code{t} values. If a vector shorter than \code{times} is supplied,
//' then this vector will be cycled. Defaults to a vector of all detected
//' years in argument \code{mpm} or argument \code{vrm}.
//' @param tweights An optional numeric vector or matrix denoting the
//' probabilities of choosing each matrix in each MPM in a stochastic
//' projection. If a matrix, then a first-order Markovian environment is
//' assumed, in which the probability of choosing a specific annual matrix
//' depends on which annual matrix is currently chosen. If an element of the
//' list is a vector, then the choice of annual matrix is assumed to be
//' independent of the current matrix. Defaults to equal weighting among
//' matrices.
//' @param density An optional data frames of class \code{lefkoDens}, which
//' provides details for density dependence in MPM elements and is created with
//' function \code{\link[lefko3]{density_input}()}. Defaults to \code{NULL}, in
//' which case no density dependence is built into matrix elements.
//' @param entry_time An optional integer vector giving the entry time for each
//' variant into each simulation. Defaults to a zero vector with length equal to
//' the number of variants to run concurrently in each simulation, as given by
//' argument \code{var_per_run}.
//' @param density_vr An optional data frame of class \code{lefkoDensVR}, which
//' provides details for density dependence in vital rate models and has been
//' created with function \code{link[lefko3]{density_vr}()}. Can only be used
//' with function-based projections. Defaults to \code{NULL}, in which case no
//' density dependence is built into vital rates.
//' @param sp_density An optional argument for use with argument \code{vrm} that
//' specifies the spatial density to be used in each time step. If used, then
//' may either be a numeric vector giving a single spatial density for each
//' time step. If vectors are shorter than specified in \code{times}, then these
//' values will be cycled.
//' @param ind_terms An optional argument providing values of individual or
//' environmental covariate values for argument \code{vrm}. Should be set to a
//' single data frame with 3 columns giving values for up to 3 covariates across
//' time (rows give the time order of these values). Unused terms within the
//' data frame must be set to \code{0} (use of \code{NA} will produce errors).
//' If the number of rows is less than \code{times}, then these values will be
//' cycled.
//' @param dev_terms An optional  data frame including 14 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival, 9: juvenile
//' observation, 10: juvenile primary size, 11: juvenile secondary size,
//' 12: juvenile tertiary size, 13: juvenile reproduction, and 14: juvenile
//' maturity transition. Unused terms must be set to \code{0} (use of \code{NA}
//' will produce errors). Single or small numbers of values per vital rate model
//' are also allowed, and if the number of rows is less than \code{times}, then
//' the terms will be cycled.
//' @param fb_sparse A logical vector indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE} for
//' each MPM.
//' @param equivalence An optional object of class \code{adaptEq} giving the
//' degree to which individuals in each stage are equivalent to one another.
//' May also be a numeric vector, in which case the vector must have the same
//' number of elements as the number of rows in the associated MPM, with each
//' element giving the effect of an individual of that age, stage, age-stage, or
//' stage-pair, depending on whether the MPM is age-based, ahistorical
//' stage-based, age-by-stage, or historical stage-based, respectively. Numeric
//' entries used in these vectors can be thought of as Lotka-Volterra
//' interaction terms, such as are used in multiple species competition models.
//' @param prebreeding An optional value stating whether the life cycle is
//' prebreeding. If no value is entered for \code{firstage}, then a value of
//' \code{TRUE} sets the minimum age to 0, and a value of \code{FALSE} sets it
//' to 1.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param theta_tol A numeric value used to indicate a maximum value to theta as
//' used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent and density
//' independent simulations. Defaults to \code{0}, which does not enforce
//' substochasticity. Alternatively, \code{1} forces all survival-transition
//' elements to range from 0.0 to 1.0, and forces fecundity to be non-negative;
//' and \code{2} forces all column rows in the survival-transition matrices to
//' total no more than 1.0, in addition to the actions outlined for option
//' \code{1}. Both settings \code{1} and \code{2} change negative fecundity
//' elements to \code{0.0}.
//' @param variant_count An integer giving the number of variants total.
//' @param var_per_run The number of variants to run in each simulation.
//' Defaults to \code{2}, resulting in pairwise invasibility analysis. See
//' \code{Notes} for details.
//' 
//' @return A list of R-defined objects, including vectors, lists, integers, and
//' data frames, for use in later stages of analysis.
//' 
//' @keywords internal
//' @noRd
Rcpp::List cleanup3_inv (Nullable<RObject> mpm = R_NilValue,
  Nullable<RObject> vrm = R_NilValue, Nullable<RObject> stageframe = R_NilValue,
  Nullable<RObject> supplement = R_NilValue, Nullable<RObject> format = R_NilValue,
  Nullable<RObject> firstage = R_NilValue, Nullable<RObject> finalage = R_NilValue,
  Nullable<RObject> fecage_min = R_NilValue, Nullable<RObject> fecage_max = R_NilValue,
  Nullable<RObject> cont = R_NilValue, Nullable<RObject> fecmod = R_NilValue,
  Nullable<RObject> start = R_NilValue, Nullable<RObject> patch = R_NilValue,
  Nullable<RObject> years = R_NilValue, Nullable<RObject> tweights = R_NilValue,
  Nullable<RObject> density = R_NilValue, Nullable<RObject> entry_time = R_NilValue,
  Nullable<RObject> density_vr = R_NilValue, Nullable<RObject> sp_density = R_NilValue,
  Nullable<RObject> ind_terms = R_NilValue, Nullable<RObject> dev_terms = R_NilValue,
  Nullable<RObject> fb_sparse = R_NilValue, Nullable<RObject> equivalence = R_NilValue,
  Nullable<RObject> prebreeding = R_NilValue, double exp_tol = 700.0,
  double theta_tol = 100000000.0, const int substoch = 0,
  const unsigned int variant_count = 1, const unsigned int var_per_run = 2) {
  
  List chosen_mpm;
  List vrm_list;
  List start_list (variant_count);
  DataFrame stageframe_df;
  DataFrame final_stageframe;
  DataFrame supplement_df;
  arma::mat final_repmatrix;
  DataFrame final_hstages;
  DataFrame final_agestages;
  DataFrame final_labels; // This and the next two might be redundant
  DataFrame labels;
  CharacterVector labels_list;
  List tweights_list (1);
  DataFrame density_df;
  DataFrame dens_index_df; // Used to be list by variant_count
  DataFrame chosen_density_vr;
  List ind_terms_num_list;
  List ind_terms_cat_list;
  List dev_terms_list;
  List sp_density_list;
  DataFrame equivalence_frame;
  DataFrame current_mpm_allstages;
  List allstages_all;
  List allmodels_all;
  
  CharacterVector patch_vec (1);
  CharacterVector year_vec; // CharacterVector of user input
  CharacterVector existing_years; // CharacterVector of years in MPM
  IntegerVector entry_time_vec (var_per_run);
  NumericVector equivalence_vec;
  
  int format_int {3};
  int stagecounts {0};
  int matrowcounts {0};
  unsigned int firstage_int {0};
  unsigned int finalage_int {0};
  unsigned int fecage_min_int {0};
  unsigned int fecage_max_int {0};
  unsigned int cont_int {0};
  double fecmod_num {1.0};
  int total_years_int {0};
  int stageframe_notNull_count {0};
  int sparse_vec_count {0};
  int start_count {0};
  int tweights_type_int {0};
  int tweights_count {0};
  int density_count {0};
  int sp_density_num_int {0};
  int entry_time_count {0};
  int density_vr_count {0};
  int equivalence_count {0};
  int dev_terms_num_int {0};
  int inda_terms_num_int {0};
  int indb_terms_num_int {0};
  int indc_terms_num_int{0};
  int inda_terms_cat_int {0};
  int indb_terms_cat_int {0};
  int indc_terms_cat_int {0}; 
  
  int preexisting_mpm_size {0};
  
  bool preexisting {false};
  bool funcbased {false};
  bool pure_leslie {false};
  bool entry_time_vec_use {false};
  bool stages_not_equal {false};
  bool trial_supp_null {false};
  bool dens_yn_bool {0};
  bool dens_vr_yn_bool {false};
  bool sparse_bool {false};
  bool prebreeding_bool {true};
  
  bool historical {false};
  
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Argument substoch must equal 0, 1, or 2.", false);
  }
  
  //Rcout << "cleanup3_inv A ";
  
  if (mpm.isNotNull()) {
    if (vrm.isNotNull() || stageframe.isNotNull()) {
      AdaptUtils::pop_error2("vrm", "stageframe", "projecting existing MPMs", 24);
    }
    
    if (is<List>(mpm)) {
      chosen_mpm = as<List>(mpm);
      
      if (chosen_mpm.hasAttribute("class")) {
        CharacterVector mpm_classes = wrap(chosen_mpm.attr("class"));
        
        bool found_lefkoMat {false};
        for (int i = 0; i < mpm_classes.length(); i++) {
          if (stringcompare_hard(String(mpm_classes(i)), "lefkoMat")) found_lefkoMat = true;
        }
        
        if (!found_lefkoMat) {
          AdaptUtils::pop_error2("mpm", "a lefkoMat object", "", 1);
        }
      } else {
        AdaptUtils::pop_error2("mpm", "a lefkoMat object", "", 1);
      }
      
    } else {
      AdaptUtils::pop_error2("mpm", "a lefkoMat object", "", 1);
    }
    
    stageframe_df = as<DataFrame>(chosen_mpm["ahstages"]);
    
    IntegerVector sf_min_age = as<IntegerVector>(stageframe_df["min_age"]);
    int min_age_length = static_cast<int>(sf_min_age.length());
    
    DataFrame sf_agestages_check = as<DataFrame>(chosen_mpm["agestages"]);
    int sfac_vars = static_cast<int>(sf_agestages_check.length());
    
    if (sfac_vars < 3) {
      //bool age_check {false};
      int orig_age {0};
      for (unsigned int sfi = 0; sfi < min_age_length; sfi++) {
        if (!IntegerVector::is_na(sf_min_age(sfi)) && !NumericVector::is_na(sf_min_age(sfi))) {
          int current_age = sf_min_age(sfi);
          if (current_age == (orig_age + 1)) {
            pure_leslie = true;
            format_int = 5;
          }
          orig_age = current_age;
        }
      }
    }
    List A_list = as<List>(chosen_mpm["A"]);
    if (is<S4>(A_list(0))) {
      sparse_bool = true;
      arma::sp_mat A_1 = as<arma::sp_mat>(A_list(0));
      preexisting_mpm_size = static_cast<int>(A_1.n_elem);
      
    } else {
      arma::mat A_1 = as<arma::mat>(A_list(0));
      preexisting_mpm_size = static_cast<int>(A_1.n_elem);
    }
    //stageframe_count = 1;
    sparse_vec_count = 1;
    preexisting = true;
  }
  
  //Rcout << "cleanup3_inv B ";
  
  if (vrm.isNotNull()) {
    if (mpm.isNotNull()) {
      throw Rcpp::exception("Function invade3 handles a single lefkoMat or a single vrm_input object only.", false);
    }
    
    if (is<List>(vrm)) {
      vrm_list = as<List>(vrm);
      
      if (format.isNotNull()) {
        if (is<NumericVector>(format) || is<IntegerVector>(format)) {
          IntegerVector format_vec = as<IntegerVector>(format);
          int format_count = static_cast<int>(format_vec.length());
          
          if (format_count != 1) {
            AdaptUtils::pop_error2("format", "a single integer", "", 1);
          }
          
          if (IntegerVector::is_na(format_vec(0))) {
            AdaptUtils::pop_error2("NA values", "format", "", 25);
          }
          
          format_int = static_cast<int>(format_vec(0));
          if (format_int != 5 && !stageframe.isNotNull()) {
            AdaptUtils::pop_error2("stageframe", "run function-based projections", "", 26);
          } else if (format_int == 5 && !stageframe.isNotNull()) {
            pure_leslie = true;
          }
        } else AdaptUtils::pop_error2("format", "a single integer", "", 1);
      } else if (funcbased) {
        if (!stageframe.isNotNull()) {
          AdaptUtils::pop_error2("stageframe", "run function-based projections", "", 26);
        }
      }
    } else AdaptUtils::pop_error2("vrm", "a vrm_input object", "", 1);
    
    if (stageframe.isNotNull()) {
      if (is<DataFrame>(stageframe)) {
        stageframe_df = as<DataFrame>(stageframe);
        //stageframe_count = 1;
        
      } else {
        if (format_int != 5) AdaptUtils::pop_error2("stageframe", "a stageframe object", "", 1);
        pure_leslie = true;
      }
      
      if (!vrm_list.hasAttribute("class")) {
        AdaptUtils::pop_error2("vrm", "a vrm_input object", "", 1);
      }
      CharacterVector chosen_vrm_class = wrap(vrm_list.attr("class"));
      
      bool found_vrmi {false};
      for (int j = 0; j < static_cast<int>(chosen_vrm_class.length()); j++) {
        if (chosen_vrm_class(j) == "vrm_input") found_vrmi = true;
      }
      
      if (!found_vrmi) AdaptUtils::pop_error2("vrm", "a vrm_input object", "", 1);
      
      if (!stageframe_df.hasAttribute("class")) {
        AdaptUtils::pop_error2("stageframe", "a stageframe object", "", 1);
      }
      CharacterVector chosen_stageframe_class = wrap(stageframe_df.attr("class"));
      
      bool found_stageframe {false};
      for (int j = 0; j < static_cast<int>(chosen_stageframe_class.length()); j++) {
        if (chosen_stageframe_class(j) == "stageframe") found_stageframe = true;
      }
      
      if (!found_stageframe) {
        AdaptUtils::pop_error2("stageframe", "a stageframe object", "", 1);
      }
      stageframe_notNull_count++;
    } else {
      if (format_int != 5) {
        throw Rcpp::exception("All non-Leslie MPMs need stageframes.", false);
      }
    }
    
    if (stageframe_notNull_count != 1 && !pure_leslie) {
      throw Rcpp::exception("Each vrm_input object must have its own stageframe.",
        false);
    }
    funcbased = true;
  }
  
  //Rcout << "cleanup3_inv C ";
  
  if (!preexisting && !funcbased) {
    throw Rcpp::exception("Cannot proceed without either argument mpms, or arguments vrms and stageframes set.",
      false);
  } else if (preexisting && funcbased) {
    throw Rcpp::exception("Cannot proceed with argument mpms, vrms, and stageframes set.",
      false);
  }
  
  LefkoInputs::RObj_DFr_input_check ("supplement", "lefkoSD", supplement_df,
    trial_supp_null, true, false, supplement);
  
  if (!funcbased && !trial_supp_null) {
    throw Rcpp::exception("Argument supplement can only be used with argument vrm.",
      false);
  }
  
  //Rcout << "cleanup3_inv D ";
  
  if (format.isNotNull()) {
    if (!funcbased) {
      AdaptUtils::pop_error2("vrm", "use argument format", "", 26);
    }
    
    if (is<IntegerVector>(format) || is<NumericVector>(format)) {
      IntegerVector format_ = as<IntegerVector>(format);
      
      format_int = format_(0);
      if (format_int < 1 || format_int > 5) {
        AdaptUtils::pop_error2("format", "an integer between 1 and 5", "", 1);
      }
    } else if (is<StringVector>(format)) {
      StringVector format_ = as<StringVector>(format);
      String format_0 = String(format_(0));
      
      if (LefkoUtils::stringcompare_simple(format_0, "hist", true)) {
        format_int = 1;
      } else if (LefkoUtils::stringcompare_simple(format_0, "ehrl", true)) {
        format_int = 1;
      } else if (LefkoUtils::stringcompare_simple(format_0, "dev", true)) {
        format_int = 2;
      } else if (LefkoUtils::stringcompare_simple(format_0, "agest", true)) {
        format_int = 4;
      } else if (LefkoUtils::stringcompare_simple(format_0, "by", true)) {
        format_int = 4;
      } else if (LefkoUtils::stringcompare_simple(format_0, "age", true)) {
        format_int = 5;
      } else if (LefkoUtils::stringcompare_simple(format_0, "stag", true)) {
        format_int = 3;
      } else {
        Rf_warningcall(R_NilValue,
          "Argument format not understood. Defaulting to stage-based ahistorical format.");
        format_int = 3;
      }
    }
    
  } else if (preexisting) {
    RObject hstages_element = as<RObject>(chosen_mpm["hstages"]);
    RObject agestages_element = as<RObject>(chosen_mpm["agestages"]);
    DataFrame core_ahstages = as<DataFrame>(chosen_mpm["ahstages"]);
    
    if (!is<LogicalVector>(hstages_element)) {
      if (is<DataFrame>(hstages_element)) {
        DataFrame hst_input = as<DataFrame>(hstages_element);
        int hst_cols = hst_input.length();
        
        if (hst_cols > 1) {
          int hst_rows = static_cast<int>(hst_input.nrows());
          int no_stages_found = static_cast<int>(core_ahstages.nrows());
          int expected_ehrlen = no_stages_found * no_stages_found;
          
          format_int = 1;
          if (hst_rows < expected_ehrlen) format_int++;
          historical = true;
        }
      }
    }
    
    if (!is<LogicalVector>(agestages_element) && format_int == 3) {
      if (is<DataFrame>(agestages_element)) {
        DataFrame ast_input = as<DataFrame>(agestages_element);
        int ast_cols = ast_input.length();
        
        if (ast_cols > 1) format_int = 4;
      }
    }
  } else if (funcbased) {
    format_int = 3;
  }
  
  //Rcout << "cleanup3_inv E ";
  
  // firstage, finalage, and cont processing for age-by-stage and Leslie MPMs
  if (firstage.isNotNull() || finalage.isNotNull() || cont.isNotNull() || 
      prebreeding.isNotNull()) {
    if (!funcbased) {
      AdaptUtils::pop_error2("vrm", "use arguments firstage, finalage, cont, and prebreeding", "", 26);
    }
    
    bool found_age_MPM {false};
    if (format_int > 3) found_age_MPM = true;
    
    if (!found_age_MPM) {
      AdaptUtils::pop_error2("firstage, finalage, and cont", "age-, function-based MPMs", "", 28);
    }
  }
  
  if (firstage.isNotNull()) {
    if (is<IntegerVector>(firstage) || is<NumericVector>(firstage)) { 
      IntegerVector firstage_vec = as<IntegerVector>(firstage);
      
      int firstage_vec_length = static_cast<int>(firstage_vec.length());
      if (firstage_vec_length != 1) {
        AdaptUtils::pop_error2("firstage", "a single integer", "", 1);
      }
      if (IntegerVector::is_na(firstage_vec(0))) {
        AdaptUtils::pop_error2("NA values", "firstage", "", 25);
      }
      
      firstage_int = static_cast<unsigned int>(firstage_vec(0));
      
      if (firstage_int < 0) {
        AdaptUtils::pop_error2("firstage", "", "", 30);
      }
      
      if (firstage_int > 1 && format_int < 4) {
        throw Rcpp::exception("Entries in argument firstage must equal 0 or 1 for MPMs without age structure.", false);
      }
      
      if (prebreeding.isNotNull()) {
        throw Rcpp::exception("Do not use both arguments prebreeding and first age.", false);
      }
    } else AdaptUtils::pop_error2("firstage", "a single integer", "", 1);
  } else {
    if (preexisting) {
      DataFrame found_agestages = as<DataFrame>(chosen_mpm["agestages"]);
      if (found_agestages.length() > 1) {
        IntegerVector found_agestages_age = found_agestages["age"];
        firstage_int = min(found_agestages_age);
      }
    } else if (prebreeding.isNotNull()) {
      if (is<LogicalVector>(prebreeding)) { 
        LogicalVector prebreeding_vec = as<LogicalVector>(prebreeding);
        
        int prebreeding_vec_length = static_cast<int>(prebreeding_vec.length());
        if (prebreeding_vec_length != 1) {
          AdaptUtils::pop_error2("prebreeding", "a single true or false value", "", 1);
        }
        if (LogicalVector::is_na(prebreeding_vec(0))) {
          AdaptUtils::pop_error2("NA values", "prebreeding", "", 25);
        }
        
        prebreeding_bool = static_cast<bool>(prebreeding_vec(0));
        
        if (prebreeding_bool) {
          firstage_int = 0;
        } else {
          firstage_int = 1;
        }
        
      } else AdaptUtils::pop_error2("prebreeding", "a single true or false value", "", 1);
    }
  }
  
  if (finalage.isNotNull()) {
    if (is<IntegerVector>(finalage) || is<NumericVector>(finalage)) { 
      IntegerVector finalage_vec = as<IntegerVector>(finalage);
      
      int finalage_vec_length = static_cast<int>(finalage_vec.length());
      if (finalage_vec_length != 1) {
        AdaptUtils::pop_error2("finalage", "a single integer", "", 1);
      }
      if (IntegerVector::is_na(finalage_vec(0))) {
        AdaptUtils::pop_error2("NA values", "finalage", "", 25);
      }
      
      finalage_int = static_cast<unsigned int>(finalage_vec(0));
      
      if (finalage_int < 0) {
        AdaptUtils::pop_error2("finalage", "", "", 30);
      }
      
      if (finalage_int > 0 && format_int < 4) {
        throw Rcpp::exception("Entries in argument finalage must equal 0 for MPMs without age structure.", false);
      }
      
      if (finalage_int < firstage_int) {
        throw Rcpp::exception("Argument finalage may not be less than value in argument firstage.", false);
      }
    } else AdaptUtils::pop_error2("finalage", "a single integer", "", 1);
  } else {
    if (preexisting) {
      DataFrame found_agestages = as<DataFrame>(chosen_mpm["agestages"]);
      if (found_agestages.length() > 1) {
        IntegerVector found_agestages_age = found_agestages["age"];
        finalage_int = max(found_agestages_age) + 1;
      }
    }
  }
  
  //Rcout << "cleanup3_inv F ";
  
  if (fecage_min.isNotNull()) {
    if (is<IntegerVector>(fecage_min) || is<NumericVector>(fecage_min)) {
      IntegerVector fecage_min_prevec = as<IntegerVector>(fecage_min);
      
      int fecage_min_length = static_cast<int>(fecage_min_prevec.length());
      if (fecage_min_length != 1) {
        AdaptUtils::pop_error2("fecage_min", "a single integer", "", 1);
      }
      if (IntegerVector::is_na(fecage_min_prevec(0))) {
        AdaptUtils::pop_error2("NA values", "fecage_min", "", 25);
      }
      
      fecage_min_int = static_cast<unsigned int>(fecage_min_prevec(0));
      
      if (fecage_min_int < 0) {
        AdaptUtils::pop_error2("fecage_min", "", "", 30);
      }
    } else AdaptUtils::pop_error2("fecage_min", "a single integer", "", 1);
  } else {
    fecage_min_int = firstage_int;
  }
  
  if (fecage_max.isNotNull()) {
    if (is<IntegerVector>(fecage_max) || is<NumericVector>(fecage_max)) {
      IntegerVector fecage_max_prevec = as<IntegerVector>(fecage_max);
      
      int fecage_max_length = static_cast<int>(fecage_max_prevec.length());
      if (fecage_max_length != 1) {
        AdaptUtils::pop_error2("fecage_max", "a single integer", "", 1);
      }
      if (IntegerVector::is_na(fecage_max_prevec(0))) {
        AdaptUtils::pop_error2("NA values", "fecage_max", "", 25);
      }
      
      fecage_max_int = static_cast<unsigned int>(fecage_max_prevec(0));
      
      if (fecage_max_int < 0) {
        AdaptUtils::pop_error2("fecage_max", "", "", 30);
      }
      
      if (fecage_max_int < fecage_min_int) {
        throw Rcpp::exception("Argument fecage_max may not be less than value in argument fecage_min.", false);
      }
    } else AdaptUtils::pop_error2("fecage_max", "a single integer", "", 1);
  } else {
    fecage_max_int = finalage_int;
  }
  
  if (fecmod.isNotNull()) {
    if (is<NumericVector>(fecmod) && funcbased) {
      NumericVector fecmod_pre = as<NumericVector>(fecmod);
      
      if (static_cast<int>(fecmod_pre.length()) != 1) {
        AdaptUtils::pop_error2("fecmod", "number of vrm inputs", "", 29);
      }
      if (NumericVector::is_na(fecmod_pre(0))) {
        AdaptUtils::pop_error2("NA values", "fecmod", "", 25);
      }
      fecmod_num = fecmod_pre(0);
      
    } else if (!funcbased) {
      AdaptUtils::pop_error2("vrm", "use argument fecmod", "", 26);
    }
  }
  
  //Rcout << "cleanup3_inv G ";
  
  if (cont.isNotNull()) {
    if (is<IntegerVector>(cont) || is<NumericVector>(cont)) { 
      IntegerVector cont_vec_temp = as<IntegerVector>(cont);
      
      int cont_vec_length = static_cast<int>(cont_vec_temp.length());
      if (cont_vec_length != 1) {
        AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
      }
      if (!IntegerVector::is_na(cont_vec_temp(0))) {
        cont_int = cont_vec_temp(0);
      }
      
      if (cont_int < 0 || cont_int > 1) {
        throw Rcpp::exception("Entries in argument cont must equal 0 or 1.",
          false);
      }
      
      if (cont_int > 0 && format_int < 4) {
        throw Rcpp::exception("Entries in argument cont must equal 0 for MPMs without age structure.", 
          false);
      }
    } if (is<LogicalVector>(cont)) {
      LogicalVector cont_vec_temp = as<LogicalVector>(cont);
      
      int cont_vec_log_length = static_cast<int>(cont_vec_temp.length());
      if (cont_vec_log_length != 1) {
        AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
      }
      
      if (cont_vec_temp(0) > 0 && !LogicalVector::is_na(cont_vec_temp[0])) cont_int = 1;
      
    } else {
      AdaptUtils::pop_error2("cont", "a single integer", "", 1);
    }
  }
  
  //Rcout << "cleanup3_inv H ";
  
  // Altered stageframe processing
  if (funcbased) {
    if (format_int < 5) {
      bool agemat = false;
      int ehrlen {1};
      //int style {0};
      //int filter {1};
      
      if (format_int == 2) ehrlen = 2;
      //if (format_int == 3) style = 1;
      if (format_int == 4) {
        agemat = true;
        //style = 2;
        //filter = 2;
      }
      if (format_int < 3) historical = true;
      
      List melchett;
      if (!trial_supp_null) {
        melchett = LefkoMats::sf_reassess_internal(stageframe_df, supplement_df,
          R_NilValue, R_NilValue, agemat, historical, ehrlen);
      } else {
        melchett = LefkoMats::sf_reassess_internal(stageframe_df, R_NilValue,
          R_NilValue, R_NilValue, agemat, historical, ehrlen);
      }
      final_stageframe = as<DataFrame>(melchett["stageframe"]);
      final_repmatrix = as<arma::mat>(melchett["repmatrix"]);
      
      if (format_int < 4) {
        DataFrame new_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
        if (new_ovtable_temp.containsElementNamed("stage3")) {
          supplement_df = new_ovtable_temp;
        } else {
          StringVector nsst3 = {};
          IntegerVector nsa2 = {};
          NumericVector nsgr = {};
          
          DataFrame intro_ovtable = DataFrame::create(_["stage3"] = nsst3,
            _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
            _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
            _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
            _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
            _["multiplier"] = clone(nsgr), _["convtype"] = clone(nsa2),
            _["convtype_t12"] = clone(nsa2), _["pop"] = clone(nsst3),
            _["patch"] = clone(nsst3), _["year2"] = clone(nsst3));
          supplement_df = intro_ovtable;
        }
      } else {
        DataFrame new_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
        if (new_ovtable_temp.containsElementNamed("stage3")) {
          supplement_df = LefkoMats::age_expanded(new_ovtable_temp,
            firstage_int, finalage_int);
        } else {
          StringVector nsst3 = {};
          IntegerVector nsa2 = {};
          NumericVector nsgr = {};
          
          DataFrame intro_ovtable = DataFrame::create(_["stage3"] = nsst3,
            _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
            _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
            _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
            _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
            _["multiplier"] = clone(nsgr), _["convtype"] = clone(nsa2),
            _["convtype_t12"] = clone(nsa2), _["pop"] = clone(nsst3),
            _["patch"] = clone(nsst3), _["year2"] = clone(nsst3));
          supplement_df = intro_ovtable;
        }
      }
      
      DataFrame chosen_stageframe_pre = clone(final_stageframe);
      
      IntegerVector removal_row = {static_cast<int>(chosen_stageframe_pre.nrows())};
      StringVector removal_var = {"stage_id"};
      DataFrame chosen_stageframe = LefkoUtils::df_remove(chosen_stageframe_pre,
        removal_row, false, true, false, false, true, as<RObject>(removal_var));
      
      if (format_int < 3) {
        hst_maker(final_hstages, chosen_stageframe, format_int);
        
        matrowcounts = static_cast<int>(final_hstages.nrows());
      } else if (format_int == 4) {
        DataFrame agestages_temp = age_maker(chosen_stageframe,
          firstage_int, finalage_int);
        final_agestages = agestages_temp;
        
        matrowcounts = static_cast<int>(agestages_temp.nrows());
      } else {
        matrowcounts = static_cast<int>(chosen_stageframe.nrows());
      }
    } else {
      bool cont_used {false};
      if (cont_int > 0) cont_used = true;
      
      DataFrame melchett = LefkoMats::sf_leslie(firstage_int, finalage_int,
        fecage_min_int, fecage_max_int, cont_used);
      DataFrame new_stageframe = melchett;
      
      DataFrame new_ovtable;
      if (!trial_supp_null) {
        new_ovtable = LefkoMats::age_expanded(supplement_df, firstage_int,
          finalage_int);
      }
      final_stageframe = new_stageframe;
      supplement_df = new_ovtable;
      
      //stageframe_count++;
    }
  } else {
    if (format_int < 3) {
      final_hstages = as<DataFrame>(chosen_mpm["hstages"]);
      
      matrowcounts = static_cast<int>(final_hstages.nrows());
    } else if (format_int == 4) {
      final_agestages = as<DataFrame>(chosen_mpm["agestages"]);
      
      matrowcounts = static_cast<int>(final_agestages.nrows());
      
      IntegerVector all_ages_agestages = final_agestages["age"];
      int min_age_agestages = min(all_ages_agestages);
      int max_age_agestages = max(all_ages_agestages);
      
      firstage_int = min_age_agestages;
      finalage_int = max_age_agestages;
      
    } else {
      DataFrame chosen_ahstages = as<DataFrame>(chosen_mpm["ahstages"]);
      
      matrowcounts = static_cast<int>(chosen_ahstages.nrows());
      final_stageframe = chosen_ahstages; // Could probably also be stageframe_df
    }
  }
  
  //Rcout << "cleanup3_inv I ";
  
  // start vector
  if (start.isNotNull()) {
    start_count = 1;
    
    if (is<DataFrame>(start)) {
      DataFrame chosen_start = as<DataFrame>(start);
      
      if (!chosen_start.hasAttribute("class")) {
        AdaptUtils::pop_error2("start", "a lefkoSV object", "", 1);
      }
      CharacterVector chosen_start_class = chosen_start.attr("class");
      
      bool found_lSt {false};
      for (int j = 0; j < static_cast<int>(chosen_start_class.length()); j++) {
        if (chosen_start_class(j) == "lefkoSV") found_lSt = true;
      }
      
      if (!found_lSt) {
        AdaptUtils::pop_error2("start", "a lefkoSV object", "", 1);
      }
      
      DataFrame chosen_stageframe;
      if (format_int < 3) {
        chosen_stageframe = final_hstages;
      } else if (format_int == 4) {
        chosen_stageframe = final_agestages;
      } else {
        chosen_stageframe = final_stageframe;
      }
      
      stagecounts = static_cast<int>(chosen_stageframe.nrows());
      if (format_int == 3 && funcbased) stagecounts--;
      
      arma::vec start_vec (stagecounts, fill::zeros);
      arma::uvec start_elems = as<arma::uvec>(chosen_start["row_num"]);
      start_elems = start_elems - 1;
      arma::vec start_values = as<arma::vec>(chosen_start["value"]);
      
      if (static_cast<int>(start_elems.max()) > (stagecounts - 1)) {
        throw Rcpp::exception("lefkoStart object includes element indices too high for associated MPM.",
          false);
      }
      
      for (int j = 0; j < static_cast<int>(start_elems.n_elem); j++) {
        start_vec(start_elems(j)) = start_values(j);
      }
      
      for (int i = 0; i < variant_count; i++) { 
        start_list(i) = start_vec;
      }
    } else if (is<List>(start)) {
      
      List start_list_pre = as<List>(start);
      int start_list_length = static_cast<int>(start_list_pre.length());
      
      if (start_list_length != 1 && start_list_length != variant_count) { 
        String eat_my_shorts = "Argument start should be composed of a single ";
        eat_my_shorts += "lefkoSV object, or a list of as many lefkoSV ";
        eat_my_shorts += "objects as rows in the data frame used in argument axis";
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
      
      for (int slelem = 0; slelem < start_list_length; slelem++) {
        if (!is<DataFrame>(start_list_pre(slelem))) {
          AdaptUtils::pop_error2("start", "a lefkoSV object or list of lefkoSV objects", "", 1);
        }
        
        DataFrame chosen_start = as<DataFrame>(start_list_pre(slelem));
        
        if (!chosen_start.hasAttribute("class")) {
          AdaptUtils::pop_error2("start", "a lefkoSV object or list of lefkoSV objects", "", 1);
        }
        CharacterVector chosen_start_class = chosen_start.attr("class");
        
        bool found_lSt {false};
        for (int j = 0; j < static_cast<int>(chosen_start_class.length()); j++) {
          if (chosen_start_class(j) == "lefkoSV") found_lSt = true;
        }
        
        if (!found_lSt) {
          AdaptUtils::pop_error2("start", "a lefkoSV object or list of lefkoSV objects", "", 1);
        }
        
        DataFrame chosen_stageframe;
        if (format_int < 3) {
          chosen_stageframe = final_hstages;
        } else if (format_int == 4) {
          chosen_stageframe = final_agestages;
        } else {
          chosen_stageframe = final_stageframe;
        }
        
        stagecounts = static_cast<int>(chosen_stageframe.nrows());
        if (format_int == 3 && funcbased) stagecounts--;
        
        arma::vec start_vec (stagecounts, fill::zeros);
        arma::uvec start_elems = as<arma::uvec>(chosen_start["row_num"]);
        start_elems = start_elems - 1;
        arma::vec start_values = as<arma::vec>(chosen_start["value"]);
        
        if (static_cast<int>(start_elems.max()) > (stagecounts - 1)) {
          throw Rcpp::exception("lefkoStart object includes element indices too high for associated MPM.",
            false);
        }
        
        for (int j = 0; j < static_cast<int>(start_elems.n_elem); j++) {
          start_vec(start_elems(j)) = start_values(j);
        }
        
        start_list(slelem) = start_vec;
      }
    } else {
      AdaptUtils::pop_error2("start", "a lefkoSV object or list of lefkoSV objects", "", 1);
    }
  } else {
    //Rcout << "FUNCTION cleanup3_inv DID NOT FIND START OPTION                  ";
    
    // Construct default list of start vectors (1 indiv / stage)
    DataFrame chosen_stageframe;
    
    if (format_int == 3 || format_int == 5) {
      chosen_stageframe = final_stageframe;
    } else if (format_int < 3) {
      chosen_stageframe = final_hstages;
    } else if (format_int == 4) {
      chosen_stageframe = final_agestages;
    }
    
    stagecounts = static_cast<int>(chosen_stageframe.nrows());
    if (format_int == 3 && funcbased) stagecounts--;
    
    arma::vec start_vec (stagecounts, fill::ones);
    
    for (int i = 0; i < variant_count; i++) {
      start_list(i) = start_vec;
    }
    start_count = variant_count;
  }
  
  //Rcout << "cleanup3_inv J ";
  
  // patches vector
  if (patch.isNotNull()) {
    if (is<CharacterVector>(patch) || is<LogicalVector>(patch)) {
      CharacterVector patch_vec_pre = as<CharacterVector>(patch);
      
      if (preexisting) {
        CharacterVector mpm_patch_vec (1);
        CharacterVector label_portion (2);
        
        DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
        final_labels = mpm_labels;
        CharacterVector mpm_labels_pop = as<CharacterVector>(mpm_labels["pop"]);
        CharacterVector mpm_labels_patch = as<CharacterVector>(mpm_labels["patch"]);
        
        if (!CharacterVector::is_na(patch_vec_pre(0))) {
          IntegerVector found_indices = index_l3(mpm_labels_patch, patch_vec_pre(0));
          if (found_indices.length() == 0) {
            throw Rcpp::exception("Some values in vector patches do not exist in entered MPMs.", 
              false);
          }
          mpm_patch_vec(0) = patch_vec_pre(0);
          
          int key_labels_index = found_indices(0);
          label_portion(0) = mpm_labels_pop(key_labels_index);
          label_portion(1) = mpm_labels_patch(key_labels_index);
        } else {
          int chosen_patch_index {0};
          for (int j = 0; j < static_cast<int>(mpm_labels_patch.length()); j++) {
            if (mpm_labels_patch(j) == "0" || CharacterVector::is_na(mpm_labels_patch(j))) {
              chosen_patch_index = j;
              break;
            }
          }
          mpm_patch_vec(0) = mpm_labels_patch(chosen_patch_index);
          
          label_portion(0) = mpm_labels_pop(chosen_patch_index);
          label_portion(1) = mpm_labels_patch(chosen_patch_index);
        }
        labels_list = label_portion;
        patch_vec = mpm_patch_vec;
        
      } else if (funcbased) {
        CharacterVector vrm_patch_vec (1);
        
        CharacterVector label_portion (2);
        label_portion(0) = "pop1";
        
        List chosen_vrm = vrm_list;
        DataFrame vrm_patchframe = as<DataFrame>(chosen_vrm["patch_frame"]);
        CharacterVector vrm_patchframe_patches = as<CharacterVector>(vrm_patchframe["patches"]);
        
        if (!CharacterVector::is_na(patch_vec_pre(0))) {
          IntegerVector found_indices = index_l3(vrm_patchframe_patches, patch_vec_pre(0));
          if (found_indices.length() == 0) {
            throw Rcpp::exception("Some values in vector patches do not exist in entered MPMs.", 
              false);
          }
          vrm_patch_vec(0) = patch_vec_pre(0);
          label_portion(1) = patch_vec_pre(0);
        } else { 
          vrm_patch_vec(0) = vrm_patchframe_patches(0);
          label_portion(1) = vrm_patchframe_patches(0);
        }
        labels_list = label_portion;
        patch_vec = vrm_patch_vec;
      }
    }
  } else {
    List labels_list_pre (1);
    
    if (preexisting) {
      CharacterVector patch_vec_pre (1);
      CharacterVector label_portion (2);
      
      DataFrame mpm_labels = as<DataFrame>(chosen_mpm["labels"]);
      final_labels = mpm_labels;
      CharacterVector mpm_labels_pop = as<CharacterVector>(mpm_labels["pop"]);
      CharacterVector mpm_labels_patch = as<CharacterVector>(mpm_labels["patch"]);
      
      int chosen_patch_index {0};
      for (int j = 0; j < static_cast<int>(mpm_labels_patch.length()); j++) {
        if (mpm_labels_patch(j) == "0" || CharacterVector::is_na(mpm_labels_patch(j))) {
          chosen_patch_index = j;
          break;
        }
      }
      
      patch_vec_pre(0) = mpm_labels_patch(chosen_patch_index);
      label_portion(0) = mpm_labels_pop(chosen_patch_index);
      label_portion(1) = mpm_labels_patch(chosen_patch_index);
      labels_list = label_portion;
      patch_vec = patch_vec_pre;
    } else if (funcbased) {
      CharacterVector patch_vec_pre (1);
      CharacterVector label_portion (2);
      label_portion(0) = "pop1";
      
      List chosen_vrm = vrm_list;
      DataFrame vrm_patchframe = as<DataFrame>(chosen_vrm["patch_frame"]);
      CharacterVector vrm_patchframe_patches = as<CharacterVector>(vrm_patchframe["patches"]);
      
      patch_vec_pre(0) = vrm_patchframe_patches(0);
      label_portion(1) = vrm_patchframe_patches(0);
      labels_list = label_portion;
      patch_vec = patch_vec_pre;
    }
  }
  
  //Rcout << "cleanup3_inv K ";
  
  // label construction
  {
    CharacterVector labels_pops (1);
    CharacterVector labels_patches (1);
    IntegerVector labels_mpms = {1};
    
    labels_pops(0) = labels_list(0);
    labels_patches(0) = labels_list(1);
    
    labels = DataFrame::create(_["mpm"] = labels_mpms, _["pop"] = labels_pops,
      _["patch"] = labels_patches);
  }
  
  //Rcout << "cleanup3_inv L ";
  
  // years vector
  if (years.isNotNull()) {
    if (is<NumericVector>(years) || is<CharacterVector>(years)) {
      CharacterVector year_vec_pre = as<CharacterVector>(years);
      
      if (preexisting) {
        CharacterVector mpm_labels_vars = final_labels.attr("names");
        IntegerVector mpm_labels_y2_var = index_l3(mpm_labels_vars, "year2");
        
        if (mpm_labels_y2_var.length() == 0) {
          throw Rcpp::exception("Input MPM appears to be mean MPM.", false);
        } else {
          CharacterVector mpm_labels_year2_full = as<CharacterVector>(final_labels["year2"]);
          CharacterVector mpm_labels_year2 = sort_unique(mpm_labels_year2_full);
          CharacterVector mly2_diffs = setdiff(year_vec, mpm_labels_year2);
          existing_years = mpm_labels_year2;
          
          total_years_int = static_cast<int>(mpm_labels_year2.length());
          
          if (mly2_diffs.length() > 0) { // Might be redundant
            throw Rcpp::exception("Some entered values in years do not exist in MPM.",
              false);
          }
          
          if (!CharacterVector::is_na(year_vec_pre(0))) {
            CharacterVector mly2_diffs = setdiff(year_vec_pre, mpm_labels_year2);
            
            year_vec = year_vec_pre;
            if (mly2_diffs.length() > 0) { 
              throw Rcpp::exception("Some entered values in years do not exist in some MPMs.",
                false);
            }
          } else {
            CharacterVector years_unique = sort_unique(mpm_labels_year2);
            year_vec = years_unique;
          }
        }
      } else if (funcbased) {
        DataFrame vrm_yearframe = as<DataFrame>(vrm_list["year_frame"]);
        
        CharacterVector vrm_yearframe_years = as<CharacterVector>(vrm_yearframe["years"]);
        if (!CharacterVector::is_na(year_vec_pre(0))) {
          CharacterVector vyy2_diffs = setdiff(year_vec_pre, vrm_yearframe_years);
          year_vec = year_vec_pre;
        
          if (vyy2_diffs.length() > 0) { 
            throw Rcpp::exception("Some entered values in years do not exist in MPM.",
              false);
          }
        } else {
          year_vec = vrm_yearframe_years;
        }
        
        CharacterVector unique_years = sort_unique(vrm_yearframe_years);
        existing_years = unique_years;
        total_years_int = static_cast<int>(unique_years.length());
      }
    } else {
      throw Rcpp::exception("Argument years is not valid.", false);
    }
  } else {
    if (preexisting) {
      CharacterVector mpm_labels_vars = final_labels.attr("names");
      IntegerVector mpm_labels_y2_var = index_l3(mpm_labels_vars, "year2");
      
      if (mpm_labels_y2_var.length() == 0) {
        total_years_int = 1;
      } else {
        CharacterVector mpm_labels_year2 = as<CharacterVector>(final_labels["year2"]);
        
        CharacterVector mly2_unique = sort_unique(mpm_labels_year2);
        int found_total_years = static_cast<int>(mly2_unique.length());
        total_years_int = found_total_years;
        
        year_vec = mly2_unique;
        existing_years = mly2_unique;
      }
    } else if (funcbased) {
      DataFrame vrm_yearframe = as<DataFrame>(vrm_list["year_frame"]);
      
      CharacterVector vrm_yearframe_years = as<CharacterVector>(vrm_yearframe["years"]);
      CharacterVector unique_years = sort_unique(vrm_yearframe_years);
      existing_years = unique_years; //vrm_yearframe_years;
      year_vec = unique_years; //vrm_yearframe_years;
      
      total_years_int = static_cast<int>(unique_years.length());
    }
  }
  
  //Rcout << "cleanup3_inv M ";
  
  // tweights list
  if (tweights.isNotNull()) {
    if (Rf_isMatrix(tweights)) {
      NumericMatrix chosen_matrix = as<NumericMatrix>(tweights);
      int mat_rows = chosen_matrix.nrow();
      int mat_cols = chosen_matrix.ncol();
      
      if (mat_rows != mat_cols) {
        throw Rcpp::exception("Matrices in argument tweights must be square.", false);
      }
      
      if (mat_rows != total_years_int) {
        throw Rcpp::exception("Matrices in argument tweights must account for all years.",
          false);
      }
      tweights_type_int = 2;
      tweights_list(0) = chosen_matrix;
      
    } else if (is<NumericVector>(tweights)) {
      NumericVector chosen_vector = as<NumericVector>(tweights);
      
      if (static_cast<int>(chosen_vector.length()) != total_years_int) {
        throw Rcpp::exception("Vectors in argument tweights must account for all years.",
          false);
      }
      tweights_type_int = 1;
      tweights_list(0) = chosen_vector;
          
    } else {
      AdaptUtils::pop_error2("tweights", "a numeric vector or matrix", "", 1);
    }
  } else {
    tweights_list(0) = R_NilValue;
    tweights_count = 1;
  }
  
  //Rcout << "cleanup3_inv N ";
  
  // density list
  if (density.isNotNull()) {
    if (is<DataFrame>(density)) {
      density_count = 1;
      DataFrame chosen_density = as<DataFrame>(density);
      
      if (chosen_density.hasAttribute("class")) {
        CharacterVector chosen_density_class = chosen_density.attr("class");
        
        for (int i = 0; i < static_cast<int>(chosen_density_class.length()); i++) {
          if (chosen_density_class(i) == "lefkoDens") dens_yn_bool = true;
        }
        if (!dens_yn_bool) AdaptUtils::pop_error2("density", "a lefkoDens object", "", 1);
        
        CharacterVector dl_stage1 = as<CharacterVector>(chosen_density["stage1"]);
        IntegerVector dl_age2 = as<IntegerVector>(chosen_density["age2"]);
        
        if (format_int < 3) {
          String eat_my_shorts = "Argument density requires real stage1 ";
          eat_my_shorts += "entries other than NA if MPMs are historical.";
              
          if (is<LogicalVector>(chosen_density["stage1"])) {
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          for (int j = 0; j < static_cast<int>(dl_stage1.length()); j++) {
            if (CharacterVector::is_na(dl_stage1(j))) {
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          }
        } else if (format_int > 3) {
          String eat_my_shorts = "Argument density requires real age2 ";
          eat_my_shorts += "entries other than NA if MPMs are age-by-stage.";
              
          if (is<LogicalVector>(chosen_density["age2"])) {
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          for (int j = 0; j < static_cast<int>(dl_age2.length()); j++) {
            if (IntegerVector::is_na(dl_age2(j)) || LogicalVector::is_na(dl_age2(j))) {
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          }
        }
      } else {
        AdaptUtils::pop_error2("density", "a lefkoDens object", "", 1);
      }
      
      Rcpp::StringVector di_stage3 = as<StringVector>(chosen_density["stage3"]);
      Rcpp::StringVector di_stage2 = as<StringVector>(chosen_density["stage2"]);
      Rcpp::StringVector di_stage1 = as<StringVector>(chosen_density["stage1"]);
      int di_size = di_stage3.length();
      
      if (format_int < 3) {
        StringVector stage3 = as<StringVector>(final_hstages["stage_2"]);
        StringVector stage2r = as<StringVector>(final_hstages["stage_1"]);
        StringVector stage2c = as<StringVector>(final_hstages["stage_2"]);
        StringVector stage1 = as<StringVector>(final_hstages["stage_1"]);
        int hst_size = stage3.length();
        
        arma::uvec hst_3(hst_size, fill::zeros);
        arma::uvec hst_2r(hst_size, fill::zeros);
        arma::uvec hst_2c(hst_size, fill::zeros);
        arma::uvec hst_1(hst_size, fill::zeros);
        
        arma::uvec di_stage32_id(di_size, fill::zeros);
        arma::uvec di_stage21_id(di_size, fill::zeros);
        arma::uvec di_index(di_size, fill::zeros);
        
        for (int j = 0; j < di_size; j++) { // Loop through each density_input line
          for (int k = 0; k < hst_size; k++) {
            if (di_stage3(j) == stage3(k)) {
              hst_3(k) = 1;
            } else {
              hst_3(k) = 0;
            }
          }
          
          for (int k = 0; k < hst_size; k++) {
            if (di_stage2(j) == stage2r(k)) {
              hst_2r(k) = 1;
            } else {
              hst_2r(k) = 0;
            }
          }
          
          for (int k = 0; k < hst_size; k++) {
            if (di_stage2(j) == stage2c(k)) {
              hst_2c(k) = 1;
            } else {
              hst_2c(k) = 0;
            }
          }
          
          for (int k = 0; k < hst_size; k++) {
            if (di_stage1(j) == stage1(k)) {
              hst_1(k) = 1;
            } else {
              hst_1(k) = 0;
            }
          }
          
          arma::uvec find_hst3 = find(hst_3);
          arma::uvec find_hst2r = find(hst_2r);
          arma::uvec find_hst2c = find(hst_2c);
          arma::uvec find_hst1 = find(hst_1);
          
          arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
          arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
          
          if (static_cast<int>(pop_32.n_elem) == 0 ||
              static_cast<int>(pop_21.n_elem) == 0) {
            throw Rcpp::exception("Some stages in argument density could not be found.", 
              false);
          }
          di_stage32_id(j) = pop_32(0);
          di_stage21_id(j) = pop_21(0);
          di_index(j) = pop_32(0) + (pop_21(0) * hst_size);
          
          hst_3.zeros();
          hst_2r.zeros();
          hst_2c.zeros();
          hst_1.zeros();
        }
        
        DataFrame dens_index_df_mpm = DataFrame::create(_["index32"] = di_stage32_id,
          _["index21"] = di_stage21_id, _["index321"] = di_index);
        dens_index_df = dens_index_df_mpm;
        density_df = chosen_density;
      } else if (format_int == 4 ) {
        IntegerVector di_age2 = as<IntegerVector>(chosen_density["age2"]);
        
        StringVector stage3 = as<StringVector>(final_agestages["stage"]);
        StringVector stage2 = as<StringVector>(final_agestages["stage"]);
        IntegerVector age2 = as<IntegerVector>(final_agestages["age"]);
        int agst_size = stage3.length();
        
        arma::uvec agst_s3(agst_size, fill::zeros);
        arma::uvec agst_a3(agst_size, fill::zeros);
        arma::uvec agst_s2(agst_size, fill::zeros);
        arma::uvec agst_a2(agst_size, fill::zeros);
        
        arma::uvec di_s3a3_id(di_size, fill::zeros);
        arma::uvec di_s2a2_id(di_size, fill::zeros);
        arma::uvec di_index(di_size, fill::zeros);
        
        for (int j = 0; j < di_size; j++) { // Loop through each density_input line
          for (int k = 0; k < agst_size; k++) {
            if (di_stage3(j) == stage3(k)) {
              agst_s3(k) = 1;
            } else {
              agst_s3(k) = 0;
            }
          }
          
          for (int k = 0; k < agst_size; k++) {
            if (di_stage2(j) == stage2(k)) {
              agst_s2(k) = 1;
            } else {
              agst_s2(k) = 0;
            }
          }
          
          for (int k = 0; k < agst_size; k++) {
            if (di_age2(j) < finalage_int) {
              if (di_age2(j) == age2(k)) {
                agst_a2(k) = 1;
                
                for (int l = 0; l < agst_size; l++) {
                  if ((di_age2(j) + 1) == age2(l)) {
                    agst_a3(l) = 1;
                  } else {
                    agst_a3(l) = 0;
                  }
                }
              } else {
                agst_a2(k) = 0;
              }
            } else {
              if (di_age2(j) == age2(k)) {
                agst_a2(k) = 1;
                agst_a3(k) = 1;
              } else {
                agst_a2(k) = 0;
                agst_a3(k) = 0;
              }
            }
          }
          
          arma::uvec find_agst_s3 = find(agst_s3);
          arma::uvec find_agst_s2 = find(agst_s2);
          arma::uvec find_agst_a3 = find(agst_a3);
          arma::uvec find_agst_a2 = find(agst_a2);
          
          arma::uvec pop_32 = intersect(find_agst_s3, find_agst_a3);
          arma::uvec pop_21 = intersect(find_agst_s2, find_agst_a2);
          
          if (static_cast<int>(pop_32.n_elem) == 0 || static_cast<int>(pop_21.n_elem) == 0) {
            throw Rcpp::exception("Some age-stages in argument density could not be found.", 
              false);
          }
          di_s3a3_id(j) = pop_32(0);
          di_s2a2_id(j) = pop_21(0);
          di_index(j) = pop_32(0) + (pop_21(0) * agst_size);
          
          agst_s3.zeros();
          agst_s2.zeros();
          agst_a3.zeros();
          agst_a2.zeros();
        }
        
        DataFrame dens_index_df_mpm = DataFrame::create(_["index32"] = di_s3a3_id,
          _["index21"] = di_s2a2_id, _["index321"] = di_index);
        dens_index_df = dens_index_df_mpm;
        density_df = chosen_density;
      } else {
        StringVector stage3 = as<StringVector>(final_stageframe["stage"]);
        StringVector stage2 = as<StringVector>(final_stageframe["stage"]);
        int ahst_size = stage3.length();
        if (funcbased) ahst_size--;
        
        arma::uvec ahst_3(ahst_size, fill::zeros);
        arma::uvec ahst_2(ahst_size, fill::zeros);
        
        arma::uvec di_stage32_id(di_size, fill::zeros);
        arma::uvec di_stage21_id(di_size, fill::zeros);
        arma::uvec di_index(di_size, fill::zeros);
        
        for (int j = 0; j < di_size; j++) { // Loop through each density_input
          for (int k = 0; k < ahst_size; k++) {
            if (di_stage3(j) == stage3(k)) {
              ahst_3(k) = 1;
            } else {
              ahst_3(k) = 0;
            }
          }
          
          for (int k = 0; k < ahst_size; k++) {
            if (di_stage2(j) == stage2(k)) {
              ahst_2(k) = 1;
            } else {
              ahst_2(k) = 0;
            }
          }
          
          arma::uvec find_ahst3 = find(ahst_3);
          arma::uvec find_ahst2 = find(ahst_2);
          di_stage32_id(j) = find_ahst3(0);
          di_stage21_id(j) = find_ahst2(0);
          di_index(j) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
          
          ahst_3.zeros();
          ahst_2.zeros();
        }
        
        DataFrame dens_index_df_mpm = DataFrame::create(_["index3"] = di_stage32_id,
          _["index2"] = di_stage21_id, _["index321"] = di_index);
        dens_index_df = dens_index_df_mpm;
        density_df = chosen_density;
      }
      
      arma::uvec dyn_style = as<arma::uvec>(chosen_density["style"]);
      arma::vec dyn_alpha = as<arma::vec>(chosen_density["alpha"]);
      arma::vec dyn_beta = as<arma::vec>(chosen_density["beta"]);
      
      for (int j = 0; j < static_cast<int>(dyn_style.n_elem); j++) {
        if (dyn_style(j) < 1 || dyn_style(j) > 4) {
          String eat_my_shorts = "Some density inputs are stated as yielding density ";
          eat_my_shorts += "dependence but not in an accepted style.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
        if (dyn_style(j) == 1) {
          if (dyn_beta(j) > exp_tol) {
            Rf_warningcall(R_NilValue,
              "Beta used in Ricker function may be too high. Results may be unpredictable.");
          } else if (dyn_beta(j) < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue,
              "Beta used in Ricker function may be too high. Results may be unpredictable.");
          }
          
        } else if (dyn_style(j) == 3) {
          double summed_stuff = dyn_alpha(j) + dyn_beta(j);
          
          if (summed_stuff > exp_tol) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high.");
            
          } else if (summed_stuff < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high.");
          }
        }
      }
      
    } else if (is<List>(density)) {
      List density_list = as<List>(density);
      density_count = static_cast<int>(density_list.length());
      if (density_count != 1) {
        throw Rcpp::exception("Please enter argument density as a single lefkoDens object.",
          false);
      }
      
      if (is<DataFrame>(density_list(0))) {
        DataFrame chosen_density = as<DataFrame>(density_list(0));
        
        if (chosen_density.hasAttribute("class")) {
          CharacterVector chosen_density_class = chosen_density.attr("class");
          bool found_lefkoDens {false};
          
          for (int j = 0; j < static_cast<int>(chosen_density_class.length()); j++) {
            if (chosen_density_class(j) == "lefkoDens") found_lefkoDens = true;
          }
          if (!found_lefkoDens) {
            AdaptUtils::pop_error2("density", "a lefkoDens object", "", 1);
          }
          
          CharacterVector dl_stage1 = as<CharacterVector>(chosen_density["stage1"]);
          IntegerVector dl_age2 = as<IntegerVector>(chosen_density["age2"]);
          
          if (format_int < 3) {
            String eat_my_shorts = "Argument density requires real stage1 ";
            eat_my_shorts += "entries other than NA if MPMs are historical.";
            
            if (is<LogicalVector>(chosen_density["stage1"])) {
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
            for (int j = 0; j < static_cast<int>(dl_stage1.length()); j++) {
              if (CharacterVector::is_na(dl_stage1(j))) {
                throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
              }
            }
          } else if (format_int > 3) {
            String eat_my_shorts = "Argument density requires real age2 ";
            eat_my_shorts += "entries other than NA if MPMs are age-by-stage.";
            
            if (is<LogicalVector>(chosen_density["age2"])) {
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
            for (int j = 0; j < static_cast<int>(dl_age2.length()); j++) {
              if (IntegerVector::is_na(dl_age2(j)) || LogicalVector::is_na(dl_age2(j))) {
                throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
              }
            }
          }
          
          dens_yn_bool = true;
        } else {
          AdaptUtils::pop_error2("density", "a lefkoDens object", "", 1);
        }
        
        Rcpp::StringVector di_stage3 = as<StringVector>(chosen_density["stage3"]);
        Rcpp::StringVector di_stage2 = as<StringVector>(chosen_density["stage2"]);
        Rcpp::StringVector di_stage1 = as<StringVector>(chosen_density["stage1"]);
        int di_size = di_stage3.length();
        
        if (format_int < 3) {
          StringVector stage3 = as<StringVector>(final_hstages["stage_2"]);
          StringVector stage2r = as<StringVector>(final_hstages["stage_1"]);
          StringVector stage2c = as<StringVector>(final_hstages["stage_2"]);
          StringVector stage1 = as<StringVector>(final_hstages["stage_1"]);
          int hst_size = stage3.length();
          
          arma::uvec hst_3(hst_size, fill::zeros);
          arma::uvec hst_2r(hst_size, fill::zeros);
          arma::uvec hst_2c(hst_size, fill::zeros);
          arma::uvec hst_1(hst_size, fill::zeros);
          
          arma::uvec di_stage32_id(di_size, fill::zeros);
          arma::uvec di_stage21_id(di_size, fill::zeros);
          arma::uvec di_index(di_size, fill::zeros);
          
          for (int j = 0; j < di_size; j++) { // Loop through each density_input line
            for (int k = 0; k < hst_size; k++) {
              if (di_stage3(j) == stage3(k)) {
                hst_3(k) = 1;
              } else {
                hst_3(k) = 0;
              }
            }
            
            for (int k = 0; k < hst_size; k++) {
              if (di_stage2(j) == stage2r(k)) {
                hst_2r(k) = 1;
              } else {
                hst_2r(k) = 0;
              }
            }
            
            for (int k = 0; k < hst_size; k++) {
              if (di_stage2(j) == stage2c(k)) {
                hst_2c(k) = 1;
              } else {
                hst_2c(k) = 0;
              }
            }
            
            for (int k = 0; k < hst_size; k++) {
              if (di_stage1(j) == stage1(k)) {
                hst_1(k) = 1;
              } else {
                hst_1(k) = 0;
              }
            }
            
            arma::uvec find_hst3 = find(hst_3);
            arma::uvec find_hst2r = find(hst_2r);
            arma::uvec find_hst2c = find(hst_2c);
            arma::uvec find_hst1 = find(hst_1);
            
            arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
            arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
            
            if (static_cast<int>(pop_32.n_elem) == 0 ||
                static_cast<int>(pop_21.n_elem) == 0) {
              throw Rcpp::exception("Some stages in argument density could not be found.", 
                false);
            }
            di_stage32_id(j) = pop_32(0);
            di_stage21_id(j) = pop_21(0);
            di_index(j) = pop_32(0) + (pop_21(0) * hst_size);
            
            hst_3.zeros();
            hst_2r.zeros();
            hst_2c.zeros();
            hst_1.zeros();
          }
          
          DataFrame dens_index_df_mpm = DataFrame::create(_["index32"] = di_stage32_id,
            _["index21"] = di_stage21_id, _["index321"] = di_index);
          dens_index_df = dens_index_df_mpm;
        } else if (format_int == 4 ) { 
          IntegerVector di_age2 = as<IntegerVector>(chosen_density["age2"]);
          
          StringVector stage3 = as<StringVector>(final_agestages["stage"]);
          StringVector stage2 = as<StringVector>(final_agestages["stage"]);
          IntegerVector age2 = as<IntegerVector>(final_agestages["age"]);
          int agst_size = stage3.length();
          
          arma::uvec agst_s3(agst_size, fill::zeros);
          arma::uvec agst_a3(agst_size, fill::zeros);
          arma::uvec agst_s2(agst_size, fill::zeros);
          arma::uvec agst_a2(agst_size, fill::zeros);
          
          arma::uvec di_s3a3_id(di_size, fill::zeros);
          arma::uvec di_s2a2_id(di_size, fill::zeros);
          arma::uvec di_index(di_size, fill::zeros);
          
          for (int j = 0; j < di_size; j++) { // Loop through each density_input line
            for (int k = 0; k < agst_size; k++) {
              if (di_stage3(j) == stage3(k)) {
                agst_s3(k) = 1;
              } else {
                agst_s3(k) = 0;
              }
            }
            
            for (int k = 0; k < agst_size; k++) {
              if (di_stage2(j) == stage2(k)) {
                agst_s2(k) = 1;
              } else {
                agst_s2(k) = 0;
              }
            }
            
            for (int k = 0; k < agst_size; k++) {
              if (di_age2(j) < finalage_int) {
                if (di_age2(j) == age2(k)) {
                  agst_a2(k) = 1;
                  
                  for (int l = 0; l < agst_size; l++) {
                    if ((di_age2(j) + 1) == age2(l)) {
                      agst_a3(l) = 1;
                    } else {
                      agst_a3(l) = 0;
                    }
                  }
                } else {
                  agst_a2(k) = 0;
                }
              } else {
                if (di_age2(j) == age2(k)) {
                  agst_a2(k) = 1;
                  agst_a3(k) = 1;
                } else {
                  agst_a2(k) = 0;
                  agst_a3(k) = 0;
                }
              }
            }
            
            arma::uvec find_agst_s3 = find(agst_s3);
            arma::uvec find_agst_s2 = find(agst_s2);
            arma::uvec find_agst_a3 = find(agst_a3);
            arma::uvec find_agst_a2 = find(agst_a2);
            
            arma::uvec pop_32 = intersect(find_agst_s3, find_agst_a3);
            arma::uvec pop_21 = intersect(find_agst_s2, find_agst_a2);
            
            if (static_cast<int>(pop_32.n_elem) == 0 || static_cast<int>(pop_21.n_elem) == 0) {
              throw Rcpp::exception("Some age-stages in argument density could not be found.", 
                false);
            }
            di_s3a3_id(j) = pop_32(0);
            di_s2a2_id(j) = pop_21(0);
            di_index(j) = pop_32(0) + (pop_21(0) * agst_size);
            
            agst_s3.zeros();
            agst_s2.zeros();
            agst_a3.zeros();
            agst_a2.zeros();
          }
          
          DataFrame dens_index_df_mpm = DataFrame::create(_["index32"] = di_s3a3_id,
            _["index21"] = di_s2a2_id, _["index321"] = di_index);
          dens_index_df = dens_index_df_mpm;
        } else {
          StringVector stage3 = as<StringVector>(final_stageframe["stage"]);
          StringVector stage2 = as<StringVector>(final_stageframe["stage"]);
          int ahst_size = stage3.length();
          if (funcbased) ahst_size--;
          
          arma::uvec ahst_3(ahst_size, fill::zeros);
          arma::uvec ahst_2(ahst_size, fill::zeros);
          
          arma::uvec di_stage32_id(di_size, fill::zeros);
          arma::uvec di_stage21_id(di_size, fill::zeros);
          arma::uvec di_index(di_size, fill::zeros);
          
          for (int j = 0; j < di_size; j++) { // Loop through each density_input
            for (int k = 0; k < ahst_size; k++) {
              if (di_stage3(j) == stage3(k)) {
                ahst_3(k) = 1;
              } else {
                ahst_3(k) = 0;
              }
            }
            
            for (int k = 0; k < ahst_size; k++) {
              if (di_stage2(j) == stage2(k)) {
                ahst_2(k) = 1;
              } else {
                ahst_2(k) = 0;
              }
            }
            
            arma::uvec find_ahst3 = find(ahst_3);
            arma::uvec find_ahst2 = find(ahst_2);
            di_stage32_id(j) = find_ahst3(0);
            di_stage21_id(j) = find_ahst2(0);
            di_index(j) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
            
            ahst_3.zeros();
            ahst_2.zeros();
          }
          DataFrame dens_index_df_mpm = DataFrame::create(_["index3"] = di_stage32_id,
            _["index2"] = di_stage21_id, _["index321"] = di_index);
          dens_index_df = dens_index_df_mpm;
        }
        
        density_df = chosen_density;
        
        arma::uvec dyn_style = as<arma::uvec>(chosen_density["style"]);
        arma::vec dyn_alpha = as<arma::vec>(chosen_density["alpha"]);
        arma::vec dyn_beta = as<arma::vec>(chosen_density["beta"]);
        
        for (int j = 0; j < static_cast<int>(dyn_style.n_elem); j++) {
          if (dyn_style(j) < 1 || dyn_style(j) > 4) {
            String eat_my_shorts = "Some density inputs are stated as yielding ";
            eat_my_shorts += "density dependence but not in an accepted style";
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          
          if (dyn_style(j) == 1) {
            if (dyn_beta(j) > exp_tol) {
              Rf_warningcall(R_NilValue,
                "Beta used in Ricker function may be too high. Results may be unpredictable");
            } else if (dyn_beta(j) < (-1.0 * exp_tol)) {
              Rf_warningcall(R_NilValue,
                "Beta used in Ricker function may be too high. Results may be unpredictable");
            }
            
          } else if (dyn_style(j) == 3) {
            double summed_stuff = dyn_alpha(j) + dyn_beta(j);
            
            if (summed_stuff > exp_tol) {
              Rf_warningcall(R_NilValue,
                "Alpha and beta used in Usher function may be too high.");
              
            } else if (summed_stuff < (-1.0 * exp_tol)) {
              Rf_warningcall(R_NilValue,
                "Alpha and beta used in Usher function may be too high.");
            }
          }
        }
      }
    } else { 
      AdaptUtils::pop_error2("density", "a lefkoDens object", "", 1);
    }
  }
  
  //Rcout << "cleanup3_inv O ";
  
  // entry time vector
  if (entry_time.isNotNull()) {
    if (is<NumericVector>(entry_time) || is<IntegerVector>(entry_time)) {
      IntegerVector entry_time_vec_pre = as<IntegerVector>(entry_time);
      entry_time_count = static_cast<int>(entry_time_vec.length());
      
      if (entry_time_count != 1 && entry_time_count != var_per_run) {
        AdaptUtils::pop_error2("entry_time", "a single integer for every variant to run concurrently", "", 1);
      }
      
      for (int i = 0; i < entry_time_count; i++) {
        if (IntegerVector::is_na(entry_time_vec_pre(i))) {
          AdaptUtils::pop_error2("NA values", "entry_time", "", 25);
        }
      }
      
      for (int i = 0; i < var_per_run; i++) {
        if (entry_time_count == 1) {
          entry_time_vec(i) = entry_time_vec_pre(0);
        } else {
          entry_time_vec(i) = entry_time_vec_pre(i);
        }
      }
      
      int entry_times_sum = sum(entry_time_vec);
      if (entry_times_sum > 0) entry_time_vec_use = true;
      
    } else {
      AdaptUtils::pop_error2("entry_time", "a single integer for every variant to run concurrently", "", 1);
    }
  }
  
  //Rcout << "cleanup3_inv P ";
  
  // vrms-only arguments
  if (!funcbased) {
    if (density_vr.isNotNull()) {
      AdaptUtils::pop_error2("vrm", "use argument density_vr", "", 26);
    }
    if (sp_density.isNotNull()) {
      AdaptUtils::pop_error2("vrm", "use argument sp_density", "", 26);
    }
    if (ind_terms.isNotNull()) {
      AdaptUtils::pop_error2("vrm", "use argument ind_terms", "", 26);
    }
    if (dev_terms.isNotNull()) {
      AdaptUtils::pop_error2("vrm", "use argument dev_terms", "", 26);
    }
    if (fb_sparse.isNotNull()) {
      AdaptUtils::pop_error2("vrm", "use argument fb_sparse", "", 26);
    }
  } else {
    // density_vr list
    IntegerVector dvr_yn_count (1);
    
    if (density_vr.isNotNull()) {
      if (is<DataFrame>(density_vr)) {
        chosen_density_vr = as<DataFrame>(density_vr);
        
        if (chosen_density_vr.hasAttribute("class")) {
          CharacterVector chosen_density_vr_class = chosen_density_vr.attr("class");
          
          for (int j = 0; j < static_cast<int>(chosen_density_vr_class.length()); j++) {
            if (chosen_density_vr_class(j) == "lefkoDensVR") dens_vr_yn_bool = true;
          }
          if (!dens_vr_yn_bool) {
            AdaptUtils::pop_error2("density_vr", "a lefkoDensVR object", "", 1);
          }
        } else {
          AdaptUtils::pop_error2("density_vr", "a lefkoDensVR object", "", 1);
        }
      } else { 
        AdaptUtils::pop_error2("density_vr", "a lefkoDensVR object", "", 1);
      }
    }
    
    // sp_density list
    if (sp_density.isNotNull()) {
      if (is<NumericVector>(sp_density)) {
        NumericVector sp_density_temp = as<NumericVector>(sp_density);
        int sp_density_count = static_cast<int>(sp_density_temp.length());
        
        if (sp_density_count != 1) {
          AdaptUtils::pop_error2("sp_density", "a single numeric value", "", 1);
        }
        
        List sp_density_initial_list (1);
        NumericVector single_value {static_cast<double>(sp_density_temp(0))};
        sp_density_initial_list(0) = single_value;
        sp_density_num_int = 1;
      } else {
        throw Rcpp::exception("Input in argument sp_density is not valid.",
          false);
      }
    } else {
      List sp_density_list_temp (1);
      
      NumericVector sp_temp {0.0};
      sp_density_list_temp(0) = sp_temp;
    }
    
    // ind_terms list
    if (ind_terms.isNotNull()) {
      if (is<DataFrame>(ind_terms)) {
        DataFrame ind_terms_df = as<DataFrame>(ind_terms);
        int idt_df_size = static_cast<int>(ind_terms_df.size());
        if (idt_df_size != 3) {
          throw Rcpp::exception("Data frame ind_terms should have exactly 3 columns.",
            false);
        }
        
        if ((!is<NumericVector>(ind_terms_df(0)) && !is<CharacterVector>(ind_terms_df(0))) ||
          (!is<NumericVector>(ind_terms_df(1)) && !is<CharacterVector>(ind_terms_df(1))) ||
          (!is<NumericVector>(ind_terms_df(2)) && !is<CharacterVector>(ind_terms_df(2)))) {
            AdaptUtils::pop_error2("ind_terms", "a data frame composed of numeric or character values", "", 1);
        }
        
        List idt_num_pre (1);
        List idt_cat_pre (1);
        
        int idt_df_nrows = static_cast<int>(ind_terms_df.nrows());
        List current_idt_cat (3);
        List current_idt_num (3);
        
        for (int j = 0; j < 3; j++) {
          if (is<CharacterVector>(ind_terms_df(j))) {
            CharacterVector current_idt_cat_col = as<CharacterVector>(ind_terms_df(j));
            NumericVector current_idt_num_col (idt_df_nrows);
            
            current_idt_cat(j) = current_idt_cat_col;
            current_idt_num(j) = current_idt_num_col;
            
            if (j == 0) {
              inda_terms_num_int = 0;
              inda_terms_cat_int = idt_df_nrows;
            } else if (j == 1) {
              indb_terms_num_int = 0;
              indb_terms_cat_int = idt_df_nrows;
            } else {
              indc_terms_num_int = 0;
              indc_terms_cat_int = idt_df_nrows;
            }
            
          } else if (is<NumericVector>(ind_terms_df(j))) {
            CharacterVector single_none {"none"};
            CharacterVector current_idt_cat_col = rep(single_none, idt_df_nrows);
            NumericVector current_idt_num_col = as<NumericVector>(ind_terms_df(j));
            
            current_idt_cat(j) = current_idt_cat_col;
            current_idt_num(j) = current_idt_num_col;
            
            if (j == 0) {
              inda_terms_num_int = idt_df_nrows;
              inda_terms_cat_int = 0;
            } else if (j == 1) {
              indb_terms_num_int = idt_df_nrows;
              indb_terms_cat_int = 0;
            } else {
              indc_terms_num_int = idt_df_nrows;
              indc_terms_cat_int = 0;
            }
          }
        }
        
        idt_num_pre(0) = current_idt_num;
        idt_cat_pre(0) = current_idt_cat;
        
        ind_terms_num_list = idt_num_pre;
        ind_terms_cat_list = idt_cat_pre;
        
      } else {
        throw Rcpp::exception("Input in argument ind_terms is not valid.",
          false);
      }
    } else {
      List ind_terms_num_list_pre (1);
      List ind_terms_cat_list_pre (1);
      
      NumericVector region_A = NumericVector::create(0.);
      CharacterVector region_B {"none"};
      
      DataFrame region_A_df = DataFrame::create(_["A"] = region_A,
        _["B"] = clone(region_A), _["C"] = clone(region_A));
      DataFrame region_B_df = DataFrame::create(_["A"] = region_B,
        _["B"] = region_B, _["C"] = region_B);
      
      ind_terms_num_list_pre(0) = region_A_df;
      ind_terms_cat_list_pre(0) = region_B_df;
      
      ind_terms_num_list = ind_terms_num_list_pre;
      ind_terms_cat_list = ind_terms_cat_list_pre;
    }
    
    // dev_terms list
    if (dev_terms.isNotNull()) {
      if (is<DataFrame>(dev_terms)) {
        DataFrame current_dev_terms = as<DataFrame>(dev_terms);
        int dvtc_df_size = static_cast<int>(current_dev_terms.size());
        int dvtc_df_nrows = static_cast<int>(current_dev_terms.nrows());
        
        if (dvtc_df_size != 14) {
          throw Rcpp::exception("Data frame in argument dev_terms must have 14 columns.",
            false);
        }
        
        for (int j = 0; j < 14; j++) {
          if (!is<NumericVector>(current_dev_terms(j))) {
            AdaptUtils::pop_error2("dev_terms", "a data frame composed of numeric values", "", 1);
          }
        }
        
        dev_terms_num_int = dvtc_df_nrows;
        
        NumericMatrix dev_matrix (14, dvtc_df_nrows); // Rows = devs / vital rate, Cols = times
        
        for (int i = 0; i < 14; i++) {
          NumericVector test_var = as<NumericVector>(current_dev_terms[i]);
          
          for (int j = 0; j < dvtc_df_nrows; j++) {
            if (!NumericVector::is_na(test_var(j))) dev_matrix(i, j) = test_var(j);
          }
        }
        
        List current_dv_dataframe_list (variant_count);
        for (int i = 0; i < variant_count; i++) {
          current_dv_dataframe_list(i) = clone(dev_matrix);
        }
        dev_terms_list = current_dv_dataframe_list;
        
      } else if (Rf_isMatrix(dev_terms)) {
        NumericMatrix input_dev_matrix = as<NumericMatrix>(dev_terms);
        int input_dv_rows = static_cast<int>(input_dev_matrix.nrow());
        int input_dv_cols = static_cast<int>(input_dev_matrix.ncol());
        
        if (input_dv_cols != 14) {
          throw Rcpp::exception("Numeric matrix in argument dev_terms must have 14 columns",
            false);
        }
        
        dev_terms_num_int = input_dv_rows;
        
        NumericMatrix dev_matrix (14, input_dv_rows); // Rows are the devs per vital rate, Cols are times
        
        for (int i = 0; i < 14; i++) {
          NumericVector test_col = input_dev_matrix(_, i);
          
          for (int j = 0; j < input_dv_rows; j++) {
            if (!NumericVector::is_na(test_col(j))) dev_matrix(i, j) = test_col(j);
          }
        }
        
        List current_dv_matrix_list (variant_count);
        for (int i = 0; i < variant_count; i++) {
          current_dv_matrix_list(i) = clone(dev_matrix);
        }
        dev_terms_list = current_dv_matrix_list;
        
      } else if (is<NumericVector>(dev_terms)) {
        NumericVector dev_terms_asnumeric = as<NumericVector>(dev_terms);
        
        int input_dv_cols = static_cast<int>(dev_terms_asnumeric.length());
        
        if (input_dv_cols != 14) {
          throw Rcpp::exception("Numeric vector in argument dev_terms must have 14 columns",
            false);
        }
        
        dev_terms_num_int = 1;
        
        NumericMatrix dev_matrix (14, 1); // Rows are the devs per vital rate, Cols are times
        
        for (int i = 0; i < 14; i++) {
          if (!NumericVector::is_na(dev_terms_asnumeric(i))) dev_matrix(i, 0) = dev_terms_asnumeric(i);
        }
        
        List current_dv_matrix_list (variant_count);
        for (int i = 0; i < variant_count; i++) {
          current_dv_matrix_list(i) = clone(dev_matrix);
        }
        dev_terms_list = current_dv_matrix_list;
        
      } else {
        AdaptUtils::pop_error2("dev_term", "a data frame or a numeric vector with 14 elements", "", 1);
      }
    }
    
    // fb_sparse
    if (fb_sparse.isNotNull()) {
      if (is<LogicalVector>(fb_sparse)) {
        LogicalVector sparse_vec = as<LogicalVector>(fb_sparse);
        sparse_vec_count = static_cast<int>(sparse_vec.length());
        
        if (sparse_vec_count != 1) {
          throw Rcpp::exception("Argument fb_sparse must be a single logical value.",
            false);
        }
        
        if (LogicalVector::is_na(sparse_vec(0))) {
          throw Rcpp::exception("No NA values are allowed in argument fb_sparse.",
            false);
        }
        
        if (sparse_vec(0) > 0) sparse_bool = true;
      }
    }
  } // End of vrm-only section
  
  //Rcout << "cleanup3_inv Q ";
  
  // equivalence interpretation
  if (equivalence.isNotNull()) {
    if (is<DataFrame>(equivalence)) {
      stages_not_equal = true;
      
      //int trial_count = 1;
      equivalence_count = 1;
      
      List equivalence_list_pre (equivalence_count);

      equivalence_frame = as<DataFrame>(equivalence);
      if (!equivalence_frame.hasAttribute("class")) {
        throw Rcpp::exception("Argument equivalence should be a data frame of class adaptEq.", 
          false);
      }
      CharacterVector eq_list_df_class = equivalence_frame.attr("class");
      bool found_adaptEq {false};
      for (int j = 0; j < static_cast<int>(eq_list_df_class.length()); j++) {
        if (eq_list_df_class(j) == "adaptEq") found_adaptEq = true;
      }
      if (!found_adaptEq) {
        throw Rcpp::exception("Argument equivalence should be a data frame of class adaptEq.",
          false);
      }
      
      IntegerVector eq_s2 = as<IntegerVector>(equivalence_frame["stage_id_2"]);
      IntegerVector eq_s1 = as<IntegerVector>(equivalence_frame["stage_id_1"]);
      IntegerVector eq_a2 = as<IntegerVector>(equivalence_frame["age2"]);
      IntegerVector eq_rn = clone(as<IntegerVector>(equivalence_frame["row_num"]));
      NumericVector eq_val = as<NumericVector>(equivalence_frame["value"]);
      
      eq_rn = eq_rn - 1;
      
      if (format_int < 3) {
        if (IntegerVector::is_na(eq_s1(0))) {
          throw Rcpp::exception("Enter stage pairs in adaptEq objects used for historical MPMs.", 
            false);
        }
        if (IntegerVector::is_na(eq_s2(0))) {
          throw Rcpp::exception("Entries in column stage2 of adaptEq objects cannot be empty except in Leslie MPMs.", false);
        }
      } else if (format_int > 3) {
        if (IntegerVector::is_na(eq_a2(0))) {
          throw Rcpp::exception("Enter ages in adaptEq objects used for age-by-stage MPMs.",
            false);
        }
        if (format_int == 4) {
          if (IntegerVector::is_na(eq_s2(0))) {
            throw Rcpp::exception("Entries in column stage2 of adaptEq objects cannot be empty except in Leslie MPMs.", false);
          }
        }
      } else {
        if (IntegerVector::is_na(eq_s2(0))) {
          throw Rcpp::exception("Entries in column stage2 of adaptEq objects cannot be empty except in Leslie MPMs.", false);
        }
      }
      
      if (max(eq_rn) > matrowcounts) {
        throw Rcpp::exception("Some row numbers in an entered adaptEq object are too high.", 
          false);
      }
      
      if (min(eq_val) < 0.0) {
        AdaptUtils::pop_error2("equivalence", "", "", 30);
      }
      
      NumericVector current_eq (matrowcounts, 1.0);
      for (int j = 0; j < static_cast<int>(eq_rn.length()); j++) {
        current_eq(eq_rn(j)) = eq_val(j);
      }
      
      equivalence_vec = current_eq;
    } else {
      throw Rcpp::exception("Argument equivalence should be a data frame of class adaptEq.", 
        false);
    }
  } else {
    equivalence_count = 1;
    NumericVector equivalance_vec_pre (equivalence_count, 1.0);
    equivalence_vec = equivalance_vec_pre;
  }
  
  //Rcout << "cleanup3_inv R ";
  
  // process stageframe, supplement, repmatrix, and allstages list for fbMPMs
  if (funcbased) {
    // Create function-based MPMs and assign them to mpm_list
    List current_vrm = vrm_list;
    int ehrlen_format {1}; // This will need to be dealt with differently later
    
    int mpm_style {1};
    int filter_style {1};
    if (format_int < 3) {
      mpm_style = 0;
      if (format_int == 2) ehrlen_format = 2;
    } else if (format_int == 4) {
      mpm_style = 2;
      filter_style = 2;
    }
    
    if (format_int < 5) {
      current_mpm_allstages = theoldpizzle(final_stageframe,
        supplement_df, final_repmatrix, firstage_int, finalage_int,
        ehrlen_format, mpm_style, cont_int, filter_style);
    } else {
      current_mpm_allstages = final_stageframe;
    }
    
    // vrm_input processing
    // Move model summaries to appropriate RObjects
    RObject current_surv_model;
    RObject current_obs_model;
    RObject current_size_model;
    RObject current_sizeb_model;
    RObject current_sizec_model;
    RObject current_repst_model;
    RObject current_fec_model;
    RObject current_jsurv_model;
    RObject current_jobs_model;
    RObject current_jsize_model;
    RObject current_jsizeb_model;
    RObject current_jsizec_model;
    RObject current_jrepst_model;
    RObject current_jmatst_model;
    DataFrame current_paramnames;
    
    DataFrame vrm_frame = as<DataFrame>(current_vrm["vrm_frame"]);
    DataFrame year_frame = as<DataFrame>(current_vrm["year_frame"]);
    DataFrame patch_frame = as<DataFrame>(current_vrm["patch_frame"]);
    DataFrame group2_frame = as<DataFrame>(current_vrm["group2_frame"]);
    DataFrame group1_frame = as<DataFrame>(current_vrm["group1_frame"]);
    DataFrame dist_frame = as<DataFrame>(current_vrm["dist_frame"]);
    NumericVector st_frame = as<NumericVector>(current_vrm["st_frame"]);
    
    CharacterVector main_effect_1 = as<CharacterVector>(vrm_frame["main_effect_1"]);
    CharacterVector effects_names = clone(main_effect_1);
    
    CharacterVector main_effect_2;
    if (main_effect_1.length() > 20) {
      main_effect_2 = as<CharacterVector>(vrm_frame["main_effect_2"]);
      
      for (int i = 0; i < main_effect_1.length(); i++) {
        if (i > 16) {
          effects_names(i) += ":";
          effects_names(i) += main_effect_2(i);
        }
      }
    }
      
    CharacterVector year_names = as<CharacterVector>(year_frame["years"]);
    CharacterVector patch_names = as<CharacterVector>(patch_frame["patches"]);
    CharacterVector group_names = as<CharacterVector>(group2_frame["groups"]);
    
    bool zi_yn = false;
    int vrm_length = vrm_frame.length();
    
    NumericVector surv_num = as<NumericVector>(vrm_frame["surv"]);
    NumericVector obs_num = as<NumericVector>(vrm_frame["obs"]);
    NumericVector sizea_num = as<NumericVector>(vrm_frame["sizea"]);
    NumericVector sizeb_num = as<NumericVector>(vrm_frame["sizeb"]);
    NumericVector sizec_num = as<NumericVector>(vrm_frame["sizec"]);
    NumericVector repst_num = as<NumericVector>(vrm_frame["repst"]);
    NumericVector fec_num = as<NumericVector>(vrm_frame["fec"]);
    NumericVector jsurv_num = as<NumericVector>(vrm_frame["jsurv"]);
    NumericVector jobs_num = as<NumericVector>(vrm_frame["jobs"]);
    NumericVector jsizea_num = as<NumericVector>(vrm_frame["jsizea"]);
    NumericVector jsizeb_num = as<NumericVector>(vrm_frame["jsizeb"]);
    NumericVector jsizec_num = as<NumericVector>(vrm_frame["jsizec"]);
    NumericVector jrepst_num = as<NumericVector>(vrm_frame["jrepst"]);
    NumericVector jmatst_num = as<NumericVector>(vrm_frame["jmatst"]);
    
    NumericVector surv_year = as<NumericVector>(year_frame["surv"]);
    NumericVector obs_year = as<NumericVector>(year_frame["obs"]);
    NumericVector sizea_year = as<NumericVector>(year_frame["sizea"]);
    NumericVector sizeb_year = as<NumericVector>(year_frame["sizeb"]);
    NumericVector sizec_year = as<NumericVector>(year_frame["sizec"]);
    NumericVector repst_year = as<NumericVector>(year_frame["repst"]);
    NumericVector fec_year = as<NumericVector>(year_frame["fec"]);
    NumericVector jsurv_year = as<NumericVector>(year_frame["jsurv"]);
    NumericVector jobs_year = as<NumericVector>(year_frame["jobs"]);
    NumericVector jsizea_year = as<NumericVector>(year_frame["jsizea"]);
    NumericVector jsizeb_year = as<NumericVector>(year_frame["jsizeb"]);
    NumericVector jsizec_year = as<NumericVector>(year_frame["jsizec"]);
    NumericVector jrepst_year = as<NumericVector>(year_frame["jrepst"]);
    NumericVector jmatst_year = as<NumericVector>(year_frame["jmatst"]);
    
    NumericVector surv_patch = as<NumericVector>(patch_frame["surv"]);
    NumericVector obs_patch = as<NumericVector>(patch_frame["obs"]);
    NumericVector sizea_patch = as<NumericVector>(patch_frame["sizea"]);
    NumericVector sizeb_patch = as<NumericVector>(patch_frame["sizeb"]);
    NumericVector sizec_patch = as<NumericVector>(patch_frame["sizec"]);
    NumericVector repst_patch = as<NumericVector>(patch_frame["repst"]);
    NumericVector fec_patch = as<NumericVector>(patch_frame["fec"]);
    NumericVector jsurv_patch = as<NumericVector>(patch_frame["jsurv"]);
    NumericVector jobs_patch = as<NumericVector>(patch_frame["jobs"]);
    NumericVector jsizea_patch = as<NumericVector>(patch_frame["jsizea"]);
    NumericVector jsizeb_patch = as<NumericVector>(patch_frame["jsizeb"]);
    NumericVector jsizec_patch = as<NumericVector>(patch_frame["jsizec"]);
    NumericVector jrepst_patch = as<NumericVector>(patch_frame["jrepst"]);
    NumericVector jmatst_patch = as<NumericVector>(patch_frame["jmatst"]);
    
    NumericVector surv_group2 = as<NumericVector>(group2_frame["surv"]);
    NumericVector obs_group2 = as<NumericVector>(group2_frame["obs"]);
    NumericVector sizea_group2 = as<NumericVector>(group2_frame["sizea"]);
    NumericVector sizeb_group2 = as<NumericVector>(group2_frame["sizeb"]);
    NumericVector sizec_group2 = as<NumericVector>(group2_frame["sizec"]);
    NumericVector repst_group2 = as<NumericVector>(group2_frame["repst"]);
    NumericVector fec_group2 = as<NumericVector>(group2_frame["fec"]);
    NumericVector jsurv_group2 = as<NumericVector>(group2_frame["jsurv"]);
    NumericVector jobs_group2 = as<NumericVector>(group2_frame["jobs"]);
    NumericVector jsizea_group2 = as<NumericVector>(group2_frame["jsizea"]);
    NumericVector jsizeb_group2 = as<NumericVector>(group2_frame["jsizeb"]);
    NumericVector jsizec_group2 = as<NumericVector>(group2_frame["jsizec"]);
    NumericVector jrepst_group2 = as<NumericVector>(group2_frame["jrepst"]);
    NumericVector jmatst_group2 = as<NumericVector>(group2_frame["jmatst"]);
    
    NumericVector surv_group1 = as<NumericVector>(group1_frame["surv"]);
    NumericVector obs_group1 = as<NumericVector>(group1_frame["obs"]);
    NumericVector sizea_group1 = as<NumericVector>(group1_frame["sizea"]);
    NumericVector sizeb_group1 = as<NumericVector>(group1_frame["sizeb"]);
    NumericVector sizec_group1 = as<NumericVector>(group1_frame["sizec"]);
    NumericVector repst_group1 = as<NumericVector>(group1_frame["repst"]);
    NumericVector fec_group1 = as<NumericVector>(group1_frame["fec"]);
    NumericVector jsurv_group1 = as<NumericVector>(group1_frame["jsurv"]);
    NumericVector jobs_group1 = as<NumericVector>(group1_frame["jobs"]);
    NumericVector jsizea_group1 = as<NumericVector>(group1_frame["jsizea"]);
    NumericVector jsizeb_group1 = as<NumericVector>(group1_frame["jsizeb"]);
    NumericVector jsizec_group1 = as<NumericVector>(group1_frame["jsizec"]);
    NumericVector jrepst_group1 = as<NumericVector>(group1_frame["jrepst"]);
    NumericVector jmatst_group1 = as<NumericVector>(group1_frame["jmatst"]);
    
    StringVector distribs = as<StringVector>(dist_frame["dist"]);
    String surv_dist = distribs(0);
    String obs_dist = distribs(1);
    String sizea_dist = distribs(2);
    String sizeb_dist = distribs(3);
    String sizec_dist = distribs(4);
    String repst_dist = distribs(5);
    String fec_dist = distribs(6);
    String jsurv_dist = distribs(7);
    String jobs_dist = distribs(8);
    String jsizea_dist = distribs(9);
    String jsizeb_dist = distribs(10);
    String jsizec_dist = distribs(11);
    String jrepst_dist = distribs(12);
    String jmatst_dist = distribs(13);
    
    double sizea_st = st_frame(2);
    double sizeb_st = st_frame(3);
    double sizec_st = st_frame(4);
    double fec_st = st_frame(6);
    double jsizea_st = st_frame(9);
    double jsizeb_st = st_frame(10);
    double jsizec_st = st_frame(11);
    
    NumericVector sizea_zi;
    NumericVector sizeb_zi;
    NumericVector sizec_zi;
    NumericVector fec_zi;
    NumericVector jsizea_zi;
    NumericVector jsizeb_zi;
    NumericVector jsizec_zi;
    
    NumericVector year_sizea_zi;
    NumericVector year_sizeb_zi;
    NumericVector year_sizec_zi;
    NumericVector year_fec_zi;
    NumericVector year_jsizea_zi;
    NumericVector year_jsizeb_zi;
    NumericVector year_jsizec_zi;
    
    NumericVector patch_sizea_zi;
    NumericVector patch_sizeb_zi;
    NumericVector patch_sizec_zi;
    NumericVector patch_fec_zi;
    NumericVector patch_jsizea_zi;
    NumericVector patch_jsizeb_zi;
    NumericVector patch_jsizec_zi;
    
    NumericVector group2_sizea_zi;
    NumericVector group2_sizeb_zi;
    NumericVector group2_sizec_zi;
    NumericVector group2_fec_zi;
    NumericVector group2_jsizea_zi;
    NumericVector group2_jsizeb_zi;
    NumericVector group2_jsizec_zi;
    
    NumericVector group1_sizea_zi;
    NumericVector group1_sizeb_zi;
    NumericVector group1_sizec_zi;
    NumericVector group1_fec_zi;
    NumericVector group1_jsizea_zi;
    NumericVector group1_jsizeb_zi;
    NumericVector group1_jsizec_zi;
    
    NumericVector dud_zi;
    
    if (vrm_length > 16) {
      zi_yn = true;
      
      sizea_zi = as<NumericVector>(vrm_frame["sizea_zi"]);
      sizeb_zi = as<NumericVector>(vrm_frame["sizeb_zi"]);
      sizec_zi = as<NumericVector>(vrm_frame["sizec_zi"]);
      fec_zi = as<NumericVector>(vrm_frame["fec_zi"]);
      jsizea_zi = as<NumericVector>(vrm_frame["jsizea_zi"]);
      jsizeb_zi = as<NumericVector>(vrm_frame["jsizeb_zi"]);
      jsizec_zi = as<NumericVector>(vrm_frame["jsizec_zi"]);
      
      year_sizea_zi = as<NumericVector>(year_frame["sizea_zi"]);
      year_sizeb_zi = as<NumericVector>(year_frame["sizeb_zi"]);
      year_sizec_zi = as<NumericVector>(year_frame["sizec_zi"]);
      year_fec_zi = as<NumericVector>(year_frame["fec_zi"]);
      year_jsizea_zi = as<NumericVector>(year_frame["jsizea_zi"]);
      year_jsizeb_zi = as<NumericVector>(year_frame["jsizeb_zi"]);
      year_jsizec_zi = as<NumericVector>(year_frame["jsizec_zi"]);
      
      patch_sizea_zi = as<NumericVector>(patch_frame["sizea_zi"]);
      patch_sizeb_zi = as<NumericVector>(patch_frame["sizeb_zi"]);
      patch_sizec_zi = as<NumericVector>(patch_frame["sizec_zi"]);
      patch_fec_zi = as<NumericVector>(patch_frame["fec_zi"]);
      patch_jsizea_zi = as<NumericVector>(patch_frame["jsizea_zi"]);
      patch_jsizeb_zi = as<NumericVector>(patch_frame["jsizeb_zi"]);
      patch_jsizec_zi = as<NumericVector>(patch_frame["jsizec_zi"]);
      
      group2_sizea_zi = as<NumericVector>(group2_frame["sizea_zi"]);
      group2_sizeb_zi = as<NumericVector>(group2_frame["sizeb_zi"]);
      group2_sizec_zi = as<NumericVector>(group2_frame["sizec_zi"]);
      group2_fec_zi = as<NumericVector>(group2_frame["fec_zi"]);
      group2_jsizea_zi = as<NumericVector>(group2_frame["jsizea_zi"]);
      group2_jsizeb_zi = as<NumericVector>(group2_frame["jsizeb_zi"]);
      group2_jsizec_zi = as<NumericVector>(group2_frame["jsizec_zi"]);
      
      group1_sizea_zi = as<NumericVector>(group1_frame["sizea_zi"]);
      group1_sizeb_zi = as<NumericVector>(group1_frame["sizeb_zi"]);
      group1_sizec_zi = as<NumericVector>(group1_frame["sizec_zi"]);
      group1_fec_zi = as<NumericVector>(group1_frame["fec_zi"]);
      group1_jsizea_zi = as<NumericVector>(group1_frame["jsizea_zi"]);
      group1_jsizeb_zi = as<NumericVector>(group1_frame["jsizeb_zi"]);
      group1_jsizec_zi = as<NumericVector>(group1_frame["jsizec_zi"]);
    }
    
    CharacterVector indcova_names;
    CharacterVector indcovb_names;
    CharacterVector indcovc_names;
    
    NumericVector surv_indcova2;
    NumericVector surv_indcovb2;
    NumericVector surv_indcovc2;
    NumericVector obs_indcova2;
    NumericVector obs_indcovb2;
    NumericVector obs_indcovc2;
    NumericVector sizea_indcova2;
    NumericVector sizea_indcovb2;
    NumericVector sizea_indcovc2;
    NumericVector sizeb_indcova2;
    NumericVector sizeb_indcovb2;
    NumericVector sizeb_indcovc2;
    NumericVector sizec_indcova2;
    NumericVector sizec_indcovb2;
    NumericVector sizec_indcovc2;
    NumericVector repst_indcova2;
    NumericVector repst_indcovb2;
    NumericVector repst_indcovc2;
    NumericVector fec_indcova2;
    NumericVector fec_indcovb2;
    NumericVector fec_indcovc2;
    NumericVector jsurv_indcova2;
    NumericVector jsurv_indcovb2;
    NumericVector jsurv_indcovc2;
    NumericVector jobs_indcova2;
    NumericVector jobs_indcovb2;
    NumericVector jobs_indcovc2;
    NumericVector jsizea_indcova2;
    NumericVector jsizea_indcovb2;
    NumericVector jsizea_indcovc2;
    NumericVector jsizeb_indcova2;
    NumericVector jsizeb_indcovb2;
    NumericVector jsizeb_indcovc2;
    NumericVector jsizec_indcova2;
    NumericVector jsizec_indcovb2;
    NumericVector jsizec_indcovc2;
    NumericVector jrepst_indcova2;
    NumericVector jrepst_indcovb2;
    NumericVector jrepst_indcovc2;
    NumericVector jmatst_indcova2;
    NumericVector jmatst_indcovb2;
    NumericVector jmatst_indcovc2;
    
    NumericVector sizea_indcova2_zi;
    NumericVector sizea_indcovb2_zi;
    NumericVector sizea_indcovc2_zi;
    NumericVector sizeb_indcova2_zi;
    NumericVector sizeb_indcovb2_zi;
    NumericVector sizeb_indcovc2_zi;
    NumericVector sizec_indcova2_zi;
    NumericVector sizec_indcovb2_zi;
    NumericVector sizec_indcovc2_zi;
    NumericVector fec_indcova2_zi;
    NumericVector fec_indcovb2_zi;
    NumericVector fec_indcovc2_zi;
    NumericVector jsizea_indcova2_zi;
    NumericVector jsizea_indcovb2_zi;
    NumericVector jsizea_indcovc2_zi;
    NumericVector jsizeb_indcova2_zi;
    NumericVector jsizeb_indcovb2_zi;
    NumericVector jsizeb_indcovc2_zi;
    NumericVector jsizec_indcova2_zi;
    NumericVector jsizec_indcovb2_zi;
    NumericVector jsizec_indcovc2_zi;
    
    NumericVector surv_indcova1;
    NumericVector surv_indcovb1;
    NumericVector surv_indcovc1;
    NumericVector obs_indcova1;
    NumericVector obs_indcovb1;
    NumericVector obs_indcovc1;
    NumericVector sizea_indcova1;
    NumericVector sizea_indcovb1;
    NumericVector sizea_indcovc1;
    NumericVector sizeb_indcova1;
    NumericVector sizeb_indcovb1;
    NumericVector sizeb_indcovc1;
    NumericVector sizec_indcova1;
    NumericVector sizec_indcovb1;
    NumericVector sizec_indcovc1;
    NumericVector repst_indcova1;
    NumericVector repst_indcovb1;
    NumericVector repst_indcovc1;
    NumericVector fec_indcova1;
    NumericVector fec_indcovb1;
    NumericVector fec_indcovc1;
    NumericVector jsurv_indcova1;
    NumericVector jsurv_indcovb1;
    NumericVector jsurv_indcovc1;
    NumericVector jobs_indcova1;
    NumericVector jobs_indcovb1;
    NumericVector jobs_indcovc1;
    NumericVector jsizea_indcova1;
    NumericVector jsizea_indcovb1;
    NumericVector jsizea_indcovc1;
    NumericVector jsizeb_indcova1;
    NumericVector jsizeb_indcovb1;
    NumericVector jsizeb_indcovc1;
    NumericVector jsizec_indcova1;
    NumericVector jsizec_indcovb1;
    NumericVector jsizec_indcovc1;
    NumericVector jrepst_indcova1;
    NumericVector jrepst_indcovb1;
    NumericVector jrepst_indcovc1;
    NumericVector jmatst_indcova1;
    NumericVector jmatst_indcovb1;
    NumericVector jmatst_indcovc1;
    
    NumericVector sizea_indcova1_zi;
    NumericVector sizea_indcovb1_zi;
    NumericVector sizea_indcovc1_zi;
    NumericVector sizeb_indcova1_zi;
    NumericVector sizeb_indcovb1_zi;
    NumericVector sizeb_indcovc1_zi;
    NumericVector sizec_indcova1_zi;
    NumericVector sizec_indcovb1_zi;
    NumericVector sizec_indcovc1_zi;
    NumericVector fec_indcova1_zi;
    NumericVector fec_indcovb1_zi;
    NumericVector fec_indcovc1_zi;
    NumericVector jsizea_indcova1_zi;
    NumericVector jsizea_indcovb1_zi;
    NumericVector jsizea_indcovc1_zi;
    NumericVector jsizeb_indcova1_zi;
    NumericVector jsizeb_indcovb1_zi;
    NumericVector jsizeb_indcovc1_zi;
    NumericVector jsizec_indcova1_zi;
    NumericVector jsizec_indcovb1_zi;
    NumericVector jsizec_indcovc1_zi;
    
    int modelsuite_length = current_vrm.length();
    CharacterVector modelsuite_names = current_vrm.attr("names");
    
    for (int i = 0; i < modelsuite_length; i++) {
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova2_frame")) {
        DataFrame indcova2_frame = as<DataFrame>(current_vrm["indcova2_frame"]);
        
        indcova_names = indcova2_frame["indcova"];
        
        surv_indcova2 = indcova2_frame["surv"];
        obs_indcova2 = indcova2_frame["obs"];
        sizea_indcova2 = indcova2_frame["sizea"];
        sizeb_indcova2 = indcova2_frame["sizeb"];
        sizec_indcova2 = indcova2_frame["sizec"];
        repst_indcova2 = indcova2_frame["repst"];
        fec_indcova2 = indcova2_frame["fec"];
        
        jsurv_indcova2 = indcova2_frame["jsurv"];
        jobs_indcova2 = indcova2_frame["jobs"];
        jsizea_indcova2 = indcova2_frame["jsizea"];
        jsizeb_indcova2 = indcova2_frame["jsizeb"];
        jsizec_indcova2 = indcova2_frame["jsizec"];
        jrepst_indcova2 = indcova2_frame["jrepst"];
        jmatst_indcova2 = indcova2_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcova2_zi = indcova2_frame["sizea_zi"];
          sizeb_indcova2_zi = indcova2_frame["sizeb_zi"];
          sizec_indcova2_zi = indcova2_frame["sizec_zi"];
          fec_indcova2_zi = indcova2_frame["fec_zi"];
          jsizea_indcova2_zi = indcova2_frame["jsizea_zi"];
          jsizeb_indcova2_zi = indcova2_frame["jsizeb_zi"];
          jsizec_indcova2_zi = indcova2_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova1_frame")) {
        DataFrame indcova1_frame = as<DataFrame>(current_vrm["indcova1_frame"]);
        
        indcova_names = indcova1_frame["indcova"];
        
        surv_indcova1 = indcova1_frame["surv"];
        obs_indcova1 = indcova1_frame["obs"];
        sizea_indcova1 = indcova1_frame["sizea"];
        sizeb_indcova1 = indcova1_frame["sizeb"];
        sizec_indcova1 = indcova1_frame["sizec"];
        repst_indcova1 = indcova1_frame["repst"];
        fec_indcova1 = indcova1_frame["fec"];
        
        jsurv_indcova1 = indcova1_frame["jsurv"];
        jobs_indcova1 = indcova1_frame["jobs"];
        jsizea_indcova1 = indcova1_frame["jsizea"];
        jsizeb_indcova1 = indcova1_frame["jsizeb"];
        jsizec_indcova1 = indcova1_frame["jsizec"];
        jrepst_indcova1 = indcova1_frame["jrepst"];
        jmatst_indcova1 = indcova1_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcova1_zi = indcova1_frame["sizea_zi"];
          sizeb_indcova1_zi = indcova1_frame["sizeb_zi"];
          sizec_indcova1_zi = indcova1_frame["sizec_zi"];
          fec_indcova1_zi = indcova1_frame["fec_zi"];
          jsizea_indcova1_zi = indcova1_frame["jsizea_zi"];
          jsizeb_indcova1_zi = indcova1_frame["jsizeb_zi"];
          jsizec_indcova1_zi = indcova1_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb2_frame")) {
        DataFrame indcovb2_frame = as<DataFrame>(current_vrm["indcovb2_frame"]);
        
        indcovb_names = indcovb2_frame["indcovb"];
        
        surv_indcovb2 = indcovb2_frame["surv"];
        obs_indcovb2 = indcovb2_frame["obs"];
        sizea_indcovb2 = indcovb2_frame["sizea"];
        sizeb_indcovb2 = indcovb2_frame["sizeb"];
        sizec_indcovb2 = indcovb2_frame["sizec"];
        repst_indcovb2 = indcovb2_frame["repst"];
        fec_indcovb2 = indcovb2_frame["fec"];
        
        jsurv_indcovb2 = indcovb2_frame["jsurv"];
        jobs_indcovb2 = indcovb2_frame["jobs"];
        jsizea_indcovb2 = indcovb2_frame["jsizea"];
        jsizeb_indcovb2 = indcovb2_frame["jsizeb"];
        jsizec_indcovb2 = indcovb2_frame["jsizec"];
        jrepst_indcovb2 = indcovb2_frame["jrepst"];
        jmatst_indcovb2 = indcovb2_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovb2_zi = indcovb2_frame["sizea_zi"];
          sizeb_indcovb2_zi = indcovb2_frame["sizeb_zi"];
          sizec_indcovb2_zi = indcovb2_frame["sizec_zi"];
          fec_indcovb2_zi = indcovb2_frame["fec_zi"];
          jsizea_indcovb2_zi = indcovb2_frame["jsizea_zi"];
          jsizeb_indcovb2_zi = indcovb2_frame["jsizeb_zi"];
          jsizec_indcovb2_zi = indcovb2_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb1_frame")) {
        DataFrame indcovb1_frame = as<DataFrame>(current_vrm["indcovb1_frame"]);
        
        indcovb_names = indcovb1_frame["indcovb"];
        
        surv_indcovb1 = indcovb1_frame["surv"];
        obs_indcovb1 = indcovb1_frame["obs"];
        sizea_indcovb1 = indcovb1_frame["sizea"];
        sizeb_indcovb1 = indcovb1_frame["sizeb"];
        sizec_indcovb1 = indcovb1_frame["sizec"];
        repst_indcovb1 = indcovb1_frame["repst"];
        fec_indcovb1 = indcovb1_frame["fec"];
        
        jsurv_indcovb1 = indcovb1_frame["jsurv"];
        jobs_indcovb1 = indcovb1_frame["jobs"];
        jsizea_indcovb1 = indcovb1_frame["jsizea"];
        jsizeb_indcovb1 = indcovb1_frame["jsizeb"];
        jsizec_indcovb1 = indcovb1_frame["jsizec"];
        jrepst_indcovb1 = indcovb1_frame["jrepst"];
        jmatst_indcovb1 = indcovb1_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovb1_zi = indcovb1_frame["sizea_zi"];
          sizeb_indcovb1_zi = indcovb1_frame["sizeb_zi"];
          sizec_indcovb1_zi = indcovb1_frame["sizec_zi"];
          fec_indcovb1_zi = indcovb1_frame["fec_zi"];
          jsizea_indcovb1_zi = indcovb1_frame["jsizea_zi"];
          jsizeb_indcovb1_zi = indcovb1_frame["jsizeb_zi"];
          jsizec_indcovb1_zi = indcovb1_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc2_frame")) {
        DataFrame indcovc2_frame = as<DataFrame>(current_vrm["indcovc2_frame"]);
        
        indcovc_names = indcovc2_frame["indcovc"];
        
        surv_indcovc2 = indcovc2_frame["surv"];
        obs_indcovc2 = indcovc2_frame["obs"];
        sizea_indcovc2 = indcovc2_frame["sizea"];
        sizeb_indcovc2 = indcovc2_frame["sizeb"];
        sizec_indcovc2 = indcovc2_frame["sizec"];
        repst_indcovc2 = indcovc2_frame["repst"];
        fec_indcovc2 = indcovc2_frame["fec"];
        
        jsurv_indcovc2 = indcovc2_frame["jsurv"];
        jobs_indcovc2 = indcovc2_frame["jobs"];
        jsizea_indcovc2 = indcovc2_frame["jsizea"];
        jsizeb_indcovc2 = indcovc2_frame["jsizeb"];
        jsizec_indcovc2 = indcovc2_frame["jsizec"];
        jrepst_indcovc2 = indcovc2_frame["jrepst"];
        jmatst_indcovc2 = indcovc2_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovc2_zi = indcovc2_frame["sizea_zi"];
          sizeb_indcovc2_zi = indcovc2_frame["sizeb_zi"];
          sizec_indcovc2_zi = indcovc2_frame["sizec_zi"];
          fec_indcovc2_zi = indcovc2_frame["fec_zi"];
          jsizea_indcovc2_zi = indcovc2_frame["jsizea_zi"];
          jsizeb_indcovc2_zi = indcovc2_frame["jsizeb_zi"];
          jsizec_indcovc2_zi = indcovc2_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc1_frame")) {
        DataFrame indcovc1_frame = as<DataFrame>(current_vrm["indcovc1_frame"]);
        
        indcovc_names = indcovc1_frame["indcovc"];
        
        surv_indcovc1 = indcovc1_frame["surv"];
        obs_indcovc1 = indcovc1_frame["obs"];
        sizea_indcovc1 = indcovc1_frame["sizea"];
        sizeb_indcovc1 = indcovc1_frame["sizeb"];
        sizec_indcovc1 = indcovc1_frame["sizec"];
        repst_indcovc1 = indcovc1_frame["repst"];
        fec_indcovc1 = indcovc1_frame["fec"];
        
        jsurv_indcovc1 = indcovc1_frame["jsurv"];
        jobs_indcovc1 = indcovc1_frame["jobs"];
        jsizea_indcovc1 = indcovc1_frame["jsizea"];
        jsizeb_indcovc1 = indcovc1_frame["jsizeb"];
        jsizec_indcovc1 = indcovc1_frame["jsizec"];
        jrepst_indcovc1 = indcovc1_frame["jrepst"];
        jmatst_indcovc1 = indcovc1_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovc1_zi = indcovc1_frame["sizea_zi"];
          sizeb_indcovc1_zi = indcovc1_frame["sizeb_zi"];
          sizec_indcovc1_zi = indcovc1_frame["sizec_zi"];
          fec_indcovc1_zi = indcovc1_frame["fec_zi"];
          jsizea_indcovc1_zi = indcovc1_frame["jsizea_zi"];
          jsizeb_indcovc1_zi = indcovc1_frame["jsizeb_zi"];
          jsizec_indcovc1_zi = indcovc1_frame["jsizec_zi"];
        }
      }
    }
    
    CharacterVector list_names = {"fixed_slopes", "year_slopes", "patch_slopes",
      "group2_slopes", "dist", "zi", "fixed_zi", "year_zi", "patch_zi",
      "group2_zi", "indcova_names", "indcova2_slopes", "indcova2_zi",
      "indcovb_names", "indcovb2_slopes", "indcovb2_zi", "indcovc_names",
      "indcovc2_slopes", "indcovc2_zi", "year_names", "patch_names",
      "group_names", "main_effect_1", "main_effect_2", "sigma_theta",
      "effects_names", "group1_slopes", "group1_zi", "indcova1_slopes",
      "indcovb1_slopes", "indcovc1_slopes", "indcova1_zi", "indcovb1_zi",
      "indcovc1_zi"};
    
    List surv_list(34);
    surv_list(0) = surv_num;
    surv_list(1) = surv_year;
    surv_list(2) = surv_patch;
    surv_list(3) = surv_group2;
    surv_list(4) = surv_dist;
    surv_list(5) = false;
    surv_list(6) = dud_zi;
    surv_list(7) = dud_zi;
    surv_list(8) = dud_zi;
    surv_list(9) = dud_zi;
    surv_list(10) = indcova_names;
    surv_list(11) = surv_indcova2;
    surv_list(12) = dud_zi;
    surv_list(13) = indcovb_names;
    surv_list(14) = surv_indcovb2;
    surv_list(15) = dud_zi;
    surv_list(16) = indcovc_names;
    surv_list(17) = surv_indcovc2;
    surv_list(18) = dud_zi;
    surv_list(19) = year_names;
    surv_list(20) = patch_names;
    surv_list(21) = group_names;
    surv_list(22) = main_effect_1;
    surv_list(23) = main_effect_2;
    surv_list(24) = 1.0;
    surv_list(25) = effects_names;
    surv_list(26) = surv_group1;
    surv_list(27) = dud_zi;
    surv_list(28) = surv_indcova1;
    surv_list(29) = surv_indcovb1;
    surv_list(30) = surv_indcovc1;
    surv_list(31) = dud_zi;
    surv_list(32) = dud_zi;
    surv_list(33) = dud_zi;
    
    List obs_list(34);
    obs_list(0) = obs_num;
    obs_list(1) = obs_year;
    obs_list(2) = obs_patch;
    obs_list(3) = obs_group2;
    obs_list(4) = obs_dist;
    obs_list(5) = false;
    obs_list(6) = dud_zi;
    obs_list(7) = dud_zi;
    obs_list(8) = dud_zi;
    obs_list(9) = dud_zi;
    obs_list(10) = indcova_names;
    obs_list(11) = obs_indcova2;
    obs_list(12) = dud_zi;
    obs_list(13) = indcovb_names;
    obs_list(14) = obs_indcovb2;
    obs_list(15) = dud_zi;
    obs_list(16) = indcovc_names;
    obs_list(17) = obs_indcovc2;
    obs_list(18) = dud_zi;
    obs_list(19) = year_names;
    obs_list(20) = patch_names;
    obs_list(21) = group_names;
    obs_list(22) = main_effect_1;
    obs_list(23) = main_effect_2;
    obs_list(24) = 1.0;
    obs_list(25) = effects_names;
    obs_list(26) = obs_group1;
    obs_list(27) = dud_zi;
    obs_list(28) = obs_indcova1;
    obs_list(29) = obs_indcovb1;
    obs_list(30) = obs_indcovc1;
    obs_list(31) = dud_zi;
    obs_list(32) = dud_zi;
    obs_list(33) = dud_zi;
    
    List sizea_list(34);
    sizea_list(0) = sizea_num;
    sizea_list(1) = sizea_year;
    sizea_list(2) = sizea_patch;
    sizea_list(3) = sizea_group2;
    sizea_list(4) = sizea_dist;
    sizea_list(5) = zi_yn;
    sizea_list(6) = sizea_zi;
    sizea_list(7) = year_sizea_zi;
    sizea_list(8) = patch_sizea_zi;
    sizea_list(9) = group2_sizea_zi;
    sizea_list(10) = indcova_names;
    sizea_list(11) = sizea_indcova2;
    sizea_list(12) = sizea_indcova2_zi;
    sizea_list(13) = indcovb_names;
    sizea_list(14) = sizea_indcovb2;
    sizea_list(15) = sizea_indcovb2_zi;
    sizea_list(16) = indcovc_names;
    sizea_list(17) = sizea_indcovc2;
    sizea_list(18) = sizea_indcovc2_zi;
    sizea_list(19) = year_names;
    sizea_list(20) = patch_names;
    sizea_list(21) = group_names;
    sizea_list(22) = main_effect_1;
    sizea_list(23) = main_effect_2;
    sizea_list(24) = sizea_st;
    sizea_list(25) = effects_names;
    sizea_list(26) = sizea_group1;
    sizea_list(27) = group1_sizea_zi;
    sizea_list(28) = sizea_indcova1;
    sizea_list(29) = sizea_indcovb1;
    sizea_list(30) = sizea_indcovc1;
    sizea_list(31) = sizea_indcova1_zi;
    sizea_list(32) = sizea_indcovb1_zi;
    sizea_list(33) = sizea_indcovc1_zi;
    
    List sizeb_list(34);
    sizeb_list(0) = sizeb_num;
    sizeb_list(1) = sizeb_year;
    sizeb_list(2) = sizeb_patch;
    sizeb_list(3) = sizeb_group2;
    sizeb_list(4) = sizeb_dist;
    sizeb_list(5) = zi_yn;
    sizeb_list(6) = sizeb_zi;
    sizeb_list(7) = year_sizeb_zi;
    sizeb_list(8) = patch_sizeb_zi;
    sizeb_list(9) = group2_sizeb_zi;
    sizeb_list(10) = indcova_names;
    sizeb_list(11) = sizeb_indcova2;
    sizeb_list(12) = sizeb_indcova2_zi;
    sizeb_list(13) = indcovb_names;
    sizeb_list(14) = sizeb_indcovb2;
    sizeb_list(15) = sizeb_indcovb2_zi;
    sizeb_list(16) = indcovc_names;
    sizeb_list(17) = sizeb_indcovc2;
    sizeb_list(18) = sizeb_indcovc2_zi;
    sizeb_list(19) = year_names;
    sizeb_list(20) = patch_names;
    sizeb_list(21) = group_names;
    sizeb_list(22) = main_effect_1;
    sizeb_list(23) = main_effect_2;
    sizeb_list(24) = sizeb_st;
    sizeb_list(25) = effects_names;
    sizeb_list(26) = sizeb_group1;
    sizeb_list(27) = group1_sizeb_zi;
    sizeb_list(28) = sizeb_indcova1;
    sizeb_list(29) = sizeb_indcovb1;
    sizeb_list(30) = sizeb_indcovc1;
    sizeb_list(31) = sizeb_indcova1_zi;
    sizeb_list(32) = sizeb_indcovb1_zi;
    sizeb_list(33) = sizeb_indcovc1_zi;
    
    List sizec_list(34);
    sizec_list(0) = sizec_num;
    sizec_list(1) = sizec_year;
    sizec_list(2) = sizec_patch;
    sizec_list(3) = sizec_group2;
    sizec_list(4) = sizec_dist;
    sizec_list(5) = zi_yn;
    sizec_list(6) = sizec_zi;
    sizec_list(7) = year_sizec_zi;
    sizec_list(8) = patch_sizec_zi;
    sizec_list(9) = group2_sizec_zi;
    sizec_list(10) = indcova_names;
    sizec_list(11) = sizec_indcova2;
    sizec_list(12) = sizec_indcova2_zi;
    sizec_list(13) = indcovb_names;
    sizec_list(14) = sizec_indcovb2;
    sizec_list(15) = sizec_indcovb2_zi;
    sizec_list(16) = indcovc_names;
    sizec_list(17) = sizec_indcovc2;
    sizec_list(18) = sizec_indcovc2_zi;
    sizec_list(19) = year_names;
    sizec_list(20) = patch_names;
    sizec_list(21) = group_names;
    sizec_list(22) = main_effect_1;
    sizec_list(23) = main_effect_2;
    sizec_list(24) = sizec_st;
    sizec_list(25) = effects_names;
    sizec_list(26) = sizec_group1;
    sizec_list(27) = group1_sizec_zi;
    sizec_list(28) = sizec_indcova1;
    sizec_list(29) = sizec_indcovb1;
    sizec_list(30) = sizec_indcovc1;
    sizec_list(31) = sizec_indcova1_zi;
    sizec_list(32) = sizec_indcovb1_zi;
    sizec_list(33) = sizec_indcovc1_zi;
    
    List repst_list(34);
    repst_list(0) = repst_num;
    repst_list(1) = repst_year;
    repst_list(2) = repst_patch;
    repst_list(3) = repst_group2;
    repst_list(4) = repst_dist;
    repst_list(5) = false;
    repst_list(6) = dud_zi;
    repst_list(7) = dud_zi;
    repst_list(8) = dud_zi;
    repst_list(9) = dud_zi;
    repst_list(10) = indcova_names;
    repst_list(11) = repst_indcova2;
    repst_list(12) = dud_zi;
    repst_list(13) = indcovb_names;
    repst_list(14) = repst_indcovb2;
    repst_list(15) = dud_zi;
    repst_list(16) = indcovc_names;
    repst_list(17) = repst_indcovc2;
    repst_list(18) = dud_zi;
    repst_list(19) = year_names;
    repst_list(20) = patch_names;
    repst_list(21) = group_names;
    repst_list(22) = main_effect_1;
    repst_list(23) = main_effect_2;
    repst_list(24) = 1.0;
    repst_list(25) = effects_names;
    repst_list(26) = repst_group1;
    repst_list(27) = dud_zi;
    repst_list(28) = repst_indcova1;
    repst_list(29) = repst_indcovb1;
    repst_list(30) = repst_indcovc1;
    repst_list(31) = dud_zi;
    repst_list(32) = dud_zi;
    repst_list(33) = dud_zi;
    
    List fec_list(34);
    fec_list(0) = fec_num;
    fec_list(1) = fec_year;
    fec_list(2) = fec_patch;
    fec_list(3) = fec_group2;
    fec_list(4) = fec_dist;
    fec_list(5) = zi_yn;
    fec_list(6) = fec_zi;
    fec_list(7) = year_fec_zi;
    fec_list(8) = patch_fec_zi;
    fec_list(9) = group2_fec_zi;
    fec_list(10) = indcova_names;
    fec_list(11) = fec_indcova2;
    fec_list(12) = fec_indcova2_zi;
    fec_list(13) = indcovb_names;
    fec_list(14) = fec_indcovb2;
    fec_list(15) = fec_indcovb2_zi;
    fec_list(16) = indcovc_names;
    fec_list(17) = fec_indcovc2;
    fec_list(18) = fec_indcovc2_zi;
    fec_list(19) = year_names;
    fec_list(20) = patch_names;
    fec_list(21) = group_names;
    fec_list(22) = main_effect_1;
    fec_list(23) = main_effect_2;
    fec_list(24) = fec_st;
    fec_list(25) = effects_names;
    fec_list(26) = fec_group1;
    fec_list(27) = group1_fec_zi;
    fec_list(28) = fec_indcova1;
    fec_list(29) = fec_indcovb1;
    fec_list(30) = fec_indcovc1;
    fec_list(31) = fec_indcova1_zi;
    fec_list(32) = fec_indcovb1_zi;
    fec_list(33) = fec_indcovc1_zi;
    
    List jsurv_list(34);
    jsurv_list(0) = jsurv_num;
    jsurv_list(1) = jsurv_year;
    jsurv_list(2) = jsurv_patch;
    jsurv_list(3) = jsurv_group2;
    jsurv_list(4) = jsurv_dist;
    jsurv_list(5) = false;
    jsurv_list(6) = dud_zi;
    jsurv_list(7) = dud_zi;
    jsurv_list(8) = dud_zi;
    jsurv_list(9) = dud_zi;
    jsurv_list(10) = indcova_names;
    jsurv_list(11) = jsurv_indcova2;
    jsurv_list(12) = dud_zi;
    jsurv_list(13) = indcovb_names;
    jsurv_list(14) = jsurv_indcovb2;
    jsurv_list(15) = dud_zi;
    jsurv_list(16) = indcovc_names;
    jsurv_list(17) = jsurv_indcovc2;
    jsurv_list(18) = dud_zi;
    jsurv_list(19) = year_names;
    jsurv_list(20) = patch_names;
    jsurv_list(21) = group_names;
    jsurv_list(22) = main_effect_1;
    jsurv_list(23) = main_effect_2;
    jsurv_list(24) = 1.0;
    jsurv_list(25) = effects_names;
    jsurv_list(26) = jsurv_group1;
    jsurv_list(27) = dud_zi;
    jsurv_list(28) = jsurv_indcova1;
    jsurv_list(29) = jsurv_indcovb1;
    jsurv_list(30) = jsurv_indcovc1;
    jsurv_list(31) = dud_zi;
    jsurv_list(32) = dud_zi;
    jsurv_list(33) = dud_zi;
    
    List jobs_list(34);
    jobs_list(0) = jobs_num;
    jobs_list(1) = jobs_year;
    jobs_list(2) = jobs_patch;
    jobs_list(3) = jobs_group2;
    jobs_list(4) = jobs_dist;
    jobs_list(5) = false;
    jobs_list(6) = dud_zi;
    jobs_list(7) = dud_zi;
    jobs_list(8) = dud_zi;
    jobs_list(9) = dud_zi;
    jobs_list(10) = indcova_names;
    jobs_list(11) = jobs_indcova2;
    jobs_list(12) = dud_zi;
    jobs_list(13) = indcovb_names;
    jobs_list(14) = jobs_indcovb2;
    jobs_list(15) = dud_zi;
    jobs_list(16) = indcovc_names;
    jobs_list(17) = jobs_indcovc2;
    jobs_list(18) = dud_zi;
    jobs_list(19) = year_names;
    jobs_list(20) = patch_names;
    jobs_list(21) = group_names;
    jobs_list(22) = main_effect_1;
    jobs_list(23) = main_effect_2;
    jobs_list(24) = 1.0;
    jobs_list(25) = effects_names;
    jobs_list(26) = jobs_group1;
    jobs_list(27) = dud_zi;
    jobs_list(28) = jobs_indcova1;
    jobs_list(29) = jobs_indcovb1;
    jobs_list(30) = jobs_indcovc1;
    jobs_list(31) = dud_zi;
    jobs_list(32) = dud_zi;
    jobs_list(33) = dud_zi;
    
    List jsizea_list(34);
    jsizea_list(0) = jsizea_num;
    jsizea_list(1) = jsizea_year;
    jsizea_list(2) = jsizea_patch;
    jsizea_list(3) = jsizea_group2;
    jsizea_list(4) = jsizea_dist;
    jsizea_list(5) = zi_yn;
    jsizea_list(6) = jsizea_zi;
    jsizea_list(7) = year_jsizea_zi;
    jsizea_list(8) = patch_jsizea_zi;
    jsizea_list(9) = group2_jsizea_zi;
    jsizea_list(10) = indcova_names;
    jsizea_list(11) = jsizea_indcova2;
    jsizea_list(12) = jsizea_indcova2_zi;
    jsizea_list(13) = indcovb_names;
    jsizea_list(14) = jsizea_indcovb2;
    jsizea_list(15) = jsizea_indcovb2_zi;
    jsizea_list(16) = indcovc_names;
    jsizea_list(17) = jsizea_indcovc2;
    jsizea_list(18) = jsizea_indcovc2_zi;
    jsizea_list(19) = year_names;
    jsizea_list(20) = patch_names;
    jsizea_list(21) = group_names;
    jsizea_list(22) = main_effect_1;
    jsizea_list(23) = main_effect_2;
    jsizea_list(24) = jsizea_st;
    jsizea_list(25) = effects_names;
    jsizea_list(26) = jsizea_group1;
    jsizea_list(27) = group1_jsizea_zi;
    jsizea_list(28) = jsizea_indcova1;
    jsizea_list(29) = jsizea_indcovb1;
    jsizea_list(30) = jsizea_indcovc1;
    jsizea_list(31) = jsizea_indcova1_zi;
    jsizea_list(32) = jsizea_indcovb1_zi;
    jsizea_list(33) = jsizea_indcovc1_zi;
    
    List jsizeb_list(34);
    jsizeb_list(0) = jsizeb_num;
    jsizeb_list(1) = jsizeb_year;
    jsizeb_list(2) = jsizeb_patch;
    jsizeb_list(3) = jsizeb_group2;
    jsizeb_list(4) = jsizeb_dist;
    jsizeb_list(5) = zi_yn;
    jsizeb_list(6) = jsizeb_zi;
    jsizeb_list(7) = year_jsizeb_zi;
    jsizeb_list(8) = patch_jsizeb_zi;
    jsizeb_list(9) = group2_jsizeb_zi;
    jsizeb_list(10) = indcova_names;
    jsizeb_list(11) = jsizeb_indcova2;
    jsizeb_list(12) = jsizeb_indcova2_zi;
    jsizeb_list(13) = indcovb_names;
    jsizeb_list(14) = jsizeb_indcovb2;
    jsizeb_list(15) = jsizeb_indcovb2_zi;
    jsizeb_list(16) = indcovc_names;
    jsizeb_list(17) = jsizeb_indcovc2;
    jsizeb_list(18) = jsizeb_indcovc2_zi;
    jsizeb_list(19) = year_names;
    jsizeb_list(20) = patch_names;
    jsizeb_list(21) = group_names;
    jsizeb_list(22) = main_effect_1;
    jsizeb_list(23) = main_effect_2;
    jsizeb_list(24) = jsizeb_st;
    jsizeb_list(25) = effects_names;
    jsizeb_list(26) = jsizeb_group1;
    jsizeb_list(27) = group1_jsizeb_zi;
    jsizeb_list(28) = jsizeb_indcova1;
    jsizeb_list(29) = jsizeb_indcovb1;
    jsizeb_list(30) = jsizeb_indcovc1;
    jsizeb_list(31) = jsizeb_indcova1_zi;
    jsizeb_list(32) = jsizeb_indcovb1_zi;
    jsizeb_list(33) = jsizeb_indcovc1_zi;
    
    List jsizec_list(34);
    jsizec_list(0) = jsizec_num;
    jsizec_list(1) = jsizec_year;
    jsizec_list(2) = jsizec_patch;
    jsizec_list(3) = jsizec_group2;
    jsizec_list(4) = jsizec_dist;
    jsizec_list(5) = zi_yn;
    jsizec_list(6) = jsizec_zi;
    jsizec_list(7) = year_jsizec_zi;
    jsizec_list(8) = patch_jsizec_zi;
    jsizec_list(9) = group2_jsizec_zi;
    jsizec_list(10) = indcova_names;
    jsizec_list(11) = jsizec_indcova2;
    jsizec_list(12) = jsizec_indcova2_zi;
    jsizec_list(13) = indcovb_names;
    jsizec_list(14) = jsizec_indcovb2;
    jsizec_list(15) = jsizec_indcovb2_zi;
    jsizec_list(16) = indcovc_names;
    jsizec_list(17) = jsizec_indcovc2;
    jsizec_list(18) = jsizec_indcovc2_zi;
    jsizec_list(19) = year_names;
    jsizec_list(20) = patch_names;
    jsizec_list(21) = group_names;
    jsizec_list(22) = main_effect_1;
    jsizec_list(23) = main_effect_2;
    jsizec_list(24) = jsizec_st;
    jsizec_list(25) = effects_names;
    jsizec_list(26) = jsizec_group1;
    jsizec_list(27) = group1_jsizec_zi;
    jsizec_list(28) = jsizec_indcova1;
    jsizec_list(29) = jsizec_indcovb1;
    jsizec_list(30) = jsizec_indcovc1;
    jsizec_list(31) = jsizec_indcova1_zi;
    jsizec_list(32) = jsizec_indcovb1_zi;
    jsizec_list(33) = jsizec_indcovc1_zi;
    
    List jrepst_list(34);
    jrepst_list(0) = jrepst_num;
    jrepst_list(1) = jrepst_year;
    jrepst_list(2) = jrepst_patch;
    jrepst_list(3) = jrepst_group2;
    jrepst_list(4) = jrepst_dist;
    jrepst_list(5) = false;
    jrepst_list(6) = dud_zi;
    jrepst_list(7) = dud_zi;
    jrepst_list(8) = dud_zi;
    jrepst_list(9) = dud_zi;
    jrepst_list(10) = indcova_names;
    jrepst_list(11) = jrepst_indcova2;
    jrepst_list(12) = dud_zi;
    jrepst_list(13) = indcovb_names;
    jrepst_list(14) = jrepst_indcovb2;
    jrepst_list(15) = dud_zi;
    jrepst_list(16) = indcovc_names;
    jrepst_list(17) = jrepst_indcovc2;
    jrepst_list(18) = dud_zi;
    jrepst_list(19) = year_names;
    jrepst_list(20) = patch_names;
    jrepst_list(21) = group_names;
    jrepst_list(22) = main_effect_1;
    jrepst_list(23) = main_effect_2;
    jrepst_list(24) = 1.0;
    jrepst_list(25) = effects_names;
    jrepst_list(26) = jrepst_group1;
    jrepst_list(27) = dud_zi;
    jrepst_list(28) = jrepst_indcova1;
    jrepst_list(29) = jrepst_indcovb1;
    jrepst_list(30) = jrepst_indcovc1;
    jrepst_list(31) = dud_zi;
    jrepst_list(32) = dud_zi;
    jrepst_list(33) = dud_zi;
    
    List jmatst_list(34);
    jmatst_list(0) = jmatst_num;
    jmatst_list(1) = jmatst_year;
    jmatst_list(2) = jmatst_patch;
    jmatst_list(3) = jmatst_group2;
    jmatst_list(4) = jmatst_dist;
    jmatst_list(5) = false;
    jmatst_list(6) = dud_zi;
    jmatst_list(7) = dud_zi;
    jmatst_list(8) = dud_zi;
    jmatst_list(9) = dud_zi;
    jmatst_list(10) = indcova_names;
    jmatst_list(11) = jmatst_indcova2;
    jmatst_list(12) = dud_zi;
    jmatst_list(13) = indcovb_names;
    jmatst_list(14) = jmatst_indcovb2;
    jmatst_list(15) = dud_zi;
    jmatst_list(16) = indcovc_names;
    jmatst_list(17) = jmatst_indcovc2;
    jmatst_list(18) = dud_zi;
    jmatst_list(19) = year_names;
    jmatst_list(20) = patch_names;
    jmatst_list(21) = group_names;
    jmatst_list(22) = main_effect_1;
    jmatst_list(23) = main_effect_2;
    jmatst_list(24) = 1.0;
    jmatst_list(25) = effects_names;
    jmatst_list(26) = jmatst_group1;
    jmatst_list(27) = dud_zi;
    jmatst_list(28) = jmatst_indcova1;
    jmatst_list(29) = jmatst_indcovb1;
    jmatst_list(30) = jmatst_indcovc1;
    jmatst_list(31) = dud_zi;
    jmatst_list(32) = dud_zi;
    jmatst_list(33) = dud_zi;
    
    current_surv_model = surv_list;
    current_obs_model = obs_list;
    current_size_model = sizea_list;
    current_sizeb_model = sizeb_list;
    current_sizec_model = sizec_list;
    current_repst_model = repst_list;
    current_fec_model = fec_list;
    
    current_jsurv_model = jsurv_list;
    current_jobs_model = jobs_list;
    current_jsize_model = jsizea_list;
    current_jsizeb_model = jsizeb_list;
    current_jsizec_model = jsizec_list;
    current_jrepst_model = jrepst_list;
    current_jmatst_model = jmatst_list;
    
    current_surv_model.attr("names") = list_names;
    current_obs_model.attr("names") = list_names;
    current_size_model.attr("names") = list_names;
    current_sizeb_model.attr("names") = list_names;
    current_sizec_model.attr("names") = list_names;
    current_repst_model.attr("names") = list_names;
    current_fec_model.attr("names") = list_names;
    current_jsurv_model.attr("names") = list_names;
    current_jobs_model.attr("names") = list_names;
    current_jsize_model.attr("names") = list_names;
    current_jsizeb_model.attr("names") = list_names;
    current_jsizec_model.attr("names") = list_names;
    current_jrepst_model.attr("names") = list_names;
    current_jmatst_model.attr("names") = list_names;
    
    DataFrame c_paramnames = paramnames_skeleton(true);
    CharacterVector modelparams = as<CharacterVector>(c_paramnames["modelparams"]);
    CharacterVector mainparams = as<CharacterVector>(c_paramnames["mainparams"]);
    CharacterVector parameter_names = as<CharacterVector>(c_paramnames["parameter_names"]);
    
    bool current_check = false;
    for (int i = 0; i < modelparams.length(); i++) {
      for (int j = 0; j < 17; j++) {
        current_check = stringcompare_hard(as<std::string>(mainparams(i)), 
          as<std::string>(main_effect_1(j)));
        if (current_check) modelparams(i) = main_effect_1(j);
      }
    }
    
    current_paramnames = DataFrame::create(_["parameter_names"] = parameter_names,
      _["mainparams"] = mainparams, _["modelparams"] = modelparams);
    
    CharacterVector current_mainyears = as<CharacterVector>(year_vec);
    //unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
    
    CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
    CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
    
    DataFrame indcova2_frame = as<DataFrame>(current_vrm["indcova2_frame"]);
    DataFrame indcovb2_frame = as<DataFrame>(current_vrm["indcovb2_frame"]);
    DataFrame indcovc2_frame = as<DataFrame>(current_vrm["indcovc2_frame"]);
    CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
    CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
    CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
    
    List surv_proxy = modelextract(current_surv_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List obs_proxy = modelextract(current_obs_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List size_proxy = modelextract(current_size_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List sizeb_proxy = modelextract(current_sizeb_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List sizec_proxy = modelextract(current_sizec_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List repst_proxy = modelextract(current_repst_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List fec_proxy = modelextract(current_fec_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    
    List jsurv_proxy = modelextract(current_jsurv_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List jobs_proxy = modelextract(current_jobs_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List jsize_proxy = modelextract(current_jsize_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List jsizeb_proxy = modelextract(current_jsizeb_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List jsizec_proxy = modelextract(current_jsizec_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List jrepst_proxy = modelextract(current_jrepst_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    List jmatst_proxy = modelextract(current_jmatst_model, current_paramnames,
      current_mainyears, current_mainpatches, current_maingroups,
      current_mainindcova, current_mainindcovb, current_mainindcovc, true);
    
    List current_vrm_extract (15);
    current_vrm_extract(0) = surv_proxy;
    current_vrm_extract(1) = obs_proxy;
    current_vrm_extract(2) = size_proxy;
    current_vrm_extract(3) = sizeb_proxy;
    current_vrm_extract(4) = sizec_proxy;
    current_vrm_extract(5) = repst_proxy;
    current_vrm_extract(6) = fec_proxy;
    current_vrm_extract(7) = jsurv_proxy;
    current_vrm_extract(8) = jobs_proxy;
    current_vrm_extract(9) = jsize_proxy;
    current_vrm_extract(10) = jsizeb_proxy;
    current_vrm_extract(11) = jsizec_proxy;
    current_vrm_extract(12) = jrepst_proxy;
    current_vrm_extract(13) = jmatst_proxy;
    current_vrm_extract(14) = current_paramnames;
    
    allmodels_all = current_vrm_extract;
  }
  
  //Rcout << "cleanup3_inv S ";
  
  // Output processing
  List out_list (68);
  
  out_list(0) = chosen_mpm;
  out_list(1) = static_cast<int>(preexisting); // Used to be mpm_count, now defunct
  out_list(2) = vrm_list;
  out_list(3) = static_cast<int>(funcbased); // Used to be vrm_count, now defunct
  out_list(4) = final_stageframe;
  out_list(5) = stageframe_df;
  out_list(6) = supplement_df;
  out_list(7) = supplement_df;
  out_list(8) = final_repmatrix;
  out_list(9) = sparse_bool;
  out_list(10) = sparse_vec_count;
  out_list(11) = format_int;
  out_list(12) = pure_leslie;
  out_list(13) = stageframe_notNull_count;
  out_list(14) = preexisting;
  out_list(15) = funcbased;
  out_list(16) = firstage_int;
  out_list(17) = finalage_int;
  out_list(18) = cont_int;
  out_list(19) = fecmod_num;
  out_list(20) = fecage_min_int;
  out_list(21) = fecage_max_int;
  out_list(22) = final_hstages;
  out_list(23) = final_agestages;
  out_list(24) = matrowcounts;
  out_list(25) = stagecounts;
  out_list(26) = start_list;
  out_list(27) = start_count;
  out_list(28) = labels_list;
  out_list(29) = labels;
  out_list(30) = patch_vec;
  out_list(31) = year_vec;
  out_list(32) = total_years_int;
  out_list(33) = tweights_list;
  out_list(34) = tweights_count;
  out_list(35) = tweights_type_int;
  out_list(36) = density_df;
  out_list(37) = dens_index_df;
  out_list(38) = dens_yn_bool;
  out_list(39) = density_count;
  out_list(40) = entry_time_vec;
  out_list(41) = entry_time_count;
  out_list(42) = entry_time_vec_use;
  out_list(43) = chosen_density_vr;
  out_list(44) = ind_terms_num_list;
  out_list(45) = ind_terms_cat_list;
  out_list(46) = dev_terms_list;
  out_list(47) = dens_vr_yn_bool;
  out_list(48) = sp_density_num_int;
  out_list(49) = dev_terms_num_int;
  out_list(50) = inda_terms_num_int;
  out_list(51) = indb_terms_num_int;
  out_list(52) = indc_terms_num_int;
  out_list(53) = inda_terms_cat_int;
  out_list(54) = indb_terms_cat_int;
  out_list(55) = indc_terms_cat_int;
  out_list(56) = density_vr_count;
  out_list(57) = sparse_vec_count;
  out_list(58) = sp_density_list;
  out_list(59) = equivalence_frame;
  out_list(60) = equivalence_vec;
  out_list(61) = equivalence_count;
  out_list(62) = stages_not_equal;
  out_list(63) = current_mpm_allstages;
  out_list(64) = allmodels_all;
  out_list(65) = historical;
  out_list(66) = preexisting_mpm_size;
  out_list(67) = prebreeding_bool;
  
  CharacterVector out_list_names = {"chosen_mpm", "mpm_count", "vrm_list",
    "vrm_count", "final_stageframe", "stageframe_df", "supplement_df", "supplement_list_fb",
    "repmatrix_list", "sparse_bool", "sparse_vec_count",
    "format_int", "pure_leslie", "stageframe_notNull_count", "preexisting",
    "funcbased", "firstage_int", "finalage_int", "cont_int",
    "fecmod_num", "fecage_min_int", "fecage_max_int", "hstages_list",
    "agestages_list", "matrowcounts", "stagecounts", "start_list",
    "start_count", "labels_list", "labels", "patch_vec", "year_list",
    "total_years_int", "tweights_list", "tweights_count", "tweights_type_int",
    "density_df", "dens_index_df", "dens_yn_bool", "density_count",
    "entry_time_vec", "entry_time_count", "entry_time_vec_use",
    "density_vr_list", "ind_terms_num_list", "ind_terms_cat_list",
    "dev_terms_list", "dens_vr_yn_bool", "sp_density_num_int",
    "dev_terms_num_int", "inda_terms_num_int", "indb_terms_num_int",
    "indc_terms_num_int", "inda_terms_cat_int", "indb_terms_cat_int",
    "indc_terms_cat_int", "density_vr_count", "sparse_vec_count",
    "sp_density_list", "equivalence_list", "equivalence_vec",
    "equivalence_count", "stages_not_equal", "allstages_all", "allmodels_all",
    "historical", "preexisting_mpm_size", "prebreeding_bool"};
  out_list.attr("names") = out_list_names;
  
  return(out_list);
}

//' Engine Running Invasion Analyses of Existing MPMs
//' 
//' Function \code{invade3_pre_core} is the main function running invasion
//' analyses with existing MPMs supplied by the user.
//' 
//' @name invade3_pre_core
//' 
//' @param Lyapunov The main data frame giving the Lyapunov coefficients
//' estimated, as well as the circumstances resulting in them. See \code{Value}
//' section below for further details.
//' @param var_run_mat A matrix giving the order of trait variants in each
//' run.
//' @param N_out The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param zero_stage_vec_list A list of vectors giving zero stage vectors for
//' each MPM, if entry times are staggered.
//' @param trait_axis A data frame of class \code{adaptAxis} holding the trait
//' data to test.
//' @param new_trait_axis A data frame giving trait axis data post-processing
//' with function \code{ta_reassess()}.
//' @param new_stageexpansion_list A list with stage expansions for all trait
//' axis data leading to matrix element changes with each list element
//' corresponding to each respective variant.
//' @param errcheck_mpms An optional list of all MPMs post-processing. Only
//' output if \code{err_check = TRUE}.
//' @param chosen_mpm An MPM in \code{lefkoMat} format.
//' @param tweights_list The tweights vector or matrix covering the MPM.
//' @param start_list The data frame or vector of starting information, ideally
//' supplied in \code{lefkoSV} format.
//' @param vrm_list The unextracted \code{vrm_input} object.
//' @param stageframe_df A stageframe object covering the MPM.
//' @param allmodels_all The list holding the extracted vrm input.
//' @param allstages_all The allstages indexing data frame used to produce MPMs.
//' @param supplement_df A supplement in \code{lefkoSD} format.
//' @param chosen_years A string vector giving the main years used.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage in each MPM relative to a reference individual. Each element of the
//' list corresponds to each respective MPM.
//' @param sp_density_num_vec A vector giving the number of
//' spatial density terms per MPM.
//' @param entry_time_vec An integer vector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param inda_terms_num_vec A vector giving the number of
//' numeric terms given in individual covariate a.
//' @param indb_terms_num_vec A vector giving the number of
//' numeric terms given in individual covariate b.
//' @param indc_terms_num_vec A vector giving the number of
//' numeric terms given in individual covariate c.
//' @param inda_terms_cat_vec A vector giving the number of
//' factor terms given in individual covariate a.
//' @param indb_terms_cat_vec A vector giving the number of
//' factor terms given in individual covariate b.
//' @param indc_terms_cat_vec A vector giving the number of
//' factor terms given in individual covariate c.
//' @param dens_vr_yn_vec A vector stating whether density dependence is used
//' in each MPM, given through \code{lefkoDensVR} objects.
//' @param tweights_type_vec A vector giving the style of \code{tweights} used
//' in each MPM.
//' @param fecmod_vec A numeric vector giving fecmod values.
//' @param patch_vec A vector giving the name of each patch used in projection.
//' @param variant_count An integer giving the number of variants to run in
//' invasion analysis.
//' @param var_per_run The number of variants to run in each projection.
//' @param nreps An integer giving the number of replicates to perform.
//' @param times An integer giving the amount of time steps to run the
//' projection for.
//' @param fitness_times An integer giving how many time steps at the end of
//' each run to use to estimate fitness.
//' @param stagecounts Integer denoting the number of stages in the MPM.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param format_int An integer giving the MPM format.
//' @param preexisting_mpm_size The number of elements per MPM matrix.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param stochastic A Boolean value indicating to perform a temporally
//' stochastic projection.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used in the MPM, given through a \code{lefkoDens} object.
//' @param entry_time_vec_use A Boolean value indicating whether entry
//' times differ among MPMs.
//' @param sparse_bool A Boolean value stating whether the MPM is in sparse
//' matrix format.
//' @param historical A Boolean value indicating whether the MPM is historical.
//' @param pure_leslie A Boolean value indicating whether the MPM is Leslie.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
void invade3_pre_core (DataFrame& Lyapunov, const arma::mat& var_run_mat,
  List& N_out, List& comm_out, List& zero_stage_vec_list, DataFrame& trait_axis,
  DataFrame& new_trait_axis, List& new_stageexpansion_list, List& errcheck_mpms,
  const List chosen_mpm, const List tweights_list, const List start_list,
  const List vrm_list, DataFrame stageframe_df, const List allmodels_all,
  const List allstages_all, const DataFrame supplement_df,
  const CharacterVector chosen_years, const List sp_density_list,
  const DataFrame density_df, const DataFrame dens_index_df,
  const List equivalence_list, const IntegerVector sp_density_num_vec,
  const IntegerVector entry_time_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_vr_yn_vec,
  const IntegerVector tweights_type_vec, const NumericVector fecmod_vec,
  const CharacterVector patch_vec, const int variant_count, const int var_per_run,
  const int nreps, const int times, const int fitness_times, const int stagecounts,
  const int substoch, const int format_int, const int preexisting_mpm_size,
  const int firstage_int, const int finalage_int, const double exp_tol,
  const double theta_tol, const bool integeronly, const bool stages_not_equal,
  const bool stochastic, const bool dens_yn_bool, const bool entry_time_vec_use,
  const bool sparse_bool, const bool historical, const bool pure_leslie,
  const bool A_only, const bool err_check, const bool err_check_extreme) {
  
  // patches?????
  
  //Rcout << "invade3_pre_core a          ";
  int var_mat_length = static_cast<int>(var_run_mat.n_rows);
  
  List A_list = as<List>(chosen_mpm["A"]);
  
  List U_list;
  List F_list;
  
  if (!A_only) {
    U_list = as<List>(chosen_mpm["U"]);
    F_list = as<List>(chosen_mpm["F"]);
  }
  
  //Rcout << "invade3_pre_core b          ";
  // Vectors of axis variants
  IntegerVector axis_variant_vec = as<IntegerVector>(trait_axis["variant"]);
  IntegerVector axis_variants_unique = sort_unique(axis_variant_vec);
  
  DataFrame trait_axis_clone = clone(trait_axis);
  new_trait_axis = AdaptUtils::ta_reassess(stageframe_df, trait_axis_clone, firstage_int,
    historical, false, pure_leslie);
  
  //Rcout << "firstage_int: " << firstage_int << endl;
  //Rcout << "finalage_int: " << finalage_int << endl;
  //Rcout << "format_int: " << format_int << endl;
  
  //Rcout << "invade3_pre_core c          ";
  // Matrix order set up and creation of zero stage vectors
  IntegerVector chosen_mats; // Eliminated matrix_choice_list
  List zvl (variant_count);
  
  IntegerVector chosen_matrix_vec;
  
  DataFrame chosen_labels = as<DataFrame>(chosen_mpm["labels"]);
  CharacterVector chosen_labels_names = as<CharacterVector>(chosen_labels.attr("names"));
  IntegerVector clm_y2 = index_l3(chosen_labels_names, "year2");
  
  CharacterVector mpm_labels_patch = as<CharacterVector>(chosen_labels["patch"]);
  IntegerVector pvis = index_l3(mpm_labels_patch, patch_vec(0)); // Is patch done properly?
  
  //Rcout << "invade3_pre_core d          ";
  if (clm_y2.length() > 0) {
    CharacterVector mpm_labels_year2 = as<CharacterVector>(chosen_labels["year2"]);
    int chosen_years_length = static_cast<int>(chosen_years.length());
    
    IntegerVector yvis = index_l3(mpm_labels_year2, chosen_years(0));
    
    if (chosen_years_length > 1) {
      for (int j = 1; j < chosen_years_length; j++) {
        IntegerVector yvis_next = index_l3(mpm_labels_year2, chosen_years(j));
        IntegerVector yvis_append = concat_int(yvis, yvis_next);
        
        yvis = sort_unique(yvis_append);
      }
    }
    
    IntegerVector chosen_mats_pre = intersect(pvis, yvis);
    chosen_mats = chosen_mats_pre;
  } else {
    IntegerVector chosen_mats_pre = sort_unique(pvis);
    chosen_mats = chosen_mats_pre;
  }
  
  //Rcout << "invade3_pre_core e          ";
  // Zero vector creation
  if (entry_time_vec_use) {
    if (!sparse_bool) {
      List chosen_A_list = as<List>(chosen_mpm["A"]);
      arma::mat chosen_A = as<arma::mat>(chosen_A_list(0));
      
      int proj_mat_rows = static_cast<int>(chosen_A.n_rows);
      
      arma::vec new_zero_vec (proj_mat_rows, fill::zeros);
      for (int i = 0; i < variant_count; i++) {
        zvl(i) =  new_zero_vec;
      }
    } else {
      List chosen_A_list = as<List>(chosen_mpm["A"]);
      arma::sp_mat chosen_A = as<arma::sp_mat>(chosen_A_list(0));
      
      int proj_mat_rows = static_cast<int>(chosen_A.n_rows);
      
      arma::vec new_zero_vec (proj_mat_rows, fill::zeros);
      for (int i = 0; i < variant_count; i++) {
        zvl(i) =  new_zero_vec;
      }
    }
  }
  
  zero_stage_vec_list = zvl;
  
  //Rcout << "invade3_pre_core f          ";
  // Year order determination
  List comm_out_pre (var_mat_length);
  List used_times (var_mat_length);
  
  // Determine order of matrices for each variant run
  for (int i = 0; i < var_mat_length; i++) {
    //Rcout << "invade3_pre_core g          ";
    int chosen_mats_length = static_cast<int>(chosen_mats.length());
    
    List all_pops_per_run (var_per_run);
    List all_used_times_per_run (var_per_run);
    for (int m = 0; m < var_per_run; m++) {
      //Rcout << "invade3_pre_core h          ";
      int current_variant_index = var_run_mat(i, m);
      
      List pop_reps (nreps);
      List used_times_reps (nreps);
      for (int j = 0; j < nreps; j++) {
        //Rcout << "invade3_pre_core i          ";
        arma::mat pops_out_pre (stagecounts, (times + 1), fill::zeros);
      
        IntegerVector years_topull;
        
        if (!stochastic) {
          //Rcout << "invade3_pre_core j          ";
          IntegerVector years_topull_pre (times);
          
          int mat_tracker {0};
          for (int k = 0; k < times; k++) {
            if (mat_tracker >= chosen_mats_length) mat_tracker = 0;
            
            years_topull_pre(k) = chosen_mats(mat_tracker);
            mat_tracker++;
          }
          years_topull = years_topull_pre;
          //Rcout << "invade3_pre_core k          ";
        } else {
          //Rcout << "invade3_pre_core l          ";
          if (tweights_type_vec(0) == 0) {
            NumericVector twinput (chosen_mats_length,
              (1.0 / static_cast<double>(chosen_mats_length)));
            years_topull = Rcpp::RcppArmadillo::sample(chosen_mats, times, true,
              twinput);
          } else if (tweights_type_vec(0) == 1) {
            NumericVector twinput = as<NumericVector>(tweights_list(0));
            NumericVector twinput_st = twinput / sum(twinput);
            
            years_topull = Rcpp::RcppArmadillo::sample(chosen_mats, times, true,
              twinput_st);
          } else if (tweights_type_vec(0) == 2) {
            arma::ivec chosen_mats_arma = as<arma::ivec>(chosen_mats);
            arma::mat twinput_mat = as<arma::mat>(tweights_list(0));
            arma::vec twinput = twinput_mat.col(0);
            twinput = twinput / sum(twinput);
            
            IntegerVector years_topull_pre (times);
            NumericVector twinput_setup (chosen_mats_length, (1.0 / 
              static_cast<double>(chosen_mats_length)));
            arma::ivec first_choice = Rcpp::RcppArmadillo::sample(chosen_mats_arma,
              times, true, twinput_setup);
            years_topull_pre(0) = chosen_mats(first_choice(0));
            
            for (int k = 1; k < times; k++) {
              arma::ivec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(chosen_mats_arma,
                1, true, twinput);
              years_topull_pre(k) = theprophecy_piecemeal(0);
              
              arma::uvec tnotb_preassigned = 
                find(chosen_mats_arma == theprophecy_piecemeal(0));
              twinput = twinput_mat.col(static_cast<int>(tnotb_preassigned(0)));
              twinput = twinput / sum(twinput);
            }
            years_topull = years_topull_pre;
            //Rcout << "invade3_pre_core m          ";
          } else {
            throw Rcpp::exception("tweights_type_vec error.", false);
          }
        }
        used_times_reps(j) = years_topull;
        pop_reps(j) = pops_out_pre;
      }
      all_used_times_per_run(m) = used_times_reps;
      all_pops_per_run(m) = pop_reps;
    }
    used_times(i) = all_used_times_per_run;
    comm_out_pre(i) = all_pops_per_run;
  }
  
  //Rcout << "invade3_pre_core o          ";
  // Here we create the modified A matrices
  List new_A_list = clone(A_list);
  CharacterVector nta_vars = new_trait_axis.names();
  
  List trait_axis_by_variant (variant_count); // Might wish to remove this later
  List stageexpansion_by_variant (variant_count); // Might wish to remove this later
  
  for (int i = 0; i < variant_count; i++) {
    IntegerVector used_i = {i + 1};
    StringVector focused_var = {"variant"};
    DataFrame current_trait_axis = LefkoUtils::df_subset(new_trait_axis, as<RObject>(used_i),
      false, true, false, false, true, as<RObject>(focused_var));
    trait_axis_by_variant(i) = current_trait_axis;
    
    int ehrlen {1};
    int style {0};
    int filter {1};
    
    if (format_int == 2) ehrlen = 2;
    if (format_int == 3 || format_int == 5) style = 1;
    if (format_int == 4) {
      //agemat = true;
      style = 2;
      filter = 2;
    }
    
    DataFrame stageexpansion = AdaptMats::thenewpizzle(stageframe_df, current_trait_axis,
      firstage_int, finalage_int, ehrlen, style, filter);
    focused_var = {"mpm_altered"};
    IntegerVector chosen_int = {1};
    DataFrame stageexpansion_reduced = LefkoUtils::df_subset(stageexpansion,
      as<RObject>(chosen_int), false, true, false, false, true, as<RObject>(focused_var));
    stageexpansion_by_variant(i) = stageexpansion_reduced;
  }
  new_stageexpansion_list = stageexpansion_by_variant;
  
  //Rcout << "invade3_pre_core p          ";
  // Main projection
  List N_out_pre (nreps);
  List errcheck_mpm_reps (nreps);
  
  for (int i = 0; i < nreps; i ++) { // 1st loop - reps i
    //Rcout << "invade3_pre_core q          ";
    List running_popvecs; //  = clone(start_list)
    List running_popvecs_startonly; //  = clone(start_list)
    arma::cube N_mpm (var_per_run, (times + 1), var_mat_length); // rows = vars, cols = times, slices = permutes 
    
    List errcheck_mpm_reps_time (times); // Could remove later
    for (int j = 0; j < times; j++) { // 2nd loop - time j
      //Rcout << "invade3_pre_core r          ";
      if (j % 10 == 0){
        Rcpp::checkUserInterrupt();
      }
      
      List errcheck_mpm_reps_time_vmt (var_mat_length); // Could remove later
      
      for (int l = 0; l < var_mat_length; l++) { // 3rd loop - permutes l
        //Rcout << "invade3_pre_core s          ";
        List errcheck_mpm_reps_time_vmt_var(var_per_run); // Could remove later
        //Rcout << "current_permutation (l): " << l << "          ";
        if (j == 0) {
          List var_popvecs_to_start (var_per_run);
          for (int n = 0; n < var_per_run; n++) {
            var_popvecs_to_start(n) = as<arma::vec>(start_list(static_cast<int>(var_run_mat(l, n))));
          }
          running_popvecs = var_popvecs_to_start;
          running_popvecs_startonly = clone(var_popvecs_to_start);
        }
        
        for (int m = 0; m < var_per_run; m++) {
          if (j == entry_time_vec(m)) {
            arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
            
            double N_current = accu(running_popvec_mpm);
            N_mpm(m, j, l) = N_current; // Used to be (k, j)
          }
        }
        
        List all_pops_per_run = as<List>(comm_out_pre(l));
        for (int m = 0; m < var_per_run; m++) { // 4th loop - var per run m
          //Rcout << "invade3_pre_core t          ";
          int current_variant_index = var_run_mat(l, m); // Equivalent to index integer k
          
          DataFrame sge_current = as<DataFrame>(new_stageexpansion_list(current_variant_index));
          
          List pop_reps = as<List>(all_pops_per_run(m));
          arma::mat pops_out = as<arma::mat>(pop_reps(i));
          
          if (j > (entry_time_vec(m) - 1)) {
            //Rcout << "invade3_pre_core u          ";
            List used_times_per_run = as<List>(used_times(l));
            List used_times_current_var = as<List>(used_times_per_run(m));
            IntegerVector current_times_vec = as<IntegerVector>(used_times_current_var(i));
            
            arma::vec running_popvec_mpm;
            if (j == entry_time_vec(m)) {
              running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
              pops_out.col(j) = running_popvec_mpm;
              
            } else {
              running_popvec_mpm = pops_out.col(j);
            }
            
            //Rcout << "invade3_pre_core v          ";
            if (!dens_yn_bool) {
              if (!sparse_bool) {
                arma::mat current_A = as<arma::mat>(A_list(current_times_vec(j)));

                if (A_only) {
                  AdaptUtils::Amat_alter(current_A, sge_current); 
                } else {
                  arma::mat current_U_unaltered = as<arma::mat>(U_list(current_times_vec(j)));
                  arma::mat current_F_unaltered = as<arma::mat>(F_list(current_times_vec(j)));
                  AdaptUtils::UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
                }
                
                if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
                running_popvec_mpm = current_A * running_popvec_mpm; 
              } else {
                arma::sp_mat current_A = as<arma::sp_mat>(A_list(current_times_vec(j)));
                
                if (A_only) {
                  AdaptUtils::sp_Amat_alter(current_A, sge_current);
                } else {
                  arma::sp_mat current_U_unaltered = as<arma::sp_mat>(U_list(current_times_vec(j)));
                  arma::sp_mat current_F_unaltered = as<arma::sp_mat>(F_list(current_times_vec(j)));
                  AdaptUtils::sp_UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
                }
                if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
                running_popvec_mpm = current_A * running_popvec_mpm;
              }
            } else {
              DataFrame used_density_input = density_df;
              DataFrame used_density_index_input = as<DataFrame>(dens_index_df);
              
              IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
              int used_delay = max(ud_delay_vec);
              
              if (j >= (used_delay - 1 )) { // Change to allow different delay Ns for different entries
                if (!stages_not_equal) {
                  arma::mat di_mat = N_mpm.slice(l);
                  arma::vec delay_issue = di_mat.col(j + 1 - used_delay);
                  
                  double delay_N_sum = arma::sum(delay_issue);
                  arma::vec new_popvec;
                  arma::mat new_projmat;
                  arma::sp_mat new_projsp;
                  
                  AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                    running_popvec_mpm, sge_current, A_list, delay_N_sum,
                    static_cast<int>(current_times_vec(j)), integeronly,
                    substoch, used_density_input, used_density_index_input,
                    false, sparse_bool, sparse_bool, false, err_check);
                  
                  running_popvec_mpm = new_popvec;
                  if (err_check_extreme) { // Could remove later
                    if (!sparse_bool) {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                    } else {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                    }
                  } // Could remove later
                } else {
                  double delay_N_sum {0.0};
                  
                  if (j > 0) {
                    for (int p = 0; p < var_per_run; p++) {
                      int current_variant_index_agg = var_run_mat(l, p);
                      
                      List current_pop_list = as<List>(comm_out_pre(l));
                      List pop_rep_list = as<List>(current_pop_list(p)); // Changed from m to p
                      arma::mat delay_pop = as<arma::mat>(pop_rep_list(i));
                      
                      arma::vec delay_pop_vec = delay_pop.col(j + 1 - used_delay);
                      arma::vec current_equiv_vec = as<arma::vec>(equivalence_list(current_variant_index_agg));
                      arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                      double delay_pop_N = arma::accu(adjusted_delay_pop_vec);
                      
                      delay_N_sum += delay_pop_N;
                    }
                  }
                  
                  arma::vec new_popvec;
                  arma::mat new_projmat;
                  arma::sp_mat new_projsp;
                  
                  AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                    running_popvec_mpm, sge_current, A_list, delay_N_sum,
                    static_cast<int>(current_times_vec(j)), integeronly,
                    substoch, used_density_input, used_density_index_input,
                    false, sparse_bool, sparse_bool, false, err_check);
                  
                  running_popvec_mpm = new_popvec;
                  if (err_check_extreme) { // Could remove later
                    if (!sparse_bool) {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                    } else {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                    }
                  } // Could remove later
                }
                //Rcout << "invade3_pre_core aa          ";
              } else {
                //Rcout << "invade3_pre_core ab          ";
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, sge_current, A_list, 0.0,
                  static_cast<int>(current_times_vec(j)), integeronly, substoch,
                  used_density_input, used_density_index_input, false,
                  sparse_bool, sparse_bool, false, err_check);
                
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) { // Could remove later
                  if (!sparse_bool) {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                  } else {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                  }
                } // Could remove later
                //Rcout << "invade3_pre_core ac          ";
              }
            }
            
            //Rcout << "invade3_pre_core ad          ";
            if (integeronly) running_popvec_mpm = floor(running_popvec_mpm);
            double N_current = arma::sum(running_popvec_mpm);
            N_mpm(m, (j + 1), l) = N_current; // Used to be (k, (j + 1))
            
            running_popvecs(m) = running_popvec_mpm;
            pops_out.col(j + 1) = running_popvec_mpm;
            //Rcout << "invade3_pre_core ae          ";
          } else {
            //Rcout << "invade3_pre_core af          ";
            arma::vec current_zero_vec = as<arma::vec>(zero_stage_vec_list(current_variant_index));
            pops_out.col(j + 1) = current_zero_vec;
          }
          pop_reps(i) = pops_out;
        } // m loop - var_per_run
        if (err_check_extreme) errcheck_mpm_reps_time_vmt(l) = errcheck_mpm_reps_time_vmt_var;
      } // l loop - var_mat_length
      if (err_check_extreme) errcheck_mpm_reps_time(j) = errcheck_mpm_reps_time_vmt;
      
    } // j loop - time
    if (err_check_extreme) errcheck_mpm_reps(i) = errcheck_mpm_reps_time; // Could remove later
    N_out_pre(i) = N_mpm;
  } // i loop - reps
  comm_out = comm_out_pre;
  if (err_check_extreme) errcheck_mpms = errcheck_mpm_reps;
  N_out = N_out_pre;
  
  //Rcout << "invade3_pre_core ag     ";
  
  IntegerVector Lyap_rows = seq(1, nreps * var_mat_length);
  CharacterVector Lyap_df_names ((3 * var_per_run) + 2);
  Lyap_df_names(0) = "simulation_num";
  Lyap_df_names(1) = "rep";
  
  String Lyap_df_names_base1 = "variant";
  String Lyap_df_names_base2 = "entrytime";
  String Lyap_df_names_base3 = "fitness_variant";
  
  List output ((3 * var_per_run) + 2);
  output(0) = Lyap_rows;
  
  //Rcout << "invade3_pre_core ah     ";
  
  int Lyap_counter {0};
  for (int m = 0; m < var_per_run; m++) {
    Lyap_counter = 0;
    
    IntegerVector Lyap_var_orig = as<IntegerVector>(Lyapunov(1+m));
    IntegerVector Lyap_rep (nreps * var_mat_length);
    IntegerVector Lyap_var (nreps * var_mat_length);
    IntegerVector Lyap_etime (nreps * var_mat_length);
    NumericVector Lyap_fitness (nreps * var_mat_length);
    
    for (int i = 0; i < nreps; i++) {
      arma::cube current_N_out_cube = as<arma::cube>(N_out(i));
      
      for (int l = 0; l < var_mat_length; l++) {
        arma::mat current_N_out = current_N_out_cube.slice(l);
        int time_length = static_cast<int>(current_N_out.n_cols);
        int start_time = time_length - fitness_times;
        
        arma::vec running_fitness(fitness_times, fill::zeros);
        arma::vec generations(fitness_times, fill::zeros);
        arma::vec intercept_ones(fitness_times, fill::ones);
        
        for (int k = 0; k < fitness_times; k++) {
          running_fitness(k) = current_N_out(m, k+start_time);
          generations(k) = static_cast<double>(k);
        }
        
        arma::mat xmat = join_rows(intercept_ones, generations);
        arma::vec Lyap_regr;
        
        AdaptUtils::fastLm_sl(Lyap_regr, running_fitness, xmat);
        
        double Lyapunov_estimate = Lyap_regr(1);
        
        Lyap_var(Lyap_counter) = Lyap_var_orig(l);
        Lyap_rep(Lyap_counter) = i+1;
        Lyap_etime(Lyap_counter) = entry_time_vec(m);
        Lyap_fitness(Lyap_counter) = Lyapunov_estimate;
        
        Lyap_counter++;
      }
    }
    output(1) = Lyap_rep;
    output(2 + m) = Lyap_var;
    output(2 + var_per_run + m) = Lyap_etime;
    output(2 + (2 * var_per_run) + m) = Lyap_fitness;
    
    String new_col_name = Lyap_df_names_base1;
    new_col_name += (m + 1);
    Lyap_df_names(2 + m) = new_col_name;
    
    String next_col_name = Lyap_df_names_base2;
    next_col_name += (m + 1);
    Lyap_df_names(2 + var_per_run + m) = next_col_name;
    
    String last_col_name = Lyap_df_names_base3;
    last_col_name += (m + 1);
    Lyap_df_names(2 + (2 * var_per_run) + m) = last_col_name;
  }
  
  output.attr("names") = Lyap_df_names;
  output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
    nreps * var_mat_length);
  StringVector needed_class {"data.frame"};
  output.attr("class") = needed_class;
    
  Lyapunov = output;
}

//' Engine Running Invasion Analyses of Function-based MPMs
//' 
//' Function \code{invade3_fb_core} is the main function running invasion
//' analyses in which matrices must be created at each time step.
//' 
//' @name invade3_fb_core
//' 
//' @param Lyapunov The main data frame giving the Lyapunov coefficients
//' estimated, as well as the circumstances resulting in them. See \code{Value}
//' section below for further details.
//' @param N_out The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param zero_stage_vec_list A list of vectors giving zero stage vectors for
//' each MPM, if entry times are staggered.
//' @param trait_axis A data frame of class \code{adaptAxis} holding the trait
//' data to test.
//' @param new_trait_axis A data frame giving trait axis data post-processing
//' with function \code{ta_reassess()}.
//' @param new_stageexpansion_list A list with stage expansions for all trait
//' axis data leading to matrix element changes with each list element
//' corresponding to each respective variant.
//' @param modified_dev_terms_list An optional list giving the vital rate
//' y-intercept deviations by variant once data from the \code{trait_axis} data
//' frame has been allocated.
//' @param errcheck_mpms An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param errcheck_mpmouts An optional list of all mpm_out matrices from MPM
//' processing. Only output if \code{err_check = "extreme"}.
//' @param tweights_list The tweights vector or matrix covering the MPM.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param vrm_list A list of \code{vrm_input} objects.
//' @param current_stageframe The main stageframe, including extra stages.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param allstages_all The allstages indexing data frame used to produce MPMs.
//' @param current_supplement A supplement in \code{lefkoSD} format.
//' @param year_vec A vector giving the main years used.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates.
//' @param dev_terms_list List of deviations for vital rate models.
//' @param density_vr_list Data frame of \code{lefkoDensVR} objects holding
//' density relationships for all 14 vital rate models.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param inda_terms_num_vec A vector giving the number of numeric terms given
//' in individual covariate a.
//' @param indb_terms_num_vec A vector giving the number of numeric terms given
//' in individual covariate b.
//' @param indc_terms_num_vec A vector giving the number of numeric terms given
//' in individual covariate c.
//' @param inda_terms_cat_vec A vector giving the number of factor terms given
//' in individual covariate a.
//' @param indb_terms_cat_vec A vector giving the number of factor terms given
//' in individual covariate b.
//' @param indc_terms_cat_vec A vector giving the number of factor terms given
//' in individual covariate c.
//' @param dens_vr_yn_vec A vector stating whether density dependence is used,
//' given through \code{lefkoDensVR} objects.
//' @param tweights_type_vec A vector giving the style of \code{tweights} used
//' in each MPM.
//' @param fecmod_vec A numeric vector giving the fecmod values.
//' @param patch_vec A vector giving the name of each patch used in projection.
//' @param variant_count An integer giving the number of variants to run in
//' invasion analysis.
//' @param var_per_run The number of variants to run in each projection.
//' @param nreps An integer giving the number of replicates to perform.
//' @param times An integer giving the amount of time steps to run the
//' projection for.
//' @param fitness_times An integer giving how many time steps at the end of
//' each run to use to estimate fitness.
//' @param stagecounts Integer denoting the number of stages in the MPM.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param format_int An integer giving the MPM format.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param dev_terms_times_int An integer giving ??? /////.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param stochastic A Boolean value indicating to perform a temporally
//' stochastic projection.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param historical A Boolean value indicating whether the MPM is historical.
//' @param pure_leslie A Boolean value indicating whether the MPM is Leslie.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
void invade3_fb_core (DataFrame& Lyapunov, const arma::mat& var_run_mat,
  List& N_out, List& comm_out, List& zero_stage_vec_list, DataFrame& trait_axis,
  DataFrame& new_trait_axis, List& new_stageexpansion_list,
  List& modified_dev_terms_list, List& errcheck_mpms, List& errcheck_mpmouts,
  const List tweights_list, const List start_list, const List vrm_list,
  DataFrame current_stageframe, const List allmodels_all, const List allstages_all,
  const DataFrame current_supplement, const CharacterVector year_vec,
  const List ind_terms_num_list, const List ind_terms_cat_list,
  const List dev_terms_list, const DataFrame density_vr_list,
  const List sp_density_list, const DataFrame density_df,
  const DataFrame dens_index_df, const List equivalence_list,
  const IntegerVector sp_density_num_vec, const IntegerVector entry_time_vec,
  const IntegerVector inda_terms_num_vec, const IntegerVector indb_terms_num_vec,
  const IntegerVector indc_terms_num_vec, const IntegerVector inda_terms_cat_vec,
  const IntegerVector indb_terms_cat_vec, const IntegerVector indc_terms_cat_vec,
  const IntegerVector dens_vr_yn_vec, const IntegerVector tweights_type_vec,
  const NumericVector fecmod_vec, const CharacterVector patch_vec,
  const int variant_count, const int var_per_run, const int nreps,
  const int times, const int fitness_times, const int stagecounts,
  const int substoch, const int format_int, const int firstage_int,
  const int finalage_int, const int dev_terms_times_int, const double exp_tol,
  const double theta_tol, const bool integeronly, const bool stochastic,
  const bool dens_yn_bool, const bool stages_not_equal, const bool sparse_bool,
  const bool historical, const bool pure_leslie, const bool A_only,
  const bool err_check, const bool err_check_extreme) {
  
  //Rcout << "invade3_fb_core:" << endl;
  
  // patch_vec???
  // density_vr_list
  // dens_vr_yn_vec
  
  // Vectors of axis variants
  IntegerVector axis_variant_vec = as<IntegerVector>(trait_axis["variant"]);
  IntegerVector axis_variants_unique = sort_unique(axis_variant_vec);
  
  int var_mat_length = static_cast<int>(var_run_mat.n_rows);
  
  DataFrame trait_axis_clone = clone(trait_axis);
  new_trait_axis = AdaptUtils::ta_reassess(current_stageframe, trait_axis_clone, firstage_int,
    historical, true, pure_leslie);
  arma::ivec variant_nta = as<arma::ivec>(new_trait_axis["variant"]);
  arma::vec surv_dev_nta = as<arma::vec>(new_trait_axis["surv_dev"]);
  arma::vec obs_dev_nta = as<arma::vec>(new_trait_axis["obs_dev"]);
  arma::vec size_dev_nta = as<arma::vec>(new_trait_axis["size_dev"]);
  arma::vec sizeb_dev_nta = as<arma::vec>(new_trait_axis["sizeb_dev"]);
  arma::vec sizec_dev_nta = as<arma::vec>(new_trait_axis["sizec_dev"]);
  arma::vec repst_dev_nta = as<arma::vec>(new_trait_axis["repst_dev"]);
  arma::vec fec_dev_nta = as<arma::vec>(new_trait_axis["fec_dev"]);
  arma::vec jsurv_dev_nta = as<arma::vec>(new_trait_axis["jsurv_dev"]);
  arma::vec jobs_dev_nta = as<arma::vec>(new_trait_axis["jobs_dev"]);
  arma::vec jsize_dev_nta = as<arma::vec>(new_trait_axis["jsize_dev"]);
  arma::vec jsizeb_dev_nta = as<arma::vec>(new_trait_axis["jsizeb_dev"]);
  arma::vec jsizec_dev_nta = as<arma::vec>(new_trait_axis["jsizec_dev"]);
  arma::vec jrepst_dev_nta = as<arma::vec>(new_trait_axis["jrepst_dev"]);
  arma::vec jmatst_dev_nta = as<arma::vec>(new_trait_axis["jmatst_dev"]);
  
  int year_counter {0};
  IntegerMatrix inda_num_terms_counter (var_mat_length, var_per_run);
  IntegerMatrix indb_num_terms_counter (var_mat_length, var_per_run);
  IntegerMatrix indc_num_terms_counter (var_mat_length, var_per_run);
  IntegerMatrix inda_cat_terms_counter (var_mat_length, var_per_run);
  IntegerMatrix indb_cat_terms_counter (var_mat_length, var_per_run);
  IntegerMatrix indc_cat_terms_counter (var_mat_length, var_per_run);
  IntegerMatrix inda_num_terms_previous (var_mat_length, var_per_run);
  IntegerMatrix indb_num_terms_previous (var_mat_length, var_per_run);
  IntegerMatrix indc_num_terms_previous (var_mat_length, var_per_run);
  IntegerMatrix inda_cat_terms_previous (var_mat_length, var_per_run);
  IntegerMatrix indb_cat_terms_previous (var_mat_length, var_per_run);
  IntegerMatrix indc_cat_terms_previous (var_mat_length, var_per_run);
  IntegerMatrix dev_num_counter (var_mat_length, var_per_run);
  IntegerMatrix sp_density_counter (var_mat_length, var_per_run);
  
  // Year order determination
  List comm_out_pre (var_mat_length);
  List used_times (var_mat_length);
  
  //Rcout << "invade3_fb_core a" << " ";
  
  
  // This next loop determines the order of matrices for each variant run
  for (int i = 0; i < var_mat_length; i++) {
    int year_length = static_cast<int>(year_vec.length());
    IntegerVector year_int_vec = seq(0, (year_length - 1));
    
    List all_pops_per_run (var_per_run);
    List all_used_times_per_run (var_per_run);
    for (int m = 0; m < var_per_run; m++) {
      int current_variant_index = var_run_mat(i, m);
      
      List pop_reps (nreps);
      List used_times_reps (nreps);
      for (int j = 0; j < nreps; j++) {
        arma::mat pops_out_pre (stagecounts, (times + 1), fill::zeros);
      
        IntegerVector years_topull;
        
        if (!stochastic) {
          IntegerVector years_topull_pre (times);
          
          int mat_tracker {0};
          for (int k = 0; k < times; k++) {
            if (mat_tracker >= year_length) mat_tracker = 0;
            
            years_topull_pre(k) = mat_tracker;
            mat_tracker++;
          }
          years_topull = years_topull_pre;
        } else {
          if (tweights_type_vec(0) == 0) {
            NumericVector twinput (year_length,
              (1.0 / static_cast<double>(year_length)));
            years_topull = Rcpp::RcppArmadillo::sample(year_int_vec, times, true,
              twinput);
          } else if (tweights_type_vec(0) == 1) {
            NumericVector twinput = as<NumericVector>(tweights_list(0));
            NumericVector twinput_st = twinput / sum(twinput);
            
            years_topull = Rcpp::RcppArmadillo::sample(year_int_vec, times, true,
              twinput_st);
          } else if (tweights_type_vec(0) == 2) {
            arma::ivec year_arma = as<arma::ivec>(year_int_vec);
            arma::mat twinput_mat = as<arma::mat>(tweights_list(0));
            arma::vec twinput = twinput_mat.col(0);
            twinput = twinput / sum(twinput);
            
            IntegerVector years_topull_pre (times);
            NumericVector twinput_setup (year_length, (1.0 / 
              static_cast<double>(year_length)));
            arma::ivec first_choice = Rcpp::RcppArmadillo::sample(year_arma,
              times, true, twinput_setup);
            years_topull_pre(0) = year_int_vec(first_choice(0));
            
            for (int k = 1; k < times; k++) {
              arma::ivec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(year_arma,
                1, true, twinput);
              years_topull_pre(k) = theprophecy_piecemeal(0);
              
              arma::uvec tnotb_preassigned = 
                find(year_arma == theprophecy_piecemeal(0));
              twinput = twinput_mat.col(static_cast<int>(tnotb_preassigned(0)));
              twinput = twinput / sum(twinput);
            }
            years_topull = years_topull_pre;
          } else {
            throw Rcpp::exception("tweights_type_vec error.", false);
          }
        }
        used_times_reps(j) = years_topull;
        pop_reps(j) = pops_out_pre;
      }
      all_used_times_per_run(m) = used_times_reps;
      all_pops_per_run(m) = pop_reps;
    }
    used_times(i) = all_used_times_per_run;
    comm_out_pre(i) = all_pops_per_run;
  }
  
  //Rcout << "invade3_fb_core b        ";
  
  List trait_axis_by_variant (variant_count); // Might wish to remove this later
  List stageexpansion_by_variant (variant_count); // Remove this later
  List stageexpansion_ta_devterms_by_variant (variant_count);
  
  for (int i = 0; i < variant_count; i++) {
    IntegerVector used_i = {i + 1};
    StringVector focused_var = {"variant"};
    DataFrame current_trait_axis = LefkoUtils::df_subset(new_trait_axis, as<RObject>(used_i),
      false, true, false, false, true, as<RObject>(focused_var));
    trait_axis_by_variant(i) = current_trait_axis;
    
    int ehrlen {1};
    int style {0};
    int filter {1};
    
    if (format_int == 2) ehrlen = 2;
    if (format_int == 3 || format_int == 5) style = 1;
    if (format_int == 4) {
      //agemat = true;
      style = 2;
      filter = 2;
    }
    
    DataFrame new_sf;
    if (format_int == 1 || format_int == 2 || format_int == 4) {
      DataFrame cloned_sf = clone (current_stageframe);
      
      StringVector csf_stages = as<StringVector>(cloned_sf["stage"]);
      IntegerVector csf_stage_id = as<IntegerVector>(cloned_sf["stage_id"]);
      
      int csf_rows = static_cast<int>(cloned_sf.nrows());
      IntegerVector chosen_rows = {(csf_rows - 1)};
      new_sf = LefkoUtils::df_shedrows(cloned_sf, chosen_rows);
    } else new_sf = current_stageframe;
    
    StringVector nsf_stages = as<StringVector>(new_sf["stage"]);
    IntegerVector nsf_stage_id = as<IntegerVector>(new_sf["stage_id"]);
    
    DataFrame stageexpansion = AdaptMats::thenewpizzle(new_sf, current_trait_axis,
      firstage_int, finalage_int, ehrlen, style, filter);
    focused_var = {"mpm_altered"};
    IntegerVector chosen_int = {1};
    DataFrame stageexpansion_reduced_mpm = LefkoUtils::df_subset(stageexpansion,
      as<RObject>(chosen_int), false, true, false, false, true, as<RObject>(focused_var));
    stageexpansion_by_variant(i) = stageexpansion_reduced_mpm;
    
    NumericVector variant_ta_devterms (14);
    int current_variant_rows = current_trait_axis.nrows();
    for (int j = 0; j < current_variant_rows; j++) {
      NumericVector cv_vt_survdev = as<NumericVector>(current_trait_axis["surv_dev"]);
      NumericVector cv_vt_obsdev = as<NumericVector>(current_trait_axis["obs_dev"]);
      NumericVector cv_vt_sizedev = as<NumericVector>(current_trait_axis["size_dev"]);
      NumericVector cv_vt_sizebdev = as<NumericVector>(current_trait_axis["sizeb_dev"]);
      NumericVector cv_vt_sizecdev = as<NumericVector>(current_trait_axis["sizec_dev"]);
      NumericVector cv_vt_repstdev = as<NumericVector>(current_trait_axis["repst_dev"]);
      NumericVector cv_vt_fecdev = as<NumericVector>(current_trait_axis["fec_dev"]);
      
      NumericVector cv_vt_jsurvdev = as<NumericVector>(current_trait_axis["jsurv_dev"]);
      NumericVector cv_vt_jobsdev = as<NumericVector>(current_trait_axis["jobs_dev"]);
      NumericVector cv_vt_jsizedev = as<NumericVector>(current_trait_axis["jsize_dev"]);
      NumericVector cv_vt_jsizebdev = as<NumericVector>(current_trait_axis["jsizeb_dev"]);
      NumericVector cv_vt_jsizecdev = as<NumericVector>(current_trait_axis["jsizec_dev"]);
      NumericVector cv_vt_jrepstdev = as<NumericVector>(current_trait_axis["jrepst_dev"]);
      NumericVector cv_vt_jmatstdev = as<NumericVector>(current_trait_axis["jmatst_dev"]);
      
      if (!NumericVector::is_na(cv_vt_survdev(j))) variant_ta_devterms(0) = variant_ta_devterms(0) + cv_vt_survdev(j);
      if (!NumericVector::is_na(cv_vt_obsdev(j))) variant_ta_devterms(1) = variant_ta_devterms(1) + cv_vt_obsdev(j);
      if (!NumericVector::is_na(cv_vt_sizedev(j))) variant_ta_devterms(2) = variant_ta_devterms(2) + cv_vt_sizedev(j);
      if (!NumericVector::is_na(cv_vt_sizebdev(j))) variant_ta_devterms(3) = variant_ta_devterms(3) + cv_vt_sizebdev(j);
      if (!NumericVector::is_na(cv_vt_sizecdev(j))) variant_ta_devterms(4) = variant_ta_devterms(4) + cv_vt_sizecdev(j);
      if (!NumericVector::is_na(cv_vt_repstdev(j))) variant_ta_devterms(5) = variant_ta_devterms(5) + cv_vt_repstdev(j);
      if (!NumericVector::is_na(cv_vt_fecdev(j))) variant_ta_devterms(6) = variant_ta_devterms(6) + cv_vt_fecdev(j);
      
      if (!NumericVector::is_na(cv_vt_jsurvdev(j))) variant_ta_devterms(7) = variant_ta_devterms(7) + cv_vt_jsurvdev(j);
      if (!NumericVector::is_na(cv_vt_jobsdev(j))) variant_ta_devterms(8) = variant_ta_devterms(8) + cv_vt_jobsdev(j);
      if (!NumericVector::is_na(cv_vt_jsizedev(j))) variant_ta_devterms(9) = variant_ta_devterms(9) + cv_vt_jsizedev(j);
      if (!NumericVector::is_na(cv_vt_jsizebdev(j))) variant_ta_devterms(10) = variant_ta_devterms(10) + cv_vt_jsizebdev(j);
      if (!NumericVector::is_na(cv_vt_jsizecdev(j))) variant_ta_devterms(11) = variant_ta_devterms(11) + cv_vt_jsizecdev(j);
      if (!NumericVector::is_na(cv_vt_jrepstdev(j))) variant_ta_devterms(12) = variant_ta_devterms(12) + cv_vt_jrepstdev(j);
      if (!NumericVector::is_na(cv_vt_jmatstdev(j))) variant_ta_devterms(13) = variant_ta_devterms(13) + cv_vt_jmatstdev(j);
    }
    stageexpansion_ta_devterms_by_variant(i) = variant_ta_devterms;
  }
  new_stageexpansion_list = stageexpansion_by_variant;
  
  //Rcout << "invade3_fb_core c        ";
  
  // Main projection
  List N_out_pre (nreps);
  List mdtl (var_mat_length);
  List errcheck_mpm_reps (nreps);
  List errcheck_mpmout_reps (nreps);
    
  for (int current_rep = 0; current_rep < nreps; current_rep++) { // 1st loop - reps current_rep
    
    //Rcout << "invade3_fb_core d current_rep: " << current_rep << "      ";
    
    List running_popvecs; // = clone(start_list);
    List running_popvecs_startonly;
    arma::cube N_mpm (var_per_run, (times + 1), var_mat_length);
    
    List errcheck_mpm_reps_time (times); 
    List errcheck_mpmout_reps_time (times); 
    
    for (int j = 0; j < times; j++) { // 2nd loop - time j
      if (j % 10 == 0){
        Rcpp::checkUserInterrupt();
      }
      
      //Rcout << "invade3_fb_core d2 j: " << j << "      ";
      List errcheck_mpm_reps_time_vmt (var_mat_length); 
      List errcheck_mpmout_reps_time_vmt (var_mat_length); 
      
      for (int i = 0; i < var_mat_length; i++) { // 3rd loop - permutes i
        List errcheck_mpm_reps_time_vmt_var(var_per_run); 
        List errcheck_mpmout_reps_time_vmt_var(var_per_run); 
        
        //Rcout << "current_permutation (i): " << i << "          ";
        if (j == 0) {
          List var_popvecs_to_start (var_per_run);
          for (int n = 0; n < var_per_run; n++) {
            var_popvecs_to_start(n) = as<arma::vec>(start_list(static_cast<int>(var_run_mat(i, n))));
          }
          running_popvecs = var_popvecs_to_start;
          running_popvecs_startonly = clone(var_popvecs_to_start);
        }
        
        for (int m = 0; m < var_per_run; m++) { // 4th loop - var per run m
          if (j == entry_time_vec(m)) {
            arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
            
            double N_current = accu(running_popvec_mpm);
            N_mpm(m, j, i) = N_current; // Used to be (k, j)
          }
        }
        
        List all_pops_per_run = as<List>(comm_out_pre(i));
        List mdtl_1 (var_per_run);
        for (int m = 0; m < var_per_run; m++) { // 4th loop - var per run m
          int current_variant_index = var_run_mat(i, m); // Equivalent to index integer k
          
          DataFrame sge_current = as<DataFrame>(new_stageexpansion_list(current_variant_index));
          
          //Rcout << "invade3_fb_core d3 i: " << i << "      ";
          List pop_reps = as<List>(all_pops_per_run(m));
          arma::mat pops_out = as<arma::mat>(pop_reps(current_rep));
          
          //Rcout << "invade3_fb_core d4 m: " << m << "      ";
          //Rcout << "entry_time_vec: " << entry_time_vec << endl;
          if (j > (entry_time_vec(m) - 1)) {
            //Rcout << "invade3_fb_core d5 current_rep: " << current_rep << "      ";
            List used_times_per_run = as<List>(used_times(i));
            List used_times_current_var = as<List>(used_times_per_run(m));
            IntegerVector current_times_vec = as<IntegerVector>(used_times_current_var(current_rep));
            
            
            
            
            
            arma::vec running_popvec_vrm;
            if (j == entry_time_vec(m)) {
              running_popvec_vrm = as<arma::vec>(running_popvecs_startonly(m));
              pops_out.col(j) = running_popvec_vrm;
            } else {
              running_popvec_vrm = pops_out.col(j);
            }
            //Rcout << "invade3_fb_core d6 current_rep: " << current_rep << "      ";
            
            List current_vrm_extract = allmodels_all; // (i)
            List current_vrm_unextract = vrm_list; // (i)
            int ehrlen_format {1}; // This will need to be dealt with differently later
            
            //Rcout << "invade3_fb_core d7 current_rep: " << current_rep << "      ";
            int mpm_style {1};
            if (format_int < 3) {
              mpm_style = 0;
              if (format_int == 2) ehrlen_format = 2;
            } else if (format_int == 4) {
              mpm_style = 2;
            }
            
            //Rcout << "invade3_fb_core d8 current_rep: " << current_rep << "      ";
            DataFrame current_mpm_allstages = allstages_all; // (i)
            
            //Rcout << "invade3_fb_core d9 current_rep: " << current_rep << "      ";
            List surv_proxy = as<List>(current_vrm_extract(0));
            List obs_proxy = as<List>(current_vrm_extract(1));
            List size_proxy = as<List>(current_vrm_extract(2));
            List sizeb_proxy = as<List>(current_vrm_extract(3));
            List sizec_proxy = as<List>(current_vrm_extract(4));
            List repst_proxy = as<List>(current_vrm_extract(5));
            List fec_proxy = as<List>(current_vrm_extract(6));
            List jsurv_proxy = as<List>(current_vrm_extract(7));
            List jobs_proxy = as<List>(current_vrm_extract(8));
            List jsize_proxy = as<List>(current_vrm_extract(9));
            List jsizeb_proxy = as<List>(current_vrm_extract(10));
            List jsizec_proxy = as<List>(current_vrm_extract(11));
            List jrepst_proxy = as<List>(current_vrm_extract(12));
            List jmatst_proxy = as<List>(current_vrm_extract(13));
            DataFrame current_paramnames = as<DataFrame>(current_vrm_extract(14));
            
            //Rcout << "invade3_fb_core d10 current_rep: " << current_rep << "      ";
            CharacterVector current_mainyears = year_vec;
            unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
            
            StringVector cveu_names = as<StringVector>(current_vrm_unextract.names()); // Remove later
            
            DataFrame group2_frame = as<DataFrame>(current_vrm_unextract["group2_frame"]);
            CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
  
            DataFrame patch_frame = as<DataFrame>(current_vrm_unextract["patch_frame"]);
            CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
            
            int patchnumber = 0;
            for (int ptl = 0; ptl < static_cast<int>(current_mainpatches.length()); ptl++) {
              if (LefkoUtils::stringcompare_simple(String(patch_vec(0)),
                  String(current_mainpatches(ptl)), false)) patchnumber = ptl;
            }
            
            //Rcout << "invade3_fb_core e        ";
            
            // Not sure if we need the next bit
            DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
            DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
            DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
            CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
            CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
            CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
            
            //Rcout << "invade3_fb_core e1        ";
            // Counter resets
            int yearnumber = current_times_vec(j); // year_counter
            //Rcout << "invade3_fb_core e2        ";
            
            //Rcout << "current_mainyears: " << current_mainyears << "               ";
            //Rcout << "yearnumber: " << yearnumber << "               ";
            CharacterVector current_year = as<CharacterVector>(current_mainyears(yearnumber));
            //Rcout << "invade3_fb_core e3        ";
            
            if (inda_num_terms_counter(i, m) >= inda_terms_num_vec(0)) inda_num_terms_counter(i, m) = 0;
            if (indb_num_terms_counter(i, m) >= indb_terms_num_vec(0)) indb_num_terms_counter(i, m) = 0;
            if (indc_num_terms_counter(i, m) >= indc_terms_num_vec(0)) indc_num_terms_counter(i, m) = 0;
            if (inda_cat_terms_counter(i, m) >= inda_terms_cat_vec(0)) inda_cat_terms_counter(i, m) = 0;
            if (indb_cat_terms_counter(i, m) >= indb_terms_cat_vec(0)) indb_cat_terms_counter(i, m) = 0;
            if (indc_cat_terms_counter(i, m) >= indc_terms_cat_vec(0)) indc_cat_terms_counter(i, m) = 0;
            
            List current_ind_terms_num = ind_terms_num_list(0);
            List current_ind_terms_cat = ind_terms_cat_list(0);
            
            //Rcout << "invade3_fb_core e4        ";
            NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
            NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
            NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
            CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
            CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
            CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
            
            //Rcout << "invade3_fb_core e5        ";
            NumericVector f2_inda = {f_inda_full(inda_num_terms_counter(i, m))}; // i
            NumericVector f1_inda = {f_inda_full(inda_num_terms_previous(i, m))};
            NumericVector f2_indb = {f_indb_full(indb_num_terms_counter(i, m))};
            NumericVector f1_indb = {f_indb_full(indb_num_terms_previous(i, m))};
            NumericVector f2_indc = {f_indc_full(indc_num_terms_counter(i, m))};
            NumericVector f1_indc = {f_indc_full(indc_num_terms_previous(i, m))};
            CharacterVector r2_inda = as<CharacterVector>(r_inda_full(inda_cat_terms_counter(i, m)));
            CharacterVector r1_inda = 
              as<CharacterVector>(r_inda_full(inda_cat_terms_previous(i, m)));
            CharacterVector r2_indb = as<CharacterVector>
              (r_indb_full(indb_cat_terms_counter(i, m)));
            CharacterVector r1_indb = as<CharacterVector>
              (r_indb_full(indb_cat_terms_previous(i, m)));
            CharacterVector r2_indc = as<CharacterVector>
              (r_indc_full(indc_cat_terms_counter(i, m)));
            CharacterVector r1_indc = 
              as<CharacterVector>(r_indc_full(indc_cat_terms_previous(i, m)));
            
            // dev_terms and vrm trait axis processing
            //Rcout << "invade3_fb_core e6        ";
            NumericVector dv_terms (14);
            if (dev_terms_times_int > 0) {
              NumericMatrix used_dv_df = as<NumericMatrix>(dev_terms_list(current_variant_index));
              if (dev_num_counter(i, m) >= dev_terms_times_int) dev_num_counter(i, m) = 0;
              dv_terms = used_dv_df(_, dev_num_counter(i, m));
              
              arma::uvec var_corresponding_elems = find(variant_nta == (i + 1));
              int vce_found = static_cast<int>(var_corresponding_elems.n_elem);
              
              //Rcout << "invade3_fb_core e11        ";
              if (vce_found > 0) {
                //Rcout << "invade3_fb_core e12        ";
                arma::vec surv_dev_nta_sub = surv_dev_nta.elem(var_corresponding_elems);
                arma::vec obs_dev_nta_sub = obs_dev_nta.elem(var_corresponding_elems);
                arma::vec size_dev_nta_sub = size_dev_nta.elem(var_corresponding_elems);
                arma::vec sizeb_dev_nta_sub = sizeb_dev_nta.elem(var_corresponding_elems);
                arma::vec sizec_dev_nta_sub = sizec_dev_nta.elem(var_corresponding_elems);
                arma::vec repst_dev_nta_sub = repst_dev_nta.elem(var_corresponding_elems);
                arma::vec fec_dev_nta_sub = fec_dev_nta.elem(var_corresponding_elems);
                arma::vec jsurv_dev_nta_sub = jsurv_dev_nta.elem(var_corresponding_elems);
                arma::vec jobs_dev_nta_sub = jobs_dev_nta.elem(var_corresponding_elems);
                arma::vec jsize_dev_nta_sub = jsize_dev_nta.elem(var_corresponding_elems);
                arma::vec jsizeb_dev_nta_sub = jsizeb_dev_nta.elem(var_corresponding_elems);
                arma::vec jsizec_dev_nta_sub = jsizec_dev_nta.elem(var_corresponding_elems);
                arma::vec jrepst_dev_nta_sub = jrepst_dev_nta.elem(var_corresponding_elems);
                arma::vec jmatst_dev_nta_sub = jmatst_dev_nta.elem(var_corresponding_elems);
                
                //Rcout << "invade3_fb_core e13        ";
                for (int vce_track = 0; vce_track < vce_found; vce_track++) {
                  if(!NumericVector::is_na(surv_dev_nta_sub(vce_track))) dv_terms(0) =
                      dv_terms(0) + surv_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(obs_dev_nta_sub(vce_track))) dv_terms(1) =
                      dv_terms(1) + obs_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(size_dev_nta_sub(vce_track))) dv_terms(2) =
                      dv_terms(2) + size_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(sizeb_dev_nta_sub(vce_track))) dv_terms(3) =
                      dv_terms(3) + sizeb_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(sizec_dev_nta_sub(vce_track))) dv_terms(4) =
                      dv_terms(4) + sizec_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(repst_dev_nta_sub(vce_track))) dv_terms(5) =
                      dv_terms(5) + repst_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(fec_dev_nta_sub(vce_track))) dv_terms(6) =
                      dv_terms(6) + fec_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jsurv_dev_nta_sub(vce_track))) dv_terms(7) =
                      dv_terms(7) + jsurv_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jobs_dev_nta_sub(vce_track))) dv_terms(8) =
                      dv_terms(8) + jobs_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jsize_dev_nta_sub(vce_track))) dv_terms(9) =
                      dv_terms(9) + jsize_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jsizeb_dev_nta_sub(vce_track))) dv_terms(10) =
                      dv_terms(10) + jsizeb_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jsizec_dev_nta_sub(vce_track))) dv_terms(11) =
                      dv_terms(11) + jsizec_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jrepst_dev_nta_sub(vce_track))) dv_terms(12) =
                      dv_terms(12) + jrepst_dev_nta_sub(vce_track);
                  if(!NumericVector::is_na(jmatst_dev_nta_sub(vce_track))) dv_terms(13) =
                      dv_terms(13) + jmatst_dev_nta_sub(vce_track);
                }
              }
            dev_num_counter(i, m) = dev_num_counter(i, m) + 1;
            }
            
            NumericVector stdbv = as<NumericVector>(stageexpansion_ta_devterms_by_variant(current_variant_index));
            
            for (int z = 0; z < 14; z++) {
              dv_terms(z) = dv_terms(z) + stdbv(z);
            }
            //Rcout << "invade3_fb_core e14        ";
            
            if (err_check) {
              mdtl(i) = mdtl_1;
            }
            bool dvr_bool {false};
            
            LogicalVector dvr_yn = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            IntegerVector dvr_style (14);
            IntegerVector dvr_time_delay (14);
            NumericVector dvr_alpha (14);
            NumericVector dvr_beta (14);
            NumericVector dens_n (14);
            
            //Rcout << "invade3_fb_core f        ";
            
            if (dens_vr_yn_vec(0) > 0) {
              dvr_bool = true;
              
              DataFrame current_dvr = density_vr_list;
              LogicalVector true_dvr_yn = as<LogicalVector>(current_dvr(1));
              IntegerVector true_dvr_style = as<IntegerVector>(current_dvr(2));
              IntegerVector true_dvr_time_delay = as<IntegerVector>(current_dvr(3));
              NumericVector true_dvr_alpha = as<NumericVector>(current_dvr(4));
              NumericVector true_dvr_beta = as<NumericVector>(current_dvr(5));
              
              dvr_yn = true_dvr_yn;
              dvr_style = true_dvr_style;
              dvr_time_delay = true_dvr_time_delay;
              dvr_alpha = true_dvr_alpha;
              dvr_beta = true_dvr_beta;
              
              int used_delay = max(true_dvr_time_delay);
                
              if (j >= (used_delay - 1)) {
                if (!stages_not_equal) {
                  arma::mat di_mat = N_mpm.slice(i);
                  arma::vec delay_issue = di_mat.col(j + 1 - used_delay);
                  double delay_N_sum = arma::sum(delay_issue);
                  
                  for (int xc = 0; xc < 14; xc++) {
                    dens_n(xc) = delay_N_sum;
                  }
                }
              }
            }
            
            double maxsize {0.0};
            double maxsizeb {0.0};
            double maxsizec {0.0};
            
            if (format_int < 5) {
              DataFrame current_allstages = allstages_all; // (i)
              
              NumericVector size3 = as<NumericVector>(current_allstages["size3"]);
              NumericVector size2n = as<NumericVector>(current_allstages["size2n"]);
              NumericVector size2o = as<NumericVector>(current_allstages["size2o"]);
              NumericVector sizeb3 = as<NumericVector>(current_allstages["sizeb3"]);
              NumericVector sizeb2n = as<NumericVector>(current_allstages["sizeb2n"]);
              NumericVector sizeb2o = as<NumericVector>(current_allstages["sizeb2o"]);
              NumericVector sizec3 = as<NumericVector>(current_allstages["sizec3"]);
              NumericVector sizec2n = as<NumericVector>(current_allstages["sizec2n"]);
              NumericVector sizec2o = as<NumericVector>(current_allstages["sizec2o"]);
              
              NumericVector maxveca = {max(size3), max(size2n), max(size2o)};
              NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
              NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
              
              maxsize = max(maxveca);
              maxsizeb = max(maxvecb);
              maxsizec = max(maxvecc);
            }
            
            //Rcout << "invade3_fb_core g        ";
            
            double dens_sp {1.0};
            
            if (sp_density_num_vec(0) > 0) {
              if (sp_density_counter(i, m) >= sp_density_num_vec(0)) sp_density_counter(i, m) = 0;
              
              NumericVector current_sp_density = as<NumericVector>(sp_density_list(0));
              dens_sp = current_sp_density(sp_density_counter(i, m));
              
              sp_density_counter(i, m) = sp_density_counter(i, m) + 1;
            }
            
            List current_mpm;
            if (format_int < 5) {
              //Rcout << "invade3_fb_core g1        ";
              current_mpm = AdaptMats::mazurekd(current_mpm_allstages,
                current_stageframe, format_int, surv_proxy, obs_proxy,
                size_proxy, sizeb_proxy, sizec_proxy, repst_proxy, fec_proxy,
                jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
                jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda, f1_inda,
                f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                r1_indb, r2_indc, r1_indc, dv_terms, dvr_bool, dvr_yn,
                dvr_style, dvr_alpha, dvr_beta, dens_n, dens_sp, fecmod_vec(0),
                maxsize, maxsizeb, maxsizec, firstage_int, finalage_int, true,
                yearnumber, patchnumber, exp_tol, theta_tol, true, err_check,
                sparse_bool, A_only);
              
              if (err_check_extreme) {
                NumericMatrix mpm_out = as<NumericMatrix>(current_mpm["out"]);
                errcheck_mpmout_reps_time_vmt_var(m) = mpm_out; 
              }
                
              //Rcout << "invade3_fb_core g2        ";
            } else {
              IntegerVector all_ages = seq(firstage_int, finalage_int);
              //DataFrame current_supplement;
              if (!(current_supplement.length() > 1)) {
                //Rcout << "invade3_fb_core g3        ";
                current_mpm = AdaptMats::mdabrowskiego(all_ages,
                  current_stageframe, surv_proxy, fec_proxy, f2_inda, f1_inda,
                  f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                  r1_indb, r2_indc, r1_indc, dv_terms(0), dv_terms(6), dens_sp,
                  fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                  dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                  exp_tol, theta_tol, sparse_bool);
                //Rcout << "invade3_fb_core g4        ";
                
              } else {
                //current_supplement = as<DataFrame>(supplement_list(0));
                
                //Rcout << "invade3_fb_core g5        ";
                current_mpm = AdaptMats::mdabrowskiego(all_ages,
                  current_stageframe, surv_proxy, fec_proxy, f2_inda, f1_inda,
                  f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                  r1_indb, r2_indc, r1_indc, dv_terms(0), dv_terms(6), dens_sp,
                  fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                  dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                  exp_tol, theta_tol, sparse_bool, current_supplement);
                //Rcout << "invade3_fb_core g6        ";
              }
            }
            //Rcout << "invade3_fb_core h        ";
            
            if (!dens_yn_bool) {
              
              if (A_only) {
                if (!sparse_bool) {
                  arma::mat current_A = as<arma::mat>(current_mpm["A"]);
                  Amat_alter(current_A, sge_current); 
                  if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; 
                  
                  running_popvec_vrm = current_A * running_popvec_vrm; 
                } else {
                  arma::sp_mat current_A = as<arma::sp_mat>(current_mpm["A"]);
                  AdaptUtils::sp_Amat_alter(current_A, sge_current);
                  if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; 
                  
                  running_popvec_vrm = current_A * running_popvec_vrm;
                }
              } else {
                if (!sparse_bool) {
                  arma::mat current_A = as<arma::mat>(current_mpm["A"]);
                  arma::mat current_U_unaltered = as<arma::mat>(current_mpm["U"]);
                  arma::mat current_F_unaltered = as<arma::mat>(current_mpm["F"]);
                  AdaptUtils::UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
                  if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
                  
                  running_popvec_vrm = current_A * running_popvec_vrm;
                } else {
                  arma::sp_mat current_A = as<arma::sp_mat>(current_mpm("A"));
                  arma::sp_mat current_U_unaltered = as<arma::sp_mat>(current_mpm["U"]);
                  arma::sp_mat current_F_unaltered = as<arma::sp_mat>(current_mpm["F"]);
                  AdaptUtils::sp_UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
                  if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
                  
                  running_popvec_vrm = current_A * running_popvec_vrm;
                }
              }
            } else {
              ///// dens_bool = true
              DataFrame used_density_input = density_df;
              DataFrame used_density_index_input = dens_index_df;
              
              IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
              int used_delay = max(ud_delay_vec);
              
              if (j >= (used_delay - 1)) {
                if (!stages_not_equal) {
                  arma::mat di_mat = N_mpm.slice(i);
                  arma::vec delay_issue = di_mat.col(j + 1 - used_delay);
                  double delay_N_sum = arma::sum(delay_issue);
                  
                  arma::vec new_popvec;
                  arma::mat new_projmat;
                  arma::sp_mat new_projsp;
                  
                  AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                    running_popvec_vrm, sge_current, current_mpm, delay_N_sum,
                    0, integeronly, substoch, used_density_input,
                    used_density_index_input, false, sparse_bool, sparse_bool,
                    false, err_check);
                  
                  running_popvec_vrm = new_popvec;
                  if (err_check_extreme) { // Could remove later
                    if (!sparse_bool) {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                    } else {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                    }
                  } // Could remove later
                } else {
                  double delay_N_sum {0.0};
                  
                  if (j > 0) {
                    for (int l = 0; l < var_per_run; l++) {
                      int current_variant_index_agg = var_run_mat(i, l);
                      List current_pop_list = as<List>(comm_out_pre(i));
                      List pop_rep_list = as<List>(current_pop_list(l));
                      arma::mat delay_pop = as<arma::mat>(current_pop_list(current_rep));
                      
                      arma::vec delay_pop_vec = delay_pop.col(j + 1 - used_delay);
                      arma::vec current_equiv_vec = as<arma::vec>(equivalence_list(current_variant_index_agg));
                      arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                      double delay_pop_N = arma::accu(adjusted_delay_pop_vec);
                      
                      delay_N_sum += delay_pop_N;
                    }
                  }
                  
                  arma::vec new_popvec;
                  arma::mat new_projmat;
                  arma::sp_mat new_projsp;
                  
                  AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                    running_popvec_vrm, sge_current, current_mpm, delay_N_sum,
                    0, integeronly, substoch, used_density_input,
                    used_density_index_input, false, sparse_bool, sparse_bool,
                    false, err_check);
                  
                  running_popvec_vrm = new_popvec;
                  if (err_check_extreme) {
                    if (!sparse_bool) {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                    } else {
                      errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                    }
                  }
                }
              } else {
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                  running_popvec_vrm, sge_current, current_mpm, 0.0, 0,
                  integeronly, substoch, used_density_input,
                  used_density_index_input, false, sparse_bool, sparse_bool,
                  false, err_check);
                
                running_popvec_vrm = new_popvec;
                if (err_check_extreme) {
                  if (!sparse_bool) {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                  } else {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                  }
                }
              }
            }
            
            //Rcout << "invade3_fb_core i        ";
            
            if (integeronly) running_popvec_vrm = floor(running_popvec_vrm);
            double N_current = arma::sum(running_popvec_vrm);
            N_mpm(m, (j + 1), i) = N_current;
            
            inda_num_terms_previous(i, m) = static_cast<int>(inda_num_terms_counter(i, m));
            indb_num_terms_previous(i, m) = static_cast<int>(indb_num_terms_counter(i, m));
            indc_num_terms_previous(i, m) = static_cast<int>(indc_num_terms_counter(i, m));
            inda_cat_terms_previous(i, m) = static_cast<int>(inda_cat_terms_counter(i, m));
            indb_cat_terms_previous(i, m) = static_cast<int>(indb_cat_terms_counter(i, m));
            indc_cat_terms_previous(i, m) = static_cast<int>(indc_cat_terms_counter(i, m));
            
            inda_num_terms_counter(i, m) = inda_num_terms_counter(i, m) + 1;
            indb_num_terms_counter(i, m) = indb_num_terms_counter(i, m) + 1;
            indc_num_terms_counter(i, m) = indc_num_terms_counter(i, m) + 1;
            inda_cat_terms_counter(i, m) = inda_cat_terms_counter(i, m) + 1;
            indb_cat_terms_counter(i, m) = indb_cat_terms_counter(i, m) + 1;
            indc_cat_terms_counter(i, m) = indc_cat_terms_counter(i, m) + 1;
            
            running_popvecs(m) = running_popvec_vrm;
            pops_out.col(j + 1) = running_popvec_vrm;
            
          } // if (j > (entry_time_vec(i) - 1))
          
          pop_reps(current_rep) = pops_out;
          // comm_out_pre(i) = reps_out;
        } // vrm loop
        if (err_check_extreme) {
          errcheck_mpmout_reps_time_vmt(i) = errcheck_mpmout_reps_time_vmt_var;
          errcheck_mpm_reps_time_vmt(i) = errcheck_mpm_reps_time_vmt_var;
        }
      }
      
      if (err_check_extreme) {
        errcheck_mpmout_reps_time(j) = errcheck_mpmout_reps_time_vmt;
        errcheck_mpm_reps_time(j) = errcheck_mpm_reps_time_vmt;
      }
      
      year_counter++;
    } // j loop
    if (err_check_extreme) {
      errcheck_mpmout_reps(current_rep) = errcheck_mpmout_reps_time;
      errcheck_mpm_reps(current_rep) = errcheck_mpm_reps_time;
    }
    N_out_pre(current_rep) = N_mpm;
    
  } // current_rep loop
  comm_out = comm_out_pre;
  N_out = N_out_pre;
  if (err_check) {
    errcheck_mpmouts = errcheck_mpmout_reps;
    errcheck_mpms = errcheck_mpm_reps;
    modified_dev_terms_list = mdtl;
  }
  
  //Rcout << "invade3_fb_core k        ";
  
  IntegerVector Lyap_rows = seq(1, nreps * var_mat_length);
  CharacterVector Lyap_df_names ((3 * var_per_run) + 2);
  Lyap_df_names(0) = "simulation_num";
  Lyap_df_names(1) = "rep";
  
  String Lyap_df_names_base1 = "variant";
  String Lyap_df_names_base2 = "entrytime";
  String Lyap_df_names_base3 = "fitness_variant";
  
  List output ((3 * var_per_run) + 2);
  output(0) = Lyap_rows;
  
  //Rcout << "invade3_fb_core l        ";
  
  int Lyap_counter {0};
  for (int m = 0; m < var_per_run; m++) {
    Lyap_counter = 0;
    
    IntegerVector Lyap_var_orig = as<IntegerVector>(Lyapunov(1+m));
    IntegerVector Lyap_rep (nreps * var_mat_length);
    IntegerVector Lyap_var (nreps * var_mat_length);
    IntegerVector Lyap_etime (nreps * var_mat_length);
    NumericVector Lyap_fitness (nreps * var_mat_length);
    
    for (int i = 0; i < nreps; i++) {
      arma::cube current_N_out_cube = as<arma::cube>(N_out(i));
      
      for (int j = 0; j < var_mat_length; j++) {
        arma::mat current_N_out = current_N_out_cube.slice(j);
        int time_length = static_cast<int>(current_N_out.n_cols);
        int start_time = time_length - fitness_times;
        
        arma::vec running_fitness(fitness_times, fill::zeros);
        arma::vec generations(fitness_times, fill::zeros);
        arma::vec intercept_ones(fitness_times, fill::ones);
        
        for (int k = 0; k < fitness_times; k++) {
          running_fitness(k) = current_N_out(m, k+start_time);
          generations(k) = static_cast<double>(k);
        }
        
        arma::mat xmat = join_rows(intercept_ones, generations);
        arma::vec Lyap_regr;
        
        AdaptUtils::fastLm_sl(Lyap_regr, running_fitness, xmat);
        
        double Lyapunov_estimate = Lyap_regr(1);
        
        Lyap_var(Lyap_counter) = Lyap_var_orig(j);
        Lyap_rep(Lyap_counter) = i+1;
        Lyap_etime(Lyap_counter) = entry_time_vec(m);
        Lyap_fitness(Lyap_counter) = Lyapunov_estimate;
        
        Lyap_counter++;
      }
    }
    output(1) = Lyap_rep;
    output(2 + m) = Lyap_var;
    output(2 + var_per_run + m) = Lyap_etime;
    output(2 + (2 * var_per_run) + m) = Lyap_fitness;
    
    String new_col_name = Lyap_df_names_base1;
    new_col_name += (m + 1);
    Lyap_df_names(2 + m) = new_col_name;
    
    String next_col_name = Lyap_df_names_base2;
    next_col_name += (m + 1);
    Lyap_df_names(2 + var_per_run + m) = next_col_name;
    
    String last_col_name = Lyap_df_names_base3;
    last_col_name += (m + 1);
    Lyap_df_names(2 + (2 * var_per_run) + m) = last_col_name;
  }
  
  output.attr("names") = Lyap_df_names;
  output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
    nreps * var_mat_length);
  StringVector needed_class {"data.frame"};
  output.attr("class") = needed_class;
    
  Lyapunov = output;
}

//' Run Pairwise and Multiple Invasion Analysis
//' 
//' Function \code{invade3} runs pairwise and multiple invasion analyses.
//' 
//' @name invade3
//' 
//' @param axis The \code{adaptAxis} object detailing all variant
//' characteristics. Essentially, a data frame giving the values of all changes
//' to vital rates and transition elements to test, where each value is change
//' is given by row.
//' @param mpm An MPM of class \code{lefkoMat}, for use if using existing MPMs.
//' @param vrm A \code{vrm_input} object corresponding to a distinct MPM that
//' will be created during analysis. Requires a stageframe, entered in argument
//' \code{stageframe}.
//' @param stageframe A stageframe defining stages and the life cycle for the
//' entered object in argument \code{vrms}. Must be of class \code{stageframe}.
//' @param supplement An optional data frame of class \code{lefkoSD} providing
//' supplemental data that should be incorporated into function-based MPMs. See
//' \code{\link[lefko3]{supplemental}()} for details. Use only with argument
//' \code{vrm}.
//' @param equivalence An optional object of class \code{adaptEq} giving the
//' degree to which individuals in each stage are equivalent to one another.
//' May also be a numeric vector, in which case the vector must have the same
//' number of elements as the number of rows in the associated MPM, with each
//' element giving the effect of an individual of that age, stage, age-stage, or
//' stage-pair, depending on whether the MPM is age-based, ahistorical
//' stage-based, age-by-stage, or historical stage-based, respectively. Numeric
//' entries used in these vectors can be thought of as Lotka-Volterra
//' interaction terms, such as are used in multiple species competition models.
//' @param starts An optional \code{lefkoSV} object, which is a data frame
//' providing the starting numbers of individuals of each stage. If not
//' provided, then all projections start with a single individual per stage.
//' @param years An optional term corresponding to a single integer vector of
//' time \code{t} values. If a vector shorter than \code{times} is supplied,
//' then this vector will be cycled. Defaults to a vector of all detected
//' years in argument \code{mpm} or argument \code{vrm}.
//' @param patches An optional single string giving a single pop-patch to be
//' used during invasion analysis. Defaults to the population-level set or the
//' first patch, depending on whether the former exists.
//' @param tweights An optional numeric vector or matrice denoting the
//' probabilities of choosing each matrix in each MPM in a stochastic
//' projection. If a matrix, then a first-order Markovian environment is
//' assumed, in which the probability of choosing a specific annual matrix
//' depends on which annual matrix is currently chosen. If an element of the
//' list is a vector, then the choice of annual matrix is assumed to be
//' independent of the current matrix. Defaults to equal weighting among
//' matrices.
//' @param format An optional integer indicating the kind of function-based MPM
//' to create, if argument \code{vrm} is provided. Possible choices include:
//' \code{1}, Ehrlen-format historical MPM; \code{2}, deVries-format historical
//' MPM; \code{3}, ahistorical MPM (default); \code{4}, age-by-stage MPM; and
//' \code{5}, Leslie (age-based) MPM. Defaults to \code{3}.
//' @param entry_time An optional integer vector giving the entry time for each
//' variant into each simulation. Defaults to a zero vector with length equal to
//' the number of variants to run concurrently in each simulation, as given by
//' argument \code{var_per_run}.
//' @param sp_density An optional argument for use with argument \code{vrm} that
//' specifies the spatial density to be used in each time step. If used, then
//' may either be a numeric vector giving a single spatial density for each
//' time step. If vectors are shorter than specified in \code{times}, then these
//' values will be cycled.
//' @param ind_terms An optional argument providing values of individual or
//' environmental covariate values for argument \code{vrm}. Should be set to a
//' single data frame with 3 columns giving values for up to 3 covariates across
//' time (rows give the time order of these values). Unused terms within the
//' data frame must be set to \code{0} (use of \code{NA} will produce errors).
//' If the number of rows is less than \code{times}, then these values will be
//' cycled.
//' @param dev_terms An optional  data frame including 14 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival, 9: juvenile
//' observation, 10: juvenile primary size, 11: juvenile secondary size,
//' 12: juvenile tertiary size, 13: juvenile reproduction, and 14: juvenile
//' maturity transition. Unused terms must be set to \code{0} (use of \code{NA}
//' will produce errors). Single or small numbers of values per vital rate model
//' are also allowed, and if the number of rows is less than \code{times}, then
//' the terms will be cycled.
//' @param fb_sparse A logical value indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE}.
//' @param firstage An optional integer used for function-based Leslie and
//' age-by-stage MPMs giving the starting age in such MPMs. Use only if the MPM
//' is both function-based and has age structure. Typically, the starting age in
//' such MPMs should be set to \code{0} if post-breeding and \code{1} if
//' pre-breeding. All other MPMs should be set to \code{0}.
//' @param finalage An optional integer used for function-based Leslie and
//' age-by-stage MPMs giving the final age in such MPMs. Use only if the MPM is
//' both function-based and has age structure.
//' @param fecage_min An optional integer used for function-based Leslie MPMs
//' giving the first age at which organisms can reproduce. Use only if the MPM
//' is both function-based and has age structure. Defaults to the value given in
//' \code{firstage}.
//' @param fecage_max An optional integer used for function-based Leslie MPMs
//' giving the final age at which organisms can reproduce. Use only if the MPM
//' is both function-based and has age structure. Defaults to the value given in
//' \code{finalage}.
//' @param cont An optional logical value for function-based Leslie and
//' age-by-stage MPMs stating whether the MPM should should include a stasis
//' transition within the final age. This should be used only when an organism
//' can maintain the demographic characteristics of the final described age
//' after reaching that age.
//' @param prebreeding An optional logical value indicating whether the life
//' cycle is prebreeding (\code{TRUE}) or postbreeding (\code{FALSE}). Defaults 
//' to \code{TRUE}.
//' @param fecmod An optional numeric value for function-based MPMs giving
//' scalar multipliers for fecundity terms, when two fecundity variables are
//' used for a collective fecundity per individual.
//' @param density An optional data frames of class \code{lefkoDens}, which
//' provides details for density dependence in MPM elements and is created with
//' function \code{\link[lefko3]{density_input}()}. Defaults to \code{NULL}, in
//' which case no density dependence is built into matrix elements.
//' @param density_vr An optional data frame of class \code{lefkoDensVR}, which
//' provides details for density dependence in vital rate models and has been
//' created with function \code{link[lefko3]{density_vr}()}. Can only be used
//' with function-based projections. Defaults to \code{NULL}, in which case no
//' density dependence is built into vital rates.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking. Can also be set to the text value
//' \code{"extreme"}, in which case all \code{err_check} output plus a multiple
//' level list with each MPM used in each time step will be output.
//' @param var_per_run The number of variants to run in each simulation.
//' Defaults to \code{2}, resulting in pairwise invasibility analysis. See
//' \code{Notes} for details.
//' @param stochastic A logical value indicating whether the projection will be
//' run as a temporally stochastic projection. Defaults to \code{FALSE}.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion down to the next lower
//' integer. Defaults to \code{FALSE}.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent and density
//' independent simulations. Defaults to \code{0}, which does not enforce
//' substochasticity. Alternatively, \code{1} forces all survival-transition
//' elements to range from 0.0 to 1.0, and forces fecundity to be non-negative;
//' and \code{2} forces all column rows in the survival-transition matrices to
//' total no more than 1.0, in addition to the actions outlined for option
//' \code{1}. Both settings \code{1} and \code{2} change negative fecundity
//' elements to \code{0.0}, while setting \code{0} does not alter fecundity.
//' @param nreps The number of replicate projections. Defaults to \code{1}.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' \code{10000}.
//' @param fitness_times An integer giving the number of time steps at the end
//' of each run to use to estimate the fitness of the respective genotype.
//' Defaults to \code{100}, but if \code{times < 100}, then is set equal to
//' \code{times}.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param theta_tol A numeric value used to indicate a maximum value to theta as
//' used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' @param A_only A logical value indicating whether to alter survival and
//' fecundity matrix elements separately prior to creating the overall \code{A}
//' matrix, or whether to alter elements directly on \code{A} matrices. Defaults
//' to \code{TRUE}, and should be kept to that setting unless some matrix
//' elements to be altered are sums of survival and fecundity transitions.
//' 
//' @return A list of class \code{adaptInv}, with the following elements:
//' \item{fitness}{A data frame giving the Lyapunov coefficients estimated for
//' each variant, per replicate.}
//' \item{variants_out}{A two-level list with the top level list having number of
//' elements equal to the number of variants, and the lower level
//' corresponding to the number of replicates. Each element of the lower level
//' list is a matrix showing the number of individuals in each stage (row) at each
//' time (column).}
//' \item{N_out}{A list with the number of elements equal to the number of
//' replicates. Each element within this list is data frame showing the number
//' of individuals of each species or genotype alive at each time. The number of
//' rows are equal to the number of MPMs used, and the columns correspond to the
//' time steps.}
//' \item{stageframe_list}{A list in which each element is the stageframe for
//' each MPM used.}
//' \item{hstages_list}{A list giving the used \code{hstages} data frames, which
//' identify the correct stage pairing for each row / column in each
//' historical MPM utilized.}
//' \item{agestages_list}{A list giving the used \code{agestages} data frames,
//' which identify the correct age-stage pairing for each row / column in each
//' age-by-stage MPM utilized.}
//' \item{labels}{A small data frame giving the the population and patch
//' identities for each MPM entered.}
//' \item{err_check}{An optional list composed of an additional six lists, each
//' of which has the number of elements equal to the number of MPMs utilized.
//' List output include \code{allstages_all}, which gives the indices of
//' estimated transitions in MPMs constructed by function \code{invade3()} from
//' input vital rate models; \code{allmodels_all}, which provides all vital rate
//' models as decomposed and interpreted by function \code{invade3()};
//' \code{equivalence_list}, which gives the stage equivalence for density
//' calculations across MPMs; \code{density_list}, which gives the
//' \code{density} inputs utilized; \code{dens_index_list}, which provides
//' indices used to identify matrix elements for density dependence; and
//' \code{density_vr_list}, which gives the \code{density_vr} inputs utilized.}
//' 
//' @section Notes:
//' 
//' The argument \code{var_per_run} establishes the style of simulation to run.
//' Entering \code{var_per_run = 1} runs each variant singly. Entering
//' \code{var_per_run = 2} runs pairwise invason analysis, trying each pair
//' permutation of variants. Greater values will lead to multiple invasion
//' analysis with different permutations of groups. For example,
//' \code{var_per_run = 3} runs each permutation of groups of three. The integer
//' set must be positive, and must not be larger than the number of variants.
//' 
//' @examples
//' library(lefko3)
//' data(cypdata)
//' 
//' sizevector <- c(0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg")
//' repvector <- c(0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.40, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, 1000, 1000),
//'   type =c(1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' 
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' cypmean <- lmean(cypmatrix2r)
//' 
//' cyp_start <- start_input(cypmean, stage2 = c("SD", "P1", "D"),
//'   value = c(1000, 200, 4))
//' 
//' c2d_4 <- density_input(cypmean, stage3 = c("P1", "P1"), stage2= c("SD", "rep"),
//'   style = 2, time_delay = 1, alpha = 0.005, beta = 0.000005, type = c(2, 2))
//' 
//' # A simple projection allows us to find a combination of density dependence
//' # and running time that produces a stable quasi-equilibrium
//' cyp_proj <- projection3(cypmean, times = 250, start_frame = cyp_start,
//'   density = c2d_4, integeronly = TRUE)
//' plot(cyp_proj)
//' 
//' cyp_ta <- trait_axis(stageframe = cypframe_raw,
//'   stage3 = rep("P1", 15),
//'   stage2 = rep("rep", 15),
//'   multiplier = seq(from = 0.1, to = 10.0, length.out = 15),
//'   type = rep(2, 15))
//' 
//' cyp_inv <- invade3(axis = cyp_ta, mpm = cypmean, density = c2d_4, times = 350,
//'   starts = cyp_start, entry_time = c(0, 250), fitness_times = 30,
//'   var_per_run = 2)
//' plot(cyp_inv)
//' 
//' @export invade3
// [[Rcpp::export(invade3)]]
List invade3 (Nullable<RObject> axis = R_NilValue, Nullable<RObject> mpm  = R_NilValue,
  Nullable<RObject> vrm = R_NilValue, Nullable<RObject> stageframe  = R_NilValue,
  Nullable<RObject> supplement = R_NilValue, Nullable<RObject> equivalence = R_NilValue,
  Nullable<RObject> starts = R_NilValue, Nullable<RObject> years = R_NilValue,
  Nullable<RObject> patches = R_NilValue, Nullable<RObject> tweights = R_NilValue,
  
  Nullable<RObject> format = R_NilValue, Nullable<RObject> entry_time = R_NilValue,
  Nullable<RObject> sp_density = R_NilValue, Nullable<RObject> ind_terms = R_NilValue,
  Nullable<RObject> dev_terms = R_NilValue, Nullable<RObject> fb_sparse = R_NilValue,
  
  Nullable<RObject> firstage = R_NilValue, Nullable<RObject> finalage = R_NilValue,
  Nullable<RObject> fecage_min = R_NilValue, Nullable<RObject> fecage_max = R_NilValue,
  Nullable<RObject> cont = R_NilValue, Nullable<RObject> prebreeding = R_NilValue,
  Nullable<RObject> fecmod = R_NilValue,
  
  Nullable<RObject> density = R_NilValue, Nullable<RObject> density_vr = R_NilValue,
  Nullable<RObject> err_check = R_NilValue,
  
  int var_per_run = 2,
  
  bool stochastic = false, bool integeronly = false, int substoch = 0,
  int nreps = 1, int times = 10000, int fitness_times = 100,
  double exp_tol = 700.0, double theta_tol = 100000000.0, bool A_only = true) {
  
  int mpm_count {0};
  int vrm_count {0};
  unsigned int variant_count {0};
  int total_mpms {0};  // This includes all MPMs and VRMs
  int stagecounts; // # stages in each MPM
  int stageframe_count {0};
  int stageframe_notNull_count {0};
  int supplement_count {0};
  int equivalence_count {0};
  int start_count {0};
  int tweights_count {0};
  int density_count {0};
  int dev_terms_times_int {0}; // number rows, corresponding to dev cycles
  int entry_time_count {0};
  int density_vr_count {0};
  int sparse_vec_count {0};
  int format_int; // MPM format (1:Ehrlen; 2:deVries; 3:ahist; 4:age-stage; 5: Leslie)
  int preexisting_mpm_size {0};
  
  bool preexisting {false}; // Are preexisting MPMs being used?
  bool funcbased {false}; // Will function-based MPMs be created?
  bool entry_time_vec_use {false}; // Are any elements in entry_time greater than 0?
  bool stages_not_equal {false}; // Are equivalence vectors supplied separating even stages?
  bool pure_leslie {false}; // Is core MPM Leslie?
  bool dens_yn_bool {false}; // density input for each MPM (0 = no, 1 = yes)
  bool prebreeding_bool {true};
  
  bool inda_char {false};
  bool indb_char {false};
  bool indc_char {false};
  bool historical {false};
  bool sparse_bool {false}; // whether MPMs are or should be in sparse format
  bool err_check_bool {false};
  bool err_check_extreme {false};
  
  // err_check processing
  LefkoInputs::RObj_TF_input_check("err_check", "extreme", err_check_bool,
    err_check_extreme, true, true, err_check);
  
  // Age-by-stage and Leslie MPM settings
  int firstage_int;
  int finalage_int;
  IntegerVector cont_vec;
  
  // Leslie MPM only
  IntegerVector fecage_min_vec;
  IntegerVector fecage_max_vec;
  
  // Main lists
  List mpm_list;
  List A_list;
  List vrm_list;
  DataFrame final_stageframe;
  DataFrame stageframe_df; // List ending in _fb are only used in function-based cases
  DataFrame supplement_df;
  DataFrame supplement_list_fb;
  List repmatrix_list;
  List equivalence_list;
  List hstages_list;
  List agestages_list;
  List start_list;
  List tweights_list;
  DataFrame density_df;
  DataFrame dens_index_df;  // Holds element index vectors for density_frames in density_df
  List sp_density_list;
  List ind_terms_num_list;
  List ind_terms_cat_list;
  List dev_terms_list;
  List modified_dev_terms_list; // Only used in fb processing
  List allstages_all; // Used in fbMPM processing
  List allmodels_all; // Used in fbMPM processirg
  List stageexpansion_list; // Get rid of this eventually
  List stageexpansion_list_vrm; // Get rid of this eventually, VRM only
  List errcheck_mpm_list;
  List errcheck_mpmout_list;
  
  DataFrame density_vr_list;
  
  IntegerVector matrowcounts; // # rows in each MPM
  NumericVector equivalence_vec; // equivalence vector if !stages_not_equal
  CharacterVector patch_vec; // choice of patch in MPM
  CharacterVector year_vec; // choice of years in MPM
  IntegerVector tweights_type_vec; // tweights input as vector (1) or matrix (2) or null (0)
  IntegerVector total_years_vec; // total # years in each MPM
  IntegerVector entry_time_vec; // times of entry for each MPM
  IntegerVector dens_vr_yn_vec; // density_vr input for each MPM (0 = no, 1 = yes)
  IntegerVector sp_density_num_vec; // # of spatial density terms per MPM
  IntegerVector inda_terms_num_vec; // # of indcova (double) times per MPM
  IntegerVector indb_terms_num_vec; // # of indcovb (double) times per MPM
  IntegerVector indc_terms_num_vec; // # of indcovc (double) times per MPM
  IntegerVector inda_terms_cat_vec; // # of indcova (cat) times per MPM
  IntegerVector indb_terms_cat_vec; // # of indcovb (cat) times per MPM
  IntegerVector indc_terms_cat_vec; // # of indcovc (cat) times per MPM
  NumericVector fecmod_vec; // fecundity multipliers for multiple offspring stages
  
  IntegerVector axis_variants_unique; // Vector giving the variant numbers
  arma::mat var_run_mat; // Holds indices for variants to run per run in invasion analysis
  List zero_stage_vec_list; // Holds zero stage vectors, if entry times are staggered
  
  DataFrame trait_axis; // DataFrame used to house trait axis
  DataFrame new_trait_axis; // Data frame housing re-assessed trait_axis
  DataFrame Lyapunov; // DataFrame showing fitness of all genotypes under all circumstances
  List comm_out; // List of mats of pop vecs (top lvl: mpm, lower lvl: reps, mats: stages x times)
  List N_out;  // List of pop size mats (top level: reps, mats: mpm by times)
  DataFrame labels;  // Data frame to provide order of MPMs
  List labels_list;  // List to hold data for data frame labels
  
  //Rcout << "invade3 a ";
  
  if (fitness_times > times) {
    Rf_warningcall(R_NilValue,
      "Argument fitness_times is too large. Reseting to value of argument times.");
    
    fitness_times = times;
  }
  
  //Rcout << "invade3 b ";
  
  if (axis.isNotNull()) {
    if (is<DataFrame>(axis)) {
      trait_axis = as<DataFrame>(axis);
      CharacterVector ta_vars = trait_axis.names();
      int ta_vars_num = static_cast<int>(ta_vars.length());
      
      if (ta_vars_num != 30) {
        throw Rcpp::exception("Argument axis is not recognized.", false);
      }
      
      //bool found_variant {false};
      //for (int i = 0; i < ta_vars_num; i++) {
      //  if (ta_vars(i) == "variant") found_variant = true;
      //}
      
      IntegerVector axis_variants = as<IntegerVector>(trait_axis["variant"]);
      axis_variants_unique = sort_unique(axis_variants);
      variant_count = static_cast<unsigned int>(axis_variants_unique.length());
      
      if (!(variant_count >= var_per_run)) {
        throw Rcpp::exception("Argument var_per_run may not be greater than the number of variants.", false);
      }
      if (var_per_run < 1) {
        throw Rcpp::exception("Argument var_per_run must be a positive integer.", false);
      }
      double exp_permutes = pow(variant_count, var_per_run);
      
      var_run_mat = AdaptUtils::exp_grid_single(variant_count, var_per_run); // Rows = runs, Cols = variants
      
    } else {
      AdaptUtils::pop_error2("axis", "a data frame of class adaptAxis", "", 1);
    }
  } else {
    throw Rcpp::exception("Argument axis is required to run invasion analyses.", false);
  }
  
  AdaptUtils::Lyapunov_df_maker(Lyapunov, var_run_mat, var_per_run, axis_variants_unique);
  
  //DataFrame final_Lyapunov = clone(Lyapunov); // Remove later
  
  List cleaned_input = cleanup3_inv(mpm, vrm, stageframe, supplement, format,
    firstage, finalage, fecage_min, fecage_max, cont, fecmod, starts, patches,
    years, tweights, density, entry_time, density_vr, sp_density, ind_terms,
    dev_terms, fb_sparse, equivalence, prebreeding, exp_tol, theta_tol,
    substoch, variant_count, var_per_run);
  
  //Rcout << "invade3 c ";
  
  mpm_list = as<List>(cleaned_input(0));
  mpm_count = static_cast<int>(cleaned_input(1));
  vrm_list = as<List>(cleaned_input(2));
  vrm_count = static_cast<int>(cleaned_input(3));
  final_stageframe = as<DataFrame>(cleaned_input(4));
  stageframe_df = as<DataFrame>(cleaned_input(5));
  supplement_df = as<List>(cleaned_input(6));
  supplement_list_fb = as<List>(cleaned_input(7));
  repmatrix_list = as<List>(cleaned_input(8));
  sparse_bool = static_cast<bool>(cleaned_input(9));
  sparse_vec_count = static_cast<int>(cleaned_input(10));
  format_int = static_cast<int>(cleaned_input(11));
  pure_leslie = static_cast<bool>(cleaned_input(12));
  stageframe_notNull_count = static_cast<int>(cleaned_input(13));
  preexisting = static_cast<bool>(cleaned_input(14));
  funcbased = static_cast<bool>(cleaned_input(15));
  firstage_int = static_cast<int>(cleaned_input(16));
  finalage_int = static_cast<int>(cleaned_input(17));
  cont_vec = as<IntegerVector>(cleaned_input(18));
  fecmod_vec = as<NumericVector>(cleaned_input(19));
  fecage_min_vec = as<IntegerVector>(cleaned_input(20));
  fecage_max_vec = as<IntegerVector>(cleaned_input(21));
  hstages_list = as<List>(cleaned_input(22));
  agestages_list = as<List>(cleaned_input(23));
  matrowcounts = as<IntegerVector>(cleaned_input(24));
  stagecounts = static_cast<int>(cleaned_input(25));
  start_list = as<List>(cleaned_input(26));
  start_count = static_cast<int>(cleaned_input(27));
  labels_list = as<List>(cleaned_input(28));
  labels = as<DataFrame>(cleaned_input(29));
  patch_vec = as<CharacterVector>(cleaned_input(30));
  year_vec = as<CharacterVector>(cleaned_input(31));
  total_years_vec = as<IntegerVector>(cleaned_input(32));
  tweights_list = as<List>(cleaned_input(33));
  tweights_count = static_cast<int>(cleaned_input(34));
  tweights_type_vec = as<IntegerVector>(cleaned_input(35));
  density_df = as<DataFrame>(cleaned_input(36));
  dens_index_df = as<DataFrame>(cleaned_input(37));
  dens_yn_bool = static_cast<bool>(cleaned_input(38));
  density_count = static_cast<int>(cleaned_input(39));
  entry_time_vec = as<IntegerVector>(cleaned_input(40));
  entry_time_count = static_cast<int>(cleaned_input(41));
  entry_time_vec_use = static_cast<bool>(cleaned_input(42));
  density_vr_list = as<DataFrame>(cleaned_input(43));
  ind_terms_num_list = as<List>(cleaned_input(44));
  ind_terms_cat_list = as<List>(cleaned_input(45));
  dev_terms_list = as<List>(cleaned_input(46));
  dens_vr_yn_vec = as<IntegerVector>(cleaned_input(47));
  sp_density_num_vec = as<IntegerVector>(cleaned_input(48));
  dev_terms_times_int = static_cast<int>(cleaned_input(49));
  inda_terms_num_vec = as<IntegerVector>(cleaned_input(50));
  indb_terms_num_vec = as<IntegerVector>(cleaned_input(51));
  indc_terms_num_vec = as<IntegerVector>(cleaned_input(52));
  inda_terms_cat_vec = as<IntegerVector>(cleaned_input(53));
  indb_terms_cat_vec = as<IntegerVector>(cleaned_input(54));
  indc_terms_cat_vec = as<IntegerVector>(cleaned_input(55)); 
  density_vr_count = static_cast<int>(cleaned_input(56));
  sparse_vec_count = static_cast<int>(cleaned_input(57));
  sp_density_list = as<List>(cleaned_input(58));
  equivalence_list = as<List>(cleaned_input(59));
  equivalence_vec = as<NumericVector>(cleaned_input(60));
  equivalence_count = static_cast<int>(cleaned_input(61));
  stages_not_equal = static_cast<bool>(cleaned_input(62));
  allstages_all = as<List>(cleaned_input(63));
  allmodels_all = as<List>(cleaned_input(64));
  historical = static_cast<bool>(cleaned_input(65));
  preexisting_mpm_size = static_cast<int>(cleaned_input(66));
  prebreeding_bool = static_cast<bool>(cleaned_input(67));
  
  total_mpms = mpm_count + vrm_count;
  
  if (total_mpms != 1) {
    throw Rcpp::exception("Function invade3 only allows a single lefkoMat or vrm_input object.",
      false);
  }
  
  //Rcout << "invade3 d ";
  
  // Projection runs
  if (preexisting) {
    //Rcout << "invade3 e ";
    invade3_pre_core (Lyapunov, var_run_mat, N_out, comm_out, zero_stage_vec_list,
      trait_axis, new_trait_axis, stageexpansion_list, errcheck_mpm_list, mpm_list,
      tweights_list, start_list, vrm_list, stageframe_df, allmodels_all,
      allstages_all, supplement_df, year_vec, sp_density_list, density_df,
      dens_index_df, equivalence_list, sp_density_num_vec, entry_time_vec,
      inda_terms_num_vec, indb_terms_num_vec, indc_terms_num_vec,
      inda_terms_cat_vec, indb_terms_cat_vec, indc_terms_cat_vec, dens_vr_yn_vec,
      tweights_type_vec, fecmod_vec, patch_vec, variant_count, var_per_run,
      nreps, times, fitness_times, stagecounts, substoch, format_int,
      preexisting_mpm_size, firstage_int, finalage_int, exp_tol, theta_tol,
      integeronly, stages_not_equal, stochastic, dens_yn_bool, entry_time_vec_use,
      sparse_bool, historical, pure_leslie, A_only, err_check_bool,
      err_check_extreme);
    
  } else if (funcbased) {
    //Rcout << "invade3 f ";
    
    invade3_fb_core (Lyapunov, var_run_mat, N_out, comm_out, zero_stage_vec_list,
      trait_axis, new_trait_axis, stageexpansion_list, modified_dev_terms_list,
      errcheck_mpm_list, errcheck_mpmout_list, tweights_list, start_list,
      vrm_list, final_stageframe, allmodels_all, allstages_all, supplement_df,
      year_vec, ind_terms_num_list, ind_terms_cat_list, dev_terms_list,
      density_vr_list, sp_density_list, density_df, dens_index_df,
      equivalence_list, sp_density_num_vec, entry_time_vec, inda_terms_num_vec,
      indb_terms_num_vec, indc_terms_num_vec, inda_terms_cat_vec,
      indb_terms_cat_vec, indc_terms_cat_vec, dens_vr_yn_vec, tweights_type_vec,
      fecmod_vec, patch_vec, variant_count, var_per_run, nreps, times,
      fitness_times, stagecounts, substoch, format_int, firstage_int,
      finalage_int, dev_terms_times_int, exp_tol, theta_tol, integeronly,
      stochastic, dens_yn_bool, stages_not_equal, sparse_bool, historical,
      pure_leslie, A_only, err_check_bool, err_check_extreme);
  } // funcbased processing
  
  //Rcout << "invade3 g ";
  
  int out_dim = 7;
  if (err_check_bool) out_dim++;
  List output (out_dim);
  
  output(0) = Lyapunov; // Needed in final output
  output(1) = comm_out; // Needed in final output
  output(2) = N_out; // Needed in final output
  output(3) = stageframe_df; // Needed in final output
  output(4) = hstages_list;
  output(5) = agestages_list;
  output(6) = labels;
  
  if (err_check_bool) {
    if (err_check_extreme) {
      List output_errcheck (15);
      
      output_errcheck(0) = allstages_all;
      output_errcheck(1) = allmodels_all;
      output_errcheck(2) = equivalence_list;
      output_errcheck(3) = density_df;
      output_errcheck(4) = dens_index_df;
      output_errcheck(5) = density_vr_list;
      output_errcheck(6) = new_trait_axis;
      output_errcheck(7) = stageexpansion_list;
      output_errcheck(8) = dev_terms_list;
      output_errcheck(9) = modified_dev_terms_list;
      output_errcheck(10) = errcheck_mpm_list;
      output_errcheck(11) = errcheck_mpmout_list;
      output_errcheck(12) = var_run_mat;
      output_errcheck(13) = start_list;
      output_errcheck(14) = final_stageframe;
      
      CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
        "equivalence_list", "density_df", "dens_index_df", "density_vr_list",
        "trait_axis_reassessed", "stageexpansion_list", "dev_terms_list",
        "modified_dev_terms_list", "modified_mpms", "fb_mpm_out_matrices",
        "var_run_mat", "start_list", "final_stageframe"};
      output_errcheck.attr("names") = output_errcheck_names;
      
      output(7) = output_errcheck;
      
      CharacterVector output_main_names = {"fitness", "variants_out", "N_out",
        "stageframe", "hstages", "agestages", "labels", "err_check"};
      output.attr("names") = output_main_names;
    } else {
      List output_errcheck (12);
      
      output_errcheck(0) = allstages_all;
      output_errcheck(1) = allmodels_all;
      output_errcheck(2) = equivalence_list;
      output_errcheck(3) = density_df;
      output_errcheck(4) = dens_index_df;
      output_errcheck(5) = density_vr_list;
      output_errcheck(6) = new_trait_axis;
      output_errcheck(7) = stageexpansion_list;
      output_errcheck(8) = dev_terms_list;
      output_errcheck(9) = modified_dev_terms_list;
      output_errcheck(10) = var_run_mat;
      output_errcheck(11) = start_list; // Remove later
      
      CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
        "equivalence_list", "density_df", "dens_index_df", "density_vr_list",
        "trait_axis_reassessed", "stageexpansion_list", "dev_terms_list",
        "modified_dev_terms_list", "var_run_mat", "start_list"};
      output_errcheck.attr("names") = output_errcheck_names;
      
      output(7) = output_errcheck;
      
      CharacterVector output_main_names = {"fitness", "variants_out", "N_out",
        "stageframe", "hstages", "agestages", "labels", "err_check"};
      output.attr("names") = output_main_names;
    }
  } else {
    CharacterVector output_main_names = {"fitness", "variants_out", "N_out",
      "stageframe", "hstages", "agestages", "labels"};
    output.attr("names") = output_main_names;
  }
  output.attr("class") = "adaptInv";
  
  return output;
}

