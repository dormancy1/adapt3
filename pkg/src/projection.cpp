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
// 5. List batch_project3  Project Multiple MPMs In Batches With Varying Starting Conditions
// 
// 6. List cleanup3_inv  Clean Up RObject Inputs for Invasion Analysis
// 7. void invpre_optim_singlerun  Core Pre-Existing MPM Projection Engine for ESS Evaluation
// 8. void ESS_optimizer_pre  Find ESS Values of Traits in Pre-Existing MPM Invasibility Analyses
// 9. void invfb_optim_singlerun  Core Function-Based Projection Engine for ESS Evaluation
// 10. void ESS_optimizer_fb  Find ESS Values of Traits in Function-Based Invasibility Analyses
// 11. void invpre_project  Core Pre-Existing MPM Projection Engine
// 12. void invpre_optim  Core Pre-Existing MPM Projection Engine for ESS Evaluation
// 13. void invade3_pre_core  Set-Up Running Invasion Analyses of Existing MPMs
// 14. void invfb_project  Core Function-Based Projection Engine
// 15. void invfb_optim  Core Function-Based Projection Engine for ESS Evaluation
// 16. void invade3_fb_core  Set-Up Function Running Invasion Analyses of Function-based MPMs
// 17. List invade3  Run Pairwise and Multiple Invasion Analysis


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
//' function-based MPMs, or for post-processing in pre-existing MPMs. If used,
//' then should be the same number of elements as the number of MPMs provided in
//' the list for argument \code{vrms}, with each element either a data frame or
//' \code{NULL}. See \code{\link[lefko3]{supplemental}()} for details.
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
//' use if no MPM has age structure. Defaults to \code{1} in age-based and age-
//' by-stage MPMs.
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
//' created with function \code{\link[lefko3]{density_input}()}, or a two-
//' layered list of such data frames, with a data frame per annual matrix per
//' MPM.
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
//' \code{vrm_input} object. Each should include 14 or 17 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival,
//' 9: juvenile observation, 10: juvenile primary size, 11: juvenile secondary
//' size, 12: juvenile tertiary size, 13: juvenile reproduction, and
//' 14: juvenile maturity transition. In addition, these may be followed by 3
//' columns designating additive deviations to: 15: individual covariate a,
//' 16: individual covariate b, and 17: individual covariate c. Unused terms
//' must be set to \code{0} (use of \code{NA} will produce errors). Single or
//' small numbers of values per vital rate model are also allowed, and if the
//' number of rows is less than \code{times}, then the terms will be cycled.
//' @param fb_sparse A logical vector indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE} for
//' each MPM.
//' @param equivalence An optional numeric vector, list of numeric vectors,
//' data frame of class \code{adaptEq} or class \code{lefkoEq}, or a list of such
//' data frames. If a numeric vector, then must have the same number of elements
//' as the number of MPMs, with each element giving the effect of an individual
//' of each MPM relative to a reference individual. If a list of vectors, then
//' the list should be composed of as many numeric vectors as MPMs, or as a two-
//' layered list of such vectors for each annual matrix per MPM, with each
//' vector giving the effect of each individual in each stage relative to a
//' reference individual. Data frames of class \code{adaptEq}, and lists of such
//' data frames, can be made with function \code{\link{equiv_input}()}. Numeric
//' entries used in these vectors can be thought of as Lotka-Volterra
//' competition terms, such as are used in multiple species competition models.
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
  List err_check_mat_indices;
  List err_check_topa3;
  
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
  IntegerVector dens_list_length_vec;
  IntegerVector eq_list_length_vec;
  
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
  
  // Rcout << "cleanup3 A    ";
  
  if (mpms.isNotNull()) {
    if (vrms.isNotNull() || stageframes.isNotNull()) {
      AdaptUtils::pop_error2("vrms", "stageframes", "projecting existing MPMms", 24);
    }
    
    if (is<List>(mpms)) { 
      List original_mpms = as<List>(mpms);
      mpm_list = clone(original_mpms);
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
  
  // Rcout << "cleanup3 B    ";
  
  if (vrms.isNotNull()) {
    // Rcout << "cleanup3 B1    ";
    
    if (is<List>(vrms)) {
      vrm_list = as<List>(vrms);
      vrm_count = static_cast<int>(vrm_list.length());
      
      if (format.isNotNull()) {
        // Rcout << "cleanup3 B3 format argument found   ";
        
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
            } else if (format_vec(i) == 5) { //&& !stageframes.isNotNull()
              found_fleslie++;
            }
          }
          if (found_fleslie == vrm_count) pure_fleslie = true;
        } else AdaptUtils::pop_error2("format", "an integer vector", "", 1);
      } else {
        IntegerVector format_vec_pre (vrm_count, 3);
        
        ///// Insert check for Leslie matrix here
        IntegerVector potential_minage_vec (vrm_count);
        IntegerVector potential_maxage_vec (vrm_count);
        
        for (int i = 0; i < vrm_count; i++) {
          
          int minage_found {0};
          int maxage_found {0};
          bool found_leslie {true};
          
          if (stageframes.isNotNull()) {
            List all_stageframes = as<List>(stageframes);
            Nullable<DataFrame> current_stageframe_test = as<Nullable<DataFrame>>(all_stageframes(i));
            
            if (current_stageframe_test.isNotNull()) {
              DataFrame cst = as<DataFrame>(current_stageframe_test);
              
              IntegerVector cst_minage = as<IntegerVector>(cst["min_age"]);
              int cst_minage_length = static_cast<int>(cst_minage.length());
              
              IntegerVector cst_minage_sorted = cst_minage.sort();
              
              for (int j = 0; j < (cst_minage_length - 1); j++) {
                if (cst_minage_sorted(j+1) != (cst_minage_sorted(j) + 1)) found_leslie = false;
                if (IntegerVector::is_na(cst_minage(j)) || IntegerVector::is_na(cst_minage(j + 1))) found_leslie = false;
                
              }
              
              if (found_leslie) {
                minage_found = min(cst_minage_sorted);
                maxage_found = max(cst_minage_sorted);
                
                if (!firstage.isNotNull()) {
                  potential_minage_vec(i) = minage_found;
                  firstage = as<Nullable<RObject>>(potential_minage_vec);
                }
                
                if (!finalage.isNotNull()) {
                  potential_maxage_vec(i) = maxage_found + 1;
                  finalage = as<Nullable<RObject>>(potential_maxage_vec);
                }
              }
            }
          } else {
            if (firstage.isNotNull() && finalage.isNotNull()) {
              IntegerVector firstage_vec_current = as<IntegerVector>(firstage);
              IntegerVector finalage_vec_current = as<IntegerVector>(finalage);
              
              if (IntegerVector::is_na(firstage_vec_current(i)) || IntegerVector::is_na(finalage_vec_current(i))) {
                found_leslie = false;
              }
            } else found_leslie = false;
          }
          if (found_leslie) format_vec_pre(i) = 5;
        }
        
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
      Rf_warningcall(R_NilValue,
        "Each vrm_input object must have its own stageframe, except for Leslie MPMs.");
    }
    
    funcbased = true;
  }
  
  if (!preexisting && !funcbased) {
    throw Rcpp::exception("Cannot proceed without either argument mpms, or arguments vrms and stageframes set.",
      false);
  } else if (preexisting && funcbased) {
    throw Rcpp::exception("Cannot proceed with argument mpms, vrms, and stageframes set.",
      false);
  }
  
  total_mpms = mpm_count + vrm_count;
  
  // Rcout << "cleanup3 C    ";
  
  if (supplements.isNotNull()) {
    if (is<List>(supplements)) {
      supplement_list_fb = as<List>(supplements);
      supplement_count = static_cast<int>(supplement_list_fb.length());
      if (funcbased && supplement_count != vrm_count) {
        AdaptUtils::pop_error2("vrms", "supplements", "lists of the same length", 27);
      } else if (!funcbased && supplement_count != mpm_count) {
        AdaptUtils::pop_error2("mpms", "supplements", "lists of the same length", 27);
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
  
  // Rcout << "cleanup3 D    ";
  
  if (format.isNotNull()) {
    if (!funcbased) {
      AdaptUtils::pop_error2("vrms", "use argument format", "", 26);
    }
    
  } else if (preexisting) {
    IntegerVector format_vec_pre (mpm_count, 3);
    
    for (int i = 0; i < mpm_count; i++) {
      List chosen_mpm = as<List>(mpm_list(i));
      
      int current_mpm_format = LefkoInputs::format_check_lM (chosen_mpm);
      format_vec_pre(i) = current_mpm_format;
    }
    format_vec = format_vec_pre;
    
  } else if (funcbased && sum(format_vec) == 0) {
    IntegerVector format_vec_pre (vrm_count, 3);
    format_vec = format_vec_pre;
  }
  
  // Rcout << "cleanup3 E    ";
  
  // firstage, finalage, and cont processing for age-by-stage and Leslie MPMs
  if (firstage.isNotNull() || finalage.isNotNull() || cont.isNotNull()) {
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
        if (firstage_vec(i) < 0 && !IntegerVector::is_na(firstage_vec(i))) {
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
    
    for (int i = 0; i < static_cast<int>(format_vec.length()); i++) {
      if (format_vec(i) > 3) firstage_vec(i) = 1;
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
        if (finalage_vec(i) < 0 && !IntegerVector::is_na(finalage_vec(i))) {
          AdaptUtils::pop_error2("finalage", "", "", 30);
        }
        
        if (finalage_vec(i) > 0 && format_vec(i) < 4) {
          throw Rcpp::exception("Entries in argument finalage must equal 0 for MPMs without age structure.", false);
        }
        
        if (finalage_vec(i) < firstage_vec(i) && !IntegerVector::is_na(finalage_vec(i))) {
          throw Rcpp::exception("Values in finalage may not be less than respective entries in firstage.", false);
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
  
  // Rcout << "cleanup3 F    ";
  
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
  
  // Rcout << "cleanup3 G    ";
  
  if (cont.isNotNull()) {
    int useable_count = vrm_count;
    if (mpm_count > vrm_count) useable_count = mpm_count;
    
    if (is<IntegerVector>(cont) || is<NumericVector>(cont)) { 
      IntegerVector cont_vec_temp = as<IntegerVector>(cont);
      
      if (static_cast<int>(cont_vec_temp.length() == 1)) {
        IntegerVector cont_maxed_out (useable_count);
        
        for (int i = 0; i < useable_count; i++) {
          cont_maxed_out(i) = cont_vec_temp(0);
        }
        cont_vec = cont_maxed_out;
      } else cont_vec = cont_vec_temp;
      
      int cont_vec_length = static_cast<int>(cont_vec.length());
      if (cont_vec_length != useable_count) {
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
        LogicalVector cont_maxed_out (useable_count);
        
        for (int i = 0; i < useable_count; i++) {
          cont_maxed_out(i) = cont_vec_temp(0);
        }
        cont_vec_log = cont_maxed_out;
      } else cont_vec_log = cont_vec_temp;
      
      int cont_vec_log_length = static_cast<int>(cont_vec_log.length());
      if (cont_vec_log_length != useable_count) {
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
  
  // Rcout << "cleanup3 H    ";
  
  // Altered stageframe processing
  if (funcbased) {
    // Rcout << "Entered funcbased stageframe processing." << endl;
    
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
        DataFrame imported_stageframe = as<DataFrame>(stageframe_list_fb(i));
        
        if (static_cast<int>(imported_stageframe.length()) > 1) {
          AdaptUtils::leslie_stageframe_updater(firstage_vec, finalage_vec,
            fecage_min_vec, fecage_max_vec, cont_vec, imported_stageframe, i);
        }
        
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
        int found_numrows = static_cast<int>(new_stageframe.nrows());
        matrows_pre(i) = found_numrows;
        
        stageframe_list_pre(i) = new_stageframe;
        supplement_list_pre(i) = new_ovtable;
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
  
  // Rcout << "cleanup3 I    ";
  
  // start vector
  if (starts.isNotNull()) {
    // Rcout << "cleanup3 I1 starts argument found    ";
    
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
    // Rcout << "cleanup3 I2 no starts argument found    ";
    
    // Construct default list of start vectors (1 indiv / stage)
    IntegerVector stagecounts_pre (stageframe_count);
    List start_list_pre (stageframe_count);
    
    for (int i = 0; i < stageframe_count; i++) {
      int scp {0};
      
      if (format_vec(i) != 5) {
        DataFrame chosen_stageframe;
        
        if (format_vec(i) == 3) {
          chosen_stageframe = as<DataFrame>(stageframe_list(i));
        } else if (format_vec(i) < 3) {
          chosen_stageframe = as<DataFrame>(hstages_list(i));
        } else if (format_vec(i) == 4) {
          chosen_stageframe = as<DataFrame>(agestages_list(i));
        }
        
        scp = static_cast<int>(chosen_stageframe.nrows());
      } else {
        Nullable<RObject> chosen_stageframe_RO = as<Nullable<RObject>>(stageframe_list(i));
        if (chosen_stageframe_RO.isNotNull()) {
          DataFrame chosen_stageframe = as<DataFrame>(stageframe_list(i));
          scp = static_cast<int>(chosen_stageframe.nrows());
        } else {
          int current_firstage {0};
          int current_finalage {0};
          
          if (IntegerVector::is_na(finalage_vec(i))) {
            throw Rcpp::exception("Function-based Leslie matrices cannot be generated without finalage provided.", false);
          } else {
            current_finalage = finalage_vec(i);
          }
          if (IntegerVector::is_na(firstage_vec(i))) {
            Rf_warningcall(R_NilValue, "Assuming that the function-based Leslie MPM life history starts with age 0.");
          } else {
            current_firstage = firstage_vec(i);
          }
          
          scp = current_finalage - current_firstage + 1;
        }
      }
 
      if (format_vec(i) == 3 && funcbased) scp--;
      stagecounts_pre(i) = scp;
      
      arma::vec start_vec (scp, fill::ones);
      start_list_pre(i) = start_vec;
    }
    
    stagecounts = stagecounts_pre;
    start_list = start_list_pre;
    start_count = stageframe_count;
  }
  
  // Rcout << "cleanup3 J    ";
  
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
  
  // Rcout << "cleanup3 K    ";
  
  // loy label df construction (no years)
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
  
  // Rcout << "cleanup3 L    ";
  
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
      int year_list_pre_length = static_cast<int>(year_list_pre.length());
      
      if (preexisting) {
        if (year_list_pre_length != mpm_count) {
          throw Rcpp::exception("List supplied in argument years must be the same length as argument mpms.", false);
        }
        
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
        if (year_list_pre_length != vrm_count) {
          throw Rcpp::exception("List supplied in argument years must be the same length as argument vrms.", false);
        }
        
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
  
  // Rcout << "cleanup3 M    ";
  
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
  
  // Rcout << "cleanup3 N    " << endl;;
  
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
      IntegerVector dens_list_length_vec_pre (density_count);
      
      for (int i = 0; i < density_count; i++) {
        if (is<DataFrame>(density_list(i))) {
          // Rcout << "cleanup3 N5 DataFrame    " << endl;;
          
          DataFrame chosen_density = as<DataFrame>(density_list(i));
          
          DataFrame hstages = as<DataFrame>(hstages_list(i));
          DataFrame agestages = as<DataFrame>(agestages_list(i));
          
          DataFrame stageframe;
          if (preexisting) {
            List chosen_mpm = as<List>(mpm_list(i));
            stageframe = as<DataFrame>(chosen_mpm["ahstages"]);
          } else {
            stageframe = as<DataFrame>(stageframe_list(i));
          }
          
          int format = format_vec(i);
          int finalage = finalage_vec(i);
          
          List dens_index;
          arma::uvec dyn_style;
          arma::vec dyn_alpha;
          arma::vec dyn_beta;
          arma::vec dyn_gamma; 
          int dens_yn_int {0};
          
          // Rcout << "cleanup3 N5h    " << endl;;
          
          AdaptUtils::density_prep (dens_index, dyn_style, dyn_alpha, dyn_beta,
            dyn_gamma, dens_yn_int, chosen_density, hstages, agestages,
            stageframe, exp_tol, format, finalage, preexisting, funcbased);
          
          // Rcout << "cleanup3 N5i    " << endl;;
          
          dens_index_list_pre(i) = dens_index;
          dens_yn_vec_temp(i) = dens_yn_int;
          
        } else if (is<List>(density_list(i))) {
          // Rcout << "cleanup3 N6 List    " << endl;;
          
          List density_second_tier = as<List>(density_list(i));
          int length_of_density_list_i = static_cast<int>(density_second_tier.length());
          dens_list_length_vec_pre(i) = length_of_density_list_i;
          List dens_index_list_pre_pre (length_of_density_list_i);
          int finalage = finalage_vec(i);
          
          int dens_yn_int {0};
          for (int j = 0; j < length_of_density_list_i; j++) {
            DataFrame chosen_density = as<DataFrame>(density_second_tier(j));
            
            DataFrame hstages = as<DataFrame>(hstages_list(i));
            DataFrame agestages = as<DataFrame>(agestages_list(i));
            
            DataFrame stageframe;
            if (preexisting) {
              List chosen_mpm = as<List>(mpm_list(i));
              stageframe = as<DataFrame>(chosen_mpm["ahstages"]);
            } else {
              stageframe = as<DataFrame>(stageframe_list(i));
            }
            
            int format = format_vec(i);
            
            List dens_index;
            arma::uvec dyn_style;
            arma::vec dyn_alpha;
            arma::vec dyn_beta;
            arma::vec dyn_gamma; 
            int dens_yn_int_i {0};
            
            AdaptUtils::density_prep (dens_index, dyn_style, dyn_alpha, dyn_beta,
              dyn_gamma, dens_yn_int_i, chosen_density, hstages, agestages,
              stageframe, exp_tol, format, finalage, preexisting, funcbased);
            
            dens_index_list_pre_pre(j) = dens_index;
            
            if (dens_yn_int_i > 0) dens_yn_int = 1;

          }
          dens_index_list_pre(i) = dens_index_list_pre_pre;
          dens_yn_vec_temp(i) = dens_yn_int;
        } else if (density_list(i) != R_NilValue) { 
          AdaptUtils::pop_error2("density", "a list of lefkoDens objects and NULL values", "", 1);
        }
      }
      dens_index_list = dens_index_list_pre;
      dens_yn_vec = dens_yn_vec_temp;
      dens_list_length_vec = dens_list_length_vec_pre;
      
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
  
  // Rcout << "cleanup3 O    ";
  
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
  
  // Rcout << "cleanup3 P    ";
  
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
          
          if (dvtc_df_size != 14 && dvtc_df_size != 17) {
            throw Rcpp::exception("Data frames in argument dev_terms must have 14 or 17 columns.",
              false);
          }
          
          for (int j = 0; j < dvtc_df_size; j++) {
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
  
  // Rcout << "cleanup3 Q    ";
  
  // equivalence interpretation
  if (equivalence.isNotNull()) {
    if (is<NumericVector>(equivalence)) {
      equivalence_vec = as<NumericVector>(equivalence);
      
      equivalence_count = static_cast<int>(equivalence_vec.length());
      IntegerVector eq_list_length_vec_pre (equivalence_count);
      eq_list_length_vec = eq_list_length_vec_pre;
      
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
      IntegerVector eq_list_length_vec_pre (equivalence_count);
      
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
            throw Rcpp::exception("Argument equivalence should include data frames of class adaptEq, or numeric vectors.", 
              false);
          }
          CharacterVector eq_list_df_class = eq_list_df.attr("class");
          bool found_adaptEq {false};
          for (int j = 0; j < static_cast<int>(eq_list_df_class.length()); j++) {
            if (eq_list_df_class(j) == "adaptEq") found_adaptEq = true;
            if (eq_list_df_class(j) == "lefkoEq") found_adaptEq = true;
          }
          if (!found_adaptEq) {
            throw Rcpp::exception("Argument equivalence should include data frames of class adaptEq, or numeric vectors.", 
              false);
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
              throw Rcpp::exception("Enter ages in adaptEq objects used for Leslie and age-by-stage MPMs.",
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
          
        } else if (is<List>(equivalence_list_temp(i))) {
          List equiv_list_current = as<List>(equivalence_list_temp(i));
          int equiv_list_current_length = static_cast<int>(equiv_list_current.length());
          eq_list_length_vec_pre(i) = equiv_list_current_length;
          List equivalence_list_pre_pre (equiv_list_current_length);
          
          for (int m = 0; m < equiv_list_current_length; m++) {
            DataFrame eq_list_df = as<DataFrame>(equiv_list_current(m));
            if (!eq_list_df.hasAttribute("class")) {
              throw Rcpp::exception("Argument equivalence should include data frames of class adaptEq, or numeric vectors.", false);
            }
            CharacterVector eq_list_df_class = eq_list_df.attr("class");
            bool found_adaptEq {false};
            for (int j = 0; j < static_cast<int>(eq_list_df_class.length()); j++) {
              if (eq_list_df_class(j) == "adaptEq") found_adaptEq = true;
              if (eq_list_df_class(j) == "lefkoEq") found_adaptEq = true;
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
                throw Rcpp::exception("Enter ages in adaptEq objects used for Leslie and age-by-stage MPMs.",
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
            
            equivalence_list_pre_pre(m) = current_eq;
          }
          equivalence_list_pre(i) = equivalence_list_pre_pre;
        }
      }
      equivalence_list = equivalence_list_pre;
      eq_list_length_vec = eq_list_length_vec_pre;
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
  
  // Rcout << "cleanup3 R    " << endl;
  
  // process stageframe, supplement, repmatrix, and allstages list for fbMPMs
  if (funcbased) {
    // Rcout << "cleanup3 R1    " << endl;
    
    // Create function-based MPMs and assign them to mpm_list
    List allstages_all_pre (vrm_count);
    List allmodels_all_pre (vrm_count);
    
    for (int i = 0; i < vrm_count; i++) {
      // Rcout << "cleanup3 R2    i: " << i << endl;;
      
      List current_vrm = as<List>(vrm_list(i));
      DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
      DataFrame current_supplement = as<DataFrame>(supplement_list(i));
      
      arma::mat current_repmatrix;
      if (format_vec(i) < 5) current_repmatrix = as<arma::mat>(repmatrix_list(i));
      
      int ehrlen_format {1}; // This will need to be dealt with differently later
      if (format_vec(i) == 2) ehrlen_format = 2;
      
      int mpm_style {1};
      int filter_style {1};
      if (format_vec(i) < 3) {
        mpm_style = 0;
      } else if (format_vec(i) == 4) {
        mpm_style = 2;
        filter_style = 2;
      }
      
      // Rcout << "cleanup3 R4    ";
      
      DataFrame current_mpm_allstages;
      if (format_vec(i) < 5) {
        current_mpm_allstages = LefkoMats::theoldpizzle(current_stageframe,
          current_supplement, current_repmatrix, firstage_vec(i), finalage_vec(i),
          ehrlen_format, mpm_style, cont_vec(i), filter_style); // Last term removes unused rows & cols
      } else {
        DataFrame leslie_allstages = clone(as<DataFrame>(stageframe_list(i)));
        current_mpm_allstages = leslie_allstages;
      }
      allstages_all_pre(i) = current_mpm_allstages;
      
      // Rcout << "cleanup3 R5    ";
      
      if (format_vec(i) < 5) {
        DataFrame chosen_stageframe_pre = clone(as<DataFrame>(stageframe_list(i)));
        
        IntegerVector removal_row = {static_cast<int>(chosen_stageframe_pre.nrows())};
        StringVector removal_var = {"stage_id"};
        DataFrame chosen_stageframe = LefkoUtils::df_remove(chosen_stageframe_pre,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        stageframe_list(i) = chosen_stageframe;
      }
      
      // Rcout << "cleanup3 R6    ";
      
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
        // Rcout << "cleanup3 R13    ";
        
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
        // Rcout << "cleanup3 R16    ";
        
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
      
      // Rcout << "cleanup3 R24    ";
      
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
      
      // Rcout << "cleanup3 R25    ";
      
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
      
      // Rcout << "cleanup3 R29    ";
      
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
  
  // Rcout << "cleanup3 S    ";
  
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
  
  // Rcout << "cleanup3 T    ";
  
  List errcheck_mpmout_vrm (vrm_count);
  
  if (fb_override) {
    // Create function-based MPMs and assign them to mpm_list
    // Rcout << "Entered fb_override    ";
    
    List mpm_list_pre (vrm_count);
    List A_list_pre (vrm_count);
    
    IntegerVector year_counter (vrm_count);
    
    for (int i = 0; i < vrm_count; i++) {
      List current_vrm_extract = as<List>(allmodels_all(i));
      List current_vrm_unextract = as<List>(vrm_list(i));
      DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
      
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
      
      // Rcout << "cleanup3 T2    ";
      
      CharacterVector current_mainyears = as<CharacterVector>(year_list(i));
      unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
      
      DataFrame group2_frame = as<DataFrame>(current_vrm_unextract["group2_frame"]);
      CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
      
      DataFrame patch_frame = as<DataFrame>(current_vrm_unextract["patch_frame"]);
      CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
      
      DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
      DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
      DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
      CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
      CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
      CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
      
      // Rcout << "cleanup3 T3    ";
      
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
      
      IntegerVector found_calls = {total_years_vec(i), sp_density_num_vec(i),
        dev_terms_num_vec(i), inda_terms_num_vec(i), indb_terms_num_vec(i),
        indc_terms_num_vec(i), inda_terms_cat_vec(i), indb_terms_cat_vec(i),
        indc_terms_cat_vec(i)};
      int needed_calls = max(found_calls);
      
      List current_building_mpm (needed_calls);
      List errcheck_mpmout_vrm_time (needed_calls);
      CharacterVector labels_year2_terms (needed_calls);
      CharacterVector labels_patch_terms (needed_calls);
      
      for (int j = 0; j < needed_calls; j++) { // time loop
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
        
        // Rcout << "cleanup3 T8    ";
        
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
        
        // Rcout << "cleanup3 T9    ";
        
        NumericVector dv_terms (17);
        if (dev_terms_num_vec(i) > 0) {
          DataFrame used_dv_df = as<DataFrame>(dev_terms_list(i));
          int used_dv_df_size = static_cast<int>(used_dv_df.size());
        
          if (dev_num_counter(i) >= dev_terms_num_vec(i)) dev_num_counter(i) = 0;
          
          for (int j = 0; j < used_dv_df_size; j++) {
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
        
        // Rcout << "cleanup3 T11    ";
        
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
          // Rcout << "cleanup3 T13    ";
          
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
          // Rcout << "cleanup3 T16    ";
          
          current_mpm = AdaptMats::mazurekd(current_mpm_allstages, // AdaptMats::mazurekd
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
            
          // Rcout << "cleanup3 T16a    ";
          
          if (err_check) {
            NumericMatrix mpm_out = as<NumericMatrix>(current_mpm["out"]);
            errcheck_mpmout_vrm_time(j) = mpm_out; 
          }
        } else {
          // Rcout << "cleanup3 T17    ";
          
          IntegerVector all_ages = seq(firstage_vec(i), finalage_vec(i));
          DataFrame current_supplement;
          if (supplement_list(i) == R_NilValue) {
            NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
            current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
              surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
              f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
              r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(i),
              finalage_vec(i), true, yearnumber, patchnumber, dvr_bool, dvr_yn,
              dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
              sparse_vec(i));
            
          } else {
            current_supplement = as<DataFrame>(supplement_list(i));
            
            NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
            current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
              surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
              f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
              r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(i),
              finalage_vec(i), true, yearnumber, patchnumber, dvr_bool, dvr_yn,
              dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
              sparse_vec(i), current_supplement);
          }
        }
        // Rcout << "cleanup3 T18    ";
        
        arma::mat pulled_matrix = as<arma::mat>(current_mpm(0));
        current_building_mpm(j) = pulled_matrix;
        
        labels_year2_terms(j) = current_year(0);
        labels_patch_terms(j) = current_mainpatches(patchnumber);
        
        year_counter(i) = year_counter(i) + 1;
      } // time loop (j)
      if (err_check) errcheck_mpmout_vrm(i) = errcheck_mpmout_vrm_time;
      
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
  
  // Rcout << "cleanup3 U    ";
  // Post-processing for preexisting MPMs with new supplements
  
  if (preexisting && supplements.isNotNull()) {
    List fortified_A_list (mpm_count);
    List err_check_indices_list (mpm_count);
    List err_check_theoldpizzle_adapt3 (mpm_count);
    
    for (int i = 0; i < mpm_count; i++) {
      // Rcout << "cleanup3 U1    " << endl;
      List err_check_indices_list_current_mpm;
      
      int ehrlen_format {1}; // This will need to be dealt with differently later
      bool historical {false};
      
      int mpm_style {1};
      int filter_style {1};
      if (format_vec(i) < 3) {
        mpm_style = 0;
      } else if (format_vec(i) == 4) {
        mpm_style = 2;
        //filter_style = 2;
      }
      if (format_vec(i) < 3) historical = true;
      if (format_vec(i) == 2) ehrlen_format = 2;
      
      DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
      DataFrame current_supplement;
      int cs_length {1};
      if (is<DataFrame>(supplement_list_fb(i))) {
        current_supplement = as<DataFrame>(supplement_list_fb(i));
        cs_length = static_cast<int>(current_supplement.length());
      }
      
      bool agemat {false};
      if (format_vec(i) == 4) agemat = true;
      
      if (cs_length > 1) {
        List melchett = LefkoMats::sf_reassess_internal(current_stageframe,
          current_supplement, R_NilValue, R_NilValue, agemat, historical,
          ehrlen_format, false);
        DataFrame new_stageframe = as<DataFrame>(melchett["stageframe"]);
        arma::mat new_repmatrix = as<arma::mat>(melchett["repmatrix"]);
        
        // Rcout << "cleanup3 U4    " << endl;
        
        DataFrame alterations = AdaptMats::theoldpizzle_adapt3(current_stageframe,
          current_supplement, new_repmatrix, firstage_vec(i), finalage_vec(i),
          ehrlen_format, mpm_style, cont_vec(i), filter_style, true);
        if (err_check) err_check_theoldpizzle_adapt3(i) = alterations;
        
        // Rcout << "cleanup3 U5    " << endl;
        
        List used_mpm = as<List>(mpm_list(i));
        List used_A_list = as<List>(used_mpm["A"]);
        List used_U_list = as<List>(used_mpm["U"]);
        List used_F_list = as<List>(used_mpm["F"]);
        
        AdaptMats::matrix_post(used_A_list, used_U_list, used_F_list,
          err_check_indices_list_current_mpm, alterations, current_stageframe,
          static_cast<int>(format_vec(i)), false, static_cast<bool>(sparse_vec(i)),
          err_check);
        
        // Rcout << "cleanup3 U10    " << endl;
        
        used_mpm["A"] = used_A_list;
        used_mpm["U"] = used_U_list;
        used_mpm["F"] = used_F_list;
        
        mpm_list(i) = used_mpm;
      }
      if (err_check) err_check_indices_list(i) = err_check_indices_list_current_mpm;
    }
    
    if (err_check) {
      err_check_mat_indices = err_check_indices_list;
      err_check_topa3 = err_check_theoldpizzle_adapt3;
    }
  }
  
  // Rcout << "cleanup3 V    ";
  // Output processing
  List out_list (76);
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
  out_list(72) = dens_list_length_vec;
  out_list(73) = eq_list_length_vec;
  out_list(74) = err_check_mat_indices;
  out_list(75) = err_check_topa3;
  
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
    "fb_override", "errcheck_mpmout_vrm", "dens_list_length_vec",
    "eq_list_length_vec", "err_check_mat_indices", "err_check_theoldpizzle_adapt3"};
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
//' @param agg_density A numeric matrix giving the aggregate community density
//' by time (column) and replicate (row).
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
//' @param dev_terms_list List of deviations for vital rate models in all MPMs,
//' plus additive deviations to inda, indb, and indc.
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
//' @param dens_list_length_vec A vector giving the number of data frames of
//' density info per MPM, if lists of such data frames are used.
//' @param eq_list_length_vec A vector giving the number of data frames of
//' equivalence info per MPM, if lists of such data frames are used.
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
void project3_pre_core (NumericMatrix& agg_density, List& N_out, List& comm_out,
  List& extreme_mpm_out, const List mpm_list, const List A_list,
  const List tweights_list, const List start_list, const List vrm_list,
  const List stageframe_list, const List allmodels_all, const List allstages_all,
  const List supplement_list, const List year_list, const List ind_terms_num_list,
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
  const IntegerVector dens_list_length_vec, const IntegerVector eq_list_length_vec,
  const IntegerVector tweights_type_vec, const NumericVector fecmod_vec,
  const LogicalVector sparse_vec, const CharacterVector patch_vec,
  const int vrm_count, const int mpm_count, const int nreps, const int times,
  const int substoch, const double exp_tol, const double theta_tol,
  const bool integeronly, const bool stages_not_equal, const bool stochastic,
  const bool entry_time_vec_use, const bool err_check,
  const bool err_check_extreme) {
  
  // Matrix order set up and creation of zero stage vectors
  // Rcout << "Entered project3_pre_core          " << endl;
  
  List matrix_choice_list (mpm_count);
  List zvl (mpm_count);
  for (int i = 0; i < mpm_count; i++) {
    // Rcout << "project3_pre_core A  i (mpm): " << i << "          " << endl;
    
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
      
    } else {
      IntegerVector chosen_mats_pre = sort_unique(pvis);
      chosen_mats = chosen_mats_pre;
    }
    
    // Rcout << "project3_pre_core B    " << endl;
    matrix_choice_list(i) = chosen_mats;
  }
  
  // Year order determination
  List comm_out_pre (mpm_count);
  List used_times (mpm_count);
  
  // Rcout << "project3_pre_core C    " << endl;
  
  for (int i = 0; i < mpm_count; i++) {
    // Rcout << "project3_pre_core C1    " << endl;
    List pop_reps (nreps);
    
    arma::mat pops_out_pre (stagecounts(i), (times + 1), fill::zeros);
    arma::vec pops_out_0 = as<arma::vec>(start_list(i));
    
    pops_out_pre.col(entry_time_vec(i)) = pops_out_0; /////
    
    IntegerVector chosen_mats = as<IntegerVector>(matrix_choice_list(i));
    int chosen_mats_length = static_cast<int>(chosen_mats.length());
    
    List used_times_mpm (nreps);
    
    for (int j = 0; j < nreps; j++) {
      // Rcout << "project3_pre_core C3    " << endl;
      
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
        // Rcout << "project3_pre_core C6    " << endl;
        
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
          // Rcout << "project3_pre_core C10    " << endl;
          
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
          // Rcout << "project3_pre_core C13    " << endl;
          
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
  
  // Rcout << "project3_pre_core D    " << endl;
  
  // Main projection
  NumericMatrix agg_density_pre (nreps, (times + 1));
  List N_out_pre (nreps);
  List extreme_mpm_reps (nreps);
  
  for (int i = 0; i < nreps; i ++) {
    List running_popvecs = clone(start_list);
    NumericMatrix N_mpm (mpm_count, (times + 1));
    List extreme_mpm_reps_times (times);
    int overall_dens_yn = sum(dens_yn_vec);
    
    for (int j = 0; j < times; j++) {
      if (j % 10 == 0){
        Rcpp::checkUserInterrupt();
      }
      
      List extreme_mpm_reps_times_mpms (mpm_count);
      
      List dense_vecs (mpm_count);
      arma::vec delay_N_sum_vec (mpm_count);
      
      if (overall_dens_yn > 0) {
        for (int k = 0; k < mpm_count; k++) {
          if (j == 0) {
            arma::vec current_popvec = as<arma::vec>(running_popvecs(k));
            dense_vecs(k) = current_popvec;
            double delay_N_sum_vec_entry = accu(current_popvec);
            delay_N_sum_vec(k) = delay_N_sum_vec_entry;
          } else {
            if (dens_yn_vec(k) > 0) {
              DataFrame used_density_input;
              DataFrame used_density_index_input;
              
              if (dens_list_length_vec(k) > 0) {
                List current_used_times_mpm = as<List>(used_times(k));
                IntegerVector current_times_vec = as<IntegerVector>(current_used_times_mpm(i));
                
                List used_density_input_list = as<List>(density_list(k));
                List used_density_index_input_list = as<List>(dens_index_list(k));
                
                used_density_input = as<DataFrame>(used_density_input_list(current_times_vec(j)));
                used_density_index_input = as<DataFrame>(used_density_index_input_list(current_times_vec(j)));
                
              } else {
                used_density_input = as<DataFrame>(density_list(k));
                used_density_index_input = as<DataFrame>(dens_index_list(k));
              }
              
              IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
              int used_delay = max(ud_delay_vec);
              
              if ((j + 1 - used_delay) > 0) {
                List reps_out = comm_out_pre(k);
                arma::mat pops_out = as<arma::mat>(reps_out(i));
                arma::vec current_popvec = pops_out.col(j);
                dense_vecs(k) = current_popvec;
                
                double delay_N_sum_vec_entry = accu(current_popvec);
                delay_N_sum_vec(k) = delay_N_sum_vec_entry;
              } else {
                List new_start_list = clone(start_list);
                arma::vec current_popvec = as<arma::vec>(new_start_list(k));
                dense_vecs(k) = current_popvec;
                
                double delay_N_sum_vec_entry = accu(current_popvec);
                delay_N_sum_vec(k) = delay_N_sum_vec_entry;
              }
            }
          }
        }
      }
      
      for (int k = 0; k < mpm_count; k++) {
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
            
            double aN_first = accu(running_popvec_mpm);
            
            if (stages_not_equal) {
              arma::vec current_equivalence_vec;
              
              if (eq_list_length_vec(k) > 0) {
                List current_equivalence_list = as<List>(equivalence_list(k));
                current_equivalence_vec = as<arma::vec>(current_equivalence_list(j));
              } else {
                current_equivalence_vec = as<arma::vec>(equivalence_list(k));
              }
              arma::vec adjusted_pop_vec = running_popvec_mpm % current_equivalence_vec;
              aN_first = accu(adjusted_pop_vec);
            }
            
            agg_density_pre(i, j) = agg_density_pre(i, j) + aN_first;
          }
          
          List chosen_As = as<List>(A_list(k));
          
          if (dens_yn_vec(k) == 0) {
            // Rcout << "project3_pre_core H no density       ";
            
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
            // Rcout << "project3_pre_core I density       ";
            
            DataFrame used_density_input;
            DataFrame used_density_index_input;
            
            if (dens_list_length_vec(k) > 0) {
              List used_density_input_list = as<List>(density_list(k));
              List used_density_index_input_list = as<List>(dens_index_list(k));
              
              used_density_input = as<DataFrame>(used_density_input_list(current_times_vec(j)));
              used_density_index_input = as<DataFrame>(used_density_index_input_list(current_times_vec(j)));
              
            } else {
              used_density_input = as<DataFrame>(density_list(k));
              used_density_index_input = as<DataFrame>(dens_index_list(k));
            }
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            if ((j + 1 - used_delay) >= 0) {
              if (!stages_not_equal) {
                // Rcout << "project3_pre_core j density stages equal       ";
                
                double delay_N_sum = accu(delay_N_sum_vec);
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, chosen_As, delay_N_sum,
                  static_cast<int>(current_times_vec(j)), integeronly,
                  substoch, used_density_input, used_density_index_input,
                  false, static_cast<bool>(sparse_vec(k)),
                  static_cast<bool>(sparse_vec(k)), false, err_check_extreme);
                
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = new_projmat;
                
              } else {
                // Rcout << "project3_pre_core k density stages not equal       ";
                
                double delay_N_sum {0.0};
                
                for (int l = 0; l < mpm_count; l++) {
                  arma::vec current_equiv_vec;
                  if (eq_list_length_vec(l) > 0) {
                    List used_equivalence_list = as<List>(equivalence_list(l));
                    current_equiv_vec = as<arma::vec>(used_equivalence_list(current_times_vec(j)));
                  
                  } else {
                    current_equiv_vec = as<arma::vec>(equivalence_list(l));
                  }
                  
                  arma::vec delay_pop_vec = as<arma::vec>(dense_vecs(l));
                  arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                  double delay_pop_N = accu(adjusted_delay_pop_vec);
                  
                  delay_N_sum += delay_pop_N;
                }
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, chosen_As, delay_N_sum,
                  static_cast<int>(current_times_vec(j)), integeronly,
                  substoch, used_density_input, used_density_index_input,
                  false, static_cast<bool>(sparse_vec(k)),
                  static_cast<bool>(sparse_vec(k)), false, err_check_extreme);
                
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_mpms(k) = new_projmat;
              }
            } else {
              // Rcout << "project3_pre_core l density initial time       ";
              
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              AdaptUtils::proj3dens_ad(new_popvec, new_projmat, new_projsp,
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
          
          double aN_current = N_current; // Aggregate density, given stage weights
          arma::vec current_equivalence_vec;
          
          if (stages_not_equal) {
            if (eq_list_length_vec(k) > 0) {
              List current_equivalence_list = as<List>(equivalence_list(k));
              int c_eq_list_length = static_cast<int>(current_equivalence_list.length());
              int current_trial_j = j % c_eq_list_length;
              
              current_equivalence_vec = as<arma::vec>(current_equivalence_list(current_trial_j));
            } else {
              current_equivalence_vec = as<arma::vec>(equivalence_list(k));
            }
            arma::vec adjusted_pop_vec = running_popvec_mpm % current_equivalence_vec;
            aN_current = accu(adjusted_pop_vec);
          }
          
          agg_density_pre(i, j + 1) = agg_density_pre(i, j + 1) + aN_current;
          
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
  
  agg_density = agg_density_pre;
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
//' @param agg_density A numeric matrix giving the aggregate community density
//' by time (column) and replicate (row).
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
//' @param dens_list_length_vec A vector giving the number of data frames of
//' density info per MPM, if lists of such data frames are used.
//' @param eq_list_length_vec A vector giving the number of data frames of
//' equivalence info per MPM, if lists of such data frames are used.
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
void project3_fb_core (NumericMatrix& agg_density, List& N_out, List& comm_out,
  List& extreme_mpm_out, List& errcheck_fb_out, const List start_list,
  const List vrm_list, const List tweights_list, const List stageframe_list,
  const List allmodels_all, const List allstages_all, const List supplement_list,
  const List year_list, const List ind_terms_num_list, const List ind_terms_cat_list,
  const List dev_terms_list, const List density_vr_list, const List sp_density_list,
  const List density_list, const List dens_index_list, const List equivalence_list,
  const IntegerVector dev_terms_num_vec, const IntegerVector sp_density_num_vec,
  const IntegerVector firstage_vec, const IntegerVector finalage_vec,
  const IntegerVector stagecounts, const IntegerVector entry_time_vec,
  const IntegerVector format_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_yn_vec,
  const IntegerVector dens_vr_yn_vec, const IntegerVector dens_list_length_vec,
  const IntegerVector eq_list_length_vec, const IntegerVector tweights_type_vec,
  const NumericVector fecmod_vec, const LogicalVector sparse_vec,
  const CharacterVector patch_vec, const int vrm_count, const int nreps,
  const int times, const int substoch, const double exp_tol,
  const double theta_tol, const bool integeronly, const bool stages_not_equal,
  const bool stochastic, const bool err_check, const bool err_check_extreme) {
  
  // start_list
  // density_vr_list
  // dens_vr_yn_vec
  
  // Rcout << "Entered project3_fb_core          " << endl;
  
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
  
  // Rcout << "project3_fb_core A          " << endl;
  
  for (int i = 0; i < vrm_count; i++) {
    List pop_reps (nreps);
    
    arma::mat pops_out_pre (stagecounts(i), (times + 1), fill::zeros);
    
    arma::vec pops_out_0 = as<arma::vec>(start_list(i));
    pops_out_pre.col(entry_time_vec(i)) = pops_out_0; /////
    
    for (int j = 0; j < nreps; j++) {
      pop_reps(j) = pops_out_pre;
    }
    
    comm_out_pre(i) = pop_reps;
  }
  
  NumericMatrix agg_density_pre (nreps, (times + 1));
  List extreme_mpm_reps (nreps);
  List errcheck_mpmout_rep (nreps);
  
  // Main projection
  // Rcout << "project3_fb_core B          " << endl;
  
  for (int current_rep = 0; current_rep < nreps; current_rep++) {
    List running_popvecs = clone(start_list);
    NumericMatrix N_vrm (vrm_count, (times + 1));
    List extreme_mpm_reps_times (times);
    List errcheck_mpmout_rep_time (times); 
    int overall_dens_yn = sum(dens_yn_vec);
    
    IntegerVector year_counter (vrm_count);
    
    for (int current_time = 0; current_time < times; current_time++) {
      // Rcout << "project3_fb_core c current_time: " << current_time << "          ";
      
      Rcpp::checkUserInterrupt();
      
      List extreme_mpm_reps_times_vrms (vrm_count);
      List errcheck_mpmout_rep_time_vrm (vrm_count); 
      
      List dense_vecs (vrm_count);
      arma::vec delay_N_sum_vec (vrm_count);
      
      if (overall_dens_yn > 0) {
        for (int k = 0; k < vrm_count; k++) {
          if (current_time == 0) {
            arma::vec current_popvec = as<arma::vec>(running_popvecs(k));
            dense_vecs(k) = current_popvec;
            
            double delay_N_sum_vec_entry = accu(current_popvec);
            delay_N_sum_vec(k) = delay_N_sum_vec_entry;
          } else {
            if (dens_yn_vec(k) > 0) {
              DataFrame used_density_input;
              DataFrame used_density_index_input;
              
              if (dens_list_length_vec(k) > 0) {
                List current_used_times_mpm = as<List>(stochastic_year_vector(k));
                IntegerVector current_times_vec = as<IntegerVector>(current_used_times_mpm(current_rep));
                
                List used_density_input_list = as<List>(density_list(k));
                List used_density_index_input_list = as<List>(dens_index_list(k));
                
                used_density_input = as<DataFrame>(used_density_input_list(current_times_vec(current_time)));
                used_density_index_input = as<DataFrame>(used_density_index_input_list(current_times_vec(current_time)));
                
              } else {
                used_density_input = as<DataFrame>(density_list(k));
                used_density_index_input = as<DataFrame>(dens_index_list(k));
              }
              
              IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
              int used_delay = max(ud_delay_vec);
              
              if ((current_time + 1 - used_delay) > 0) {
                List reps_out = comm_out_pre(k);
                arma::mat pops_out = as<arma::mat>(reps_out(current_rep));
                arma::vec current_popvec = pops_out.col(current_time);
                dense_vecs(k) = current_popvec;
                
                double delay_N_sum_vec_entry = accu(current_popvec);
                delay_N_sum_vec(k) = delay_N_sum_vec_entry;
              } else {
                List new_start_list = clone(start_list);
                arma::vec current_popvec = as<arma::vec>(new_start_list(k));
                dense_vecs(k) = current_popvec;
                
                double delay_N_sum_vec_entry = accu(current_popvec);
                delay_N_sum_vec(k) = delay_N_sum_vec_entry;
              }
            }
          }
        }
      }
      
      for (int i = 0; i < vrm_count; i++) {
        List reps_out = comm_out_pre(i);
        arma::mat pops_out = as<arma::mat>(reps_out(current_rep));
        
        if (current_time > (entry_time_vec(i) - 1)) {
          arma::vec running_popvec_vrm = as<arma::vec>(running_popvecs(i));
          
          if (current_time == entry_time_vec(i)) {
            pops_out.col(current_time) = running_popvec_vrm;
            
            double N_current = accu(running_popvec_vrm);
            N_vrm(i, current_time) = N_current;
            
            double aN_first = accu(running_popvec_vrm);
            
            if (stages_not_equal) {
              arma::vec current_equivalence_vec;
              
              if (eq_list_length_vec(i) > 0) {
                List current_equivalence_list = as<List>(equivalence_list(i));
                current_equivalence_vec = as<arma::vec>(current_equivalence_list(current_time));
              } else {
                current_equivalence_vec = as<arma::vec>(equivalence_list(i));
              }
              arma::vec adjusted_pop_vec = running_popvec_vrm % current_equivalence_vec;
              aN_first = accu(adjusted_pop_vec);
            }
            
            agg_density_pre(current_rep, current_time) = agg_density_pre(current_rep, current_time) + aN_first;
          }
          
          List current_vrm_extract = as<List>(allmodels_all(i));
          List current_vrm_unextract = as<List>(vrm_list(i));
          DataFrame current_stageframe = as<DataFrame>(stageframe_list(i));
          
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
          
          // Rcout << "project3_fb_core B3          " << endl;
          
          CharacterVector current_mainyears = as<CharacterVector>(year_list(i));
          unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
          
          DataFrame group2_frame = as<DataFrame>(current_vrm_unextract["group2_frame"]);
          CharacterVector current_maingroups = as<CharacterVector>(group2_frame["groups"]);
          
          DataFrame patch_frame = as<DataFrame>(current_vrm_unextract["patch_frame"]);
          CharacterVector current_mainpatches = as<CharacterVector>(patch_frame["patches"]);
          
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
          IntegerVector pvis = index_l3(current_mainpatches, patch_vec(i));
          if (static_cast<int>(pvis.length()) == 0) {
            throw Rcpp::exception("Value for argument patch not found.", false);
          }
          int patchnumber = pvis(0);
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
          
          NumericVector dv_terms (17);
          if (dev_terms_num_vec(i) > 0) {
            DataFrame used_dv_df = as<DataFrame>(dev_terms_list(i));
            int used_dv_df_size = static_cast<int>(used_dv_df.size());
          
            if (dev_num_counter(i) >= dev_terms_num_vec(i)) dev_num_counter(i) = 0;
            
            for (int j = 0; j < used_dv_df_size; j++) {
              dv_terms(j) = used_dv_df(dev_num_counter(i), j);
            }
            dev_num_counter(i) = dev_num_counter(i) + 1;
          }
          
          bool dvr_bool {false};
          
          // Rcout << "project3_fb_core B5          " << endl;
          
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
          
          // Rcout << "project3_fb_core B6          " << endl;
          
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
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe, surv_proxy,
                fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc,
                r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, dvt3, dv_terms(0),
                dv_terms(6), dens_sp, fecmod_vec(i), finalage_vec(i), true,
                yearnumber, patchnumber, dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
                dens_n, exp_tol, theta_tol, sparse_vec(i));
              
            } else {
              current_supplement = as<DataFrame>(supplement_list(i));
              
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe, surv_proxy,
                fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc,
                r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, dvt3, dv_terms(0),
                dv_terms(6), dens_sp, fecmod_vec(i), finalage_vec(i), true,
                yearnumber, patchnumber, dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
                dens_n, exp_tol, theta_tol, sparse_vec(i), current_supplement);
            }
          }
          
          if (dens_yn_vec(i) == 0) {
            // Rcout << "project3_fb_core C no density       " << endl;
            
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
            // Rcout << "project3_fb_core D density       " << endl;
            
            DataFrame used_density_input;
            DataFrame used_density_index_input;
            
            if (dens_list_length_vec(i) > 0) {
              List used_density_input_list = as<List>(density_list(i));
              List used_density_index_input_list = as<List>(dens_index_list(i));
              
              int current_used_index = yearnumber;
              if (current_used_index >= static_cast<int>(dens_list_length_vec(i))) current_used_index = 0;
              used_density_input = as<DataFrame>(used_density_input_list(current_used_index));
              used_density_index_input = as<DataFrame>(used_density_index_input_list(current_used_index));
              
            } else {
              used_density_input = as<DataFrame>(density_list(i));
              used_density_index_input = as<DataFrame>(dens_index_list(i));
            }
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            if ((current_time + 1 - used_delay) >= 0) { // Changed to allow different delay Ns
              if (!stages_not_equal) {
                double delay_N_sum = accu(delay_N_sum_vec);
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_ad(new_popvec, new_projmat, new_projsp,
                  running_popvec_vrm, current_mpm, delay_N_sum, 0, integeronly,
                  substoch, used_density_input, used_density_index_input, false,
                  static_cast<bool>(sparse_vec(i)),
                  static_cast<bool>(sparse_vec(i)), false, err_check_extreme);
                
                running_popvec_vrm = new_popvec;
                if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = new_projmat;
              } else {
                double delay_N_sum {0.0};
                
                for (int l = 0; l < vrm_count; l++) {
                  arma::vec delay_pop_vec = as<arma::vec>(dense_vecs(l));
                  
                  arma::vec current_equiv_vec;
                  if (eq_list_length_vec(l) > 0) {
                    List used_equivalence_list = as<List>(equivalence_list(l));
                    int current_used_index = yearnumber;
                    if (current_used_index >= static_cast<int>(eq_list_length_vec(l))) current_used_index = 0;
                    current_equiv_vec = as<arma::vec>(used_equivalence_list(current_used_index));
                  
                  } else {
                    current_equiv_vec = as<arma::vec>(equivalence_list(l));
                  }
                  
                  arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                  double delay_pop_N = accu(adjusted_delay_pop_vec);
                  
                  delay_N_sum += delay_pop_N;
                }
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_ad(new_popvec, new_projmat, new_projsp,
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
              
              AdaptUtils::proj3dens_ad(new_popvec, new_projmat, new_projsp,
                running_popvec_vrm, current_mpm, 0.0, 0, integeronly, substoch,
                used_density_input, used_density_index_input, false,
                static_cast<bool>(sparse_vec(i)), static_cast<bool>(sparse_vec(i)),
                false, err_check_extreme);
              
              running_popvec_vrm = new_popvec;
              if (err_check_extreme) extreme_mpm_reps_times_vrms(i) = new_projmat;
            }
          }
          // Rcout << "project3_fb_core F       " << endl;
          
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
          
          double aN_current = N_current; // Aggregate density, given stage weights
          arma::vec current_equivalence_vec;
          
          if (stages_not_equal) {
            if (eq_list_length_vec(i) > 0) {
              List current_equivalence_list = as<List>(equivalence_list(i));
              current_equivalence_vec = as<arma::vec>(current_equivalence_list(current_time));
            } else {
              current_equivalence_vec = as<arma::vec>(equivalence_list(i));
            }
            arma::vec adjusted_pop_vec = running_popvec_vrm % current_equivalence_vec;
            aN_current = accu(adjusted_pop_vec);
          }
          
          agg_density_pre(current_rep, current_time + 1) = agg_density_pre(current_rep, current_time + 1) + aN_current;
        }
        year_counter(i) = year_counter(i) + 1;
        
        reps_out(current_rep) = pops_out;
        comm_out_pre(i) = reps_out;
      } // vrm loop
      
      // Rcout << "project3_fb_core G       " << endl;
      
      if (err_check_extreme) extreme_mpm_reps_times(current_time) = extreme_mpm_reps_times_vrms;
      if (err_check) errcheck_mpmout_rep_time(current_time) = errcheck_mpmout_rep_time_vrm;
    } // current_time loop
    comm_out = comm_out_pre;
    N_out_pre(current_rep) = N_vrm;
    if (err_check_extreme) extreme_mpm_reps(current_rep) = extreme_mpm_reps_times;
    if (err_check) errcheck_mpmout_rep(current_rep) = errcheck_mpmout_rep_time;
    
  } // current_rep loop
  
  agg_density = agg_density_pre;
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
//' that provide supplemental data that should be incorporated into MPMs. If
//' used, then should have the same number of elements as the number of MPMs
//' provided in the list for argument \code{vrms}, with MPMs that do not need
//' supplemental data entered as \code{NULL} elements in this list. See
//' \code{\link[lefko3]{supplemental}()} for details of how these data frames
//' are used in function-based MPM construction.
//' @param equivalence An optional numeric vector, list of numeric vectors,
//' data frame of class \code{adaptEq} or class \code{lefkoEq}, or a list of such
//' data frames. If a numeric vector, then must have the same number of elements
//' as the number of MPMs, with each element giving the effect of an individual
//' of each MPM relative to a reference individual. If a list of vectors, then
//' the list should be composed of as many numeric vectors as MPMs, or as a two-
//' layered list of such vectors for each annual matrix per MPM, with each
//' vector giving the effect of each individual in each stage relative to a
//' reference individual. Data frames of class \code{adaptEq}, and lists of such
//' data frames, can be made with function \code{\link{equiv_input}()}. Numeric
//' entries used in these vectors can be thought of as Lotka-Volterra
//' competition terms, such as are used in multiple species competition models.
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
//' \code{vrm_input} object. Each should include either 14 or 17 columns, and up
//' to \code{times} rows showing the values of the deviation terms to be added
//' to each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival, 9: juvenile
//' observation, 10: juvenile primary size, 11: juvenile secondary size,
//' 12: juvenile tertiary size, 13: juvenile reproduction, and 14: juvenile
//' maturity transition. In addition, users may supply 3 more columns
//' designating additive deviations to: 15: individual covariate a,
//' 16: individual covariate b, and 17: individual covariate c. Unused terms
//' must be set to \code{0} (use of \code{NA} will produce errors). Single or
//' small numbers of values per vital rate model are also allowed, and if the
//' number of rows is less than \code{times}, then the terms will be cycled.
//' @param fb_sparse A logical vector indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE} for
//' each MPM.
//' @param firstage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the starting ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Typically,
//' the starting age in such MPMs should be set to \code{0} if post-breeding and
//' \code{1} if pre-breeding. All other MPMs should be set to \code{0}. Do not
//' use if no MPM has age structure. Defaults to \code{1} in Leslie and
//' age-by-stage MPMs.
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
//' created with function \code{\link[lefko3]{density_input}()}, or a two-
//' layered list of such data frames, with a data frame per annual matrix per
//' MPM.
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
//' \item{agg_density}{A numeric matrix giving the aggregate community density
//' by time (column) and replicate (row). Here, the aggregate density is
//' calculated using any stage weights supplied, and is exactly equal to the
//' aggegate density used in density dependence calculations by this function.}
//' \item{N_out}{A list with the number of elements equal to the number of
//' replicates. Each element within this list is a numeric matrix showing the number
//' of individuals of each species or genotype alive at each time. The number of
//' rows are equal to the number of MPMs used, and the columns correspond to the
//' time steps.}
//' \item{structure}{A two-level list with the top level list having number of
//' elements equal to the number of replicates, and the lower level
//' corresponding to the number of MPMs. Each element of the lower level list is
//' a numeric matrix showing the number of individuals in each stage at each
//' time. Rows and columns in the data frames correspond to stages and time
//' steps, respectively.}
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
//' \item{err_check}{An optional list composed of six additional lists, each
//' of which has the number of elements equal to the number of MPMs utilized.
//' List output include \code{allstages_all}, which gives the indices of
//' estimated transitions in MPMs constructed by function \code{project3()} from
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
//' cyparaw_v1 <- verticalize3(data = cypa_data, noyears = 18, firstyear = 1994,
//'   individcol = "plant_id", blocksize = 2, sizeacol = "Inf.94",
//'   sizebcol = "Veg.94", repstracol = "Inf.94", fecacol = "Inf.94",
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
//'   year = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
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
//' summary(cyp_comm_proj)
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
  //int total_mpms {0};  // This includes all MPMs and VRMs
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
  
  //bool inda_char {false};
  //bool indb_char {false};
  //bool indc_char {false};
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
  List err_check_mat_indices; // matrix post-processing vectors giving matrix element numbers associated with changes
  List err_check_theoldpizzle_adapt3; // matrix post-processing indexing data frame
  
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
  IntegerVector dens_list_length_vec;
  IntegerVector eq_list_length_vec;
  List fb_mpmout;
  
  NumericMatrix agg_density; // A numeric matrix giving the aggregate community density by time (column) and replicate (row)
  List comm_out; // List of matrices of population vectors (top lvl: mpm, lower lvl: reps, mats: stages x times)
  List N_out;  // List of pop size matrices (top level: reps, mats: mpm by times)
  DataFrame labels;  // Data frame to provide order of MPMs
  List labels_list;  // List to hold data for data frame labels
  
  // Rcout << "project3 A " << endl;
  
  List cleaned_input = cleanup3(mpms, vrms, stageframes, supplements, format,
    firstage, finalage, fecage_min, fecage_max, cont, fecmod, starts, patches,
    years, tweights, density, entry_time, density_vr, sp_density, ind_terms,
    dev_terms, fb_sparse, equivalence, exp_tol, theta_tol, prep_mats, substoch,
    force_fb, stochastic, err_check_bool);
  
  // Rcout << "project3 B " << endl;
  
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
  dens_list_length_vec = as<IntegerVector>(cleaned_input(72));
  eq_list_length_vec = as<IntegerVector>(cleaned_input(73));
  err_check_mat_indices = as<List>(cleaned_input(74));
  err_check_theoldpizzle_adapt3 = as<List>(cleaned_input(75));
  
  // Rcout << "project3 C " << endl;
  
  List final_out_matrices;
  //total_mpms = mpm_count + vrm_count;
  
  // Projection runs
  if (preexisting || fb_override) { // Preexisting MPMs
    // Rcout << "project3 D pre_core" << endl;
    
    project3_pre_core (agg_density, N_out, comm_out, extreme_mpm_out, mpm_list,
      A_list, tweights_list, start_list, vrm_list, stageframe_list, allmodels_all,
      allstages_all, supplement_list, year_list, ind_terms_num_list,
      ind_terms_cat_list, dev_terms_list, density_vr_list, sp_density_list,
      density_list, dens_index_list, equivalence_list, dev_terms_num_vec,
      sp_density_num_vec, firstage_vec, finalage_vec, stagecounts,
      entry_time_vec, format_vec, inda_terms_num_vec, indb_terms_num_vec,
      indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
      indc_terms_cat_vec, dens_yn_vec, dens_vr_yn_vec, dens_list_length_vec,
      eq_list_length_vec, tweights_type_vec, fecmod_vec, sparse_vec, patch_vec,
      vrm_count, mpm_count, nreps, times, substoch, exp_tol, theta_tol,
      integeronly, stages_not_equal, stochastic,  entry_time_vec_use,
      err_check_bool, err_check_extreme);
    final_out_matrices = fb_mpmout;
    
  } else if (funcbased && !fb_override) { // Function-based MPMs
    // Rcout << "project3 E fb_core" << endl;
    
    project3_fb_core (agg_density, N_out, comm_out, extreme_mpm_out, err_check_fb_out,
      start_list, vrm_list, tweights_list, stageframe_list, allmodels_all,
      allstages_all, supplement_list, year_list, ind_terms_num_list,
      ind_terms_cat_list, dev_terms_list, density_vr_list, sp_density_list,
      density_list, dens_index_list, equivalence_list, dev_terms_num_vec,
      sp_density_num_vec, firstage_vec, finalage_vec, stagecounts,
      entry_time_vec, format_vec, inda_terms_num_vec, indb_terms_num_vec,
      indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
      indc_terms_cat_vec, dens_yn_vec, dens_vr_yn_vec, dens_list_length_vec,
      eq_list_length_vec, tweights_type_vec, fecmod_vec, sparse_vec, patch_vec,
      vrm_count, nreps, times, substoch, exp_tol, theta_tol, integeronly,
      stages_not_equal, stochastic, err_check_bool, err_check_extreme);
    final_out_matrices = err_check_fb_out;
  }
  
  // Restructure comm_out to structure
  List new_structure (nreps);
  int num_mpms = mpm_count;
  if (num_mpms < vrm_count) num_mpms = vrm_count;
  
  for (int i = 0; i < nreps; i++) {
    List lovely_mpms (num_mpms);
    
    for (int j = 0; j < num_mpms; j++) {
      List this_comm_out_list = as<List>(comm_out(j));
      NumericMatrix this_comm_out_matrix = as<NumericMatrix>(this_comm_out_list(i));
      
      lovely_mpms(j) = this_comm_out_matrix;
    }
    new_structure(i) = lovely_mpms;
  }
  
  // Rcout << "project3 F " << endl;
  
  int out_dim = 7;
  if (err_check_bool) out_dim++;
  
  List output (out_dim);
  
  output(0) = agg_density;
  output(1) = N_out; // Needed in final output
  output(2) = new_structure; // Needed in final output
  output(3) = stageframe_list; // Needed in final output
  output(4) = hstages_list;
  output(5) = agestages_list;
  output(6) = labels;
  
  if (err_check_extreme) {
    
    List output_errcheck (10);
    output_errcheck(0) = allstages_all;
    output_errcheck(1) = allmodels_all;
    output_errcheck(2) = equivalence_list;
    output_errcheck(3) = density_list;
    output_errcheck(4) = dens_index_list;
    output_errcheck(5) = density_vr_list;
    output_errcheck(6) = final_out_matrices;
    output_errcheck(7) = extreme_mpm_out;
    output_errcheck(8) = err_check_mat_indices;
    output_errcheck(9) = err_check_theoldpizzle_adapt3;
    
    CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
      "equivalence_list", "density_list", "dens_index_list", "density_vr_list",
      "fb_mpm_out_matrices", "modified_mpms", "post_processing_index_vectors",
      "post_processing_supp_frame"};
    output_errcheck.attr("names") = output_errcheck_names;
    
    output(7) = output_errcheck;
    
    CharacterVector output_main_names = {"agg_density", "N_out", "structure",
      "stageframe_list", "hstages_list", "agestages_list", "labels",
      "err_check"};
    output.attr("names") = output_main_names;
    
  } else if (err_check_bool) {
    List output_errcheck (9);
    output_errcheck(0) = allstages_all;
    output_errcheck(1) = allmodels_all;
    output_errcheck(2) = equivalence_list;
    output_errcheck(3) = density_list;
    output_errcheck(4) = dens_index_list;
    output_errcheck(5) = density_vr_list;
    output_errcheck(6) = final_out_matrices;
    output_errcheck(7) = err_check_mat_indices;
    output_errcheck(8) = err_check_theoldpizzle_adapt3;
    
    CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
      "equivalence_list", "density_list", "dens_index_list", "density_vr_list",
      "fb_mpm_out_matrices", "post_processing_index_vectors",
      "post_processing_supp_frame"};
    output_errcheck.attr("names") = output_errcheck_names;
    
    output(7) = output_errcheck;
    
    CharacterVector output_main_names = {"agg_density", "N_out", "structure",
      "stageframe_list", "hstages_list", "agestages_list", "labels",
      "err_check"};
    output.attr("names") = output_main_names;
    
  } else {
    CharacterVector output_main_names = {"agg_density", "N_out", "structure",
      "stageframe_list", "hstages_list", "agestages_list", "labels"};
    output.attr("names") = output_main_names;
  }
  output.attr("class") = "adaptProj";
  
  return output;
}

//' Project Multiple MPMs In Batches With Varying Starting Conditions
//' 
//' Function \code{batch_project3} runs function \code{project3()} across series
//' of differing beginning conditions. Currently used to run sensitivity
//' analyses.
//' 
//' @name batch_project3
//' 
//' @param used_mpms An integer vector detailing the MPMs entered in argument
//' \code{mpms} to run the batch projection on. Although all MPMs entered in
//' argument \code{mpms} will be projected, the alterations given in
//' \code{givenrate}, \code{offset}, or \code{multiplier} will only be performed
//' with those specified here. Defaults to \code{"all"}, which is equivalent to
//' an integer vector holding the sequence of 1 through the number of MPMs
//' entered.
//' @param givenrate A numeric vector giving values to replace targeted matrix
//' elements with. If used, then arguments \code{offset} and \code{multiplier}
//' cannot be used. Defaults to \code{NULL}.
//' @param offset A numeric vector giving the values to offset matrix elements
//' by. If used, then arguments \code{givenrate} and \code{multiplier} cannot be
//' used. Defaults to \code{c(0.005, 0.010, 0.015, 0.020, 0.025)}.
//' @param multiplier A numeric vector giving values to multiply targeted matrix
//' elements by. If used, then arguments \code{givenrate} and \code{offset}
//' cannot be used. Defaults to \code{NULL}.
//' @param all_elems A logical value indicating whether to use the alterations
//' specified in argument \code{givenrate}, \code{offset}, or \code{multiplier}
//' on all matrix elements. Defaults to \code{FALSE}.
//' @param quiet A logical value indicating whether to block messages during
//' the batch projection informing the user of which projection is currently
//' being performed, and how many total projections are planned. Defaults to
//' \code{TRUE}, which silences all such messages.
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
//' that provide supplemental data that should be incorporated into MPMs. If
//' used, then should have the same number of elements as the number of MPMs
//' provided in the list for argument \code{vrms}, with MPMs that do not need
//' supplemental data entered as \code{NULL} elements in this list. See
//' \code{\link[lefko3]{supplemental}()} for details of how these data frames
//' are used in function-based MPM construction.
//' @param equivalence An optional numeric vector, list of numeric vectors,
//' data frame of class \code{adaptEq} or class \code{lefkoEq}, or a list of such
//' data frames. If a numeric vector, then must have the same number of elements
//' as the number of MPMs, with each element giving the effect of an individual
//' of each MPM relative to a reference individual. If a list of vectors, then
//' the list should be composed of as many numeric vectors as MPMs, or as a two-
//' layered list of such vectors for each annual matrix per MPM, with each
//' vector giving the effect of each individual in each stage relative to a
//' reference individual. Data frames of class \code{adaptEq}, and lists of such
//' data frames, can be made with function \code{\link{equiv_input}()}. Numeric
//' entries used in these vectors can be thought of as Lotka-Volterra
//' competition terms, such as are used in multiple species competition models.
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
//' \code{vrm_input} object. Each should include either 14 or 17 columns, and up
//' to \code{times} rows showing the values of the deviation terms to be added
//' to each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival, 9: juvenile
//' observation, 10: juvenile primary size, 11: juvenile secondary size,
//' 12: juvenile tertiary size, 13: juvenile reproduction, and 14: juvenile
//' maturity transition. In addition, users may supply 3 more columns
//' designating additive deviations to: 15: individual covariate a,
//' 16: individual covariate b, and 17: individual covariate c. Unused terms
//' must be set to \code{0} (use of \code{NA} will produce errors). Single or
//' small numbers of values per vital rate model are also allowed, and if the
//' number of rows is less than \code{times}, then the terms will be cycled.
//' @param fb_sparse A logical vector indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE} for
//' each MPM.
//' @param firstage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the starting ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Typically,
//' the starting age in such MPMs should be set to \code{0} if post-breeding and
//' \code{1} if pre-breeding. All other MPMs should be set to \code{0}. Do not
//' use if no MPM has age structure. Defaults to \code{1} in Leslie and
//' age-by-stage MPMs.
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
//' created with function \code{\link[lefko3]{density_input}()}, or a two-
//' layered list of such data frames, with a data frame per annual matrix per
//' MPM.
//' @param density_vr An optional list of data frames of class
//' \code{lefkoDensVR}, which provide details for density dependence in vital
//' rate models and have been created with function
//' \code{link[lefko3]{density_vr}()}. If used, then one such data frame per MPM
//' is required. MPMs to be run without vital describing density dependence
//' relationships in vital rates should be set to \code{NULL}. Can only be used
//' with function-based projections.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking with each projection. Can also be set
//' to the text value \code{"extreme"}, in which case all \code{err_check}
//' output plus a multiple level list with each MPM used in each time step will
//' be output.
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
//' @return A list of class \code{adaptProjBatch}, which is composed of three
//' elements:
//' \item{proj_out}{A multi-layered list, with the top number of elements equal
//' to the number of \code{used_mpms}. Each element is a list, with the number
//' of elements equal to the number of alterations tried. Each of these elements
//' is, in turn, an \code{adaptProj} object.}
//' \item{ref}{A list with the same structure as \code{proj_out}. This list is
//' composed of supplements (\code{lefkoSD} objects) giving the exact
//' alterations performed in the associated MCM.}
//' \item{control}{An integer vector giving the identities of the MPM altered
//' in each top list.}
//' 
//' @section Notes:
//' 
//' This function is currently used to run sensitivity analyses, in which the
//' slope of the final population size is regressed against the the starting
//' offset values for each element. Each regression is run on the offset values
//' within each element.
//' 
//' The only arguments unique to this function are arguments \code{used_mpms},
//' \code{givenrate}, \code{offset}, \code{multiplier}, \code{all_elems}, and
//' \code{quiet}. All others are passed to function \code{project3}.
//' 
//' Setting \code{all_elems = TRUE} will lead to very time-intensive analysis.
//' We encourage users to break down their analyses into smaller batches, in 
//' order to make them more tractable.
//' 
//' @examples
//' library(lefko3)
//' data(cypdata)
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
//' sizevector <- c(0, 0, 3.0, 15)
//' stagevector <- c("P1", "D", "Sm", "Lg")
//' repvector <- c(0, 0, 1, 1)
//' obsvector <- c(0, 0, 1, 1)
//' matvector <- c(0, 1, 1, 1)
//' immvector <- c(1, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1)
//' binvec <- c(0, 0.5, 2.5, 9.5)
//' 
//' cypframe_small_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypraw_v2 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_small_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypraw_v3 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   NAas0 = TRUE, NRasRep = TRUE)
//' 
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 1500, 500),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' cypmean <- lmean(cypmatrix2r)
//' 
//' cypsupp2r_small <- supplemental(stage3 = c("D", "Sm", "Lg", "P1"),
//'   stage2 = c("P1", "P1", "P1", "rep"), eststage3 = c(NA, "Sm", "Lg", NA),
//'   eststage2 = c(NA, "D", "D", NA), givenrate = c(0.05, NA, NA, NA),
//'   offset = c(NA, NA, -0.1, NA), multiplier = c(NA, NA, NA, 0.5),
//'   type =c(1, 1, 1, 3), stageframe = cypframe_small_raw, historical = FALSE)
//' cypmatrix2r_small <- rlefko2(data = cypraw_v2, stageframe = cypframe_small_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r_small,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' cypmean_small <- lmean(cypmatrix2r_small)
//' 
//' cypmatrixL_small <- rleslie(data = cypraw_v3, start_age = 1, last_age = 4,
//'   continue = TRUE, fecage_min = 3, year = "all", pop = NA, patch = "all",
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cyp_mpms1 <- list(cypmatrix2r, cypmatrix2r_small, cypmatrixL_small)
//' 
//' c2d_4 <- density_input(cypmean, stage3 = c("P1", "P1"), stage2= c("SD", "rep"),
//'   style = 1, time_delay = 1, alpha = 1, beta = 0.0005, type = c(2, 2))
//' c2d_4a <- density_input(cypmean_small, stage3 = c("P1", "P1"), stage2= c("P1", "rep"),
//'   style = 1, time_delay = 1, alpha = 1, beta = 0.0005, type = c(2, 2))
//' cypL_dv <- density_input(cypmatrixL_small, stage3 = c("Age1"), stage2 = c("rep"),
//'   style = c(1), alpha = c(0.5), beta = c(1.0), type = c(2))
//' cyp_density <- list(c2d_4, c2d_4a, cypL_dv)
//' 
//' cyp_start1 <- start_input(cypmatrix2r, stage2 = c("SD", "P1", "D"),
//'   value = c(100, 200, 4))
//' cyp_start2 <- start_input(cypmatrix2r_small, stage2 = c("P1", "D"),
//'   value = c(10, 2000))
//' cypL_start_1 <- start_input(cypmatrixL_small, stage2 = c("Age1"),
//'   value = c(200))
//' cyp_start <- list(cyp_start1, cyp_start2, cypL_start_1)
//' 
//' new_supplement_cyp2_small <- sup_skeleton(2)
//' new_supplement_cyp2_small$stage3 <- c("D", "Sm")
//' new_supplement_cyp2_small$stage2 <- c("Lg", "Lg")
//' new_supplement_cyp2_small$convtype <- c(1, 1)
//' used_supplements <- list(new_supplement_cyp2_small,
//'   new_supplement_cyp2_small, NULL)
//' 
//' aaa1_prj_batch2 <- batch_project3(used_mpms = "all", all_elems = FALSE,
//'   mpms =  cyp_mpms1, entry_time = c(0, 5, 8), times = 15, nreps = 3,
//'   supplement = used_supplements, integeronly = TRUE, density = cyp_density)
//' 
//' @export batch_project3
// [[Rcpp::export(batch_project3)]]
Rcpp::List batch_project3 (Nullable<RObject> used_mpms = R_NilValue,
  Nullable<RObject> givenrate = R_NilValue, Nullable<RObject> offset = R_NilValue,
  Nullable<RObject> multiplier = R_NilValue, Nullable<RObject> all_elems = R_NilValue,
  Nullable<RObject> quiet = R_NilValue, Nullable<RObject> mpms = R_NilValue,
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
  Nullable<RObject> err_check = R_NilValue, bool stochastic = false,
  bool integeronly = false, int substoch = 0, int nreps = 1, int times = 10000,
  int prep_mats = 20, bool force_fb = false, double exp_tol = 700.0,
  double theta_tol = 100000000.0) {
  
  // Rcout << "Entered batch_project3" << endl;
  
  int mpms_entered {0};
  int trial_alterations {0};
  bool sensitivity_bool {false};
  bool all_elems_bool {false};
  bool quiet_bool {true};
  bool using_vrms {false};
  bool using_supplements {false};
  
  bool givenrate_used {false};
  bool offset_used {false};
  bool multiplier_used {false};
  bool err_check_bool {true};
  
  IntegerVector used_mpm_vec;
  IntegerVector format_vec;
  NumericVector givenrate_vec;
  NumericVector offset_vec;
  NumericVector multiplier_vec;
  
  if (mpms.isNotNull()) {
    RObject mpms_true = RObject(mpms);
    
    if (is<List>(mpms_true)) {
      List mpms_list = as<List>(mpms_true);
      mpms_entered = static_cast<int>(mpms_list.length());
      
      for (int i = 0; i < mpms_entered; i++) {
        if (!is<List>(mpms_list(i))) throw Rcpp::exception("Argument mpms should be a list of lefkoMat objects.", false);
        
        List current_mpm = as<List>(mpms_list(i));
        CharacterVector current_mpm_class = as<CharacterVector>(current_mpm.attr("class"));
        int num_of_classes = static_cast<int>(current_mpm_class.length());
        
        bool found_lefkoMat {false};
        for (int j = 0; j < num_of_classes; j++) {
          if (current_mpm_class(j) == "lefkoMat") found_lefkoMat = true;
        }
        
        if (!found_lefkoMat) throw Rcpp::exception("Argument mpms should be a list of lefkoMat objects.", false);
      }
    }
  } else if (vrms.isNotNull()) {
    using_vrms = true;
    
    RObject vrms_true = RObject(vrms);
    
    if (is<List>(vrms_true)) {
      List vrms_list = as<List>(vrms_true);
      mpms_entered = static_cast<int>(vrms_list.length());
      
      for (int i = 0; i < mpms_entered; i++) {
        if (!is<List>(vrms_list(i))) throw Rcpp::exception("Argument vrms should be a list of vrm_input objects.", false);
        
        List current_vrm = as<List>(vrms_list(i));
        CharacterVector current_vrm_class = as<CharacterVector>(current_vrm.attr("class"));
        int num_of_classes = static_cast<int>(current_vrm_class.length());
        
        bool found_vrm_input {false};
        for (int j = 0; j < num_of_classes; j++) {
          if (current_vrm_class(j) == "vrm_input") found_vrm_input = true;
        }
        
        if (!found_vrm_input) throw Rcpp::exception("Argument vrms should be a list of vrm_input objects.", false);
      }
    }
  } else {
    throw Rcpp::exception("Batch projection requires MPMs supplied through either argument mpms or argument vrms.", false);
  }
  
  if (supplements.isNotNull()) {
    RObject supplements_RO = RObject(supplements);
    List supplements_list = as<List>(supplements_RO);
    int supplements_list_length = static_cast<int>(supplements_list.length());
    
    for (int i = 0; i < supplements_list_length; i++) {
      Nullable<RObject> current_supplement = as<Nullable<RObject>>(supplements_list(i));
      
      if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
        List current_supplement_list = as<List>(current_supplement);
        CharacterVector current_supplement_class = as<CharacterVector>(current_supplement_list.attr("class"));
        int current_supplement_class_length = static_cast<int>(current_supplement_class.length());
        
        for (int j = 0; j < current_supplement_class_length; j++) {
          if (current_supplement_class(j) == "lefkoSD") using_supplements = true;
        }
        
      } else if (current_supplement.isNotNull()) {
        throw Rcpp::exception("Some supplements are not in a recognized format.", false);
      } //else throw Rcpp::exception("Some input supplements are not in a recognized format.", false);
    }
  }
  
  if (used_mpms.isNotNull()) {
    RObject used_mpms_RO = RObject(used_mpms);
    
    if (is<NumericVector>(used_mpms_RO) || is<IntegerVector>(used_mpms_RO)) {
      IntegerVector used_mpms_iv = as<IntegerVector>(used_mpms_RO);
      int used_mpms_iv_length = static_cast<int>(used_mpms_iv.length());
      
      for (int i = 0; i < used_mpms_iv_length; i++) {
        if (used_mpms_iv(i) < 1 || used_mpms_iv(i) > mpms_entered) {
          throw Rcpp::exception("Use only integers from 1 to the number of MPMs entered in vector used_mpms.", false);
        }
      }
      
      used_mpm_vec = used_mpms_iv;
      
    } else if (is<CharacterVector>(used_mpms_RO)) {
      CharacterVector used_mpms_cv = as<CharacterVector>(used_mpms_RO);
      
      if (LefkoInputs::stringcompare_input(String(used_mpms_cv(0)), "all")) {
        used_mpm_vec = seq(1, mpms_entered);
      } else {
        throw Rcpp::exception("Batch projection requires MPMs supplied through either argument mpms or argument vrms.", false);
      }
    } else {
      throw Rcpp::exception("Value entered in argument used_mpms is not recognized.", false);
    }
  } else {
    used_mpm_vec = seq(1, mpms_entered);
  }
  
  if (givenrate.isNotNull()) {
    RObject givenrate_RO = RObject(givenrate);
    
    if (is<NumericVector>(givenrate_RO) || is<IntegerVector>(givenrate_RO)) {
      givenrate_vec = as<NumericVector>(givenrate_RO);
      givenrate_used = true;
      trial_alterations = static_cast<int>(givenrate_vec.length());
    } else {
      throw Rcpp::exception("Argument givenrate must be a numeric vector.", false);
    }
  }  
  
  if (multiplier.isNotNull()) {
    RObject multiplier_RO = RObject(multiplier);
    
    if (is<NumericVector>(multiplier_RO) || is<IntegerVector>(multiplier_RO)) {
      multiplier_vec = as<NumericVector>(multiplier_RO);
      multiplier_used = true;
      trial_alterations = static_cast<int>(multiplier_vec.length());
    } else {
      throw Rcpp::exception("Argument multiplier must be a numeric vector.", false);
    }
  }  
  
  if (offset.isNotNull()) {
    RObject offset_RO = RObject(offset);
    
    if (is<NumericVector>(offset_RO) || is<IntegerVector>(offset_RO)) {
      offset_vec = as<NumericVector>(offset_RO);
      offset_used = true;
      trial_alterations = static_cast<int>(offset_vec.length());
    } else {
      throw Rcpp::exception("Argument offset must be a numeric vector.", false);
    }
  } else {
    if (!multiplier_used && !givenrate_used) {
      offset_vec = {0.005, 0.010, 0.015, 0.020, 0.025};
      offset_used = true;
      trial_alterations = 5;
    }
  }
  
  if (!multiplier_used && !givenrate_used && !offset_used) {
    throw Rcpp::exception ("No batch alterations input.", false);
  }
  
  if ((multiplier_used && givenrate_used) || (offset_used && givenrate_used) || (multiplier_used && offset_used)) {
    throw Rcpp::exception ("Batch projection currently handles multipliers, offset, or givenrates, but not in combination.",     
      false);
  }
  
  // Probably need to get rid of this one
  if (all_elems.isNotNull()) {
    RObject all_elems_true = RObject(all_elems);
    all_elems_bool = LefkoInputs::yesno_to_logic (all_elems_true, "all_elems");
  }
  
  // Probably need to get rid of this one
  if (quiet.isNotNull()) {
    RObject quiet_true = RObject(quiet);
    quiet_bool = LefkoInputs::yesno_to_logic (quiet_true, "quiet");
  }
  
  if (format.isNotNull()) {
    RObject format_RO = RObject(format);
    
    if (is<NumericVector>(format_RO) || is<IntegerVector>(format_RO)) {
      format_vec = as<IntegerVector>(format_RO);
      
      int format_vec_length = static_cast<int>(format_vec.length());
      if (format_vec_length != mpms_entered) {
        throw Rcpp::exception("Argument format must be an integer vector with as many elements as MPMs used.", false);
      }
      
      for (int i = 0; i < format_vec_length; i++) {
        if (format_vec(i) < 1 || format_vec(i) > 5) {
          throw Rcpp::exception("Vector format must be composed of integers from 1 through 5.", false);
        }
        
        if (!using_vrms) {
          RObject mpms_true = RObject(mpms);
          List mpms_list = as<List>(mpms_true);
          
          for (int i = 0; i < mpms_entered; i++) {
            List current_mpm = as<List>(mpms_list(i));
            int current_format = LefkoInputs::format_check_lM (current_mpm);
            
            if (current_format != format_vec(i)) {
              throw Rcpp::exception("Formats entered in argument format do not match entered MPMs.", false);
            }
          }
        }
      }
    } else {
      throw Rcpp::exception("Argument format must be an integer vector with as many elements as MPMs used.", false);
    }
  } else {
    if (!using_vrms) {
      RObject mpms_true = RObject(mpms);
      List mpms_list = as<List>(mpms_true);
      
      IntegerVector format_vec_pre (mpms_entered);
      
      for (int i = 0; i < mpms_entered; i++) {
        List current_mpm = as<List>(mpms_list(i));
        int current_format = LefkoInputs::format_check_lM (current_mpm);
        format_vec_pre(i) = current_format;
      }
      
      format_vec = format_vec_pre;
    } else {
      IntegerVector format_vec_pre (mpms_entered, 3);
      format_vec = format_vec_pre;
    }
  }
  
  if (!all_elems_bool && !using_supplements) {
    throw Rcpp::exception("Argument supplements is required to hold at least 1 lefkoSD object if all_elems = FALSE.", false);
  }
  // Check inputs
  // Rcout << "mpms_entered: " << mpms_entered << endl;
  // Rcout << "sensitivity_bool: " << sensitivity_bool << endl;
  // Rcout << "all_elems_bool: " << all_elems_bool << endl;
  // Rcout << "quiet_bool: " << quiet_bool << endl;
  // Rcout << "using_vrms: " << using_vrms << endl;
  // Rcout << "using_supplements: " << using_supplements << endl;
  // Rcout << "givenrate_used: " << givenrate_used << endl;
  // Rcout << "offset_used: " << offset_used << endl;
  // Rcout << "multiplier_used: " << multiplier_used << endl;
  // Rcout << "used_mpm_vec: " << used_mpm_vec << endl;
  // Rcout << "offset_vec: " << offset_vec << endl;
  // Rcout << "format_vec: " << format_vec << endl;
  
  List projections_list (mpms_entered);
  List err_check_list (mpms_entered);
  int overall_projection_tracker {0};
  int planned_projections {0};
  
  // Calculate the total number of projections
  for (int i = 0; i < mpms_entered; i++) {
    // Rcout << "batch_project3 A mpm" << endl;
    
    int current_mpm_int = (used_mpm_vec(i) - 1);
    int numstages {0};
    DataFrame current_stageframe;
    CharacterVector current_stageframe_stage;
    
    if (using_vrms) {
      if (stageframes.isNotNull()) {
        RObject stageframes_true = RObject(stageframes);
        List stageframes_list = as<List>(stageframes);
        current_stageframe = as<DataFrame>(stageframes_list(current_mpm_int));
      }
      
    } else {
      RObject mpms_true = RObject(mpms);
      List mpms_list = as<List>(mpms_true);
      List current_mpm = as<List>(mpms_list(current_mpm_int));
      current_stageframe = as<DataFrame>(current_mpm["ahstages"]);
    }
    
    // Rcout << "batch_project3 B" << endl;
    
    if (static_cast<int>(current_stageframe.nrows()) > 1) {
      current_stageframe_stage = as<CharacterVector>(current_stageframe["stage"]);
      numstages = static_cast<int>(current_stageframe_stage.length());
      
      if (all_elems_bool) {
        planned_projections = planned_projections + (numstages * numstages * trial_alterations);
      } else if (using_supplements) {
        RObject supplements_RO = RObject(supplements);
        List supplements_list = as<List>(supplements_RO);
        DataFrame current_supplement = as<DataFrame>(supplements_list(i));
        
        if (static_cast<int>(current_supplement.length()) > 1) planned_projections = planned_projections + trial_alterations;
      }
    }
  }
  
  // Rcout << "batch_project3 C" << endl;
  
  // Main projection routine
  for (int i = 0; i < mpms_entered; i++) {
    // Rcout << "batch_project3 D mpm " << i << " started." << endl;
    
    int current_mpm_int = (used_mpm_vec(i) - 1);
    int numstages {0};
    
    DataFrame current_stageframe;
    CharacterVector current_stageframe_stage;
    
    if (using_vrms) {
      if (stageframes.isNotNull()) {
        RObject stageframes_true = RObject(stageframes);
        List stageframes_list = as<List>(stageframes);
        current_stageframe = clone(as<DataFrame>(stageframes_list(current_mpm_int)));
      }
      
    } else {
      RObject mpms_true = RObject(mpms);
      List mpms_list = as<List>(mpms_true);
      List current_mpm = as<List>(mpms_list(current_mpm_int));
      current_stageframe = clone(as<DataFrame>(current_mpm["ahstages"]));
    }
    
    if (static_cast<int>(current_stageframe.length()) > 1) {
      current_stageframe_stage = as<CharacterVector>(current_stageframe["stage"]);
      numstages = static_cast<int>(current_stageframe_stage.length());
      
      if (format_vec(current_mpm_int) == 2) {
        bool found_AlmostBorn {false};
        
        for (int j = 0; j < numstages; j++) {
          if (current_stageframe_stage(j) == "AlmostBorn") found_AlmostBorn = true;
        }
      
        if (!found_AlmostBorn) {
          current_stageframe_stage.push_back("AlmostBorn");
          numstages++;
        }
      }
    }
    
    Nullable<RObject> current_supplement;
    if (using_supplements) {
      RObject supplements_RO = RObject(supplements);
      List supplements_list = as<List>(supplements_RO);
      current_supplement = as<Nullable<RObject>>(supplements_list(current_mpm_int));
    } else current_supplement = R_NilValue;
    
    if (all_elems_bool) {
      if (format_vec(current_mpm_int) == 1) {
        // Rcout << "batch_project3 F historical Ehrlen" << endl;
        
        int current_planned_projections = {numstages * numstages * numstages * trial_alterations};
        List projections_list_current_mpm_supp (current_planned_projections);
        List err_check_list_current_mpm_supp (current_planned_projections);
          
        int current_tracker {0};
        for (int stage1 = 0; stage1 < numstages; stage1++) {
          for (int stage2 = 0; stage2 < numstages; stage2++) {
            for (int stage3 = 0; stage3 < numstages; stage3++) {
              for (int try1 = 0; try1 < trial_alterations; try1++) {
                Rcpp::checkUserInterrupt();
                if (!quiet_bool) {
                  Rcout << "Currently on projection " << (current_tracker + 1) <<
                    " out of " << current_planned_projections << " for Ehrlen-formatted historical mpm " <<
                    (current_mpm_int + 1) << "." << endl;
                }
                
                DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
                CharacterVector new_stage3 = {current_stageframe_stage(stage3)};
                CharacterVector new_stage2 = {current_stageframe_stage(stage2)};
                CharacterVector new_stage1 = {current_stageframe_stage(stage1)};
                IntegerVector new_convtype = {1};
                
                current_skeleton["stage3"] = new_stage3;
                current_skeleton["stage2"] = new_stage2;
                current_skeleton["stage1"] = new_stage1;
                current_skeleton["convtype"] = new_convtype;
                
                if (givenrate_used) {
                  NumericVector new_altered = {givenrate_vec(try1)};
                  current_skeleton["givenrate"] = new_altered;
                } else if (offset_used) {
                  NumericVector new_altered = {offset_vec(try1)};
                  current_skeleton["offset"] = new_altered;
                } else if (multiplier_used) {
                  NumericVector new_altered = {multiplier_vec(try1)};
                  current_skeleton["multiplier"] = new_altered;
                }
                
                DataFrame new_supplement;
                if (current_supplement.isNotNull()) {
                  if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                    DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                    DataFrame used_supplement = clone(current_supplement_base);
                    
                    new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                  } else new_supplement = clone(current_skeleton);
                } else new_supplement = clone(current_skeleton);
                
                if (err_check_bool) {
                  //DataFrame new_skeleton = clone(new_supplement);
                  DataFrame err_supplement = clone(new_supplement);
                  IntegerVector mpm_id = {used_mpm_vec(i)};
                  err_supplement.push_front(mpm_id);
                  
                  CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                  err_supp_names(0) = "mpm_id";
                  err_supplement.attr("names") = err_supp_names;
                  err_check_list_current_mpm_supp(current_tracker) = err_supplement;
                }
                
                Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                  new_supplement, current_mpm_int, mpms_entered);
                
                List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                  equivalence, starts, years, patches, tweights, format, entry_time,
                  sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                  fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                  stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                  exp_tol, theta_tol);
                
                projections_list_current_mpm_supp(current_tracker) = current_projection;
                
                current_tracker++;
                overall_projection_tracker++;
              }
            }
          }
        }
        projections_list(i) = projections_list_current_mpm_supp;
        err_check_list(i) = err_check_list_current_mpm_supp;
        
      } else if (format_vec(current_mpm_int) == 2) {
        // Rcout << "batch_project3 G historical deVries" << endl;
        
        int current_planned_projections = {(((numstages - 1) * (numstages - 1) * numstages) + 
          ((numstages - 1) * (numstages - 1) * (numstages - 1))) * trial_alterations};
        List projections_list_current_mpm_supp (current_planned_projections);
        List err_check_list_current_mpm_supp (current_planned_projections);
          
        int current_tracker {0};
        
        // First all survival-focused runs
        for (int stage1 = 0; stage1 < numstages; stage1++) {
          for (int stage2 = 0; stage2 < (numstages - 1); stage2++) {
            for (int stage3 = 0; stage3 < (numstages - 1); stage3++) {
              for (int try1 = 0; try1 < trial_alterations; try1++) {
                Rcpp::checkUserInterrupt();
                if (!quiet_bool) {
                  Rcout << "Currently on projection " << (current_tracker + 1) <<
                    " out of " << current_planned_projections << " for deVries-formatted historical mpm " <<
                    (current_mpm_int + 1) << "." << endl;
                }
                
                DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
                CharacterVector new_stage3 = {current_stageframe_stage(stage3)};
                CharacterVector new_stage2 = {current_stageframe_stage(stage2)};
                CharacterVector new_stage1 = {current_stageframe_stage(stage1)};
                IntegerVector new_convtype = {1};
                
                current_skeleton["stage3"] = new_stage3;
                current_skeleton["stage2"] = new_stage2;
                current_skeleton["stage1"] = new_stage1;
                current_skeleton["convtype"] = new_convtype;
                
                if (givenrate_used) {
                  NumericVector new_altered = {givenrate_vec(try1)};
                  current_skeleton["givenrate"] = new_altered;
                } else if (offset_used) {
                  NumericVector new_altered = {offset_vec(try1)};
                  current_skeleton["offset"] = new_altered;
                } else if (multiplier_used) {
                  NumericVector new_altered = {multiplier_vec(try1)};
                  current_skeleton["multiplier"] = new_altered;
                }
                
                DataFrame new_supplement;
                if (current_supplement.isNotNull()) {
                  if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                    DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                    DataFrame used_supplement = clone(current_supplement_base);
                    
                    new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                  } else new_supplement = clone(current_skeleton);
                } else new_supplement = clone(current_skeleton);
                
                if (err_check_bool) {
                  //DataFrame new_skeleton = clone(new_supplement);
                  DataFrame err_supplement = clone(new_supplement);
                  IntegerVector mpm_id = {used_mpm_vec(i)};
                  err_supplement.push_front(mpm_id);
                  
                  CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                  err_supp_names(0) = "mpm_id";
                  err_supplement.attr("names") = err_supp_names;
                  err_check_list_current_mpm_supp(current_tracker) = err_supplement;
                }
                
                Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                  new_supplement, current_mpm_int, mpms_entered);
                
                List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                  equivalence, starts, years, patches, tweights, format, entry_time,
                  sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                  fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                  stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                  exp_tol, theta_tol);
                
                projections_list_current_mpm_supp(current_tracker) = current_projection;
                
                current_tracker++;
                overall_projection_tracker++;
              }
            }
          }
        }
        
        // Now all reproductive runs
        for (int stage1 = 0; stage1 < (numstages - 1); stage1++) {
          for (int stage2 = 0; stage2 < (numstages - 1); stage2++) {
            for (int stage3 = 0; stage3 < (numstages - 1); stage3++) {
              for (int try1 = 0; try1 < trial_alterations; try1++) {
                Rcpp::checkUserInterrupt();
                if (!quiet_bool) {
                  Rcout << "Currently on projection " << (current_tracker + 1) <<
                    " out of " << current_planned_projections << " for deVries-formatted historical mpm " <<
                    (current_mpm_int + 1) << "." << endl;
                }
                
                DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
                CharacterVector new_stage3 = {current_stageframe_stage(stage3)};
                CharacterVector new_stage2 = {current_stageframe_stage(stage2)};
                CharacterVector new_stage1 = {current_stageframe_stage(stage1)};
                IntegerVector new_convtype = {2};
                
                current_skeleton["stage3"] = new_stage3;
                current_skeleton["stage2"] = new_stage2;
                current_skeleton["stage1"] = new_stage1;
                current_skeleton["convtype"] = new_convtype;
                
                if (givenrate_used) {
                  NumericVector new_altered = {givenrate_vec(try1)};
                  current_skeleton["givenrate"] = new_altered;
                } else if (offset_used) {
                  NumericVector new_altered = {offset_vec(try1)};
                  current_skeleton["offset"] = new_altered;
                } else if (multiplier_used) {
                  NumericVector new_altered = {multiplier_vec(try1)};
                  current_skeleton["multiplier"] = new_altered;
                }
                
                DataFrame new_supplement;
                if (current_supplement.isNotNull()) {
                  if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                    DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                    DataFrame used_supplement = clone(current_supplement_base);
                    
                    new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                  } else new_supplement = clone(current_skeleton);
                } else new_supplement = clone(current_skeleton);
                
                if (err_check_bool) {
                  //DataFrame new_skeleton = clone(new_supplement);
                  DataFrame err_supplement = clone(new_supplement);
                  IntegerVector mpm_id = {used_mpm_vec(i)};
                  err_supplement.push_front(mpm_id);
                  
                  CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                  err_supp_names(0) = "mpm_id";
                  err_supplement.attr("names") = err_supp_names;
                  err_check_list_current_mpm_supp(current_tracker) = err_supplement;
                }
                
                Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                  new_supplement, current_mpm_int, mpms_entered);
                
                List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                  equivalence, starts, years, patches, tweights, format, entry_time,
                  sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                  fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                  stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                  exp_tol, theta_tol);
                
                projections_list_current_mpm_supp(current_tracker) = current_projection;
                
                current_tracker++;
                overall_projection_tracker++;
              }
            }
          }
        }
        projections_list(i) = projections_list_current_mpm_supp;
        err_check_list(i) = err_check_list_current_mpm_supp;
        
      } else if (format_vec(current_mpm_int) == 3) {
        // Rcout << "batch_project3 H ahistorical" << endl;
        
        int current_planned_projections = {numstages * numstages * trial_alterations};
        List projections_list_current_mpm_supp (current_planned_projections);
        List err_check_list_current_mpm_supp (current_planned_projections);
          
        int current_tracker {0};
        for (int stage2 = 0; stage2 < numstages; stage2++) {
          for (int stage3 = 0; stage3 < numstages; stage3++) {
            for (int try1 = 0; try1 < trial_alterations; try1++) {
              Rcpp::checkUserInterrupt();
              if (!quiet_bool) {
                Rcout << "Currently on projection " << (current_tracker + 1) <<
                  " out of " << current_planned_projections << " for ahistorical stage-based mpm " <<
                  (current_mpm_int + 1) << "." << endl;
              }
              
              DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
              CharacterVector new_stage3 = {current_stageframe_stage(stage3)};
              CharacterVector new_stage2 = {current_stageframe_stage(stage2)};
              IntegerVector new_convtype = {1};
              
              current_skeleton["stage3"] = new_stage3;
              current_skeleton["stage2"] = new_stage2;
              current_skeleton["convtype"] = new_convtype;
              
              if (givenrate_used) {
                NumericVector new_altered = {givenrate_vec(try1)};
                current_skeleton["givenrate"] = new_altered;
              } else if (offset_used) {
                NumericVector new_altered = {offset_vec(try1)};
                current_skeleton["offset"] = new_altered;
              } else if (multiplier_used) {
                NumericVector new_altered = {multiplier_vec(try1)};
                current_skeleton["multiplier"] = new_altered;
              }
              
              DataFrame new_supplement;
              if (current_supplement.isNotNull()) {
                if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                  DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                  DataFrame used_supplement = clone(current_supplement_base);
                  
                  new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                } else new_supplement = clone(current_skeleton);
              } else new_supplement = clone(current_skeleton);
              
              if (err_check_bool) {
                //DataFrame new_skeleton = clone(new_supplement);
                DataFrame err_supplement = clone(new_supplement);
                IntegerVector mpm_id = {used_mpm_vec(i)};
                err_supplement.push_front(mpm_id);
                
                CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                err_supp_names(0) = "mpm_id";
                err_supplement.attr("names") = err_supp_names;
                err_check_list_current_mpm_supp(current_tracker) = err_supplement;
              }
              
              Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                new_supplement, current_mpm_int, mpms_entered);
              
              List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                equivalence, starts, years, patches, tweights, format, entry_time,
                sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                exp_tol, theta_tol);
              
              projections_list_current_mpm_supp(current_tracker) = current_projection;
              
              current_tracker++;
              overall_projection_tracker++;
            }
          }
        }
        projections_list(i) = projections_list_current_mpm_supp;
        err_check_list(i) = err_check_list_current_mpm_supp;
        
        
      } else if (format_vec(current_mpm_int) == 4) {
        // Rcout << "batch_project3 I age-by-stage" << endl;
        
        // Need to set up a test if ahstages is NULL, to create a pseudo-ahstages
        int firstage_int {0};
        int finalage_int {0};
        int fecage_min_int {0};
        int fecage_max_int {0};
        bool cont_bool {false};
        
        arma::ivec current_stageframe_repstatus_arma = as<arma::ivec>(current_stageframe["repstatus"]);
        arma::ivec current_stageframe_entry_arma = as<arma::ivec>(current_stageframe["entrystage"]);
        arma::ivec current_stageframe_min_age_arma = as<arma::ivec>(current_stageframe["min_age"]);
        
        arma::uvec rep_stages = find(current_stageframe_repstatus_arma);
        arma::uvec entry_stages = find(current_stageframe_entry_arma);
        arma::ivec rep_min_ages = current_stageframe_min_age_arma.elem(rep_stages);
        
        int total_rep_stages = static_cast<int>(rep_stages.n_elem);
        int total_entry_stages = static_cast<int>(entry_stages.n_elem);
        fecage_min_int = min(rep_min_ages);
        
        DataFrame current_agestages;
        if (!using_vrms) {
          RObject mpms_true = RObject(mpms);
          List mpms_list = as<List>(mpms_true);
          List current_mpm = as<List>(mpms_list(current_mpm_int));
          current_agestages = clone(as<DataFrame>(current_mpm["agestages"]));
          
          IntegerVector current_agestages_age = current_agestages["age"];
          firstage_int = min(current_agestages_age);
          finalage_int = max(current_agestages_age);
          
          if (cont.isNotNull()) {
            if (is<IntegerVector>(cont) || is<NumericVector>(cont)) { 
              IntegerVector cont_vec_temp = as<IntegerVector>(cont);
              IntegerVector cont_vec;
              
              if (static_cast<int>(cont_vec_temp.length() == 1)) {
                IntegerVector cont_maxed_out (mpms_entered);
                
                for (int i = 0; i < mpms_entered; i++) {
                  cont_maxed_out(i) = cont_vec_temp(0);
                }
                cont_vec = cont_maxed_out;
              } else cont_vec = cont_vec_temp;
              
              int cont_vec_length = static_cast<int>(cont_vec.length());
              if (cont_vec_length != mpms_entered) {
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
              if (cont_vec(current_mpm_int) > 1) cont_bool = true;
            } if (is<LogicalVector>(cont)) {
              LogicalVector cont_vec_temp = as<LogicalVector>(cont);
              LogicalVector cont_vec_log;
              
              if (static_cast<int>(cont_vec_temp.length() == 1)) {
                LogicalVector cont_maxed_out (mpms_entered);
                
                for (int i = 0; i < mpms_entered; i++) {
                  cont_maxed_out(i) = cont_vec_temp(0);
                }
                cont_vec_log = cont_maxed_out;
              } else cont_vec_log = cont_vec_temp;
              
              int cont_vec_log_length = static_cast<int>(cont_vec_log.length());
              if (cont_vec_log_length != mpms_entered) {
                AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < cont_vec_log_length; i++) {
                if (LogicalVector::is_na(cont_vec_log(i))) {
                  cont_bool = false;
                } else if (cont_vec_log(i)) {
                  cont_bool = true;
                } else if (!cont_vec_log(i)) {
                  cont_bool = false;
                }
              }
            } else {
              AdaptUtils::pop_error2("cont", "an integer vector", "", 1);
            }
          }
          fecage_max_int = finalage_int;
          
        } else {
          if (firstage.isNotNull()) {
            if (is<IntegerVector>(firstage) || is<NumericVector>(firstage)) { 
              IntegerVector firstage_vec = as<IntegerVector>(firstage);
              
              int firstage_vec_length = static_cast<int>(firstage_vec.length());
              if (firstage_vec_length != mpms_entered) {
                AdaptUtils::pop_error2("firstage", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < firstage_vec_length; i++) {
                if (firstage_vec(i) < 0 && !IntegerVector::is_na(firstage_vec(i))) {
                  AdaptUtils::pop_error2("firstage", "", "", 30);
                }
                
                if (firstage_vec(i) > 0 && format_vec(i) < 4) {
                  throw Rcpp::exception("Entries in argument firstage must equal 0 for MPMs without age structure.", false);
                }
                
                if (firstage_vec(i) < firstage_vec(i) && !IntegerVector::is_na(firstage_vec(i))) {
                  throw Rcpp::exception("Values in firstage may not be less than respective entries in firstage.", false);
                }
              }
              
              firstage_int = firstage_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("firstage", "an integer vector", "", 1);
          }
          
          if (finalage.isNotNull()) {
            if (is<IntegerVector>(finalage) || is<NumericVector>(finalage)) { 
              IntegerVector finalage_vec = as<IntegerVector>(finalage);
              
              int finalage_vec_length = static_cast<int>(finalage_vec.length());
              if (finalage_vec_length != mpms_entered) {
                AdaptUtils::pop_error2("finalage", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < finalage_vec_length; i++) {
                if (finalage_vec(i) < 0 && !IntegerVector::is_na(finalage_vec(i))) {
                  AdaptUtils::pop_error2("finalage", "", "", 30);
                }
                
                if (finalage_vec(i) > 0 && format_vec(i) < 4) {
                  throw Rcpp::exception("Entries in argument finalage must equal 0 for MPMs without age structure.", false);
                }
              }
              
              finalage_int = finalage_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("finalage", "an integer vector", "", 1);
          } else throw Rcpp::exception("No information provided on max age in age-by-stage MPMs.", false);
          
          if (fecage_min.isNotNull()) {
            if (is<IntegerVector>(fecage_min) || is<NumericVector>(fecage_min)) {
              IntegerVector fecage_min_prevec = as<IntegerVector>(fecage_min);
              IntegerVector fecage_min_vec;
              
              int fecage_min_length = static_cast<int>(fecage_min_prevec.length());
              
              if (fecage_min_length == mpms_entered) {
                fecage_min_vec = fecage_min_prevec;
              } else if (fecage_min_length == 1) {
                fecage_min_vec = rep(fecage_min_prevec(0), mpms_entered);
              } else AdaptUtils::pop_error2("fecage_min", "vrm_input objects", "vrms", 2);
              
              fecage_min_int = fecage_min_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("fecage_min", "an integer vector", "", 1);
          }
          
          if (fecage_max.isNotNull()) {
            if (is<IntegerVector>(fecage_max) || is<NumericVector>(fecage_max)) {
              IntegerVector fecage_max_prevec = as<IntegerVector>(fecage_max);
              IntegerVector fecage_max_vec;
              
              int fecage_max_length = static_cast<int>(fecage_max_prevec.length());
              
              if (fecage_max_length == mpms_entered) {
                fecage_max_vec = fecage_max_prevec;
              } else if (fecage_max_length == 1) {
                fecage_max_vec = rep(fecage_max_prevec(0), mpms_entered);
              } else AdaptUtils::pop_error2("fecage_max", "vrm_input objects", "vrms", 2);
              
              fecage_max_int = fecage_max_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("fecage_max", "an integer vector", "", 1);
          } else fecage_max_int = finalage_int;
          
          if (cont.isNotNull()) {
            if (is<IntegerVector>(cont) || is<NumericVector>(cont)) { 
              IntegerVector cont_vec_temp = as<IntegerVector>(cont);
              IntegerVector cont_vec;
              
              if (static_cast<int>(cont_vec_temp.length() == 1)) {
                IntegerVector cont_maxed_out (mpms_entered);
                
                for (int i = 0; i < mpms_entered; i++) {
                  cont_maxed_out(i) = cont_vec_temp(0);
                }
                cont_vec = cont_maxed_out;
              } else cont_vec = cont_vec_temp;
              
              int cont_vec_length = static_cast<int>(cont_vec.length());
              if (cont_vec_length != mpms_entered) {
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
              if (cont_vec(current_mpm_int) > 1) cont_bool = true;
            } if (is<LogicalVector>(cont)) {
              LogicalVector cont_vec_temp = as<LogicalVector>(cont);
              LogicalVector cont_vec_log;
              
              if (static_cast<int>(cont_vec_temp.length() == 1)) {
                LogicalVector cont_maxed_out (mpms_entered);
                
                for (int i = 0; i < mpms_entered; i++) {
                  cont_maxed_out(i) = cont_vec_temp(0);
                }
                cont_vec_log = cont_maxed_out;
              } else cont_vec_log = cont_vec_temp;
              
              int cont_vec_log_length = static_cast<int>(cont_vec_log.length());
              if (cont_vec_log_length != mpms_entered) {
                AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < cont_vec_log_length; i++) {
                if (LogicalVector::is_na(cont_vec_log(i))) {
                  cont_bool = false;
                } else if (cont_vec_log(i)) {
                  cont_bool = true;
                } else if (!cont_vec_log(i)) {
                  cont_bool = false;
                }
              }
            } else {
              AdaptUtils::pop_error2("cont", "an integer vector", "", 1);
            }
          }
          
        }
        
        IntegerVector age_sequence = seq(firstage_int, finalage_int);
        IntegerVector rep_age_sequence = seq(fecage_min_int, fecage_max_int);
        
        int total_ages = (finalage_int - firstage_int); // + 1
        int total_rep_ages = (fecage_max_int - fecage_min_int); // + 1
        
        if (cont_bool) {
          age_sequence.push_back(finalage_int);
          rep_age_sequence.push_back(fecage_max_int);
          
          total_ages++;
          total_rep_ages++;
        }
        
        int current_planned_projections = {((total_ages * numstages * numstages) + 
          (total_rep_ages * total_rep_stages * total_entry_stages)) * trial_alterations};
        List projections_list_current_mpm_supp (current_planned_projections);
        List err_check_list_current_mpm_supp (current_planned_projections);
        
        int current_tracker {0};
        
        // Rcout << "Start of main survival batch." << endl;
        // First all survival-focused runs
        for (int age2 = 0; age2 < total_ages; age2++) {
          for (int stage2 = 0; stage2 < numstages; stage2++) {
            for (int stage3 = 0; stage3 < numstages; stage3++) {
              for (int try1 = 0; try1 < trial_alterations; try1++) {
                Rcpp::checkUserInterrupt();
                if (!quiet_bool) {
                  Rcout << "Currently on projection " << (current_tracker + 1) <<
                    " out of " << current_planned_projections << " for age-by-stage mpm " <<
                    (current_mpm_int + 1) << "." << endl;
                }
                
                DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
                CharacterVector new_stage3 = {current_stageframe_stage(stage3)};
                CharacterVector new_stage2 = {current_stageframe_stage(stage2)};
                IntegerVector new_age2 = {age_sequence(age2)};
                IntegerVector new_convtype = {1};
                
                current_skeleton["stage3"] = new_stage3;
                current_skeleton["stage2"] = new_stage2;
                current_skeleton["age2"] = new_age2;
                current_skeleton["convtype"] = new_convtype;
                
                if (givenrate_used) {
                  NumericVector new_altered = {givenrate_vec(try1)};
                  current_skeleton["givenrate"] = new_altered;
                } else if (offset_used) {
                  NumericVector new_altered = {offset_vec(try1)};
                  current_skeleton["offset"] = new_altered;
                } else if (multiplier_used) {
                  NumericVector new_altered = {multiplier_vec(try1)};
                  current_skeleton["multiplier"] = new_altered;
                }
                
                DataFrame new_supplement;
                if (current_supplement.isNotNull()) {
                  if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                    DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                    DataFrame used_supplement = clone(current_supplement_base);
                    
                    new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                  } else new_supplement = clone(current_skeleton);
                } else new_supplement = clone(current_skeleton);
                
                if (err_check_bool) {
                  //DataFrame new_skeleton = clone(new_supplement);
                  DataFrame err_supplement = clone(new_supplement);
                  IntegerVector mpm_id = {used_mpm_vec(i)};
                  err_supplement.push_front(mpm_id);
                  
                  CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                  err_supp_names(0) = "mpm_id";
                  err_supplement.attr("names") = err_supp_names;
                  err_check_list_current_mpm_supp(current_tracker) = err_supplement;
                }
                
                Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                  new_supplement, current_mpm_int, mpms_entered);
                
                List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                  equivalence, starts, years, patches, tweights, format, entry_time,
                  sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                  fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                  stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                  exp_tol, theta_tol);
                
                projections_list_current_mpm_supp(current_tracker) = current_projection;
                
                current_tracker++;
                overall_projection_tracker++;
              }
            }
          }
        }
        
        // Rcout << "Start of main reproduction batch." << endl;
        // Now all reproductive runs
        for (int age2 = 0; age2 < total_rep_ages; age2++) {
          for (int stage2 = 0; stage2 < total_rep_stages; stage2++) {
            for (int stage3 = 0; stage3 < total_entry_stages; stage3++) {
              for (int try1 = 0; try1 < trial_alterations; try1++) {
                Rcpp::checkUserInterrupt();
                if (!quiet_bool) {
                  Rcout << "Currently on projection " << (current_tracker + 1) <<
                    " out of " << current_planned_projections << " for age-by-stage mpm " <<
                    (current_mpm_int + 1) << "." << endl;
                }
                
                DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
                CharacterVector new_stage3 = {current_stageframe_stage(rep_stages(stage3))};
                CharacterVector new_stage2 = {current_stageframe_stage(rep_stages(stage2))};
                IntegerVector new_age2 = {rep_age_sequence(age2)};
                IntegerVector new_convtype = {2};
                
                current_skeleton["stage3"] = new_stage3;
                current_skeleton["stage2"] = new_stage2;
                current_skeleton["age2"] = new_age2;
                current_skeleton["convtype"] = new_convtype;
                
                if (givenrate_used) {
                  NumericVector new_altered = {givenrate_vec(try1)};
                  current_skeleton["givenrate"] = new_altered;
                } else if (offset_used) {
                  NumericVector new_altered = {offset_vec(try1)};
                  current_skeleton["offset"] = new_altered;
                } else if (multiplier_used) {
                  NumericVector new_altered = {multiplier_vec(try1)};
                  current_skeleton["multiplier"] = new_altered;
                }
                
                DataFrame new_supplement;
                if (current_supplement.isNotNull()) {
                  if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                    DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                    DataFrame used_supplement = clone(current_supplement_base);
                    
                    new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                  } else new_supplement = clone(current_skeleton);
                } else new_supplement = clone(current_skeleton);
                
                if (err_check_bool) {
                  //DataFrame new_skeleton = clone(new_supplement);
                  DataFrame err_supplement = clone(new_supplement);
                  IntegerVector mpm_id = {used_mpm_vec(i)};
                  err_supplement.push_front(mpm_id);
                  
                  CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                  err_supp_names(0) = "mpm_id";
                  err_supplement.attr("names") = err_supp_names;
                  err_check_list_current_mpm_supp(current_tracker) = err_supplement;
                }
                
                Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                  new_supplement, current_mpm_int, mpms_entered);
                
                List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                  equivalence, starts, years, patches, tweights, format, entry_time,
                  sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                  fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                  stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                  exp_tol, theta_tol);
                
                projections_list_current_mpm_supp(current_tracker) = current_projection;
                
                current_tracker++;
                overall_projection_tracker++;
              }
            }
          }
        }
        projections_list(i) = projections_list_current_mpm_supp;
        err_check_list(i) = err_check_list_current_mpm_supp;
        
      } else if (format_vec(current_mpm_int) == 5) {
        // Rcout << "batch_project3 J Leslie" << endl;
        
        // Need to set up a test if ahstages is NULL, to create a pseudo-ahstages
        int firstage_int {0};
        int finalage_int {0};
        int fecage_min_int {0};
        int fecage_max_int {0};
        bool cont_bool {false};
        
        if (!(static_cast<int>(current_stageframe.length()) > 1)) {
          if (firstage.isNotNull()) {
            if (is<IntegerVector>(firstage) || is<NumericVector>(firstage)) { 
              IntegerVector firstage_vec = as<IntegerVector>(firstage);
              
              int firstage_vec_length = static_cast<int>(firstage_vec.length());
              if (firstage_vec_length != mpms_entered) {
                AdaptUtils::pop_error2("firstage", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < firstage_vec_length; i++) {
                if (firstage_vec(i) < 0 && !IntegerVector::is_na(firstage_vec(i))) {
                  AdaptUtils::pop_error2("firstage", "", "", 30);
                }
                
                if (firstage_vec(i) > 0 && format_vec(i) < 4) {
                  throw Rcpp::exception("Entries in argument firstage must equal 0 for MPMs without age structure.", false);
                }
                
                if (firstage_vec(i) < firstage_vec(i) && !IntegerVector::is_na(firstage_vec(i))) {
                  throw Rcpp::exception("Values in firstage may not be less than respective entries in firstage.", false);
                }
              }
              
              firstage_int = firstage_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("firstage", "an integer vector", "", 1);
          }
          
          if (finalage.isNotNull()) {
            if (is<IntegerVector>(finalage) || is<NumericVector>(finalage)) { 
              IntegerVector finalage_vec = as<IntegerVector>(finalage);
              
              int finalage_vec_length = static_cast<int>(finalage_vec.length());
              if (finalage_vec_length != mpms_entered) {
                AdaptUtils::pop_error2("finalage", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < finalage_vec_length; i++) {
                if (finalage_vec(i) < 0 && !IntegerVector::is_na(finalage_vec(i))) {
                  AdaptUtils::pop_error2("finalage", "", "", 30);
                }
                
                if (finalage_vec(i) > 0 && format_vec(i) < 4) {
                  throw Rcpp::exception("Entries in argument finalage must equal 0 for MPMs without age structure.", false);
                }
              }
              
              finalage_int = finalage_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("finalage", "an integer vector", "", 1);
          } else throw Rcpp::exception("No information provided on max age in Leslie MPMs.", false);
          
          if (fecage_min.isNotNull()) {
            if (is<IntegerVector>(fecage_min) || is<NumericVector>(fecage_min)) {
              IntegerVector fecage_min_prevec = as<IntegerVector>(fecage_min);
              IntegerVector fecage_min_vec;
              
              int fecage_min_length = static_cast<int>(fecage_min_prevec.length());
              
              if (fecage_min_length == mpms_entered) {
                fecage_min_vec = fecage_min_prevec;
              } else if (fecage_min_length == 1) {
                fecage_min_vec = rep(fecage_min_prevec(0), mpms_entered);
              } else AdaptUtils::pop_error2("fecage_min", "vrm_input objects", "vrms", 2);
              
              fecage_min_int = fecage_min_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("fecage_min", "an integer vector", "", 1);
          }
          
          if (fecage_max.isNotNull()) {
            if (is<IntegerVector>(fecage_max) || is<NumericVector>(fecage_max)) {
              IntegerVector fecage_max_prevec = as<IntegerVector>(fecage_max);
              IntegerVector fecage_max_vec;
              
              int fecage_max_length = static_cast<int>(fecage_max_prevec.length());
              
              if (fecage_max_length == mpms_entered) {
                fecage_max_vec = fecage_max_prevec;
              } else if (fecage_max_length == 1) {
                fecage_max_vec = rep(fecage_max_prevec(0), mpms_entered);
              } else AdaptUtils::pop_error2("fecage_max", "vrm_input objects", "vrms", 2);
              
              fecage_max_int = fecage_max_vec(current_mpm_int);
            } else AdaptUtils::pop_error2("fecage_max", "an integer vector", "", 1);
          } else fecage_max_int = finalage_int;
          
          if (cont.isNotNull()) {
            if (is<IntegerVector>(cont) || is<NumericVector>(cont)) { 
              IntegerVector cont_vec_temp = as<IntegerVector>(cont);
              IntegerVector cont_vec;
              
              if (static_cast<int>(cont_vec_temp.length() == 1)) {
                IntegerVector cont_maxed_out (mpms_entered);
                
                for (int i = 0; i < mpms_entered; i++) {
                  cont_maxed_out(i) = cont_vec_temp(0);
                }
                cont_vec = cont_maxed_out;
              } else cont_vec = cont_vec_temp;
              
              int cont_vec_length = static_cast<int>(cont_vec.length());
              if (cont_vec_length != mpms_entered) {
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
              if (cont_vec(current_mpm_int) > 1) cont_bool = true;
            } if (is<LogicalVector>(cont)) {
              LogicalVector cont_vec_temp = as<LogicalVector>(cont);
              LogicalVector cont_vec_log;
              
              if (static_cast<int>(cont_vec_temp.length() == 1)) {
                LogicalVector cont_maxed_out (mpms_entered);
                
                for (int i = 0; i < mpms_entered; i++) {
                  cont_maxed_out(i) = cont_vec_temp(0);
                }
                cont_vec_log = cont_maxed_out;
              } else cont_vec_log = cont_vec_temp;
              
              int cont_vec_log_length = static_cast<int>(cont_vec_log.length());
              if (cont_vec_log_length != mpms_entered) {
                AdaptUtils::pop_error2("cont", "number of MPMs to project", "", 29);
              }
              
              for (int i = 0; i < cont_vec_log_length; i++) {
                if (LogicalVector::is_na(cont_vec_log(i))) {
                  cont_bool = false;
                } else if (cont_vec_log(i)) {
                  cont_bool = true;
                } else if (!cont_vec_log(i)) {
                  cont_bool = false;
                }
              }
            } else {
              AdaptUtils::pop_error2("cont", "an integer vector", "", 1);
            }
          }
          
          CharacterVector new_current_stageframe_stage (finalage_int - firstage_int + 1);
          
          for (int j = 0; j < (finalage_int - firstage_int + 1); j++) {
            String new_age = "Age";
            new_age += (firstage_int + j);
            
            new_current_stageframe_stage(j) = new_age;
          }
          current_stageframe_stage = new_current_stageframe_stage;
          
        } else {
          IntegerVector current_stageframe_minage = as<IntegerVector>(current_stageframe["min_age"]);
          IntegerVector current_stageframe_maxage = as<IntegerVector>(current_stageframe["max_age"]);
          
          arma::ivec current_stageframe_repstatus_arma = as<arma::ivec>(current_stageframe["repstatus"]);
          arma::uvec repstages_arma = find(current_stageframe_repstatus_arma);
          
          int current_stageframe_minage_length = static_cast<int>(current_stageframe_minage.length());
          
          int minage_NAs {0};
          int maxage_NAs {0};
          for (int j = 0; j < current_stageframe_minage_length; j++) {
            if (IntegerVector::is_na(current_stageframe_minage(j))) minage_NAs++;
            if (IntegerVector::is_na(current_stageframe_maxage(j))) {
              if (j == (current_stageframe_minage_length - 1)) cont_bool = true;
              maxage_NAs++;
            }
          }
          
          if (minage_NAs == 0) {
            int minage_min = min(current_stageframe_minage);
            int minage_max = max(current_stageframe_minage);
            
            firstage_int = minage_min;
            finalage_int = minage_max;
          } else {
            throw Rcpp::exception("Information on first age and final age are required for all Leslie MPMs.", false);
          }
          
          IntegerVector trial_all_ages = seq(firstage_int, finalage_int);
          fecage_min_int = trial_all_ages(min(repstages_arma));
          fecage_max_int = trial_all_ages(max(repstages_arma));
        }
        
        if (fecage_min_int < firstage_int) fecage_min_int = firstage_int;
        if (fecage_max_int < fecage_min_int) fecage_max_int = finalage_int;
        
        IntegerVector allstages = seq(firstage_int, finalage_int);
        IntegerVector rep_allstages = seq(fecage_min_int, fecage_max_int);
        
        numstages = (finalage_int - firstage_int + 1);
        int rep_numstages = (fecage_max_int - fecage_min_int + 1);
        
        int current_planned_projections = {(numstages + rep_numstages) * trial_alterations};
        List projections_list_current_mpm_supp (current_planned_projections);
        List err_check_list_current_mpm_supp (current_planned_projections);
        
        int current_tracker {0};
        
        // Rcout << "Survival projections starting." << endl;
        // Survival projections
        for (int age2 = 0; age2 < numstages; age2++) {
          for (int try1 = 0; try1 < trial_alterations; try1++) {
            if ((!cont_bool && age2 < (numstages - 1)) || (cont_bool)) {
              Rcpp::checkUserInterrupt();
              if (!quiet_bool) {
                Rcout << "Currently on projection " << (current_tracker + 1) <<
                  " out of " << current_planned_projections << " for Leslie mpm " <<
                  (current_mpm_int + 1) << "." << endl;
              }
              
              DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
              CharacterVector new_stage2 = {current_stageframe_stage(age2)};
              IntegerVector new_age2 = {allstages(age2)};
              CharacterVector new_stage3;
              if (age2 < (numstages - 1)) {
                new_stage3 = {current_stageframe_stage(0)};
              } else if (cont_bool) {
                new_stage3 = {current_stageframe_stage(age2)};
              }
              IntegerVector new_convtype = {1};
              
              current_skeleton["age2"] = new_age2;
              current_skeleton["stage3"] = new_stage3;
              current_skeleton["stage2"] = new_stage2;
              current_skeleton["convtype"] = new_convtype;
              
              if (givenrate_used) {
                NumericVector new_altered = {givenrate_vec(try1)};
                current_skeleton["givenrate"] = new_altered;
              } else if (offset_used) {
                NumericVector new_altered = {offset_vec(try1)};
                current_skeleton["offset"] = new_altered;
              } else if (multiplier_used) {
                NumericVector new_altered = {multiplier_vec(try1)};
                current_skeleton["multiplier"] = new_altered;
              }
              
              DataFrame new_supplement;
              if (current_supplement.isNotNull()) {
                if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                  DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                  DataFrame used_supplement = clone(current_supplement_base);
                  
                  new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                } else new_supplement = clone(current_skeleton);
              } else new_supplement = clone(current_skeleton);
              
              if (err_check_bool) {
                //DataFrame new_skeleton = clone(new_supplement);
                DataFrame err_supplement = clone(new_supplement);
                IntegerVector mpm_id = {used_mpm_vec(i)};
                err_supplement.push_front(mpm_id);
                
                CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                err_supp_names(0) = "mpm_id";
                err_supplement.attr("names") = err_supp_names;
                err_check_list_current_mpm_supp(current_tracker) = err_supplement;
              }
              
              Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                new_supplement, current_mpm_int, mpms_entered);
              
              List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                equivalence, starts, years, patches, tweights, format, entry_time,
                sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                exp_tol, theta_tol);
              
              projections_list_current_mpm_supp(current_tracker) = current_projection;
              
              current_tracker++;
              overall_projection_tracker++;
            }
          }
        }
        
        // Rcout << "Fecundity projections starting." << endl;
        // Fecundity projections
        for (int age2 = 0; age2 < rep_numstages; age2++) {
          for (int try1 = 0; try1 < trial_alterations; try1++) {
            if ((!cont_bool && age2 < (numstages - 1)) || (cont_bool)) {
              Rcpp::checkUserInterrupt();
              if (!quiet_bool) {
                Rcout << "Currently on projection " << (current_tracker + 1) <<
                  " out of " << current_planned_projections << " for Leslie mpm " <<
                  (current_mpm_int + 1) << "." << endl;
              }
              
              DataFrame current_skeleton = LefkoInputs::sp_skeleton(1);
              CharacterVector new_stage2 = {current_stageframe_stage(rep_allstages(age2) - 1)};
              IntegerVector new_age2 = {rep_allstages(age2)};
              CharacterVector new_stage3 = {current_stageframe_stage(0)};
              IntegerVector new_convtype = {2};
              
              current_skeleton["age2"] = new_age2;
              current_skeleton["stage3"] = new_stage3;
              current_skeleton["stage2"] = new_stage2;
              current_skeleton["convtype"] = new_convtype;
              
              if (givenrate_used) {
                NumericVector new_altered = {givenrate_vec(try1)};
                current_skeleton["givenrate"] = new_altered;
              } else if (offset_used) {
                NumericVector new_altered = {offset_vec(try1)};
                current_skeleton["offset"] = new_altered;
              } else if (multiplier_used) {
                NumericVector new_altered = {multiplier_vec(try1)};
                current_skeleton["multiplier"] = new_altered;
              }
              
              DataFrame new_supplement;
              if (current_supplement.isNotNull()) {
                if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
                  DataFrame current_supplement_base = as<DataFrame>(current_supplement);
                  DataFrame used_supplement = clone(current_supplement_base);
                  
                  new_supplement = LefkoUtils::df_rbind(used_supplement, current_skeleton);
                } else new_supplement = clone(current_skeleton);
              } else new_supplement = clone(current_skeleton);
              
              if (err_check_bool) {
                //DataFrame new_skeleton = clone(new_supplement);
                DataFrame err_supplement = clone(new_supplement);
                IntegerVector mpm_id = {used_mpm_vec(i)};
                err_supplement.push_front(mpm_id);
                
                CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
                err_supp_names(0) = "mpm_id";
                err_supplement.attr("names") = err_supp_names;
                err_check_list_current_mpm_supp(current_tracker) = err_supplement;
              }
              
              Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
                new_supplement, current_mpm_int, mpms_entered);
              
              List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
                equivalence, starts, years, patches, tweights, format, entry_time,
                sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
                fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
                stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
                exp_tol, theta_tol);
              
              projections_list_current_mpm_supp(current_tracker) = current_projection;
              
              current_tracker++;
              overall_projection_tracker++;
            }
          }
        }
        projections_list(i) = projections_list_current_mpm_supp;
        err_check_list(i) = err_check_list_current_mpm_supp;
        
      }
    } else {
      List projections_list_current_mpm_supp (trial_alterations);
      List err_check_list_current_mpm_supp (trial_alterations);
      
      if (is<DataFrame>(current_supplement) || is<List>(current_supplement)) {
        DataFrame current_supplement_df = as<DataFrame>(current_supplement);
        int current_supplement_rows = static_cast<int>(current_supplement_df.nrows());
        
        for (int try1 = 0; try1 < trial_alterations; try1++) {
          Rcpp::checkUserInterrupt();
          if (!quiet_bool) {
            Rcout << "Currently on projection " << (overall_projection_tracker + 1) <<
              " out of " << planned_projections << " for mpm " << (current_mpm_int + 1) << "." << endl;
          }
          
          DataFrame used_supplement = clone(current_supplement_df);
          
          if (givenrate_used) {
            NumericVector givenrate_new = rep(givenrate_vec(try1), current_supplement_rows);
            used_supplement["givenrate"] = givenrate_new;
          } else if (offset_used) {
            NumericVector offset_new = rep(offset_vec(try1), current_supplement_rows);
            used_supplement["offset"] = offset_new;
          } else if (multiplier_used) {
            NumericVector multiplier_new = rep(multiplier_vec(try1), current_supplement_rows);
            used_supplement["multiplier"] = multiplier_new;
          }
          
          arma::ivec used_convtype = as<arma::ivec>(used_supplement["convtype"]);
          arma::uvec used_convtype_toosmall = find(used_convtype < 1);
          arma::uvec used_convtype_toolarge = find(used_convtype > 3);
          int ucts = used_convtype_toosmall.n_elem;
          int uctl = used_convtype_toolarge.n_elem;
          if ((ucts + uctl) > 0) {
            throw Rcpp::exception("Argument convtype in each supplement must be set to either 1, 2, or 3.", false);
          }
          
          if (err_check_bool && current_supplement.isNotNull()) {
            DataFrame err_supplement = clone(used_supplement);
            IntegerVector mpm_id = {used_mpm_vec(i)};
            err_supplement.push_front(mpm_id);
                  
            CharacterVector err_supp_names = as<CharacterVector>(err_supplement.attr("names"));
            err_supp_names(0) = "mpm_id";
            err_supplement.attr("names") = err_supp_names;
            
            err_check_list_current_mpm_supp(try1) = err_supplement;
          }
          
          Nullable<RObject> new_supplements = AdaptUtils::supplements_replacer (supplements,
            used_supplement, current_mpm_int, mpms_entered);
          
          List current_projection = project3 (mpms, vrms, stageframes, new_supplements,
            equivalence, starts, years, patches, tweights, format, entry_time,
            sp_density, ind_terms, dev_terms, fb_sparse, firstage, finalage,
            fecage_min, fecage_max, cont, fecmod, density, density_vr, err_check,
            stochastic, integeronly, substoch, nreps, times, prep_mats, force_fb,
            exp_tol, theta_tol);
          
          projections_list_current_mpm_supp(try1) = current_projection;
          
          overall_projection_tracker++;
        }
        projections_list(i) = projections_list_current_mpm_supp;
        err_check_list(i) = err_check_list_current_mpm_supp;
        
      } else if (!current_supplement.isNotNull()) {
        Rf_warningcall(R_NilValue, "If all_elems = FALSE, then all MPMs require supplements defining elements to alter.");
      }
    }
  }
  
  int final_list_length {3};
  CharacterVector output_main_names;
  
  List output (final_list_length);
  output(0) = projections_list;
  output(1) = err_check_list;
  output(2) = used_mpm_vec;
  
  output_main_names = {"proj_out", "ref", "control"};
  
  output.attr("names") = output_main_names;
  output.attr("class") = "adaptProjBatch";
  
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
//' @param dev_terms An optional list of data frames, one for each
//' \code{vrm_input} object. Each should include 14 or 17 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival,
//' 9: juvenile observation, 10: juvenile primary size, 11: juvenile secondary
//' size, 12: juvenile tertiary size, 13: juvenile reproduction, and
//' 14: juvenile maturity transition. In addition, these may be followed by 3
//' columns designating additive deviations to: 15: individual covariate a,
//' 16: individual covariate b, and 17: individual covariate c. Unused terms
//' must be set to \code{0} (use of \code{NA} will produce errors). Single or
//' small numbers of values per vital rate model are also allowed, and if the
//' number of rows is less than \code{times}, then the terms will be cycled.
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
  
  // Rcout << "cleanup3_inv A ";
  
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
  
  // Rcout << "cleanup3_inv B ";
  
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
  
  // Rcout << "cleanup3_inv C ";
  
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
  
  // Rcout << "cleanup3_inv D ";
  
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
  
  // Rcout << "cleanup3_inv E ";
  
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
  
  // Rcout << "cleanup3_inv F ";
  
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
  
  // Rcout << "cleanup3_inv G ";
  
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
  
  // Rcout << "cleanup3_inv H ";
  
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
  
  // Rcout << "cleanup3_inv I ";
  
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
    // Rcout << "FUNCTION cleanup3_inv DID NOT FIND START OPTION                  ";
    
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
  
  // Rcout << "cleanup3_inv J ";
  
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
  
  // Rcout << "cleanup3_inv K ";
  
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
  
  // Rcout << "cleanup3_inv L ";
  
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
  
  // Rcout << "cleanup3_inv M ";
  
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
  
  // Rcout << "cleanup3_inv N ";
  
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
      arma::vec dyn_gamma = as<arma::vec>(chosen_density["gamma"]);
      
      for (int j = 0; j < static_cast<int>(dyn_style.n_elem); j++) {
        if (dyn_style(j) < 1 || dyn_style(j) > 6) {
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
        arma::vec dyn_gamma = as<arma::vec>(chosen_density["gamma"]);
        
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
  
  // Rcout << "cleanup3_inv O ";
  
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
  
  // Rcout << "cleanup3_inv P ";
  
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
        
        if (dvtc_df_size != 14 && dvtc_df_size != 17) {
          throw Rcpp::exception("Data frame in argument dev_terms must have 14 or 17 columns.",
            false);
        }
        
        for (int j = 0; j < dvtc_df_size; j++) {
          if (!is<NumericVector>(current_dev_terms(j))) {
            AdaptUtils::pop_error2("dev_terms", "a data frame composed of numeric values", "", 1);
          }
        }
        
        dev_terms_num_int = dvtc_df_nrows;
        
        NumericMatrix dev_matrix (17, dvtc_df_nrows); // Rows = devs / vital rate, Cols = times
        
        for (int i = 0; i < dvtc_df_size; i++) {
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
        
        if (input_dv_cols != 14 && input_dv_cols != 17) {
          throw Rcpp::exception("Numeric matrix in argument dev_terms must have 14 or 17 columns",
            false);
        }
        
        dev_terms_num_int = input_dv_rows;
        
        NumericMatrix dev_matrix (17, input_dv_rows); // Rows are the devs per vital rate, Cols are times
        
        for (int i = 0; i < input_dv_cols; i++) {
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
        
        if (input_dv_cols != 14 && input_dv_cols != 17) {
          throw Rcpp::exception("Numeric vector in argument dev_terms must have 14 or 17 columns",
            false);
        }
        
        dev_terms_num_int = 1;
        
        NumericMatrix dev_matrix (17, 1); // Rows are the devs per vital rate, Cols are times
        
        for (int i = 0; i < input_dv_cols; i++) {
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
  
  // Rcout << "cleanup3_inv Q ";
  
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
  
  // Rcout << "cleanup3_inv R ";
  
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
  
  // Rcout << "cleanup3_inv S ";
  
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

//' Core Pre-Existing MPM Projection Engine for ESS Evaluation
//' 
//' Function \code{invpre_optim_single} runs single simulation instances of the
//' optimization projections in function \code{invade3_pre_core}, and is used to
//' estimate ESS values.
//' 
//' @name invpre_optim_singlerun
//' 
//' @param out_df A data frame created by reference to hold the fitness values
//' produced by the run.
//' @param base_trait_axis The currently used trait axis.
//' @param N_out A matrix giving the final population sizes of the run. Supplied
//' as a reference.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function. Needs to be essentially
//' blank to work properly. Needs to me a lst with only one level, with elements
//' as matrics.
//' @param new_stageexpansion_list A list with stage expansions only for the
//' variants to run, in order.
//' @param used_times A list of year numbers for each time per run.
//' @param zero_stage_vec_list A list of population stage vectors full of zeros.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param A_list A list of all A matrices.
//' @param U_list A list of all U matrices.
//' @param F_list A list of all F matrices.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param times An integer giving the number of occasions to project.
//' @param fitness_times The number of occasions at the end of each run to use
//' to estimate the Lyapunov coefficient.
//' @param format_int An integer giving the MPM format.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' @param conv_threshold The lower limit for the absolute value of fitness,
//' below which fitness is rounded to 0. Defaults to 0.00000001.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' @param zap_min A Boolean value describing whether to round fitness values
//' below the value given in \code{threshold}.
//' 
//' @return An arma style matrix is produced and returned via reference.
//' 
//' @keywords internal
//' @noRd
inline void invpre_optim_singlerun (DataFrame& out_df,
  DataFrame& base_trait_axis, arma::mat& N_out,
  List& comm_out, List& new_stageexpansion_list, List& used_times,
  List& zero_stage_vec_list, const List start_list, const List equivalence_list,
  const List A_list, const List U_list, const List F_list,
  const DataFrame density_df, const DataFrame dens_index_df,
  const IntegerVector entry_time_vec, const int times, const int fitness_times,
  const int format_int, const int firstage_int, const int finalage_int,
  const int substoch, const double exp_tol, const double theta_tol,
  const double conv_threshold, const bool sparse_bool, const bool A_only,
  const bool stages_not_equal, const bool integeronly, const bool dens_yn_bool,
  const bool zap_min) {
  
  // Rcout << "iposr A" << endl;
  
  int var_to_run {2};
  List running_popvecs; //  = clone(start_list)
  List running_popvecs_startonly; //  = clone(start_list)
  arma::mat N_mpm (2, (times + 1)); // rows = vars, cols = times
  
  CharacterVector stage3_nta = as<CharacterVector>(base_trait_axis["stage3"]);
  CharacterVector stage2_nta = as<CharacterVector>(base_trait_axis["stage2"]);
  CharacterVector stage1_nta = as<CharacterVector>(base_trait_axis["stage1"]);
  IntegerVector age3_nta = as<IntegerVector>(base_trait_axis["age3"]);
  IntegerVector age2_nta = as<IntegerVector>(base_trait_axis["age2"]);
  CharacterVector eststage3_nta = as<CharacterVector>(base_trait_axis["eststage3"]);
  CharacterVector eststage2_nta = as<CharacterVector>(base_trait_axis["eststage2"]);
  CharacterVector eststage1_nta = as<CharacterVector>(base_trait_axis["eststage1"]);
  IntegerVector estage3_nta = as<IntegerVector>(base_trait_axis["estage3"]);
  IntegerVector estage2_nta = as<IntegerVector>(base_trait_axis["estage2"]);
  NumericVector givenrate_nta = as<NumericVector>(base_trait_axis["givenrate"]);
  NumericVector offset_nta = as<NumericVector>(base_trait_axis["offset"]);
  NumericVector multiplier_nta = as<NumericVector>(base_trait_axis["multiplier"]);
  IntegerVector convtype_nta = as<IntegerVector>(base_trait_axis["convtype"]);
  IntegerVector convtype_t12_nta = as<IntegerVector>(base_trait_axis["convtype_t12"]);
  CharacterVector year2_nta = as<CharacterVector>(base_trait_axis["year2"]);
  IntegerVector mpm_altered_nta = as<IntegerVector>(base_trait_axis["mpm_altered"]);
  IntegerVector vrm_altered_nta = as<IntegerVector>(base_trait_axis["vrm_altered"]);
  
  arma::ivec variant_nta = as<arma::ivec>(base_trait_axis["variant"]);
  arma::vec surv_dev_nta = as<arma::vec>(base_trait_axis["surv_dev"]);
  arma::vec obs_dev_nta = as<arma::vec>(base_trait_axis["obs_dev"]);
  arma::vec size_dev_nta = as<arma::vec>(base_trait_axis["size_dev"]);
  arma::vec sizeb_dev_nta = as<arma::vec>(base_trait_axis["sizeb_dev"]);
  arma::vec sizec_dev_nta = as<arma::vec>(base_trait_axis["sizec_dev"]);
  arma::vec repst_dev_nta = as<arma::vec>(base_trait_axis["repst_dev"]);
  arma::vec fec_dev_nta = as<arma::vec>(base_trait_axis["fec_dev"]);
  arma::vec jsurv_dev_nta = as<arma::vec>(base_trait_axis["jsurv_dev"]);
  arma::vec jobs_dev_nta = as<arma::vec>(base_trait_axis["jobs_dev"]);
  arma::vec jsize_dev_nta = as<arma::vec>(base_trait_axis["jsize_dev"]);
  arma::vec jsizeb_dev_nta = as<arma::vec>(base_trait_axis["jsizeb_dev"]);
  arma::vec jsizec_dev_nta = as<arma::vec>(base_trait_axis["jsizec_dev"]);
  arma::vec jrepst_dev_nta = as<arma::vec>(base_trait_axis["jrepst_dev"]);
  arma::vec jmatst_dev_nta = as<arma::vec>(base_trait_axis["jmatst_dev"]);
  arma::vec indcova_nta = as<arma::vec>(base_trait_axis["indcova"]);
  arma::vec indcovb_nta = as<arma::vec>(base_trait_axis["indcovb"]);
  arma::vec indcovc_nta = as<arma::vec>(base_trait_axis["indcovc"]);
  
  for (int j = 0; j < times; j++) { // 2nd loop - time j
    // Rcout << "iposr A1          ";
    if (j % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    if (j == 0) {
      List var_popvecs_to_start (var_to_run);
      for (int n = 0; n < var_to_run; n++) {
        var_popvecs_to_start(n) = as<arma::vec>(start_list(static_cast<int>(n))); // May need to change (n) to (0)
      }
      running_popvecs = var_popvecs_to_start;
      running_popvecs_startonly = clone(var_popvecs_to_start);
    }
    
    for (int m = 0; m < var_to_run; m++) {
      if (j == entry_time_vec(m)) {
        arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(0)); // the 0 was originally m
        
        double N_current = accu(running_popvec_mpm);
        N_mpm(m, j) = N_current;
      }
    }
    
    for (int m = 0; m < var_to_run; m++) { // main variant run
      // Rcout << "iposr A3          ";
      DataFrame sge_current = as<DataFrame>(new_stageexpansion_list(m)); // new_stageexpansion_list to be subset & ordered
      arma::mat pops_out = as<arma::mat>(comm_out(m));
      
      if (j > (entry_time_vec(m) - 1)) {
        List used_times_higher = as<List>(used_times(0));
        List used_times_lower = as<List>(used_times_higher(0));
        IntegerVector current_times_vec = as<IntegerVector>(used_times_lower(0)); // Maybe m instead of 0?
        
        arma::vec running_popvec_mpm;
        if (j == entry_time_vec(m)) {
          running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
          pops_out.col(j) = running_popvec_mpm;
          
        } else {
          running_popvec_mpm = pops_out.col(j);
        }
        
        // Rcout << "iposr A5          ";
        if (!dens_yn_bool) {
          if (!sparse_bool) {
            // Rcout << "iposr A6          ";
            arma::mat current_A = as<arma::mat>(A_list(current_times_vec(j)));
            
            if (A_only) {
              AdaptUtils::Amat_alter(current_A, sge_current); 
            } else {
              arma::mat current_U_unaltered = as<arma::mat>(U_list(current_times_vec(j)));
              arma::mat current_F_unaltered = as<arma::mat>(F_list(current_times_vec(j)));
              AdaptUtils::UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
            }
            
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
            // Rcout << "iposr A11          ";
            running_popvec_mpm = current_A * running_popvec_mpm;
          }
        } else {
          // Rcout << "iposr A13          ";
          DataFrame used_density_input = density_df;
          DataFrame used_density_index_input = as<DataFrame>(dens_index_df);
          
          IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
          int used_delay = max(ud_delay_vec);
          
          if (j >= (used_delay - 1 )) { // Change to allow different delay Ns for different entries
            if (!stages_not_equal) {
              arma::vec delay_issue = N_mpm.col(j + 1 - used_delay);
              
              double delay_N_sum = arma::sum(delay_issue);
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                running_popvec_mpm, sge_current, A_list, delay_N_sum,
                static_cast<int>(current_times_vec(j)), integeronly,
                substoch, used_density_input, used_density_index_input,
                false, sparse_bool, sparse_bool, false, false);
              
              // Rcout << "iposr A18          ";
              running_popvec_mpm = new_popvec;
            } else {
              // Rcout << "iposr A20          ";
              double delay_N_sum {0.0};
              
              if (j > 0) {
                for (int p = 0; p < var_to_run; p++) {
                  arma::vec delay_pop_vec = N_mpm.col(j + 1 - used_delay);
                  arma::vec current_equiv_vec = as<arma::vec>(equivalence_list(0));
                  arma::vec adjusted_delay_pop_vec = delay_pop_vec % current_equiv_vec;
                  double delay_pop_N = arma::accu(adjusted_delay_pop_vec);
                  
                  delay_N_sum += delay_pop_N;
                }
              }
              
              // Rcout << "iposr A21          ";
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                running_popvec_mpm, sge_current, A_list, delay_N_sum,
                static_cast<int>(current_times_vec(j)), integeronly,
                substoch, used_density_input, used_density_index_input,
                false, sparse_bool, sparse_bool, false, false);
              
              // Rcout << "iposr A22          ";
              running_popvec_mpm = new_popvec;
            }
          } else {
            // Rcout << "iposr A24          ";
            arma::vec new_popvec;
            arma::mat new_projmat;
            arma::sp_mat new_projsp;
            
            AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
              running_popvec_mpm, sge_current, A_list, 0.0,
              static_cast<int>(current_times_vec(j)), integeronly, substoch,
              used_density_input, used_density_index_input, false,
              sparse_bool, sparse_bool, false, false);
            
            // Rcout << "iposr A25          ";
            running_popvec_mpm = new_popvec;
          }
        }
        
        // Rcout << "iposr A27          ";
        if (integeronly) running_popvec_mpm = floor(running_popvec_mpm);
        double N_current = arma::sum(running_popvec_mpm);
        N_mpm(m, (j + 1)) = N_current;
        
        running_popvecs(m) = running_popvec_mpm;
        pops_out.col(j + 1) = running_popvec_mpm;
      } else {
        arma::vec current_zero_vec = as<arma::vec>(zero_stage_vec_list(0)); // The 0 was originally m
        pops_out.col(j + 1) = current_zero_vec;
      }
    } // m loop - var_per_run
  } // j loop - time
  N_out = N_mpm;
  
  arma::rowvec N_out_row0 = N_out.row(0);
  arma::rowvec N_out_row1 = N_out.row(1);
  
  List N_out_list (1);
  N_out_list(0) = N_out;
  
  // Rcout << "ifosr D got to Lyapunov_creator" << endl;
  AdaptUtils::Lyapunov_creator (out_df, N_out_list, entry_time_vec, 1,
    var_to_run, 2, times, fitness_times, conv_threshold, 1, true, zap_min);
  
  // Rcout << "ifosr E" << endl;
  // This bit needs to be redone to allow multiple variants, as right now it really only allows 2
  NumericVector out_fitness = as<NumericVector>(out_df(6));
  NumericVector out_fitness_e995 = as<NumericVector>(out_df(7));
  
  // Rcout << "ifosr F" << endl;
  CharacterVector out_names = as<CharacterVector>(out_df.names());
  
  NumericVector final_fitness (var_to_run);
  final_fitness(0) = out_fitness(0);
  final_fitness(1) = out_fitness_e995(0);
  
  LogicalVector converged (2);
  
  List processed_fitness_out (38);
  
  processed_fitness_out(0) = as<IntegerVector>(wrap(variant_nta));
  processed_fitness_out(1) = stage3_nta;
  processed_fitness_out(2) = stage2_nta;
  processed_fitness_out(3) = stage1_nta;
  processed_fitness_out(4) = age3_nta;
  processed_fitness_out(5) = age2_nta;
  processed_fitness_out(6) = eststage3_nta;
  processed_fitness_out(7) = eststage2_nta;
  processed_fitness_out(8) = eststage1_nta;
  processed_fitness_out(9) = estage3_nta;
  processed_fitness_out(10) = estage2_nta;
  processed_fitness_out(11) = givenrate_nta;
  processed_fitness_out(12) = offset_nta;
  processed_fitness_out(13) = multiplier_nta;
  processed_fitness_out(14) = convtype_nta;
  processed_fitness_out(15) = convtype_t12_nta;
  processed_fitness_out(16) = wrap(surv_dev_nta);
  processed_fitness_out(17) = wrap(obs_dev_nta);
  processed_fitness_out(18) = wrap(size_dev_nta);
  processed_fitness_out(19) = wrap(sizeb_dev_nta);
  processed_fitness_out(20) = wrap(sizec_dev_nta);
  processed_fitness_out(21) = wrap(repst_dev_nta);
  processed_fitness_out(22) = wrap(fec_dev_nta);
  processed_fitness_out(23) = wrap(jsurv_dev_nta);
  processed_fitness_out(24) = wrap(jobs_dev_nta);
  processed_fitness_out(25) = wrap(jsize_dev_nta);
  processed_fitness_out(26) = wrap(jsizeb_dev_nta);
  processed_fitness_out(27) = wrap(jsizec_dev_nta);
  processed_fitness_out(28) = wrap(jrepst_dev_nta);
  processed_fitness_out(29) = wrap(jmatst_dev_nta);
  processed_fitness_out(30) = wrap(indcova_nta);
  processed_fitness_out(31) = wrap(indcovb_nta);
  processed_fitness_out(32) = wrap(indcovc_nta);
  processed_fitness_out(33) = year2_nta;
  processed_fitness_out(34) = mpm_altered_nta;
  processed_fitness_out(35) = vrm_altered_nta;
  processed_fitness_out(36) = final_fitness;
  processed_fitness_out(37) = converged;
  
  CharacterVector processed_out_names = {"variant", "stage3", "stage2", "stage1",
    "age3", "age2", "eststage3", "eststage2", "eststage1", "estage3", "estage2",
    "givenrate", "offset", "multiplier", "convtype", "convtype_t12", "surv_dev",
    "obs_dev", "size_dev", "sizeb_dev", "sizec_dev", "repst_dev", "fec_dev",
    "jsurv_dev", "jobs_dev", "jsize_dev", "jsizeb_dev", "jsizec_dev",
    "jrepst_dev", "jmatst_dev", "indcova", "indcovb", "indcovc", "year2",
    "mpm_altered", "vrm_altered", "fitness", "converged"};
  CharacterVector processed_out_class = {"data.frame"};
  processed_fitness_out.attr("class") = processed_out_class;
  processed_fitness_out.attr("names") = processed_out_names;
  
  StringVector row_names(static_cast<int>(final_fitness.length()));
  for (int i = 0; i < static_cast<int>(final_fitness.length()); i++) {
    row_names(i) = std::to_string(i+1);
  }
  processed_fitness_out.attr("row.names") = row_names;
  
  out_df = processed_fitness_out;
}

//' Find ESS Values of Traits in Pre-Existing MPM Invasibility Analyses
//' 
//' This function attempts to estimate any ESS trait values in the trait axis.
//' These ESS trait values are initially identified as possible given flipped
//' signs of elasticity fitness at points along the input trait axis to use in
//' optimization.
//' 
//' @name ESS_optimizer_pre
//' 
//' @param ESS_Lyapunov A data frame to hold the ESS trait values, built by
//' reference.
//' @param ESS_trait_axis A data frame giving the minimum and maximum values of
//' variable traits, and the constant values of all others, in the trait axis.
//' @param Lyapunov_optim Main data frame giving Lyapunov coefficients for all
//' trait combinations developed for the ESS optima table. Holds elasticity
//' fitness values.
//' @param optim_trait_axis Main trait axis data frame corresponding to
//' \code{Lyapunov_optim}.
//' @param ESS_var_traits An integer vector modified by this function by
//' reference, indicating the actual traits that vary. The element order is:
//' 1, givenrate; 2, offset; 3, multiplier; 4, surv_dev; 5, obs_dev;
//' 6, size_dev; 7, sizeb_dev; 8, sizec_dev; 9, repst_dev; 10, fec_dev;
//' 11, jsurv_dev; 12, jobs_dev; 13, jsize_dev; 14, jsizeb_dev; 15, jsizec_dev;
//' 16, jrepst_dev; 17, jmatst_dev; 18, indcova; 19, indcovb; and 20, indcovc.
//' @param flipped_traits An integer vector with 17 elements, giving the
//' identities of variables with both positive and negative values.
//' @param new_stageexpansion_list A list with stage expansions for all
//' variant data used in ESS evaluation. This list includes an extra layer of
//' list elements, corresponding to the optim_ta and optim_ta_995 data.
//' @param used_times A list of year numbers for each time per run.
//' @param zero_stage_vec_list A list of population stage vectors full of zeros.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param A_list A list of all A matrices.
//' @param U_list A list of all U matrices.
//' @param F_list A list of all F matrices.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param stageframe_df A stageframe object covering the MPM.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param times An integer giving the number of occasions to project.
//' @param fitness_times AN integer giving the number of occasions at the end of
//' the run to use to estimate Lyapunov coefficients.
//' @param format_int An integer giving the MPM format.
//' @param stagecounts An integer giving the number of stages in the life
//' history.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' vital rate models (particulaly negative binomial models).
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' @param conv_threshold The convergence threshold value for Lyapunov
//' coefficients estimated in ESS optimization.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis. A relic value that equals the number
//' of variants in the entered trait axis.
//' @param opt_res_orig The original value of \code{opt_res}, prior to the
//' determination of the number of variable traits. Equal to \code{opt_res} if
//' the number of variable traits is 1, and to the square root of \code{opt_res}
//' if the number of variable traits is 2.
//' @param ehrlen An integer stating if historical MPMs should be in Ehrlen
//' format.
//' @param style An integer giving the style (e.g. ahistorical vs. historical)
//' of MPM.
//' @param loop_max The maximum number of times to attempt optimization, if
//' convergence does not occur.
//' @param filter An integer giving cleanup options for MPMs.
//' @param elast_mult A multiplier for traits to assess the elasticity of
//' fitness in trait optimization. Defaults to 0.995.
//' @param zap_min A Boolean value describing whether to round fitness values
//' below the value given in \code{threshold}.
//' 
//' @return A final data frame giving the zeros associated with variable traits.
//' 
//' @keywords internal
//' @noRd
inline void ESS_optimizer_pre (DataFrame& ESS_Lyapunov,
  DataFrame& ESS_trait_axis, DataFrame& Lyapunov_optim,
  DataFrame& optim_trait_axis, IntegerVector& ESS_var_traits,
  IntegerVector& flipped_traits, List& new_stageexpansion_list, List& used_times,
  List& zero_stage_vec_list, const List start_list, const List equivalence_list,
  const List A_list, const List U_list, const List F_list,
  const DataFrame density_df, const DataFrame dens_index_df,
  const DataFrame stageframe_df, const IntegerVector entry_time_vec,
  const int times, const int fitness_times, const int format_int,
  const int stagecounts, const int firstage_int, const int finalage_int,
  const int substoch, const double exp_tol, const double theta_tol,
  const bool sparse_bool, const bool A_only, const bool stages_not_equal,
  const bool integeronly, const bool dens_yn_bool, const double conv_threshold,
  const int opt_res, const int opt_res_orig, const int ehrlen, const int style,
  const int loop_max, const int filter, double elast_mult, const bool zap_min) {
  
  // Rcout << "ESS_optimizer_pre A" << endl;
  
  arma::uvec potential_optima (opt_res, fill::zeros);
  arma::vec inv_fitness = as<NumericVector>(Lyapunov_optim["fitness_variant2_e995"]);
  IntegerVector ESS_var_traits_corresponding_indices = {11, 12, 13, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
  int main_loop_breakpoint = static_cast<int>(inv_fitness.n_elem);
  
  // Rcout << "ESS_optimizer_pre B" << endl;
  for (int j = 0; j < (main_loop_breakpoint - 1); j++) {
    double base_inv_fit = inv_fitness(j);
    double next_inv_fit = inv_fitness(j + 1);
    
    if (base_inv_fit < 0. && next_inv_fit >= 0.) potential_optima(j) = 1;
    if (base_inv_fit >= 0. && next_inv_fit < 0.) potential_optima(j) = 1;
  }
  
  // Rcout << "ESS_optimizer_pre C" << endl;
  unsigned int total_optima = arma::sum(potential_optima);
  
  arma::uvec potential_optima_indices = find(potential_optima);
  
  arma::uvec ESS_var_traits_arma = as<arma::uvec>(ESS_var_traits);
  arma::uvec ESS_vta_indices = find(ESS_var_traits_arma);
  int vars_to_alter = static_cast<int>(ESS_vta_indices.n_elem);
  
  if (vars_to_alter > 1) {
    Rf_warningcall(R_NilValue, "Trait optimization may not run properly with more than 1 tested trait.");
  }
  // Rcout << "total_optima (potential): " << total_optima << endl;
  List final_output (total_optima);
  // Rcout << "ESS_optimizer_pre D" << endl;
  
  for (int i = 0; i < total_optima; i++) {
    // Rcout << "ESS_optimizer_pre E entered optimization trait loop i: " << i << endl;
    bool found_converged {false};
    
    arma::mat N_out;
    List comm_out (2);
    
    for (int m = 0; m < 2; m++) {
      arma::mat pops_out (stagecounts, (times + 1), fill::zeros);
      comm_out(m) = pops_out;
    }
    
    IntegerVector core_index_list = {(static_cast<int>(potential_optima_indices(i))), 
      (static_cast<int>(potential_optima_indices(i) + 1))};
    DataFrame core_trait_axis_instance = AdaptUtils::df_indices(optim_trait_axis, core_index_list);
    
    DataFrame reference_variants = clone(core_trait_axis_instance);
    NumericVector ref_fitness = {static_cast<double>(inv_fitness(potential_optima_indices(i))), 
      static_cast<double>(inv_fitness(potential_optima_indices(i) + 1))};
    // Rcout << "ref_fitness (lower bound invasion fitness, higher bound invasion fitness): " << ref_fitness << endl;
    reference_variants["fitness"] = ref_fitness;

    //DataFrame old_reference_variants = clone(reference_variants); // This should be remarked out
    bool opt_needed {false};
    
    if (abs(static_cast<double>(inv_fitness(potential_optima_indices(i)))) < conv_threshold) {
      // Rcout << "abs of lower bound invader fitness < threshold, so no further optimization needed" << endl;
      IntegerVector chosen_one = {0};
      DataFrame optim_point = AdaptUtils::df_indices(reference_variants, chosen_one);
      
      found_converged = true;
      LogicalVector found_converged_vec = {found_converged};
      optim_point["converged"] = found_converged_vec;
      
      final_output(i) = optim_point;
      //converged(i) = 1;
    } else if (abs(static_cast<double>(inv_fitness(potential_optima_indices(i) + 1))) < conv_threshold) {
      // Rcout << "abs of higher bound invader fitness < threshold, so no further optimization needed" << endl;
      IntegerVector chosen_one = {1};
      DataFrame optim_point = AdaptUtils::df_indices(reference_variants, chosen_one);
      
      found_converged = true;
      LogicalVector found_converged_vec = {found_converged};
      optim_point["converged"] = found_converged_vec;
      
      final_output(i) = optim_point;
      //converged(i) = 1;
    } else {
      // Rcout << "lower and upper bound invader fitness above threshold, so further optimization needed" << endl;
      opt_needed = true;
    }
    
    DataFrame variants_to_test = clone(reference_variants); // These values will be overwritten in the cloned data frame
    
    DataFrame variant_L_bound;
    DataFrame variant_H_bound;
    
    if (inv_fitness(potential_optima_indices(i)) < inv_fitness(potential_optima_indices(i) + 1)) {
      IntegerVector chosen_low = {0};
      IntegerVector chosen_high = {1};
      DataFrame variant_L_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_low);
      DataFrame variant_H_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_high);
      
      variant_L_bound = variant_L_bound_pre;
      variant_H_bound = variant_H_bound_pre;
    } else {
      IntegerVector chosen_low = {1};
      IntegerVector chosen_high = {0};
      DataFrame variant_L_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_low);
      DataFrame variant_H_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_high);
      
      variant_L_bound = variant_L_bound_pre;
      variant_H_bound = variant_H_bound_pre;
    }
    
    // Rcout << "ESS_optimizer_pre H" << endl;
    
    // Initial point set-up and fitness check
    DataFrame point_1; // Will be set to 1/3 distance from first bound to second bound
    DataFrame point_2; // Will be set to 2/3 distance from first bound to second bound
    
    { // Point 1 initialization
      // Rcout << "ESS_optimizer_pre I point_1 initialization " << endl;
      
      arma::uvec variant_nta_new = {1, 2};
      arma::vec givenrate_nta_new;
      arma::vec offset_nta_new;
      arma::vec multiplier_nta_new;
      
      NumericVector variants_vars_to_test;
      
      for (int j = 0; j < vars_to_alter; j++) {
        // Rcout << "ESS_optimizer_pre J" << endl;
        double base_mean = {0.};
        
        // New optim_ta
        NumericVector core_var = as<NumericVector>(core_trait_axis_instance(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        int flipped_trait_status = flipped_traits(ESS_vta_indices(j));
        double high_value = static_cast<double>(core_var(0));
        double low_value = static_cast<double>(core_var(1));
        
        // Rcout << "ESS_optimizer_pre K" << endl;
        NumericVector base_values = {high_value, low_value};
        double current_mult = 1./3.;
        base_mean = low_value + current_mult * (high_value - low_value);
        
        variants_vars_to_test = as<NumericVector>(variants_to_test(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        variants_vars_to_test(0) = base_mean;
        
        if (base_mean < 0. && flipped_trait_status == 1) {
          variants_vars_to_test(1) = base_mean + (base_mean - (base_mean * elast_mult));
        } else {
          variants_vars_to_test(1) = base_mean * elast_mult;
        }
      }
      // Rcout << "ESS_optimizer_pre L" << endl;
      
      givenrate_nta_new = as<arma::vec>(variants_to_test(11));
      offset_nta_new = as<arma::vec>(variants_to_test(12));
      multiplier_nta_new = as<arma::vec>(variants_to_test(13));
      
      // New stageexpansion
      IntegerVector chosen_int = {0};
      DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
      DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
      
      DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
        variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
        filter);
      DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
        variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
        filter);
      
      chosen_int = 1;
      StringVector focused_var = {"mpm_altered"};
      DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "ESS_optimizer_pre M" << endl;
      
      DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      
      List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
        _["e995"] = stageexpansion_995_reduced);
      
      // Rcout << "ESS_optimizer_pre N" << endl;
      DataFrame ESS_out_values;
      
      invpre_optim_singlerun (point_1, variants_to_test, N_out, comm_out, sge_to_test, used_times,
        zero_stage_vec_list, start_list, equivalence_list, A_list, U_list,
        F_list, density_df, dens_index_df, entry_time_vec, times, fitness_times, format_int,
        firstage_int, finalage_int, substoch, exp_tol, theta_tol, conv_threshold,
        sparse_bool, A_only, stages_not_equal, integeronly, dens_yn_bool, zap_min);
      
    }
    
    { // Point 2 initialization
      // Rcout << "ESS_optimizer_pre O point_2 initialization " << endl;
      
      arma::uvec variant_nta_new = {1, 2};
      arma::vec givenrate_nta_new;
      arma::vec offset_nta_new;
      arma::vec multiplier_nta_new;
      
      NumericVector variants_vars_to_test;
      
      for (int j = 0; j < vars_to_alter; j++) {
        // Rcout << "ESS_optimizer_pre P" << endl;
        double base_mean = {0.};
        
        // New optim_ta
        NumericVector core_var = as<NumericVector>(core_trait_axis_instance(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        int flipped_trait_status = flipped_traits(ESS_vta_indices(j));
        double high_value = static_cast<double>(core_var(0));
        double low_value = static_cast<double>(core_var(1));
        
        // Rcout << "ESS_optimizer_pre Q" << endl;
        NumericVector base_values = {high_value, low_value};
        double current_mult = 2./3.;
        base_mean = low_value + current_mult * (high_value - low_value);
        
        variants_vars_to_test = as<NumericVector>(variants_to_test(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        variants_vars_to_test(0) = base_mean;
        
        if (base_mean < 0. && flipped_trait_status == 1) {
          variants_vars_to_test(1) = base_mean + (base_mean - (base_mean * elast_mult));
        } else {
          variants_vars_to_test(1) = base_mean * elast_mult;
        }
      }
      // Rcout << "ESS_optimizer_pre R" << endl;
      
      givenrate_nta_new = as<arma::vec>(variants_to_test(11));
      offset_nta_new = as<arma::vec>(variants_to_test(12));
      multiplier_nta_new = as<arma::vec>(variants_to_test(13));
      
      // New stageexpansion
      IntegerVector chosen_int = {0};
      DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
      DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
      
      DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
        variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
        filter);
      DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
        variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
        filter);
      
      chosen_int = 1;
      StringVector focused_var = {"mpm_altered"};
      DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "ESS_optimizer_pre S" << endl;
      
      DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      
      List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
        _["e995"] = stageexpansion_995_reduced);
      
      // Rcout << "ESS_optimizer_pre T" << endl;
      DataFrame ESS_out_values;
      
      invpre_optim_singlerun (point_2, variants_to_test, N_out, comm_out, sge_to_test, used_times,
        zero_stage_vec_list, start_list, equivalence_list, A_list, U_list,
        F_list, density_df, dens_index_df, entry_time_vec, times, fitness_times, format_int,
        firstage_int, finalage_int, substoch, exp_tol, theta_tol, conv_threshold,
        sparse_bool, A_only, stages_not_equal, integeronly, dens_yn_bool, zap_min);
    }
    
    int search_mode {0};
    int loop_tracker {0};
    
    DataFrame variant_L; // This single row dataframe will mark the test point with the lower abs fitness
    DataFrame variant_H; // This single row dataframe will mark the test point with the higher abs fitness
    DataFrame variant_C; // Centroid
    
    int worst_point {2};
    double x1 {0.};
    double x2 {0.};
    double conv_test {0.};
    double xR {0.}; // test value to keep in memory from previous reflection round
    double xE {0.}; // test value to keep in memory from previous expansion round
    double xC_new {0.}; // test value to keep in memory from previous contraction round
    
    double xR_fitness {0.}; // fitness associated with test value to keep in memory from previous reflection round
    double xE_fitness {0.}; // fitness associated with test value to keep in memory from previous reflection round
    double xC_new_fitness {0.}; // fitness associated with test value to keep in memory from previous reflection round
    
    double e995_fitness_1 {0.};
    double e995_fitness_2 {0.};
    double e995_fitness_L {0.};
    double e995_fitness_H {0.};
    
    DataFrame variant_R;
    DataFrame variant_E;
    DataFrame variant_C_new;
    
    DataFrame ESS_test_values;
    DataFrame ESS_out_values;
    
    while (opt_needed && loop_tracker < loop_max) {
      // Rcout << "\n\n\nESS_optimizer_pre U entered optimization main loop, with loop_tracker = " << loop_tracker << endl;
      
      arma::uvec variant_nta_new = {1, 2};
      arma::vec givenrate_nta_new;
      arma::vec offset_nta_new;
      arma::vec multiplier_nta_new;
      
      NumericVector variants_vars_to_test;
      
      for (int j = 0; j < vars_to_alter; j++) {
        // Rcout << "ESS_optimizer_pre V j: " << j << endl;
        double base_mean = {0.};
        
        NumericVector point1_fitness_values = as<NumericVector>(point_1["fitness"]);
        NumericVector point2_fitness_values = as<NumericVector>(point_2["fitness"]);
        e995_fitness_1 = point1_fitness_values(1); 
        e995_fitness_2 = point2_fitness_values(1); 
        
        // Rcout << "ESS_optimizer_pre W" << endl;
        if (abs(e995_fitness_1) < abs(e995_fitness_2)) {
          variant_L = clone(point_1);
          variant_H = clone(point_2);
          e995_fitness_L = e995_fitness_1;
          e995_fitness_H = e995_fitness_2;
          worst_point = 2;
        } else {
          variant_H = clone(point_1);
          variant_L = clone(point_2);
          e995_fitness_L = e995_fitness_2;
          e995_fitness_H = e995_fitness_1;
          worst_point = 1;
        }
        
        // Rcout << "ESS_optimizer_pre X" << endl;
        variant_C = clone(variant_L);
        double base_test_value = {0.};
        
        // Rcout << "ESS_optimizer_pre X1" << endl;
        NumericVector core_var_C = as<NumericVector>(variant_C(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_pre X2" << endl;
        NumericVector core_var_H = as<NumericVector>(variant_H(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_pre X3" << endl;
        NumericVector core_var_L_bound = as<NumericVector>(variant_L_bound(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_pre X4" << endl;
        NumericVector core_var_H_bound = as<NumericVector>(variant_H_bound(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_pre X5" << endl;
        int flipped_trait_status = flipped_traits(ESS_vta_indices(j));
        // Rcout << "ESS_optimizer_pre X6" << endl;
        NumericVector point_1_trait = as<NumericVector>(point_1(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j)))); // Not necessary - used for testing
        // Rcout << "ESS_optimizer_pre X7" << endl;
        NumericVector point_2_trait = as<NumericVector>(point_2(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j)))); // Not necessary - used for testing
        // Rcout << "ESS_optimizer_pre X8" << endl;
        
        double xC = static_cast<double>(core_var_C(0));
        double xL = xC;
        double xH = static_cast<double>(core_var_H(0)); // Used to be index 1
        double low_value_bound = static_cast<double>(core_var_L_bound(0));
        double high_value_bound = static_cast<double>(core_var_H_bound(0)); // Used to be index 1
        
        x1 = xL;
        x2 = xH;
        conv_test = abs(x1 - x2);
        
        // Rcout << "\nESS_optimizer_pre Y search mode choice" << endl;
        // Rcout << "original low_value_bound: " << low_value_bound << endl;
        // Rcout << "original high_value_bound: " << high_value_bound << endl;
        // Rcout << "original point_1 trait value: " << point_1_trait << endl;
        // Rcout << "original point_2 trait value: " << point_2_trait << endl;
        // Rcout << "worst_point: " << worst_point << endl;
        // Rcout << "original xC: " << xC << endl;
        // Rcout << "original xL: " << xL << endl;
        // Rcout << "original xH: " << xH << endl;
        // Rcout << "original xR: " << xR << endl;
        // Rcout << "original xE: " << xE << endl;
        // Rcout << "original xC_new: " << xC_new << endl;
        // Rcout << "original point_1 invasion fitness: " << e995_fitness_1 << endl;
        // Rcout << "original point_2 invasion fitness: " << e995_fitness_2 << endl;
        // Rcout << "original e995_fitness_L: " << e995_fitness_L << endl;
        // Rcout << "original e995_fitness_H: " << e995_fitness_H << endl;
        // Rcout << "original xR_fitness: " << xR_fitness << endl;
        // Rcout << "original xE_fitness: " << xE_fitness << endl;
        // Rcout << "original xC_new_fitness: " << xC_new_fitness << endl;
        
        // Rcout << "ESS_optimizer_pre Z" << endl;
        //NumericVector base_values = {high_value, low_value};
        if (search_mode == 0) { // Reflection
          // Rcout << "\nReflection mode (search mode = 0)" << endl;
          
          double alpha = 0.8;
          base_test_value = xC + alpha * (xC - xH);
          
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          xR = base_test_value;
          
          // Rcout << "value of xC (used in base_test_value calculation): " << xC << endl;
          // Rcout << "value of xH (used in base_test_value calculation): " << xH << endl;
          // Rcout << "value of xR (calculated during reflection): " << xR << endl;
          
        } else if (search_mode == 1) { // Expansion
          // Rcout << "\nExpansion mode (search mode = 1)" << endl;
          
          double gamma = 2.0;
          base_test_value = xC + gamma * (xR - xC);
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          xE = base_test_value;
          
          // Rcout << "value of xC (used in base_test_value calculation): " << xC << endl;
          // Rcout << "value of xR (used in base_test_value calculation): " << xR << endl;
          // Rcout << "value of xE (calculated during expansion): " << xE << endl;
          
        } else if (search_mode == 2) { // Contraction
          // Rcout << "\nContraction mode (search mode = 2)" << endl;
          
          double rho = 0.5;
          
          if (worst_point == 1) {
            if (abs(xR_fitness) < abs(e995_fitness_1)) {
              base_test_value = xC + rho * (xR - xC);
            } else if (abs(xR_fitness) >= abs(e995_fitness_1)) {
              base_test_value = xC + rho * (xH - xC);
            }
          } else if (worst_point == 2) {
            if (abs(xR_fitness) < abs(e995_fitness_2)) {
              base_test_value = xC + rho * (xR - xC);
            } else if (abs(xR_fitness) >= abs(e995_fitness_2)) {
              base_test_value = xC + rho * (xH - xC);
            }
          }
          
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          xC_new = base_test_value;
          
          // Rcout << "value of xC (used in base_test_value calculation): " << xC << endl;
          // Rcout << "value of xH (used in base_test_value calculation): " << xH << endl;
          // Rcout << "value of xC_new (calculated during contraction): " << xC_new << endl;
          
        } else if (search_mode == 3) { // Shrinking
          // Rcout << "\nShrinking mode (search mode = 3)" << endl;
          
          double sigma = 0.5;
          
          NumericVector core_var_p1 = as<NumericVector>(point_1(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          NumericVector core_var_p2 = as<NumericVector>(point_2(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          
          double p1 = static_cast<double>(core_var_p1(0));
          double p2 = static_cast<double>(core_var_p2(0));
          
          base_test_value = xL + sigma * (p1 - xL);
          double base_test_value2 = xL + sigma * (p2 - xL);
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          
          NumericVector maximal_test2 = {base_test_value2, low_value_bound};
          NumericVector minimal_test2 = {maximal_test2(0), high_value_bound};
          
          base_test_value2 = minimal_test2(0);
          
          // Rcout << "value of xL (used in base_test_value calculation and base_test_value2): " << xL << endl;
          // Rcout << "value of p1 (used in base_test_value calculation and base_test_value2): " << p1 << endl;
          // Rcout << "value of p2 (used in base_test_value calculation and base_test_value2): " << p2 << endl;
          // Rcout << "value of base_test_value (calculated during shrinking): " << base_test_value << endl;
          // Rcout << "value of base_test_value2 (calculated during shrinking): " << base_test_value2 << endl;
          
          variants_vars_to_test = as<NumericVector>(variants_to_test(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          variants_vars_to_test(0) = base_test_value;
          
          if (base_test_value < 0. && flipped_trait_status == 1) {
            variants_vars_to_test(1) = base_test_value + (base_test_value - (base_test_value * elast_mult));
          } else {
            variants_vars_to_test(1) = base_test_value * elast_mult;
          }
          
          givenrate_nta_new = as<arma::vec>(variants_to_test(11));
          offset_nta_new = as<arma::vec>(variants_to_test(12));
          multiplier_nta_new = as<arma::vec>(variants_to_test(13));
          
          // New stageexpansion
          IntegerVector chosen_int = {0};
          DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
          DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
          
          DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
            filter);
          DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
            filter);
          
          chosen_int = 1;
          StringVector focused_var = {"mpm_altered"};
          DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_pre Z1" << endl;
          
          DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          
          List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
            _["e995"] = stageexpansion_995_reduced);
          
          // Rcout << "ESS_optimizer_pre Z2" << endl;
          DataFrame ESS_out_values;
          
          invpre_optim_singlerun (point_1, variants_to_test, N_out, comm_out, sge_to_test, used_times,
            zero_stage_vec_list, start_list, equivalence_list, A_list, U_list,
            F_list, density_df, dens_index_df, entry_time_vec, times, fitness_times, format_int,
            firstage_int, finalage_int, substoch, exp_tol, theta_tol, conv_threshold,
            sparse_bool, A_only, stages_not_equal, integeronly, dens_yn_bool, zap_min);
          
          variants_vars_to_test = as<NumericVector>(variants_to_test(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          variants_vars_to_test(0) = base_test_value2;
          
          if (base_test_value2 < 0. && flipped_trait_status == 1) {
            variants_vars_to_test(1) = base_test_value2 + (base_test_value2 - (base_test_value2 * elast_mult));
          } else {
            variants_vars_to_test(1) = base_test_value2 * elast_mult;
          }
          
          givenrate_nta_new = as<arma::vec>(variants_to_test(11));
          offset_nta_new = as<arma::vec>(variants_to_test(12));
          multiplier_nta_new = as<arma::vec>(variants_to_test(13));
          
          // New stageexpansion
          chosen_int = {0};
          variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
          variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
          
          stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
            filter);
          stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
            filter);
          
          chosen_int = 1;
          focused_var = {"mpm_altered"};
          stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_pre Z3" << endl;
          
          stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          
          sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
            _["e995"] = stageexpansion_995_reduced);
          
          // Rcout << "ESS_optimizer_pre Z4" << endl;
          invpre_optim_singlerun (point_2, variants_to_test, N_out, comm_out, sge_to_test, used_times,
            zero_stage_vec_list, start_list, equivalence_list, A_list, U_list,
            F_list, density_df, dens_index_df, entry_time_vec, times, fitness_times, format_int,
            firstage_int, finalage_int, substoch, exp_tol, theta_tol, conv_threshold,
            sparse_bool, A_only, stages_not_equal, integeronly, dens_yn_bool, zap_min);
          
        }
        
        variants_vars_to_test = as<NumericVector>(variants_to_test(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        variants_vars_to_test(0) = base_mean;
        
        if (base_mean < 0. && flipped_trait_status == 1) {
          variants_vars_to_test(1) = base_mean + (base_mean - (base_mean * elast_mult));
        } else {
          variants_vars_to_test(1) = base_mean * elast_mult;
        }
      //}
        // Rcout << "ESS_optimizer_pre AA" << endl;
        
        if (search_mode < 3) {
          // Rcout << "ESS_optimizer_pre AB" << endl;
          givenrate_nta_new = as<arma::vec>(variants_to_test(11));
          offset_nta_new = as<arma::vec>(variants_to_test(12));
          multiplier_nta_new = as<arma::vec>(variants_to_test(13));
          
          givenrate_nta_new = as<arma::vec>(variants_to_test(11));
          offset_nta_new = as<arma::vec>(variants_to_test(12));
          multiplier_nta_new = as<arma::vec>(variants_to_test(13));
          
          // New stageexpansion
          IntegerVector chosen_int = {0};
          DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
          DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
          
          DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
            filter);
          DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
            filter);
          
          chosen_int = 1;
          StringVector focused_var = {"mpm_altered"};
          DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_pre AC" << endl;
          
          DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          
          List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
            _["e995"] = stageexpansion_995_reduced);
          
          // Rcout << "ESS_optimizer_pre AD" << endl;
          DataFrame ESS_out_values;
          
          invpre_optim_singlerun (ESS_out_values, variants_to_test, N_out, comm_out, sge_to_test, used_times,
            zero_stage_vec_list, start_list, equivalence_list, A_list, U_list,
            F_list, density_df, dens_index_df, entry_time_vec, times, fitness_times, format_int,
            firstage_int, finalage_int, substoch, exp_tol, theta_tol, conv_threshold,
            sparse_bool, A_only, stages_not_equal, integeronly, dens_yn_bool, zap_min);
          
          // Rcout << "ESS_optimizer_fb T" << endl;
          NumericVector current_round_fitness_values = as<NumericVector>(ESS_test_values["fitness"]);
          // Rcout << "current_round_fitness_values: " << current_round_fitness_values << endl;
          double res_fitness_test = current_round_fitness_values(0);
          double e995_fitness_test = current_round_fitness_values(1); 
          // Rcout << "new invader res_fitness_test: " << res_fitness_test << endl;
          // Rcout << "new invader e995_fitness_test: " << e995_fitness_test << endl;
          
          NumericVector ref_fitness_values = as<NumericVector>(reference_variants["fitness"]);
          NumericVector abs_fitness_values = abs(ref_fitness_values);
          
          LogicalVector converged = as<LogicalVector>(ESS_test_values["converged"]);
          // Rcout << "ref_fitness_values: " << ref_fitness_values << endl;
          // Rcout << "abs_fitness_values: " << abs_fitness_values << endl;
          
          // Rcout << "Entered search mode determination section" << endl;
          if (search_mode == 0) { // Reflection
            variant_R = clone(ESS_test_values);
            xR_fitness = e995_fitness_test;
            // Rcout << "xR_fitness: " << xR_fitness << endl;
            // Rcout << "e995_fitness_L: " << e995_fitness_L << endl;
            // Rcout << "e995_fitness_H: " << e995_fitness_H << endl;
            if (abs(xR_fitness) < abs(e995_fitness_L)) {
              // Rcout << "ESS_optimizer_pre AD1 search mode 0, moving to search mode 1" << endl;
              search_mode = 1;
            } else if (abs(xR_fitness) < abs(e995_fitness_H)) {
              // Rcout << "ESS_optimizer_pre AD2 search mode 0, staying in search mode 0" << endl;
              variant_H = clone(ESS_test_values);
              e995_fitness_H = xR_fitness;
              
              if (worst_point == 2) {
                point_2 = clone(variant_H);
                e995_fitness_2 = e995_fitness_H;
              } else if (worst_point == 1) {
                point_1 = clone(variant_H);
                e995_fitness_1 = e995_fitness_H;
              }
            } else {
              // Rcout << "ESS_optimizer_pre AD3 search mode 0, moving to search mode 2" << endl;
              search_mode = 2;
            }
            
          } else if (search_mode == 1) { // Expansion
            // Rcout << "ESS_optimizer_pre AD4 search mode 1, moving to search mode 0" << endl;
            variant_E = clone(ESS_test_values);
            xE_fitness = e995_fitness_test;
            // Rcout << "xE_fitness: " << xE_fitness << endl;
            // Rcout << "xR_fitness: " << xR_fitness << endl;
            // Rcout << "e995_fitness_L: " << e995_fitness_L << endl;
            // Rcout << "e995_fitness_H: " << e995_fitness_H << endl;
            if (abs(xE_fitness) < abs(xR_fitness)) {
              // Rcout << "abs(xE_fitness) < abs(xR_fitness)" << endl;
              if (worst_point == 2) {
                // Rcout << "Replacing point_2 with xE" << endl;
                point_2 = clone(variant_E);
                e995_fitness_2 = xE_fitness;
              } else if (worst_point == 1) {
                // Rcout << "Replacing point_1 with xE" << endl;
                point_1 = clone(variant_E);
                e995_fitness_1 = xE_fitness;
              }
              
              //variant_H = clone(variant_E);
              //e995_fitness_H = xE_fitness;
            } else {
              // Rcout << "abs(xE_fitness) >= abs(xR_fitness)" << endl;
              if (worst_point == 2) {
                // Rcout << "Replacing point_2 with xR" << endl;
                point_2 = clone(variant_R);
                e995_fitness_2 = xR_fitness;
              } else if (worst_point == 1) {
                // Rcout << "Replacing point_1 with xR" << endl;
                point_1 = clone(variant_R);
                e995_fitness_1 = xR_fitness;
              }
              
              //variant_H = clone(variant_R);
              //e995_fitness_H = xR_fitness;
            }
            search_mode = 0;
            
          } else if (search_mode == 2) { // Contraction
            variant_C_new = clone(ESS_test_values);
            xC_new_fitness = e995_fitness_test;
            // Rcout << "xC_new_fitness: " << xC_new_fitness << endl;
            // Rcout << "e995_fitness_L: " << e995_fitness_L << endl;
            // Rcout << "e995_fitness_H: " << e995_fitness_H << endl;
            
            if (worst_point == 1) {
              if ((abs(xC_new_fitness) < abs(xR_fitness)) && (abs(xR_fitness) < abs(e995_fitness_1))) {
                // Rcout << "ESS_optimizer_pre AD5 search mode 2, moving to search mode 0" << endl;
                point_1 = clone(variant_C_new);
                e995_fitness_1 = xC_new_fitness;
                search_mode = 0;
              } else if ((abs(xC_new_fitness) < abs(e995_fitness_1)) && (abs(xR_fitness) >= abs(e995_fitness_1))) {
                // Rcout << "ESS_optimizer_pre AD6 search mode 2, moving to search mode 0" << endl;
                point_1 = clone(variant_C_new);
                e995_fitness_1 = xC_new_fitness;
                search_mode = 0;
              } else search_mode = 3;
            } else if (worst_point == 2) {
              if ((abs(xC_new_fitness) < abs(xR_fitness)) && (abs(xR_fitness) < abs(e995_fitness_2))) {
                // Rcout << "ESS_optimizer_pre AD7 search mode 2, moving to search mode 0" << endl;
                point_2 = clone(variant_C_new);
                e995_fitness_2 = xC_new_fitness;
                search_mode = 0;
              } else if ((abs(xC_new_fitness) < abs(e995_fitness_2)) && (abs(xR_fitness) >= abs(e995_fitness_2))) {
                // Rcout << "ESS_optimizer_pre AD8 search mode 2, moving to search mode 0" << endl;
                point_2 = clone(variant_C_new);
                e995_fitness_2 = xC_new_fitness;
                search_mode = 0;
              } else search_mode = 3;
            }
          }
          
        } else if (search_mode == 3) {
          // Rcout << "ESS_optimizer_pre AD9 search mode 3, moving to search mode 0" << endl;
          search_mode = 0;
          
          NumericVector point1_fitness_values = as<NumericVector>(point_1["fitness"]);
          NumericVector point2_fitness_values = as<NumericVector>(point_2["fitness"]);
          e995_fitness_1 = point1_fitness_values(1); 
          e995_fitness_2 = point2_fitness_values(1); 
          
          if (abs(e995_fitness_1) < abs(e995_fitness_2)) {
            variant_L = clone(point_1);
            variant_H = clone(point_2);
            e995_fitness_L = e995_fitness_1;
            e995_fitness_H = e995_fitness_2;
          } else {
            variant_H = clone(point_1);
            variant_L = clone(point_2);
            e995_fitness_L = e995_fitness_2;
            e995_fitness_H = e995_fitness_1;
          }
            
          variant_C = clone(variant_L);
        }
        
        // Rcout << "ESS_optimizer_pre AD10" << endl;
      }
      
      loop_tracker++;
      // Rcout << "ESS_optimizer_pre AE " << endl;
      
      if (abs(e995_fitness_1) <= conv_threshold || abs(e995_fitness_2) <= conv_threshold ||
          abs(e995_fitness_1 - e995_fitness_2) <= conv_threshold) {
        //bool found_optimum {true};
        if (abs(e995_fitness_1) <= conv_threshold || abs(e995_fitness_1 - e995_fitness_2) <= conv_threshold) {
          // Rcout << "Exporting point 1 on convergence" << endl;
          ESS_test_values = clone(point_1);
        } else {
          // Rcout << "Exporting point 2 on convergence" << endl;
          ESS_test_values = clone(point_2);
        }
        
        IntegerVector next_chosen_one = {0};
        ESS_out_values = AdaptUtils::df_indices(ESS_test_values, next_chosen_one);
        IntegerVector further_one = {1};
        DataFrame NeededFitnessRow = AdaptUtils::df_indices(ESS_test_values, further_one);
        
        NumericVector true_fitness_vec = as<NumericVector>(NeededFitnessRow["fitness"]);
        ESS_out_values["fitness"] = true_fitness_vec;
        
        found_converged = true;
        LogicalVector found_converged_vec = {found_converged};
        ESS_out_values["converged"] = found_converged_vec;
        
        opt_needed = false;
        break;
      }
    }
    // Rcout << "ESS_optimizer_pre AF" << endl;
    //IntegerVector next_chosen_one = {0};
    //ESS_out_values = AdaptUtils::df_indices(ESS_test_values, next_chosen_one);
    
    final_output(i) = ESS_out_values;
  }
  
  // Rcout << "ESS_optimizer_pre AI" << endl;
  DataFrame final_output_df;
  
  if (total_optima > 1) {
    DataFrame df1 = as<DataFrame>(final_output(0));
    DataFrame df2 = as<DataFrame>(final_output(1));
    CharacterVector df1_names = as<CharacterVector>(df1.names());
    CharacterVector df2_names = as<CharacterVector>(df2.names());
    
    //int df1_rows = static_cast<int>(df1.nrows());
    //int df2_rows = static_cast<int>(df2.nrows());
    
    DataFrame final_output_df_pre = LefkoUtils::df_rbind(as<DataFrame>(final_output(0)), as<DataFrame>(final_output(1)));
    
    // Rcout << "ESS_optimizer_pre AJ" << endl;
    if (total_optima > 2) {
      for (int i = 2; i < total_optima; i++) {
        final_output_df_pre = LefkoUtils::df_rbind(final_output_df_pre, as<DataFrame>(final_output(i)));
      }
    }
    
    // Rcout << "ESS_optimizer_pre AK" << endl;
    int fodfp_size = final_output_df_pre.nrows();
    IntegerVector new_variant_index = seq(1, fodfp_size);
    final_output_df_pre["variant"] = new_variant_index;
    
    final_output_df = final_output_df_pre;
  } else if (total_optima == 1) {
    DataFrame final_output_df_pre = as<DataFrame>(final_output(0));
    
    int fodfp_size = final_output_df_pre.nrows();
    IntegerVector new_variant_index = seq(1, fodfp_size);
    final_output_df_pre["variant"] = new_variant_index;
    
    final_output_df = final_output_df_pre;
  } else {
    final_output_df = R_NilValue;
  }
  
  int optim_rows = static_cast<int>(final_output_df.nrows());
  
  ESS_Lyapunov = final_output_df;
}

//' Core Function-Based Projection Engine for ESS Evaluation
//' 
//' Function \code{invfb_optim_singlerun} runs single simulation instances of the
//' optimization projections in function \code{invade3_fb_core}, and is used to
//' estimate ESS values.
//' 
//' @name invfb_optim_singlerun
//' 
//' @param out_df A data frame created by reference to hold the fitness values
//' produced by the run.
//' @param surv_dev_nta The survival column in the reassessed trait axis.
//' @param obs_dev_nta The observation status column in the reassessed trait
//' axis.
//' @param size_dev_nta The primary size column in the reassessed trait axis.
//' @param sizeb_dev_nta The secondary size column in the reassessed trait axis.
//' @param sizec_dev_nta The tertiary size column in the reassessed trait axis.
//' @param repst_dev_nta The reproductive status column in the reassessed trait
//' axis.
//' @param fec_dev_nta The fecundity column in the reassessed trait axis.
//' @param jsurv_dev_nta The juvenile survival column in the reassessed trait
//' axis.
//' @param jobs_dev_nta The juvenile observation status column in the reassessed
//' trait axis.
//' @param jsize_dev_nta The juvenile primary size column in the reassessed
//' trait axis.
//' @param jsizeb_dev_nta The juvenile secondary size column in the reassessed
//' trait axis.
//' @param jsizec_dev_nta The juvenile tertiary size column in the reassessed
//' trait axis.
//' @param jrepst_dev_nta The juvenile reproductive status column in the
//' reassessed trait axis.
//' @param jmatst_dev_nta The juvenile maturity status column in the reassessed
//' trait axis.
//' @param indcova_nta The individual covariate a column in the reassessed
//' trait axis.
//' @param indcovb_nta The individual covariate b column in the reassessed
//' trait axis.
//' @param indcovc_nta The individual covariate c column in the reassessed
//' trait axis.
//' @param variant_nta The variant column in the reassessed 995 trait axis.
//' @param base_trait_axis The currently used trait axis.
//' @param N_out A matrix giving the final population sizes of the run. Supplied
//' as a reference.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function. Needs to be essentially
//' blank to work properly. Needs to me a lst with only one level, with elements
//' as matrices.
//' @param errcheck_mpm A list of all MPMs post-processing. Only output if
//' \code{err_check = "extreme"}.
//' @param errcheck_mpmout A list of all mpm_out matrices from MPM processing. 
//' Only output if \code{err_check = "extreme"}.
//' @param new_stageexpansion_list A list with stage expansions for all
//' variant data used in ESS evaluation. This list includes an extra layer of
//' list elements, corresponding to the optim_ta and optim_ta_995 data.
//' @param used_times A list of year numbers for each time per run.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param vrm_list A list of \code{vrm_input} objects.
//' @param allstages_all The allstages indexing data frame used to produce MPMs.
//' @param dev_terms_list List of deviations for vital rate models.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates.
//' @param stageexpansion_ta_devterms_by_variant A list giving trait axis info
//' by variant, with each variant given a list element.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param density_vr_list Data frame of \code{lefkoDensVR} objects holding
//' density relationships for all 14 vital rate models.
//' @param current_stageframe The main stageframe, including extra stages.
//' @param current_supplement A supplement in \code{lefkoSD} format.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms.
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
//' @param fecmod_vec A numeric vector giving the fecmod values.
//' @param year_vec A vector giving the main years used.
//' @param patch_vec A vector giving the name of each patch used in projection.
//' @param times An integer giving the number of occasions to project.
//' @param fitness_times An integer giving the number of occasions at the end of
//' each run to use to estimate Lyapunov coefficients.
//' @param format_int An integer giving the MPM format.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param dev_terms_times_int A vector giving the number of occasions over
//' which vital rate y-intercept deviations cycle.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis. Note that the version used in this
//' function is actually equivalent to \code{opt_res_true}.
//' @param opt_res_orig The original value of \code{opt_res}, prior to the
//' determination of the number of variable traits. Equal to \code{opt_res} if
//' the number of variable traits is 1, and to the square root of \code{opt_res}
//' if the number of variable traits is 2.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param conv_threshold The lower limit for the absolute value of fitness,
//' below which fitness is rounded to 0. Defaults to 0.00000001.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param zap_min A Boolean value describing whether to round fitness values
//' below the value given in \code{threshold}.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
inline void invfb_optim_singlerun (DataFrame& out_df, arma::vec& surv_dev_nta,
  arma::vec& obs_dev_nta, arma::vec& size_dev_nta, arma::vec& sizeb_dev_nta,
  arma::vec& sizec_dev_nta, arma::vec& repst_dev_nta, arma::vec& fec_dev_nta,
  arma::vec& jsurv_dev_nta, arma::vec& jobs_dev_nta, arma::vec& jsize_dev_nta,
  arma::vec& jsizeb_dev_nta, arma::vec& jsizec_dev_nta,
  arma::vec& jrepst_dev_nta, arma::vec& jmatst_dev_nta, arma::vec& indcova_nta,
  arma::vec& indcovb_nta, arma::vec& indcovc_nta, arma::uvec& variant_nta,
  DataFrame& base_trait_axis, arma::mat& N_out, List& comm_out,
  List& errcheck_mpm, List& errcheck_mpmout, List& new_stageexpansion_list,
  List& used_times, const List allmodels_all, const List vrm_list,
  const List allstages_all, const List dev_terms_list,
  const List ind_terms_num_list, const List ind_terms_cat_list,
  List& stageexpansion_ta_devterms_by_variant, const List sp_density_list,
  const List start_list, const List equivalence_list,
  const DataFrame density_vr_list, DataFrame& current_stageframe,
  const DataFrame current_supplement, const DataFrame density_df,
  const DataFrame dens_index_df, const IntegerVector entry_time_vec,
  const IntegerVector sp_density_num_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_vr_yn_vec,
  const NumericVector fecmod_vec, const CharacterVector year_vec,
  const CharacterVector patch_vec, const int times, const int fitness_times,
  const int format_int, const int firstage_int, const int finalage_int,
  const int dev_terms_times_int, const int substoch, const int opt_res,
  const int opt_res_orig, const double exp_tol, const double theta_tol,
  const double conv_threshold, const bool sparse_bool, const bool A_only,
  const bool stages_not_equal, const bool integeronly, const bool dens_yn_bool,
  const bool err_check, const bool zap_min) {
  
  // Rcout << "ifosr A" << endl;
  
  int var_to_run {2};
  List running_popvecs; //  = clone(start_list)
  List running_popvecs_startonly; //  = clone(start_list)
  arma::mat N_mpm (2, (times + 1)); // rows = vars, cols = times
  
  //int year_counter {0};
  CharacterVector stage3_nta = as<CharacterVector>(base_trait_axis["stage3"]);
  CharacterVector stage2_nta = as<CharacterVector>(base_trait_axis["stage2"]);
  CharacterVector stage1_nta = as<CharacterVector>(base_trait_axis["stage1"]);
  IntegerVector age3_nta = as<IntegerVector>(base_trait_axis["age3"]);
  IntegerVector age2_nta = as<IntegerVector>(base_trait_axis["age2"]);
  CharacterVector eststage3_nta = as<CharacterVector>(base_trait_axis["eststage3"]);
  CharacterVector eststage2_nta = as<CharacterVector>(base_trait_axis["eststage2"]);
  CharacterVector eststage1_nta = as<CharacterVector>(base_trait_axis["eststage1"]);
  IntegerVector estage3_nta = as<IntegerVector>(base_trait_axis["estage3"]);
  IntegerVector estage2_nta = as<IntegerVector>(base_trait_axis["estage2"]);
  NumericVector givenrate_nta = as<NumericVector>(base_trait_axis["givenrate"]);
  NumericVector offset_nta = as<NumericVector>(base_trait_axis["offset"]);
  NumericVector multiplier_nta = as<NumericVector>(base_trait_axis["multiplier"]);
  IntegerVector convtype_nta = as<IntegerVector>(base_trait_axis["convtype"]);
  IntegerVector convtype_t12_nta = as<IntegerVector>(base_trait_axis["convtype_t12"]);
  CharacterVector year2_nta = as<CharacterVector>(base_trait_axis["year2"]);
  IntegerVector mpm_altered_nta = as<IntegerVector>(base_trait_axis["mpm_altered"]);
  IntegerVector vrm_altered_nta = as<IntegerVector>(base_trait_axis["vrm_altered"]);
  
  IntegerVector inda_num_terms_counter (var_to_run);
  IntegerVector indb_num_terms_counter (var_to_run);
  IntegerVector indc_num_terms_counter (var_to_run);
  IntegerVector inda_cat_terms_counter (var_to_run);
  IntegerVector indb_cat_terms_counter (var_to_run);
  IntegerVector indc_cat_terms_counter (var_to_run);
  IntegerVector inda_num_terms_previous (var_to_run);
  IntegerVector indb_num_terms_previous (var_to_run);
  IntegerVector indc_num_terms_previous (var_to_run);
  IntegerVector inda_cat_terms_previous (var_to_run);
  IntegerVector indb_cat_terms_previous (var_to_run);
  IntegerVector indc_cat_terms_previous (var_to_run);
  IntegerVector dev_num_counter (var_to_run);
  IntegerVector sp_density_counter (var_to_run);
  
  List errcheck_mpm_times (times);
  List errcheck_mpmout_times (times);
  
  for (int j = 0; j < times; j++) { // 2nd loop - time j
    // Rcout << "ifosr A1          ";
      if (j % 10 == 0){
        Rcpp::checkUserInterrupt();
      }
      
      if (j == 0) {
        List var_popvecs_to_start (var_to_run);
        for (int n = 0; n < var_to_run; n++) {
          var_popvecs_to_start(n) = as<arma::vec>(start_list(0));
        }
        running_popvecs = var_popvecs_to_start;
        running_popvecs_startonly = clone(var_popvecs_to_start);
      }
      
      for (int m = 0; m < var_to_run; m++) { // 4th loop - var per run m
        if (j == entry_time_vec(m)) {
          arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(0));
          
          double N_current = accu(running_popvec_mpm);
          N_mpm(m, j) = N_current;
        }
      }
      
      List errcheck_mpm_times_vtr (var_to_run);
      List errcheck_mpmout_times_vtr (var_to_run);
      
      for (int m = 0; m < var_to_run; m++) { // Main variant section
        // Rcout << "ifosr A2 m: " << m << "          ";
        DataFrame sge_current = as<DataFrame>(new_stageexpansion_list(m)); // new_stageexpansion_list to be subset & ordered
        arma::mat pops_out = as<arma::mat>(comm_out(m));
        
        // Rcout << "ifosr A3          ";
        
        if (j > (entry_time_vec(m) - 1)) {
          List used_times_higher = as<List>(used_times(0));
          List used_times_lower = as<List>(used_times_higher(0));
          IntegerVector current_times_vec = as<IntegerVector>(used_times_lower(0)); // Maybe m instead of 0?

          arma::vec running_popvec_vrm;
          if (j == entry_time_vec(m)) {
            running_popvec_vrm = as<arma::vec>(running_popvecs_startonly(m));
            pops_out.col(j) = running_popvec_vrm;
          } else {
            running_popvec_vrm = pops_out.col(j);
          }

          List current_vrm_extract = allmodels_all; // (i)
          List current_vrm_unextract = vrm_list; // (i)
          //int ehrlen_format {1}; // This will need to be dealt with differently later
          
          // Rcout << "ifosr A7       ";
          //int mpm_style {1};
          //if (format_int < 3) {
          //  mpm_style = 0;
          //  if (format_int == 2) ehrlen_format = 2;
          //} else if (format_int == 4) {
          //  mpm_style = 2;
          //}
          
          // Rcout << "ifosr A8       ";
          DataFrame current_mpm_allstages = allstages_all; // (i)
          
          // Rcout << "ifosr A9       ";
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
          
          // Rcout << "ifosr A10       ";
          CharacterVector current_mainyears = year_vec;
          //unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
          
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
          
          // Not sure if we need the next bit
          DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
          DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
          DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
          CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
          CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
          CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
          
          // Rcout << "ifosr A12        ";
          // Counter resets
          int yearnumber = current_times_vec(j); // year_counter
          CharacterVector current_year = as<CharacterVector>(current_mainyears(yearnumber));
          
          if (inda_num_terms_counter(m) >= inda_terms_num_vec(0)) inda_num_terms_counter(m) = 0;
          if (indb_num_terms_counter(m) >= indb_terms_num_vec(0)) indb_num_terms_counter(m) = 0;
          if (indc_num_terms_counter(m) >= indc_terms_num_vec(0)) indc_num_terms_counter(m) = 0;
          if (inda_cat_terms_counter(m) >= inda_terms_cat_vec(0)) inda_cat_terms_counter(m) = 0;
          if (indb_cat_terms_counter(m) >= indb_terms_cat_vec(0)) indb_cat_terms_counter(m) = 0;
          if (indc_cat_terms_counter(m) >= indc_terms_cat_vec(0)) indc_cat_terms_counter(m) = 0;
          
          List current_ind_terms_num = ind_terms_num_list(0);
          List current_ind_terms_cat = ind_terms_cat_list(0);
          
          NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
          NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
          NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
          CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
          CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
          CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
          
          NumericVector f2_inda = {f_inda_full(inda_num_terms_counter(m))}; // i
          NumericVector f1_inda = {f_inda_full(inda_num_terms_previous(m))};
          NumericVector f2_indb = {f_indb_full(indb_num_terms_counter(m))};
          NumericVector f1_indb = {f_indb_full(indb_num_terms_previous(m))};
          NumericVector f2_indc = {f_indc_full(indc_num_terms_counter(m))};
          NumericVector f1_indc = {f_indc_full(indc_num_terms_previous(m))};
          CharacterVector r2_inda = as<CharacterVector>(r_inda_full(inda_cat_terms_counter(m)));
          CharacterVector r1_inda = 
            as<CharacterVector>(r_inda_full(inda_cat_terms_previous(m)));
          CharacterVector r2_indb = as<CharacterVector>
            (r_indb_full(indb_cat_terms_counter(m)));
          CharacterVector r1_indb = as<CharacterVector>
            (r_indb_full(indb_cat_terms_previous(m)));
          CharacterVector r2_indc = as<CharacterVector>
            (r_indc_full(indc_cat_terms_counter(m)));
          CharacterVector r1_indc = 
            as<CharacterVector>(r_indc_full(indc_cat_terms_previous(m)));
          
          // dev_terms and vrm trait axis processing
          NumericVector dv_terms (17);
          arma::uvec var_corresponding_elems;
          int vce_found {0};
          if (dev_terms_times_int > 0) {
            NumericMatrix used_dv_df = clone(as<NumericMatrix>(dev_terms_list(0)));
            if (dev_num_counter(m) >= dev_terms_times_int) dev_num_counter(m) = 0;
            dv_terms = used_dv_df(_, dev_num_counter(m));
          }
          
          var_corresponding_elems = find(variant_nta == (m + 1));
          vce_found = static_cast<int>(var_corresponding_elems.n_elem);
          
          if (vce_found > 0) {
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
            arma::vec indcova_nta_sub = indcova_nta.elem(var_corresponding_elems);
            arma::vec indcovb_nta_sub = indcovb_nta.elem(var_corresponding_elems);
            arma::vec indcovc_nta_sub = indcovc_nta.elem(var_corresponding_elems);
            
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
              if(!NumericVector::is_na(indcova_nta_sub(vce_track))) dv_terms(14) =
                  dv_terms(14) + indcova_nta_sub(vce_track);
              if(!NumericVector::is_na(indcovb_nta_sub(vce_track))) dv_terms(15) =
                  dv_terms(15) + indcovb_nta_sub(vce_track);
              if(!NumericVector::is_na(indcovc_nta_sub(vce_track))) dv_terms(16) =
                  dv_terms(16) + indcovc_nta_sub(vce_track);
            }
          }
          dev_num_counter(m) = dev_num_counter(m) + 1;
          
          bool dvr_bool {false};
          
          LogicalVector dvr_yn = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          IntegerVector dvr_style (14);
          IntegerVector dvr_time_delay (14);
          NumericVector dvr_alpha (14);
          NumericVector dvr_beta (14);
          NumericVector dens_n (14);
          
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
                arma::vec delay_issue = N_mpm.col(j + 1 - used_delay);
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
          
          double dens_sp {1.0};
          
          if (sp_density_num_vec(0) > 0) {
            if (sp_density_counter(m) >= sp_density_num_vec(0)) sp_density_counter(m) = 0;
            
            NumericVector current_sp_density = as<NumericVector>(sp_density_list(0));
            dens_sp = current_sp_density(sp_density_counter(m));
            
            sp_density_counter(m) = sp_density_counter(m) + 1;
          }
          
          List current_mpm;
          if (format_int < 5) {
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
            
            if (err_check) errcheck_mpm_times_vtr(m) = current_mpm;
            // Rcout << "ifosr A29        ";
          } else {
            IntegerVector all_ages = seq(firstage_int, finalage_int);
            if (!(current_supplement.length() > 1)) {
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages,
                current_stageframe, surv_proxy, fec_proxy, f2_inda, f1_inda,
                f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                r1_indb, r2_indc, r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp,
                fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                exp_tol, theta_tol, sparse_bool);
            if (err_check) errcheck_mpm_times_vtr(m) = current_mpm;
              // Rcout << "ifosr A30        ";
              
            } else {
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages,
                current_stageframe, surv_proxy, fec_proxy, f2_inda, f1_inda,
                f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                r1_indb, r2_indc, r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp,
                fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                exp_tol, theta_tol, sparse_bool, current_supplement);
              if (err_check) errcheck_mpm_times_vtr(m) = current_mpm;
              // Rcout << "ifosr A32        ";
            }
          }
          // Rcout << "ifosr A33        ";
          
          if (!dens_yn_bool) {
            if (A_only) {
              if (!sparse_bool) {
                arma::mat current_A = as<arma::mat>(current_mpm["A"]);
                Amat_alter(current_A, sge_current); 
                
                if (err_check) errcheck_mpmout_times_vtr(m) = current_A;
                running_popvec_vrm = current_A * running_popvec_vrm; 
              } else {
                arma::sp_mat current_A = as<arma::sp_mat>(current_mpm["A"]);
                AdaptUtils::sp_Amat_alter(current_A, sge_current);
                
                if (err_check) errcheck_mpmout_times_vtr(m) = current_A;
                running_popvec_vrm = current_A * running_popvec_vrm;
              }
            } else {
              if (!sparse_bool) {
                arma::mat current_A = as<arma::mat>(current_mpm["A"]);
                arma::mat current_U_unaltered = as<arma::mat>(current_mpm["U"]);
                arma::mat current_F_unaltered = as<arma::mat>(current_mpm["F"]);
                AdaptUtils::UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
                
                if (err_check) errcheck_mpmout_times_vtr(m) = current_A;
                running_popvec_vrm = current_A * running_popvec_vrm;
              } else {
                arma::sp_mat current_A = as<arma::sp_mat>(current_mpm("A"));
                arma::sp_mat current_U_unaltered = as<arma::sp_mat>(current_mpm["U"]);
                arma::sp_mat current_F_unaltered = as<arma::sp_mat>(current_mpm["F"]);
                AdaptUtils::sp_UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
                
                if (err_check) errcheck_mpmout_times_vtr(m) = current_A;
                running_popvec_vrm = current_A * running_popvec_vrm;
              }
            }
          } else {
            // Rcout << "ifosr A35 (mat mult with density)        ";
            DataFrame used_density_input = density_df;
            DataFrame used_density_index_input = dens_index_df;
            
            NumericVector udi_alpha = as<NumericVector>(used_density_input["alpha"]);
            
            IntegerVector udii_index321 = as<IntegerVector>(used_density_index_input["index321"]);
            // Rcout << "udii_index321: " << udii_index321 << endl;
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            if (j >= (used_delay - 1)) {
              if (!stages_not_equal) {
                arma::vec delay_issue = N_mpm.col(j + 1 - used_delay);
                double delay_N_sum = arma::sum(delay_issue);
                
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                  running_popvec_vrm, sge_current, current_mpm, delay_N_sum,
                  0, integeronly, substoch, used_density_input,
                  used_density_index_input, false, sparse_bool, sparse_bool,
                  false, err_check);
                
                if (err_check) errcheck_mpmout_times_vtr(m) = new_projmat;
                running_popvec_vrm = new_popvec;
              } else {
                double delay_N_sum {0.0};
                
                if (j > 0) {
                  for (int l = 0; l < 2; l++) {
                    arma::mat delay_pop = as<arma::mat>(comm_out(m));
                    
                    arma::vec delay_pop_vec = delay_pop.col(j + 1 - used_delay);
                    arma::vec current_equiv_vec = as<arma::vec>(equivalence_list(0));
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
                
                if (err_check) errcheck_mpmout_times_vtr(m) = new_projmat;
                running_popvec_vrm = new_popvec;
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
              
              if (err_check) errcheck_mpmout_times_vtr(m) = new_projmat;
              running_popvec_vrm = new_popvec;
            }
          }
          
          // Rcout << "ifosr A36        ";
          if (integeronly) running_popvec_vrm = floor(running_popvec_vrm);
          double N_current = arma::sum(running_popvec_vrm);
          N_mpm(m, (j + 1)) = N_current;
          
          inda_num_terms_previous(m) = static_cast<int>(inda_num_terms_counter(m));
          indb_num_terms_previous(m) = static_cast<int>(indb_num_terms_counter(m));
          indc_num_terms_previous(m) = static_cast<int>(indc_num_terms_counter(m));
          inda_cat_terms_previous(m) = static_cast<int>(inda_cat_terms_counter(m));
          indb_cat_terms_previous(m) = static_cast<int>(indb_cat_terms_counter(m));
          indc_cat_terms_previous(m) = static_cast<int>(indc_cat_terms_counter(m));
          
          inda_num_terms_counter(m) = inda_num_terms_counter(m) + 1;
          indb_num_terms_counter(m) = indb_num_terms_counter(m) + 1;
          indc_num_terms_counter(m) = indc_num_terms_counter(m) + 1;
          inda_cat_terms_counter(m) = inda_cat_terms_counter(m) + 1;
          indb_cat_terms_counter(m) = indb_cat_terms_counter(m) + 1;
          indc_cat_terms_counter(m) = indc_cat_terms_counter(m) + 1;
          
          running_popvecs(m) = running_popvec_vrm;
          pops_out.col(j + 1) = running_popvec_vrm;
          
        } // if (j > (entry_time_vec(i) - 1))
        
        // Rcout << "ifosr A37        ";
        comm_out(m) = pops_out;
      } // main variant section
      
    // Rcout << "ifosr C3              ";
    errcheck_mpm_times(j) = errcheck_mpm_times_vtr;
    errcheck_mpmout_times(j) = errcheck_mpmout_times_vtr;
    //year_counter++;
  } // j loop
  // Rcout << "ifosr C5              ";
  N_out = N_mpm;
  
  if (err_check) {
    errcheck_mpm = errcheck_mpm_times;
    errcheck_mpmout = errcheck_mpmout_times;
  }
  
  arma::rowvec N_out_row0 = N_out.row(0);
  arma::rowvec N_out_row1 = N_out.row(1);
  
  List N_out_list (1);
  N_out_list(0) = N_out;
  
  // Rcout << "ifosr D got to Lyapunov_creator" << endl;
  AdaptUtils::Lyapunov_creator (out_df, N_out_list, entry_time_vec, 1, var_to_run, 2,
    times, fitness_times, conv_threshold, 1, true, zap_min);
  
  // Rcout << "ifosr E" << endl;
  // This bit needs to be redone to allow multiple variants, as right now it really only allows 2
  NumericVector out_fitness = as<NumericVector>(out_df(6));
  NumericVector out_fitness_e995 = as<NumericVector>(out_df(7));
  
  // Rcout << "ifosr F" << endl;
  CharacterVector out_names = as<CharacterVector>(out_df.names());
  
  NumericVector final_fitness (var_to_run);
  final_fitness(0) = out_fitness(0);
  final_fitness(1) = out_fitness_e995(0);
  
  LogicalVector converged (2);
  
  List processed_fitness_out (38);
  
  processed_fitness_out(0) = as<IntegerVector>(wrap(variant_nta));
  processed_fitness_out(1) = stage3_nta;
  processed_fitness_out(2) = stage2_nta;
  processed_fitness_out(3) = stage1_nta;
  processed_fitness_out(4) = age3_nta;
  processed_fitness_out(5) = age2_nta;
  processed_fitness_out(6) = eststage3_nta;
  processed_fitness_out(7) = eststage2_nta;
  processed_fitness_out(8) = eststage1_nta;
  processed_fitness_out(9) = estage3_nta;
  processed_fitness_out(10) = estage2_nta;
  processed_fitness_out(11) = givenrate_nta;
  processed_fitness_out(12) = offset_nta;
  processed_fitness_out(13) = multiplier_nta;
  processed_fitness_out(14) = convtype_nta;
  processed_fitness_out(15) = convtype_t12_nta;
  processed_fitness_out(16) = wrap(surv_dev_nta);
  processed_fitness_out(17) = wrap(obs_dev_nta);
  processed_fitness_out(18) = wrap(size_dev_nta);
  processed_fitness_out(19) = wrap(sizeb_dev_nta);
  processed_fitness_out(20) = wrap(sizec_dev_nta);
  processed_fitness_out(21) = wrap(repst_dev_nta);
  processed_fitness_out(22) = wrap(fec_dev_nta);
  processed_fitness_out(23) = wrap(jsurv_dev_nta);
  processed_fitness_out(24) = wrap(jobs_dev_nta);
  processed_fitness_out(25) = wrap(jsize_dev_nta);
  processed_fitness_out(26) = wrap(jsizeb_dev_nta);
  processed_fitness_out(27) = wrap(jsizec_dev_nta);
  processed_fitness_out(28) = wrap(jrepst_dev_nta);
  processed_fitness_out(29) = wrap(jmatst_dev_nta);
  processed_fitness_out(30) = wrap(indcova_nta);
  processed_fitness_out(31) = wrap(indcovb_nta);
  processed_fitness_out(32) = wrap(indcovc_nta);
  processed_fitness_out(33) = year2_nta;
  processed_fitness_out(34) = mpm_altered_nta;
  processed_fitness_out(35) = vrm_altered_nta;
  processed_fitness_out(36) = final_fitness;
  processed_fitness_out(37) = converged;
  
  CharacterVector processed_out_names = {"variant", "stage3", "stage2", "stage1",
    "age3", "age2", "eststage3", "eststage2", "eststage1", "estage3", "estage2",
    "givenrate", "offset", "multiplier", "convtype", "convtype_t12", "surv_dev",
    "obs_dev", "size_dev", "sizeb_dev", "sizec_dev", "repst_dev", "fec_dev",
    "jsurv_dev", "jobs_dev", "jsize_dev", "jsizeb_dev", "jsizec_dev",
    "jrepst_dev", "jmatst_dev", "indcova", "indcovb", "indcovc", "year2",
    "mpm_altered", "vrm_altered", "fitness", "converged"};
  CharacterVector processed_out_class = {"data.frame"};
  processed_fitness_out.attr("class") = processed_out_class;
  processed_fitness_out.attr("names") = processed_out_names;
  
  StringVector row_names(static_cast<int>(final_fitness.length()));
  for (int i = 0; i < static_cast<int>(final_fitness.length()); i++) {
    row_names(i) = std::to_string(i+1);
  }
  processed_fitness_out.attr("row.names") = row_names;
  
  out_df = processed_fitness_out;
}

//' Find ESS Values of Traits in Function-Based MPM Invasibility Analyses
//' 
//' @name ESS_optimizer_fb
//' 
//' @param ESS_Lyapunov A data frame to hold the ESS trait values, built by
//' reference.
//' @param ESS_trait_axis A data frame giving the minimum and maximum values of
//' variable traits, and the constant values of all others, in the trait axis.
//' @param Lyapunov_optim Main data frame giving Lyapunov coefficients for all
//' trait combinations developed for the ESS optima table. Holds elasticity
//' fitness values.
//' @param optim_trait_axis Main trait axis data frame corresponding to
//' \code{Lyapunov_optim}.
//' @param ESS_var_traits An integer vector modifed by this function by
//' reference, indicating the actual traits that vary. The element order is:
//' 1, givenrate; 2, offset; 3, multiplier; 4, surv_dev; 5, obs_dev;
//' 6, size_dev; 7, sizeb_dev; 8, sizec_dev; 9, repst_dev; 10, fec_dev;
//' 11, jsurv_dev; 12, jobs_dev; 13, jsize_dev; 14, jsizeb_dev; 15, jsizec_dev;
//' 16, jrepst_dev; 17, jmatst_dev; 18, indcova; 19, indcovb; and 20, indcovc.
//' @param flipped_traits An integer vector with 17 elements, giving the
//' identities of variables with both positive and negative values.
//' @param surv_dev_nta The survival column in the reassessed trait axis.
//' @param obs_dev_nta The observation status column in the reassessed trait
//' axis.
//' @param size_dev_nta The primary size column in the reassessed trait axis.
//' @param sizeb_dev_nta The secondary size column in the reassessed trait axis.
//' @param sizec_dev_nta The tertiary size column in the reassessed trait axis.
//' @param repst_dev_nta The reproductive status column in the reassessed trait
//' axis.
//' @param fec_dev_nta The fecundity column in the reassessed trait axis.
//' @param jsurv_dev_nta The juvenile survival column in the reassessed trait
//' axis.
//' @param jobs_dev_nta The juvenile observation status column in the reassessed
//' trait axis.
//' @param jsize_dev_nta The juvenile primary size column in the reassessed
//' trait axis.
//' @param jsizeb_dev_nta The juvenile secondary size column in the reassessed
//' trait axis.
//' @param jsizec_dev_nta The juvenile tertiary size column in the reassessed
//' trait axis.
//' @param jrepst_dev_nta The juvenile reproductive status column in the
//' reassessed trait axis.
//' @param jmatst_dev_nta The juvenile maturity status column in the reassessed
//' trait axis.
//' @param indcova_nta The individual covariate a column in the reassessed
//' trait axis.
//' @param indcovb_nta The individual covariate b column in the reassessed
//' trait axis.
//' @param indcovc_nta The individual covariate c column in the reassessed
//' trait axis.
//' @param variant_nta The variant column in the reassessed 995 trait axis.
//' @param new_stageexpansion_list A list with stage expansions for all
//' variant data used in ESS evaluation. This list includes an extra layer of
//' list elements, corresponding to the optim_ta and optim_ta_995 data.
//' @param used_times A list of year numbers for each time per run.
//' @param errcheck_mpm An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param errcheck_mpmout An optional list of all mpm_out matrices from MPM
//' processing. Only output if \code{err_check = "extreme"}.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param vrm_list A list of \code{vrm_input} objects.
//' @param allstages_all The allstages indexing data frame used to produce MPMs.
//' @param dev_terms_list List of deviations for vital rate models.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates.
//' @param stageexpansion_ta_devterms_by_variant A list giving trait axis info
//' by variant, with each variant given a list element.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param density_vr_list Data frame of \code{lefkoDensVR} objects holding
//' density relationships for all 14 vital rate models.
//' @param current_stageframe The main stageframe, including extra stages.
//' @param current_supplement A supplement in \code{lefkoSD} format.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param stageframe_df A stageframe object covering the MPM.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms.
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
//' @param fecmod_vec A numeric vector giving the fecmod values.
//' @param year_vec A vector giving the main years used.
//' @param patch_vec A vector giving the name of each patch used in projection.
//' @param times An integer giving the number of occasions to project.
//' @param fitness_times An integer giving the number of occasions at the end of
//' each run to use to estimate Lyapunov coefficients.
//' @param format_int An integer giving the MPM format.
//' @param stagecounts An integer giving the number of stages in the life
//' history.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param dev_terms_times_int A vector giving the number of occasions over
//' which vital rate y-intercept deviations cycle.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' @param conv_threshold The convergence threshold for Lyapunov coefficients
//' used in ESS optimization.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis. A relic value currently equal to the
//' number of variants in the entered trait axis.
//' @param opt_res_orig The original value of \code{opt_res}, prior to the
//' determination of the number of variable traits. Equal to \code{opt_res} if
//' the number of variable traits is 1, and to the square root of \code{opt_res}
//' if the number of variable traits is 2.
//' @param ehrlen An integer stating if historical MPMs should be in Ehrlen
//' format.
//' @param style An integer giving the style (e.g. ahistorical vs. historical)
//' of MPM.
//' @param loop_max The maximum number of times to attempt optimization, if
//' convergence does not occur.
//' @param filter An integer giving cleanup options for MPMs.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param elast_mult A multiplier for traits to assess the elasticity of
//' fitness in trait optimization. Defaults to 0.995.
//' @param zap_min A Boolean value describing whether to round fitness values
//' below the value given in \code{threshold}.
//' 
//' @return A final data frame giving the zeros associated with variable traits.
//' 
//' @keywords internal
//' @noRd
inline void ESS_optimizer_fb (DataFrame& ESS_Lyapunov, DataFrame& ESS_trait_axis,
  DataFrame& Lyapunov_optim, DataFrame& optim_trait_axis,
  IntegerVector& ESS_var_traits, IntegerVector& flipped_traits,
  arma::vec& surv_dev_nta, arma::vec& obs_dev_nta, arma::vec& size_dev_nta,
  arma::vec& sizeb_dev_nta, arma::vec& sizec_dev_nta, arma::vec& repst_dev_nta,
  arma::vec& fec_dev_nta, arma::vec& jsurv_dev_nta, arma::vec& jobs_dev_nta,
  arma::vec& jsize_dev_nta, arma::vec& jsizeb_dev_nta,
  arma::vec& jsizec_dev_nta, arma::vec& jrepst_dev_nta,
  arma::vec& jmatst_dev_nta, arma::vec& indcova_nta, arma::vec& indcovb_nta,
  arma::vec& indcovc_nta, arma::ivec& variant_nta,
  List& new_stageexpansion_list, List& used_times, List& errcheck_mpm,
  List& errcheck_mpmout, const List allmodels_all, const List vrm_list,
  const List allstages_all, const List dev_terms_list,
  const List ind_terms_num_list, const List ind_terms_cat_list,
  List& stageexpansion_ta_devterms_by_variant, const List sp_density_list,
  const List start_list, const List equivalence_list,
  const DataFrame density_vr_list, DataFrame& current_stageframe,
  const DataFrame current_supplement, const DataFrame density_df,
  const DataFrame dens_index_df, const DataFrame stageframe_df,
  const IntegerVector entry_time_vec, const IntegerVector sp_density_num_vec,
  const IntegerVector inda_terms_num_vec, const IntegerVector indb_terms_num_vec,
  const IntegerVector indc_terms_num_vec, const IntegerVector inda_terms_cat_vec,
  const IntegerVector indb_terms_cat_vec, const IntegerVector indc_terms_cat_vec,
  const IntegerVector dens_vr_yn_vec, const NumericVector fecmod_vec,
  const CharacterVector year_vec, const CharacterVector patch_vec,
  const int times, const int fitness_times, const int format_int,
  const int stagecounts, const int firstage_int, const int finalage_int,
  const int dev_terms_times_int, const int substoch,  const double exp_tol,
  const double theta_tol, const bool sparse_bool, const bool A_only,
  const bool stages_not_equal, const bool integeronly, const bool dens_yn_bool,
  const double conv_threshold, const int opt_res, const int opt_res_orig,
  const int ehrlen, const int style, const int loop_max, const int filter,
  const bool errcheck, double elast_mult, const bool zap_min) {
  
  // Rcout << "Entered ESS_optimizer_fb" << endl;
  // Rcout << "flipped_traits: " << flipped_traits << endl;
  
  arma::uvec potential_optima (opt_res, fill::zeros);
  arma::vec inv_fitness = as<NumericVector>(Lyapunov_optim["fitness_variant2_e995"]);
  IntegerVector ESS_var_traits_corresponding_indices = {11, 12, 13, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
  
  //int Lyapunov_optim_rows = Lyapunov_optim.nrows();
  int main_loop_breakpoint = static_cast<int>(inv_fitness.n_elem);
  if (opt_res != opt_res_orig) main_loop_breakpoint = opt_res_orig;
  
  // Rcout << "ESS_optimizer_fb A" << endl;
  for (int j = 0; j < (main_loop_breakpoint - 1); j++) {
    double base_inv_fit = inv_fitness(j);
    double next_inv_fit = inv_fitness(j + 1);
    
    if (base_inv_fit < 0. && next_inv_fit >= 0.) potential_optima(j) = 1;
    if (base_inv_fit >= 0. && next_inv_fit < 0.) potential_optima(j) = 1;
  }
  
  unsigned int total_optima = arma::sum(potential_optima);
  arma::uvec potential_optima_indices = find(potential_optima);
  arma::uvec ESS_var_traits_arma = as<arma::uvec>(ESS_var_traits);
  arma::uvec ESS_vta_indices = find(ESS_var_traits_arma);
  int vars_to_alter = static_cast<int>(ESS_vta_indices.n_elem);
  
  if (vars_to_alter > 1) {
    Rf_warningcall(R_NilValue, "Trait optimization may not run properly with more than 1 tested trait.");
  }
  // Rcout << "ESS_optimizer_fb B" << endl;
  // Rcout << "total_optima (potential): " << total_optima << endl;
  
  List final_output (total_optima);
  for (int i = 0; i < total_optima; i++) {
    // Rcout << "ESS_optimizer_fb C entered optimization trait loop" << endl;
    arma::mat N_out;
    List comm_out (2);
    bool found_converged {false};
    
    for (int m = 0; m < 2; m++) {
      arma::mat pops_out (stagecounts, (times + 1), fill::zeros);
      comm_out(m) = pops_out;
    }
    
    IntegerVector core_index_list = {(static_cast<int>(potential_optima_indices(i))), 
      (static_cast<int>(potential_optima_indices(i) + 1))};
    DataFrame core_trait_axis_instance = AdaptUtils::df_indices(optim_trait_axis, core_index_list);
    
    DataFrame reference_variants = clone(core_trait_axis_instance);
    NumericVector ref_fitness = {static_cast<double>(inv_fitness(potential_optima_indices(i))), 
      static_cast<double>(inv_fitness(potential_optima_indices(i) + 1))};
    // Rcout << "ref_fitness (lower bound invasion fitness, higher bound invasion fitness): " << ref_fitness << endl;
    
    reference_variants["fitness"] = ref_fitness;
    //DataFrame old_reference_variants = clone(reference_variants);
    bool opt_needed {false};
    
    if (abs(static_cast<double>(inv_fitness(potential_optima_indices(i)))) < conv_threshold) {
      // Rcout << "abs of lower bound invader fitness < threshold, so no further optimization needed" << endl;
      IntegerVector chosen_one = {0};
      DataFrame optim_point = AdaptUtils::df_indices(reference_variants, chosen_one);
      
      found_converged = true;
      LogicalVector found_converged_vec = {found_converged};
      optim_point["converged"] = found_converged_vec;
      
      final_output(i) = optim_point;
    } else if (abs(static_cast<double>(inv_fitness(potential_optima_indices(i) + 1))) < conv_threshold) {
      // Rcout << "abs of higher bound invader fitness < threshold, so no further optimization needed" << endl;
      IntegerVector chosen_one = {1};
      DataFrame optim_point = AdaptUtils::df_indices(reference_variants, chosen_one);
      
      found_converged = true;
      LogicalVector found_converged_vec = {found_converged};
      optim_point["converged"] = found_converged_vec;
      
      final_output(i) = optim_point;
    } else {
      // Rcout << "lower and upper bound invader fitness above threshold, so further optimization needed" << endl;
      opt_needed = true;
    }
    
    int loop_tracker {0};
    DataFrame variants_to_test = clone(reference_variants); // These values will be overwritten in the cloned data frame
    
    DataFrame variant_L_bound; // These dataframes are single row
    DataFrame variant_H_bound; // They mark the boundaries of the search
    
    if (inv_fitness(potential_optima_indices(i)) < inv_fitness(potential_optima_indices(i) + 1)) {
      IntegerVector chosen_low = {0};
      IntegerVector chosen_high = {1};
      DataFrame variant_L_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_low);
      DataFrame variant_H_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_high);
      
      variant_L_bound = variant_L_bound_pre;
      variant_H_bound = variant_H_bound_pre;
    } else {
      IntegerVector chosen_low = {1};
      IntegerVector chosen_high = {0};
      DataFrame variant_L_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_low);
      DataFrame variant_H_bound_pre = AdaptUtils::df_indices(variants_to_test, chosen_high);
      
      variant_L_bound = variant_L_bound_pre;
      variant_H_bound = variant_H_bound_pre;
    }
    
    // Rcout << "ESS_optimizer_fb D" << endl;
    DataFrame variant_L; // Test point with lower abs fitness
    DataFrame variant_H; // Test point with higher abs fitness
    DataFrame variant_C; // Centroid
    
    // Initial point set-up and fitness check
    DataFrame point_1; // 1/3 distance from low bound to high bound
    DataFrame point_2; // 2/3 distance from low bound to high bound
    
    {
      { // Point 1 initialization
        // Rcout << "ESS_optimizer_fb E point_1 initialization " << endl;
        arma::uvec variant_nta_new = {1, 2};
        arma::vec surv_dev_nta_new;
        arma::vec obs_dev_nta_new;
        arma::vec size_dev_nta_new;
        arma::vec sizeb_dev_nta_new;
        arma::vec sizec_dev_nta_new;
        arma::vec repst_dev_nta_new;
        arma::vec fec_dev_nta_new;
        arma::vec jsurv_dev_nta_new;
        arma::vec jobs_dev_nta_new;
        arma::vec jsize_dev_nta_new;
        arma::vec jsizeb_dev_nta_new;
        arma::vec jsizec_dev_nta_new;
        arma::vec jrepst_dev_nta_new;
        arma::vec jmatst_dev_nta_new;
        arma::vec indcova_nta_new;
        arma::vec indcovb_nta_new;
        arma::vec indcovc_nta_new;
        
        NumericVector variants_vars_to_test;
        
        for (int j = 0; j < vars_to_alter; j++) {
          double base_mean {0.};
          
          // Rcout << "ESS_optimizer_fb E1 j: " << j << endl;
          NumericVector core_var = as<NumericVector>(reference_variants(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          int flipped_trait_status = flipped_traits(ESS_vta_indices(j));
          // Rcout << "ESS_optimizer_fb E2" << endl;
          double low_value = static_cast<double>(core_var(0));
          double high_value = static_cast<double>(core_var(1));
          
          // Rcout << "ESS_optimizer_fb E3" << endl;
          NumericVector base_values = {high_value, low_value};
          double current_mult {1./3.};
          base_mean = low_value + current_mult * (high_value - low_value);
          
          variants_vars_to_test = as<NumericVector>(variants_to_test(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          variants_vars_to_test(0) = base_mean;
          
          if (base_mean < 0. && flipped_trait_status == 1) {
            variants_vars_to_test(1) = base_mean + (base_mean - (base_mean * elast_mult));
          } else {
            variants_vars_to_test(1) = base_mean * elast_mult;
          }
        }
        
        // Rcout << "ESS_optimizer_fb E4" << endl;
        surv_dev_nta_new = as<arma::vec>(variants_to_test(16));
        obs_dev_nta_new = as<arma::vec>(variants_to_test(17));
        size_dev_nta_new = as<arma::vec>(variants_to_test(18));
        sizeb_dev_nta_new = as<arma::vec>(variants_to_test(19));
        sizec_dev_nta_new = as<arma::vec>(variants_to_test(20));
        repst_dev_nta_new = as<arma::vec>(variants_to_test(21));
        fec_dev_nta_new = as<arma::vec>(variants_to_test(22));
        jsurv_dev_nta_new = as<arma::vec>(variants_to_test(23));
        jobs_dev_nta_new = as<arma::vec>(variants_to_test(24));
        jsize_dev_nta_new = as<arma::vec>(variants_to_test(25));
        jsizeb_dev_nta_new = as<arma::vec>(variants_to_test(26));
        jsizec_dev_nta_new = as<arma::vec>(variants_to_test(27));
        jrepst_dev_nta_new = as<arma::vec>(variants_to_test(28));
        jmatst_dev_nta_new = as<arma::vec>(variants_to_test(29));
        indcova_nta_new = as<arma::vec>(variants_to_test(30));
        indcovb_nta_new = as<arma::vec>(variants_to_test(31));
        indcovc_nta_new = as<arma::vec>(variants_to_test(32));
        
        // Rcout << "ESS_optimizer_fb E5" << endl;
        
        // New stageexpansion
        IntegerVector chosen_int = {0};
        DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
        // Rcout << "ESS_optimizer_fb E6" << endl;
        DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
        
        DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
          variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
          filter);
        // Rcout << "ESS_optimizer_fb E7" << endl;
        DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
          variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
          filter);
        
        // Rcout << "ESS_optimizer_fb E8" << endl;
        
        chosen_int = 1;
        StringVector focused_var = {"mpm_altered"};
        DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
          as<RObject>(chosen_int), false, true, false, false, true,
          as<RObject>(focused_var));
        // Rcout << "ESS_optimizer_fb E9" << endl;
        
        DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
          as<RObject>(chosen_int), false, true, false, false, true,
          as<RObject>(focused_var));
        // Rcout << "ESS_optimizer_fb E10" << endl;
        
        List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
          _["e995"] = stageexpansion_995_reduced);
        // Rcout << "ESS_optimizer_fb E11" << endl;
        
        invfb_optim_singlerun (point_1, surv_dev_nta_new, obs_dev_nta_new,
          size_dev_nta_new, sizeb_dev_nta_new, sizec_dev_nta_new, repst_dev_nta_new, 
          fec_dev_nta_new, jsurv_dev_nta_new, jobs_dev_nta_new, jsize_dev_nta_new, 
          jsizeb_dev_nta_new, jsizec_dev_nta_new, jrepst_dev_nta_new,
          jmatst_dev_nta_new, indcova_nta_new, indcovb_nta_new, indcovc_nta_new,
          variant_nta_new, variants_to_test, N_out, comm_out, errcheck_mpm,
          errcheck_mpmout, sge_to_test, used_times, allmodels_all, vrm_list,
          allstages_all, dev_terms_list, ind_terms_num_list, ind_terms_cat_list,
          stageexpansion_ta_devterms_by_variant, sp_density_list,
          start_list, equivalence_list, density_vr_list, current_stageframe,
          current_supplement, density_df, dens_index_df, entry_time_vec,
          sp_density_num_vec, inda_terms_num_vec, indb_terms_num_vec,
          indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
          indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec, year_vec, patch_vec,
          times, fitness_times, format_int, firstage_int, finalage_int,
          dev_terms_times_int, substoch, opt_res, opt_res_orig, exp_tol,
          theta_tol, conv_threshold, sparse_bool, A_only, stages_not_equal,
          integeronly, dens_yn_bool, errcheck, zap_min);
        
        // Rcout << "ESS_optimizer_fb E12" << endl;
      }
        
      { // Point 2 initialization
        // Rcout << "ESS_optimizer_fb F point_2 initialization" << endl;
        arma::uvec variant_nta_new = {1, 2};
        arma::vec surv_dev_nta_new;
        arma::vec obs_dev_nta_new;
        arma::vec size_dev_nta_new;
        arma::vec sizeb_dev_nta_new;
        arma::vec sizec_dev_nta_new;
        arma::vec repst_dev_nta_new;
        arma::vec fec_dev_nta_new;
        arma::vec jsurv_dev_nta_new;
        arma::vec jobs_dev_nta_new;
        arma::vec jsize_dev_nta_new;
        arma::vec jsizeb_dev_nta_new;
        arma::vec jsizec_dev_nta_new;
        arma::vec jrepst_dev_nta_new;
        arma::vec jmatst_dev_nta_new;
        arma::vec indcova_nta_new;
        arma::vec indcovb_nta_new;
        arma::vec indcovc_nta_new;
        
        NumericVector variants_vars_to_test;
        
        for (int j = 0; j < vars_to_alter; j++) {
          double base_mean {0.};
          
          NumericVector core_var = as<NumericVector>(reference_variants(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          int flipped_trait_status = flipped_traits(ESS_vta_indices(j));
          // Rcout << "ESS_optimizer_fb F1" << endl;
          double low_value = static_cast<double>(core_var(0));
          double high_value = static_cast<double>(core_var(1));
          
          // Rcout << "ESS_optimizer_fb F2" << endl;
          NumericVector base_values = {high_value, low_value};
          double current_mult {2./3.};
          base_mean = low_value + current_mult * (high_value - low_value);
          
          variants_vars_to_test = as<NumericVector>(variants_to_test(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          variants_vars_to_test(0) = base_mean;
          
          if (base_mean < 0. && flipped_trait_status == 1) {
            variants_vars_to_test(1) = base_mean + (base_mean - (base_mean * elast_mult));
          } else {
            variants_vars_to_test(1) = base_mean * elast_mult;
          }
        }
        
        // Rcout << "ESS_optimizer_fb F3" << endl;
        surv_dev_nta_new = as<arma::vec>(variants_to_test(16));
        obs_dev_nta_new = as<arma::vec>(variants_to_test(17));
        size_dev_nta_new = as<arma::vec>(variants_to_test(18));
        sizeb_dev_nta_new = as<arma::vec>(variants_to_test(19));
        sizec_dev_nta_new = as<arma::vec>(variants_to_test(20));
        repst_dev_nta_new = as<arma::vec>(variants_to_test(21));
        fec_dev_nta_new = as<arma::vec>(variants_to_test(22));
        jsurv_dev_nta_new = as<arma::vec>(variants_to_test(23));
        jobs_dev_nta_new = as<arma::vec>(variants_to_test(24));
        jsize_dev_nta_new = as<arma::vec>(variants_to_test(25));
        jsizeb_dev_nta_new = as<arma::vec>(variants_to_test(26));
        jsizec_dev_nta_new = as<arma::vec>(variants_to_test(27));
        jrepst_dev_nta_new = as<arma::vec>(variants_to_test(28));
        jmatst_dev_nta_new = as<arma::vec>(variants_to_test(29));
        indcova_nta_new = as<arma::vec>(variants_to_test(30));
        indcovb_nta_new = as<arma::vec>(variants_to_test(31));
        indcovc_nta_new = as<arma::vec>(variants_to_test(32));
        
        // Rcout << "ESS_optimizer_fb F4" << endl;
        
        // New stageexpansion
        IntegerVector chosen_int = {0};
        DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
        // Rcout << "ESS_optimizer_fb F5" << endl;
        DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
        
        DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
          variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
          filter);
        // Rcout << "ESS_optimizer_fb F6" << endl;
        DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
          variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
          filter);
        
        // Rcout << "ESS_optimizer_fb F7" << endl;
        
        chosen_int = 1;
        StringVector focused_var = {"mpm_altered"};
        DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
          as<RObject>(chosen_int), false, true, false, false, true,
          as<RObject>(focused_var));
        // Rcout << "ESS_optimizer_fb F8" << endl;
        
        DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
          as<RObject>(chosen_int), false, true, false, false, true,
          as<RObject>(focused_var));
        // Rcout << "ESS_optimizer_fb F9" << endl;
        
        List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
          _["e995"] = stageexpansion_995_reduced);
        // Rcout << "ESS_optimizer_fb F10" << endl;
        
        invfb_optim_singlerun (point_2, surv_dev_nta_new, obs_dev_nta_new,
          size_dev_nta_new, sizeb_dev_nta_new, sizec_dev_nta_new, repst_dev_nta_new, 
          fec_dev_nta_new, jsurv_dev_nta_new, jobs_dev_nta_new, jsize_dev_nta_new, 
          jsizeb_dev_nta_new, jsizec_dev_nta_new, jrepst_dev_nta_new,
          jmatst_dev_nta_new, indcova_nta_new, indcovb_nta_new, indcovc_nta_new,
          variant_nta_new, variants_to_test, N_out, comm_out, errcheck_mpm,
          errcheck_mpmout, sge_to_test, used_times, allmodels_all, vrm_list,
          allstages_all, dev_terms_list, ind_terms_num_list, ind_terms_cat_list,
          stageexpansion_ta_devterms_by_variant, sp_density_list,
          start_list, equivalence_list, density_vr_list, current_stageframe,
          current_supplement, density_df, dens_index_df, entry_time_vec,
          sp_density_num_vec, inda_terms_num_vec, indb_terms_num_vec,
          indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
          indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec, year_vec, patch_vec,
          times, fitness_times, format_int, firstage_int, finalage_int,
          dev_terms_times_int, substoch, opt_res, opt_res_orig, exp_tol,
          theta_tol, conv_threshold, sparse_bool, A_only, stages_not_equal,
          integeronly, dens_yn_bool, errcheck, zap_min);
        
        // Rcout << "ESS_optimizer_fb F11" << endl;
      }
    }
    
    // Rcout << "ESS_optimizer_fb G" << endl;
    int search_mode {0};
    int worst_point {2};
    double x1 {0.};
    double x2 {0.};
    double conv_test {0.};
    double xR {0.}; // test value to keep in memory from previous reflection round
    double xE {0.}; // test value to keep in memory from previous expansion round
    double xC_new {0.}; // test value to keep in memory from previous contraction round
    
    double xR_fitness {0.}; // fitness associated with test value to keep in memory from previous reflection round
    double xE_fitness {0.}; // fitness associated with test value to keep in memory from previous reflection round
    double xC_new_fitness {0.}; // fitness associated with test value to keep in memory from previous reflection round
    
    double e995_fitness_1 {0.};
    double e995_fitness_2 {0.};
    double e995_fitness_L {0.};
    double e995_fitness_H {0.};
    
    DataFrame variant_R;
    DataFrame variant_E;
    DataFrame variant_C_new;
    
    DataFrame ESS_test_values;
    DataFrame ESS_out_values;
    
    while (opt_needed && loop_tracker < loop_max) {
      // Rcout << "\n\n\nESS_optimizer_fb H entered optimization main loop, with loop_tracker = " << loop_tracker << endl;
      
      arma::uvec variant_nta_new = {1, 2};
      arma::vec surv_dev_nta_new;
      arma::vec obs_dev_nta_new;
      arma::vec size_dev_nta_new;
      arma::vec sizeb_dev_nta_new;
      arma::vec sizec_dev_nta_new;
      arma::vec repst_dev_nta_new;
      arma::vec fec_dev_nta_new;
      arma::vec jsurv_dev_nta_new;
      arma::vec jobs_dev_nta_new;
      arma::vec jsize_dev_nta_new;
      arma::vec jsizeb_dev_nta_new;
      arma::vec jsizec_dev_nta_new;
      arma::vec jrepst_dev_nta_new;
      arma::vec jmatst_dev_nta_new;
      arma::vec indcova_nta_new;
      arma::vec indcovb_nta_new;
      arma::vec indcovc_nta_new;
      
      NumericVector variants_vars_to_test;
      
      for (int j = 0; j < vars_to_alter; j++) {
        // Rcout << "ESS_optimizer_fb I j: " << j << endl;
        
        NumericVector point1_fitness_values = as<NumericVector>(point_1["fitness"]);
        NumericVector point2_fitness_values = as<NumericVector>(point_2["fitness"]);
        e995_fitness_1 = point1_fitness_values(1); 
        e995_fitness_2 = point2_fitness_values(1); 
        
        // Rcout << "ESS_optimizer_fb J" << endl;
        if (abs(e995_fitness_1) < abs(e995_fitness_2)) {
          variant_L = clone(point_1);
          variant_H = clone(point_2);
          e995_fitness_L = e995_fitness_1;
          e995_fitness_H = e995_fitness_2;
          worst_point = 2;
        } else {
          variant_H = clone(point_1);
          variant_L = clone(point_2);
          e995_fitness_L = e995_fitness_2;
          e995_fitness_H = e995_fitness_1;
          worst_point = 1;
        }
        
        // Rcout << "ESS_optimizer_fb J1" << endl;
        variant_C = clone(variant_L);
        double base_test_value = {0.};
        
        // Rcout << "ESS_optimizer_fb J2" << endl;
        NumericVector core_var_C = as<NumericVector>(variant_C(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_fb J3" << endl;
        NumericVector core_var_H = as<NumericVector>(variant_H(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_fb J4" << endl;
        NumericVector core_var_L_bound = as<NumericVector>(variant_L_bound(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_fb J5" << endl;
        NumericVector core_var_H_bound = as<NumericVector>(variant_H_bound(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        // Rcout << "ESS_optimizer_fb J6" << endl;
        int flipped_trait_status = flipped_traits(ESS_vta_indices(j));
        // Rcout << "ESS_optimizer_fb J7" << endl;
        NumericVector point_1_trait = as<NumericVector>(point_1(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j)))); // Not necessary - used for testing
        // Rcout << "ESS_optimizer_fb J8" << endl;
        NumericVector point_2_trait = as<NumericVector>(point_2(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j)))); // Not necessary - used for testing
        // Rcout << "ESS_optimizer_fb J9" << endl;
        
        double xC = static_cast<double>(core_var_C(0));
        double xL = xC;
        double xH = static_cast<double>(core_var_H(0)); // Used to be index 1
        double low_value_bound = static_cast<double>(core_var_L_bound(0));
        double high_value_bound = static_cast<double>(core_var_H_bound(0)); // Used to be index 1
        
        x1 = xL;
        x2 = xH;
        conv_test = abs(x1 - x2);
        
        // Rcout << "\nESS_optimizer_fb K search mode choice" << endl;
        // Rcout << "original low_value_bound: " << low_value_bound << endl;
        // Rcout << "original high_value_bound: " << high_value_bound << endl;
        // Rcout << "original point_1 trait value: " << point_1_trait << endl;
        // Rcout << "original point_2 trait value: " << point_2_trait << endl;
        // Rcout << "worst_point: " << worst_point << endl;
        // Rcout << "original xC: " << xC << endl;
        // Rcout << "original xL: " << xL << endl;
        // Rcout << "original xH: " << xH << endl;
        // Rcout << "original xR: " << xR << endl;
        // Rcout << "original xE: " << xE << endl;
        // Rcout << "original xC_new: " << xC_new << endl;
        // Rcout << "original point_1 invasion fitness: " << e995_fitness_1 << endl;
        // Rcout << "original point_2 invasion fitness: " << e995_fitness_2 << endl;
        // Rcout << "original e995_fitness_L: " << e995_fitness_L << endl;
        // Rcout << "original e995_fitness_H: " << e995_fitness_H << endl;
        // Rcout << "original xR_fitness: " << xR_fitness << endl;
        // Rcout << "original xE_fitness: " << xE_fitness << endl;
        // Rcout << "original xC_new_fitness: " << xC_new_fitness << endl;
        
        //NumericVector base_values = {high_value, low_value};
        if (search_mode == 0) { // Reflection
          // Rcout << "\nReflection mode (search mode = 0)" << endl;
          
          double alpha = 0.8;
          base_test_value = xC + alpha * (xC - xH);
          
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          xR = base_test_value;
          
          // Rcout << "value of xC (used in base_test_value calculation): " << xC << endl;
          // Rcout << "value of xH (used in base_test_value calculation): " << xH << endl;
          // Rcout << "value of xR (calculated during reflection): " << xR << endl;
          
        } else if (search_mode == 1) { // Expansion
          // Rcout << "\nExpansion mode (search mode = 1)" << endl;
          
          double gamma = 2.0;
          base_test_value = xC + gamma * (xR - xC);
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          xE = base_test_value;
          
          // Rcout << "value of xC (used in base_test_value calculation): " << xC << endl;
          // Rcout << "value of xR (used in base_test_value calculation): " << xR << endl;
          // Rcout << "value of xE (calculated during expansion): " << xE << endl;
          
        } else if (search_mode == 2) { // Contraction
          // Rcout << "\nContraction mode (search mode = 2)" << endl;
          
          double rho = 0.5;
          
          if (worst_point == 1) {
            if (abs(xR_fitness) < abs(e995_fitness_1)) {
              base_test_value = xC + rho * (xR - xC);
            } else if (abs(xR_fitness) >= abs(e995_fitness_1)) {
              base_test_value = xC + rho * (xH - xC);
            }
          } else if (worst_point == 2) {
            if (abs(xR_fitness) < abs(e995_fitness_2)) {
              base_test_value = xC + rho * (xR - xC);
            } else if (abs(xR_fitness) >= abs(e995_fitness_2)) {
              base_test_value = xC + rho * (xH - xC);
            }
          }
          
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          xC_new = base_test_value;
          
          // Rcout << "value of xC (used in base_test_value calculation): " << xC << endl;
          // Rcout << "value of xH (used in base_test_value calculation): " << xH << endl;
          // Rcout << "value of xC_new (calculated during contraction): " << xC_new << endl;
          
        } else if (search_mode == 3) { // Shrinking
          // Rcout << "\nShrinking mode (search mode = 3)" << endl;
          
          double sigma = 0.5;
          
          NumericVector core_var_p1 = as<NumericVector>(point_1(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          NumericVector core_var_p2 = as<NumericVector>(point_2(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          
          double p1 = static_cast<double>(core_var_p1(0));
          double p2 = static_cast<double>(core_var_p2(0));
          
          base_test_value = xL + sigma * (p1 - xL);
          double base_test_value2 = xL + sigma * (p2 - xL);
          NumericVector maximal_test_vec = {base_test_value, low_value_bound};
          double maximal_test = max(maximal_test_vec);
          NumericVector minimal_test_vec = {maximal_test, high_value_bound};
          double minimal_test = min(minimal_test_vec);
          
          base_test_value = minimal_test;
          
          NumericVector maximal_test2 = {base_test_value2, low_value_bound};
          NumericVector minimal_test2 = {maximal_test2(0), high_value_bound};
          
          base_test_value2 = minimal_test2(0);
          
          // Rcout << "value of xL (used in base_test_value calculation and base_test_value2): " << xL << endl;
          // Rcout << "value of p1 (used in base_test_value calculation and base_test_value2): " << p1 << endl;
          // Rcout << "value of p2 (used in base_test_value calculation and base_test_value2): " << p2 << endl;
          // Rcout << "value of base_test_value (calculated during shrinking): " << base_test_value << endl;
          // Rcout << "value of base_test_value2 (calculated during shrinking): " << base_test_value2 << endl;
          
          // Rcout << "ESS_optimizer_fb K1" << endl;
          
          variants_vars_to_test = as<NumericVector>(variants_to_test(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          variants_vars_to_test(0) = base_test_value;
          
          if (base_test_value < 0. && flipped_trait_status == 1) {
            variants_vars_to_test(1) = base_test_value + (base_test_value - (base_test_value * elast_mult));
          } else {
            variants_vars_to_test(1) = base_test_value * elast_mult;
          }
          
          // Rcout << "ESS_optimizer_fb K2" << endl;
          
          surv_dev_nta_new = as<arma::vec>(variants_to_test(16));
          obs_dev_nta_new = as<arma::vec>(variants_to_test(17));
          size_dev_nta_new = as<arma::vec>(variants_to_test(18));
          sizeb_dev_nta_new = as<arma::vec>(variants_to_test(19));
          sizec_dev_nta_new = as<arma::vec>(variants_to_test(20));
          repst_dev_nta_new = as<arma::vec>(variants_to_test(21));
          fec_dev_nta_new = as<arma::vec>(variants_to_test(22));
          jsurv_dev_nta_new = as<arma::vec>(variants_to_test(23));
          jobs_dev_nta_new = as<arma::vec>(variants_to_test(24));
          jsize_dev_nta_new = as<arma::vec>(variants_to_test(25));
          jsizeb_dev_nta_new = as<arma::vec>(variants_to_test(26));
          jsizec_dev_nta_new = as<arma::vec>(variants_to_test(27));
          jrepst_dev_nta_new = as<arma::vec>(variants_to_test(28));
          jmatst_dev_nta_new = as<arma::vec>(variants_to_test(29));
          indcova_nta_new = as<arma::vec>(variants_to_test(30));
          indcovb_nta_new = as<arma::vec>(variants_to_test(31));
          indcovc_nta_new = as<arma::vec>(variants_to_test(32));
          
          // Rcout << "ESS_optimizer_fb K3" << endl;
          
          // New stageexpansion
          IntegerVector chosen_int = {0};
          DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
          // Rcout << "ESS_optimizer_fb K4" << endl;
          DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
          
          DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
            filter);
          // Rcout << "ESS_optimizer_fb K5" << endl;
          DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
            filter);
          
          // Rcout << "ESS_optimizer_fb K6" << endl;
          
          chosen_int = 1;
          StringVector focused_var = {"mpm_altered"};
          DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_fb K7" << endl;
          
          DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_fb K8" << endl;
          
          List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
            _["e995"] = stageexpansion_995_reduced);
          // Rcout << "ESS_optimizer_fb K9" << endl;
          
          invfb_optim_singlerun (point_1, surv_dev_nta_new, obs_dev_nta_new,
            size_dev_nta_new, sizeb_dev_nta_new, sizec_dev_nta_new, repst_dev_nta_new, 
            fec_dev_nta_new, jsurv_dev_nta_new, jobs_dev_nta_new, jsize_dev_nta_new, 
            jsizeb_dev_nta_new, jsizec_dev_nta_new, jrepst_dev_nta_new,
            jmatst_dev_nta_new, indcova_nta_new, indcovb_nta_new,
            indcovc_nta_new, variant_nta_new, variants_to_test, N_out, comm_out,
            errcheck_mpm, errcheck_mpmout, sge_to_test, used_times, allmodels_all,
            vrm_list, allstages_all, dev_terms_list, ind_terms_num_list,
            ind_terms_cat_list, stageexpansion_ta_devterms_by_variant, sp_density_list,
            start_list, equivalence_list, density_vr_list, current_stageframe,
            current_supplement, density_df, dens_index_df, entry_time_vec,
            sp_density_num_vec, inda_terms_num_vec, indb_terms_num_vec,
            indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
            indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec, year_vec, patch_vec,
            times, fitness_times, format_int, firstage_int, finalage_int,
            dev_terms_times_int, substoch, opt_res, opt_res_orig, exp_tol,
            theta_tol, conv_threshold, sparse_bool, A_only, stages_not_equal,
            integeronly, dens_yn_bool, errcheck, zap_min);
          
          // Rcout << "ESS_optimizer_fb K10" << endl;
          variants_vars_to_test = as<NumericVector>(variants_to_test(
              ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
          variants_vars_to_test(0) = base_test_value2;
          
          if (base_test_value2 < 0. && flipped_trait_status == 1) {
            variants_vars_to_test(1) = base_test_value2 + (base_test_value2 - (base_test_value2 * elast_mult));
          } else {
            variants_vars_to_test(1) = base_test_value2 * elast_mult;
          }
          
          // Rcout << "ESS_optimizer_fb K11" << endl;
          surv_dev_nta_new = as<arma::vec>(variants_to_test(16));
          obs_dev_nta_new = as<arma::vec>(variants_to_test(17));
          size_dev_nta_new = as<arma::vec>(variants_to_test(18));
          sizeb_dev_nta_new = as<arma::vec>(variants_to_test(19));
          sizec_dev_nta_new = as<arma::vec>(variants_to_test(20));
          repst_dev_nta_new = as<arma::vec>(variants_to_test(21));
          fec_dev_nta_new = as<arma::vec>(variants_to_test(22));
          jsurv_dev_nta_new = as<arma::vec>(variants_to_test(23));
          jobs_dev_nta_new = as<arma::vec>(variants_to_test(24));
          jsize_dev_nta_new = as<arma::vec>(variants_to_test(25));
          jsizeb_dev_nta_new = as<arma::vec>(variants_to_test(26));
          jsizec_dev_nta_new = as<arma::vec>(variants_to_test(27));
          jrepst_dev_nta_new = as<arma::vec>(variants_to_test(28));
          jmatst_dev_nta_new = as<arma::vec>(variants_to_test(29));
          indcova_nta_new = as<arma::vec>(variants_to_test(30));
          indcovb_nta_new = as<arma::vec>(variants_to_test(31));
          indcovc_nta_new = as<arma::vec>(variants_to_test(32));
          
          // Rcout << "ESS_optimizer_fb K12" << endl;
          // New stageexpansion
          chosen_int = {0};
          variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
          // Rcout << "ESS_optimizer_fb K13" << endl;
          variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
          
          stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
            filter);
          // Rcout << "ESS_optimizer_fb K14" << endl;
          stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
            filter);
          
          // Rcout << "ESS_optimizer_fb K15" << endl;
          chosen_int = 1;
          focused_var = {"mpm_altered"};
          stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_fb K16" << endl;
          
          stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_fb K17" << endl;
          
          sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
            _["e995"] = stageexpansion_995_reduced);
          // Rcout << "ESS_optimizer_fb K18" << endl;
          
          invfb_optim_singlerun (point_2, surv_dev_nta_new, obs_dev_nta_new,
            size_dev_nta_new, sizeb_dev_nta_new, sizec_dev_nta_new, repst_dev_nta_new, 
            fec_dev_nta_new, jsurv_dev_nta_new, jobs_dev_nta_new, jsize_dev_nta_new, 
            jsizeb_dev_nta_new, jsizec_dev_nta_new, jrepst_dev_nta_new,
            jmatst_dev_nta_new, indcova_nta_new, indcovb_nta_new,
            indcovc_nta_new,variant_nta_new, variants_to_test, N_out, comm_out,
            errcheck_mpm, errcheck_mpmout, sge_to_test, used_times, allmodels_all,
            vrm_list, allstages_all, dev_terms_list, ind_terms_num_list,
            ind_terms_cat_list, stageexpansion_ta_devterms_by_variant, sp_density_list,
            start_list, equivalence_list, density_vr_list, current_stageframe,
            current_supplement, density_df, dens_index_df, entry_time_vec,
            sp_density_num_vec, inda_terms_num_vec, indb_terms_num_vec,
            indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
            indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec, year_vec, patch_vec,
            times, fitness_times, format_int, firstage_int, finalage_int,
            dev_terms_times_int, substoch, opt_res, opt_res_orig, exp_tol,
            theta_tol, conv_threshold, sparse_bool, A_only, stages_not_equal,
            integeronly, dens_yn_bool, errcheck, zap_min);
          
          // Rcout << "ESS_optimizer_fb K19" << endl;
        }
        
        variants_vars_to_test = as<NumericVector>(variants_to_test(
            ESS_var_traits_corresponding_indices(ESS_vta_indices(j))));
        variants_vars_to_test(0) = base_test_value;
        
        if (base_test_value < 0. && flipped_trait_status == 1) {
          variants_vars_to_test(1) = base_test_value + (base_test_value - (base_test_value * elast_mult));
        } else {
          variants_vars_to_test(1) = base_test_value * elast_mult;
        }
        // Rcout << "ESS_optimizer_fb K20" << endl;
        
      //}
        
        if (search_mode < 3) {
          // Rcout << "ESS_optimizer_fb L" << endl;
          surv_dev_nta_new = as<arma::vec>(variants_to_test(16));
          obs_dev_nta_new = as<arma::vec>(variants_to_test(17));
          size_dev_nta_new = as<arma::vec>(variants_to_test(18));
          sizeb_dev_nta_new = as<arma::vec>(variants_to_test(19));
          sizec_dev_nta_new = as<arma::vec>(variants_to_test(20));
          repst_dev_nta_new = as<arma::vec>(variants_to_test(21));
          fec_dev_nta_new = as<arma::vec>(variants_to_test(22));
          jsurv_dev_nta_new = as<arma::vec>(variants_to_test(23));
          jobs_dev_nta_new = as<arma::vec>(variants_to_test(24));
          jsize_dev_nta_new = as<arma::vec>(variants_to_test(25));
          jsizeb_dev_nta_new = as<arma::vec>(variants_to_test(26));
          jsizec_dev_nta_new = as<arma::vec>(variants_to_test(27));
          jrepst_dev_nta_new = as<arma::vec>(variants_to_test(28));
          jmatst_dev_nta_new = as<arma::vec>(variants_to_test(29));
          indcova_nta_new = as<arma::vec>(variants_to_test(30));
          indcovb_nta_new = as<arma::vec>(variants_to_test(31));
          indcovc_nta_new = as<arma::vec>(variants_to_test(32));
          
          // Rcout << "ESS_optimizer_fb M" << endl;
          // New stageexpansion
          IntegerVector chosen_int = {0};
          DataFrame variants_to_test_main = AdaptUtils::df_indices(variants_to_test, 0);
          // Rcout << "ESS_optimizer_fb N" << endl;
          DataFrame variants_to_test_995 = AdaptUtils::df_indices(variants_to_test, 1);
          
          DataFrame stageexpansion_main = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_main, firstage_int, finalage_int, ehrlen, style,
            filter);
          // Rcout << "ESS_optimizer_fb O" << endl;
          DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
            variants_to_test_995, firstage_int, finalage_int, ehrlen, style,
            filter);
          
          // Rcout << "ESS_optimizer_fb P" << endl;
          chosen_int = 1;
          StringVector focused_var = {"mpm_altered"};
          DataFrame stageexpansion_main_reduced = LefkoUtils::df_subset(stageexpansion_main,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_fb Q" << endl;
          
          DataFrame stageexpansion_995_reduced = LefkoUtils::df_subset(stageexpansion_995,
            as<RObject>(chosen_int), false, true, false, false, true,
            as<RObject>(focused_var));
          // Rcout << "ESS_optimizer_fb R" << endl;
          
          List sge_to_test = Rcpp::List::create(_["main"] = stageexpansion_main_reduced,
            _["e995"] = stageexpansion_995_reduced);
          // Rcout << "ESS_optimizer_fb S" << endl;
          
          DataFrame ESS_test_values;
          
          invfb_optim_singlerun (ESS_test_values, surv_dev_nta_new, obs_dev_nta_new,
            size_dev_nta_new, sizeb_dev_nta_new, sizec_dev_nta_new, repst_dev_nta_new, 
            fec_dev_nta_new, jsurv_dev_nta_new, jobs_dev_nta_new, jsize_dev_nta_new, 
            jsizeb_dev_nta_new, jsizec_dev_nta_new, jrepst_dev_nta_new,
            jmatst_dev_nta_new, indcova_nta_new, indcovb_nta_new,
            indcovc_nta_new, variant_nta_new, variants_to_test, N_out, comm_out,
            errcheck_mpm, errcheck_mpmout, sge_to_test, used_times, allmodels_all,
            vrm_list, allstages_all, dev_terms_list, ind_terms_num_list,
            ind_terms_cat_list, stageexpansion_ta_devterms_by_variant, sp_density_list,
            start_list, equivalence_list, density_vr_list, current_stageframe,
            current_supplement, density_df, dens_index_df, entry_time_vec,
            sp_density_num_vec, inda_terms_num_vec, indb_terms_num_vec,
            indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
            indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec, year_vec, patch_vec,
            times, fitness_times, format_int, firstage_int, finalage_int,
            dev_terms_times_int, substoch, opt_res, opt_res_orig, exp_tol,
            theta_tol, conv_threshold, sparse_bool, A_only, stages_not_equal,
            integeronly, dens_yn_bool, errcheck, zap_min);
          
          // Rcout << "ESS_optimizer_fb T" << endl;
          NumericVector current_round_fitness_values = as<NumericVector>(ESS_test_values["fitness"]);
          // Rcout << "current_round_fitness_values: " << current_round_fitness_values << endl;
          double res_fitness_test = current_round_fitness_values(0);
          double e995_fitness_test = current_round_fitness_values(1); 
          // Rcout << "new invader res_fitness_test: " << res_fitness_test << endl;
          // Rcout << "new invader e995_fitness_test: " << e995_fitness_test << endl;
          
          NumericVector ref_fitness_values = as<NumericVector>(reference_variants["fitness"]);
          NumericVector abs_fitness_values = abs(ref_fitness_values);
          
          LogicalVector converged = as<LogicalVector>(ESS_test_values["converged"]);
          // Rcout << "ref_fitness_values: " << ref_fitness_values << endl;
          // Rcout << "abs_fitness_values: " << abs_fitness_values << endl;
          
          // Rcout << "Entered search mode determination section" << endl;
          if (search_mode == 0) { // Reflection
            variant_R = clone(ESS_test_values);
            xR_fitness = e995_fitness_test;
            // Rcout << "xR_fitness: " << xR_fitness << endl;
            // Rcout << "e995_fitness_L: " << e995_fitness_L << endl;
            // Rcout << "e995_fitness_H: " << e995_fitness_H << endl;
            if (abs(xR_fitness) < abs(e995_fitness_L)) {
              // Rcout << "ESS_optimizer_fb T1 search mode 0, moving to search mode 1" << endl;
              search_mode = 1;
            } else if (abs(xR_fitness) < abs(e995_fitness_H)) {
              // Rcout << "ESS_optimizer_fb T2 search mode 0, staying in search mode 0" << endl;
              variant_H = clone(ESS_test_values);
              e995_fitness_H = xR_fitness;
              
              if (worst_point == 2) {
                point_2 = clone(variant_H);
                e995_fitness_2 = e995_fitness_H;
              } else if (worst_point == 1) {
                point_1 = clone(variant_H);
                e995_fitness_1 = e995_fitness_H;
              }
            } else {
              // Rcout << "ESS_optimizer_fb T3 search mode 0, moving to search mode 2" << endl;
              search_mode = 2;
            }
            
          } else if (search_mode == 1) { // Expansion
            // Rcout << "ESS_optimizer_fb T4 search mode 1, moving to search mode 0" << endl;
            variant_E = clone(ESS_test_values);
            xE_fitness = e995_fitness_test;
            // Rcout << "xE_fitness: " << xE_fitness << endl;
            // Rcout << "xR_fitness: " << xR_fitness << endl;
            // Rcout << "e995_fitness_L: " << e995_fitness_L << endl;
            // Rcout << "e995_fitness_H: " << e995_fitness_H << endl;
            if (abs(xE_fitness) < abs(xR_fitness)) {
              // Rcout << "abs(xE_fitness) < abs(xR_fitness)" << endl;
              if (worst_point == 2) {
                // Rcout << "Replacing point_2 with xE" << endl;
                point_2 = clone(variant_E);
                e995_fitness_2 = xE_fitness;
              } else if (worst_point == 1) {
                // Rcout << "Replacing point_1 with xE" << endl;
                point_1 = clone(variant_E);
                e995_fitness_1 = xE_fitness;
              }
              
              //variant_H = clone(variant_E);
              //e995_fitness_H = xE_fitness;
            } else {
              // Rcout << "abs(xE_fitness) >= abs(xR_fitness)" << endl;
              if (worst_point == 2) {
                // Rcout << "Replacing point_2 with xR" << endl;
                point_2 = clone(variant_R);
                e995_fitness_2 = xR_fitness;
              } else if (worst_point == 1) {
                // Rcout << "Replacing point_1 with xR" << endl;
                point_1 = clone(variant_R);
                e995_fitness_1 = xR_fitness;
              }
              
              //variant_H = clone(variant_R);
              //e995_fitness_H = xR_fitness;
            }
            search_mode = 0;
            
          } else if (search_mode == 2) { // Contraction
            variant_C_new = clone(ESS_test_values);
            xC_new_fitness = e995_fitness_test;
            // Rcout << "xC_new_fitness: " << xC_new_fitness << endl;
            // Rcout << "e995_fitness_L: " << e995_fitness_L << endl;
            // Rcout << "e995_fitness_H: " << e995_fitness_H << endl;
            
            if (worst_point == 1) {
              if ((abs(xC_new_fitness) < abs(xR_fitness)) && (abs(xR_fitness) < abs(e995_fitness_1))) {
                // Rcout << "ESS_optimizer_fb T5 search mode 2, moving to search mode 0" << endl;
                point_1 = clone(variant_C_new);
                e995_fitness_1 = xC_new_fitness;
                search_mode = 0;
              } else if ((abs(xC_new_fitness) < abs(e995_fitness_1)) && (abs(xR_fitness) >= abs(e995_fitness_1))) {
                // Rcout << "ESS_optimizer_fb T6 search mode 2, moving to search mode 0" << endl;
                point_1 = clone(variant_C_new);
                e995_fitness_1 = xC_new_fitness;
                search_mode = 0;
              } else search_mode = 3;
            } else if (worst_point == 2) {
              if ((abs(xC_new_fitness) < abs(xR_fitness)) && (abs(xR_fitness) < abs(e995_fitness_2))) {
                // Rcout << "ESS_optimizer_fb T7 search mode 2, moving to search mode 0" << endl;
                point_2 = clone(variant_C_new);
                e995_fitness_2 = xC_new_fitness;
                search_mode = 0;
              } else if ((abs(xC_new_fitness) < abs(e995_fitness_2)) && (abs(xR_fitness) >= abs(e995_fitness_2))) {
                // Rcout << "ESS_optimizer_fb T8 search mode 2, moving to search mode 0" << endl;
                point_2 = clone(variant_C_new);
                e995_fitness_2 = xC_new_fitness;
                search_mode = 0;
              } else search_mode = 3;
            }
          }
          
        } else if (search_mode == 3) {
          // Rcout << "ESS_optimizer_fb T7 search mode 3, moving to search mode 0" << endl;
          search_mode = 0;
          
          NumericVector point1_fitness_values = as<NumericVector>(point_1["fitness"]);
          NumericVector point2_fitness_values = as<NumericVector>(point_2["fitness"]);
          e995_fitness_1 = point1_fitness_values(1); 
          e995_fitness_2 = point2_fitness_values(1); 
          
          if (abs(e995_fitness_1) < abs(e995_fitness_2)) {
            variant_L = clone(point_1);
            variant_H = clone(point_2);
            e995_fitness_L = e995_fitness_1;
            e995_fitness_H = e995_fitness_2;
          } else {
            variant_H = clone(point_1);
            variant_L = clone(point_2);
            e995_fitness_L = e995_fitness_2;
            e995_fitness_H = e995_fitness_1;
          }
            
          variant_C = clone(variant_L);
        }
        // Rcout << "ESS_optimizer_fb T8" << endl;
        
      }
      loop_tracker++;
      // Rcout << "ESS_optimizer_fb W" << endl;
      
      if (abs(e995_fitness_1) <= conv_threshold || abs(e995_fitness_2) <= conv_threshold ||
          abs(e995_fitness_1 - e995_fitness_2) <= conv_threshold) {
        //bool found_optimum {true};
        if (abs(e995_fitness_1) <= conv_threshold || abs(e995_fitness_1 - e995_fitness_2) <= conv_threshold) {
          // Rcout << "Exporting point 1 on convergence" << endl;
          ESS_test_values = clone(point_1);
        } else {
          // Rcout << "Exporting point 2 on convergence" << endl;
          ESS_test_values = clone(point_2);
        }
        
        IntegerVector next_chosen_one = {0};
        ESS_out_values = AdaptUtils::df_indices(ESS_test_values, next_chosen_one);
        IntegerVector further_one = {1};
        DataFrame NeededFitnessRow = AdaptUtils::df_indices(ESS_test_values, further_one);
        
        NumericVector true_fitness_vec = as<NumericVector>(NeededFitnessRow["fitness"]);
        ESS_out_values["fitness"] = true_fitness_vec;
        
        found_converged = true;
        LogicalVector found_converged_vec = {found_converged};
        ESS_out_values["converged"] = found_converged_vec;
        
        opt_needed = false;
        break;
      }
    }
    // Rcout << "ESS_optimizer_fb X" << endl;
    //IntegerVector next_chosen_one = {0};
    //ESS_out_values = AdaptUtils::df_indices(ESS_test_values, next_chosen_one);
    
    final_output(i) = ESS_out_values;
  }
  
  // Rcout << "ESS_optimizer_fb Y finished all optimization loops" << endl;
  DataFrame final_output_df;
  
  if (total_optima > 1) {
    DataFrame df1 = as<DataFrame>(final_output(0));
    DataFrame df2 = as<DataFrame>(final_output(1));
    CharacterVector df1_names = as<CharacterVector>(df1.names());
    CharacterVector df2_names = as<CharacterVector>(df2.names());
    
    //int df1_rows = static_cast<int>(df1.nrows());
    //int df2_rows = static_cast<int>(df2.nrows());
    
    DataFrame final_output_df_pre = LefkoUtils::df_rbind(as<DataFrame>(final_output(0)), as<DataFrame>(final_output(1)));
    
    // Rcout << "ESS_optimizer_fb Z" << endl;
    if (total_optima > 2) {
      for (int i = 2; i < total_optima; i++) {
        final_output_df_pre = LefkoUtils::df_rbind(final_output_df_pre, as<DataFrame>(final_output(i)));
      }
    }
    
    // Rcout << "ESS_optimizer_fb AA" << endl;
    int fodfp_size = final_output_df_pre.nrows();
    IntegerVector new_variant_index = seq(1, fodfp_size);
    final_output_df_pre["variant"] = new_variant_index;
    
    final_output_df = final_output_df_pre;
  } else if (total_optima == 1) {
    // Rcout << "ESS_optimizer_fb AB" << endl;
    DataFrame final_output_df_pre = as<DataFrame>(final_output(0));
    
    int fodfp_size = final_output_df_pre.nrows();

    IntegerVector new_variant_index = seq(1, fodfp_size);
    final_output_df_pre["variant"] = new_variant_index;
    
    final_output_df = final_output_df_pre;
  } else {
    final_output_df = R_NilValue;
  }
  
  // Rcout << "ESS_optimizer_fb AC" << endl;
  int optim_rows = static_cast<int>(final_output_df.nrows());
  
  ESS_Lyapunov = final_output_df;
}

//' Core Pre-Existing MPM Projection Engine
//' 
//' Function \code{invpre_project} runs the projections in function
//' \code{invade3_pre_core}.
//' 
//' @name invpre_project
//' 
//' @param var_run_mat A matrix giving the the variants to be run in each
//' projection, with rows giving the projections and columns giving the
//' variants.
//' @param N_out_pre The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out_pre The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param new_stageexpansion_list A list with stage expansions for all trait
//' axis data leading to matrix element changes with each list element
//' corresponding to each respective variant.
//' @param errcheck_mpm_reps An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param used_times A list of year numbers for each time per run.
//' @param zero_stage_vec_list A list of population stage vectors full of zeros.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param A_list A list of allA matrices.
//' @param U_list A list of all U matrices.
//' @param F_list A list of all F matrices.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param var_per_run An integer giving the number of variants per run.
//' @param times An integer giving the number of occasions to project.
//' @param var_mat_length An integer giving the number of rows in the variant
//' matrix.
//' @param format_int An integer giving the MPM format.
//' @param current_rep The integer giving the current replicate.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' 
//' @return Arguments 2 through 7 are directly manipulated without any values
//' returned.
//' 
//' @keywords internal
//' @noRd
inline void invpre_project (const arma::mat var_run_mat, List& N_out_pre,
  List& comm_out_pre, List& new_stageexpansion_list, List& errcheck_mpm_reps,
  List& used_times, List& zero_stage_vec_list, const List start_list,
  const List equivalence_list, const List A_list, const List U_list,
  const List F_list, const DataFrame density_df, const DataFrame dens_index_df,
  const IntegerVector entry_time_vec, const int var_per_run, const int times,
  const int var_mat_length, const int format_int, const int current_rep,
  const int firstage_int, const int finalage_int, const int substoch,
  const double exp_tol, const double theta_tol, const bool err_check,
  const bool err_check_extreme, const bool sparse_bool, const bool A_only,
  const bool stages_not_equal, const bool integeronly, const bool dens_yn_bool) {
    
  // Rcout << "invpre_project A" << endl;
  int i = current_rep;
  
  List running_popvecs; //  = clone(start_list)
  List running_popvecs_startonly; //  = clone(start_list)
  arma::cube N_mpm (var_per_run, (times + 1), var_mat_length); // rows = vars, cols = times, slices = permutes 
  
  List errcheck_mpm_reps_time (times); // Could remove later
  for (int j = 0; j < times; j++) { // 2nd loop - time j
    // Rcout << "invpre_project r          ";
    if (j % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    List errcheck_mpm_reps_time_vmt (var_mat_length); // Could remove later
    
    for (int l = 0; l < var_mat_length; l++) { // 3rd loop - permutes l
      // Rcout << "invpre_project s          ";
      List errcheck_mpm_reps_time_vmt_var(var_per_run); // Could remove later
      // Rcout << "current_permutation (l): " << l << "          ";
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
        // Rcout << "invpre_project t          ";
        int current_variant_index = var_run_mat(l, m); // Equivalent to index integer k
        
        DataFrame sge_current = as<DataFrame>(new_stageexpansion_list(current_variant_index));
        
        List pop_reps = as<List>(all_pops_per_run(m));
        arma::mat pops_out = as<arma::mat>(pop_reps(i));
        
        if (j > (entry_time_vec(m) - 1)) {
          // Rcout << "invpre_project u          ";
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
          
          // Rcout << "invpre_project v          ";
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
              // Rcout << "invpre_project aa          ";
            } else {
              // Rcout << "invpre_project ab          ";
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
              // Rcout << "invpre_project ac          ";
            }
          }
          
          // Rcout << "invpre_project ad          ";
          if (integeronly) running_popvec_mpm = floor(running_popvec_mpm);
          double N_current = arma::sum(running_popvec_mpm);
          N_mpm(m, (j + 1), l) = N_current; // Used to be (k, (j + 1))
          
          running_popvecs(m) = running_popvec_mpm;
          pops_out.col(j + 1) = running_popvec_mpm;
          // Rcout << "invpre_project ae          ";
        } else {
          // Rcout << "invpre_project af          ";
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
  
  // Rcout << "invpre_project B" << endl;
}

//' Core Pre-Existing MPM Projection Engine for ESS Evaluation
//' 
//' Function \code{invpre_optim} runs the optimization projections in function
//' \code{invade3_pre_core} used to estimate ESS values.
//' 
//' @name invpre_optim
//' 
//' @param N_out_pre The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out_pre The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param new_stageexpansion_list A list with stage expansions for all
//' variant data used in ESS evaluation. This list includes an extra layer of
//' list elements, corresponding to the optim_ta and optim_ta_995 data.
//' @param errcheck_mpm_reps An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param used_times A list of year numbers for each time per run.
//' @param zero_stage_vec_list A list of population stage vectors full of zeros.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param A_list A list of all A matrices.
//' @param U_list A list of all U matrices.
//' @param F_list A list of all F matrices.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param var_per_run An integer giving the number of variants per run.
//' @param times An integer giving the number of occasions to project.
//' @param var_mat_length An integer giving the number of rows in the variant
//' matrix.
//' @param format_int An integer giving the MPM format.
//' @param current_rep The integer giving the current replicate.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis. Note that the version used in this
//' function is actually equivalent to \code{opt_res_true}.
//' @param opt_res_orig The original value of \code{opt_res}, prior to the
//' determination of the number of variable traits. Equal to \code{opt_res} if
//' the number of variable traits is 1, and to the square root of \code{opt_res}
//' if the number of variable traits is 2.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' @param threshold The lower limit for the absolute value of fitness, below
//' which fitness is rounded to 0. Defaults to 0.00000001.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' 
//' @return Arguments 2 through 7 are directly manipulated without any values
//' returned.
//' 
//' @keywords internal
//' @noRd
inline void invpre_optim (List& N_out_pre, List& comm_out_pre,
  List& new_stageexpansion_list, List& errcheck_mpm_reps, List& used_times,
  List& zero_stage_vec_list, const List start_list, const List equivalence_list,
  const List A_list, const List U_list, const List F_list,
  const DataFrame density_df, const DataFrame dens_index_df,
  const IntegerVector entry_time_vec, const int var_per_run, const int times,
  const int var_mat_length, const int format_int, const int current_rep,
  const int firstage_int, const int finalage_int, const int substoch,
  const int opt_res, const int opt_res_orig, const double exp_tol,
  const double theta_tol, const double threshold, const bool err_check,
  const bool err_check_extreme, const bool sparse_bool, const bool A_only,
  const bool stages_not_equal, const bool integeronly, const bool dens_yn_bool) {
  
  // Rcout << "invpre_optim A" << endl;
  int i = current_rep;
  
  List running_popvecs; //  = clone(start_list)
  List running_popvecs_startonly; //  = clone(start_list)
  arma::cube N_mpm (2, (times + 1), opt_res); // rows = vars, cols = times, slices = permutes 
  
  List errcheck_mpm_reps_time (times); // Could remove later
  
  for (int j = 0; j < times; j++) { // 2nd loop - time j
    // Rcout << "invpre_optim A1          ";
    if (j % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    List errcheck_mpm_reps_time_vmt (opt_res); // Could remove later
    
    for (int l = 0; l < opt_res; l++) { // 3rd loop - permutes l
      // Rcout << "invpre_optim A2          ";
      List errcheck_mpm_reps_time_vmt_var(2); // Could remove later
      // Rcout << "current_permutation (l): " << l << "          ";
      if (j == 0) {
        List var_popvecs_to_start (2);
        for (int n = 0; n < var_per_run; n++) {
          var_popvecs_to_start(n) = as<arma::vec>(start_list(static_cast<int>(0)));
        }
        running_popvecs = var_popvecs_to_start;
        running_popvecs_startonly = clone(var_popvecs_to_start);
      }
      
      for (int m = 0; m < 2; m++) {
        if (j == entry_time_vec(m)) {
          arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
          
          double N_current = accu(running_popvec_mpm);
          N_mpm(m, j, l) = N_current; // Used to be (k, j)
        }
      }
      
      List all_pops_per_run = as<List>(comm_out_pre(l));
      
      { // main variant run
        int m = 0;
        // Rcout << "invpre_optim A3          ";
        int current_variant_index = l; // Equivalent to index integer k
        
        List sge_list = as<List>(new_stageexpansion_list(current_variant_index)); // Does this lose year differences?
        DataFrame sge_current = as<DataFrame>(sge_list(m));
        
        List pop_reps = as<List>(all_pops_per_run(m));
        arma::mat pops_out = as<arma::mat>(pop_reps(i));
          
        if (j > (entry_time_vec(m) - 1)) {
          // Rcout << "invpre_optim A4          ";
          List used_times_per_run = as<List>(used_times(0));
          List used_times_current_var = as<List>(used_times_per_run(m));
          IntegerVector current_times_vec = as<IntegerVector>(used_times_current_var(i));
          
          arma::vec running_popvec_mpm;
          if (j == entry_time_vec(m)) {
            running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
            pops_out.col(j) = running_popvec_mpm;
            
          } else {
            running_popvec_mpm = pops_out.col(j);
          }
          
          // Rcout << "invpre_optim A5          ";
          if (!dens_yn_bool) {
            if (!sparse_bool) {
              // Rcout << "invpre_optim A6          ";
              arma::mat current_A = as<arma::mat>(A_list(current_times_vec(j)));
              
              // Rcout << "invpre_optim A7          ";
              if (A_only) {
                AdaptUtils::Amat_alter(current_A, sge_current); 
              } else {
                arma::mat current_U_unaltered = as<arma::mat>(U_list(current_times_vec(j)));
                arma::mat current_F_unaltered = as<arma::mat>(F_list(current_times_vec(j)));
                AdaptUtils::UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
              }
              
              // Rcout << "invpre_optim A8          ";
              if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
              running_popvec_mpm = current_A * running_popvec_mpm; 
            } else {
              // Rcout << "invpre_optim A9          ";
              arma::sp_mat current_A = as<arma::sp_mat>(A_list(current_times_vec(j)));
              
              // Rcout << "invpre_optim A10          ";
              if (A_only) {
                AdaptUtils::sp_Amat_alter(current_A, sge_current);
              } else {
                arma::sp_mat current_U_unaltered = as<arma::sp_mat>(U_list(current_times_vec(j)));
                arma::sp_mat current_F_unaltered = as<arma::sp_mat>(F_list(current_times_vec(j)));
                AdaptUtils::sp_UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
              }
              // Rcout << "invpre_optim A11          ";
              if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
              running_popvec_mpm = current_A * running_popvec_mpm;
              // Rcout << "invpre_optim A12          ";
            }
          } else {
            // Rcout << "invpre_optim A13          ";
            DataFrame used_density_input = density_df;
            DataFrame used_density_index_input = as<DataFrame>(dens_index_df);
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            // Rcout << "invpre_optim A14          ";
            if (j >= (used_delay - 1 )) { // Change to allow different delay Ns for different entries
              if (!stages_not_equal) {
                // Rcout << "invpre_optim A15          ";
                arma::mat di_mat = N_mpm.slice(l);
                arma::vec delay_issue = di_mat.col(j + 1 - used_delay);
                
                // Rcout << "invpre_optim A16          ";
                double delay_N_sum = arma::sum(delay_issue);
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                // Rcout << "invpre_optim A17          ";
                AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, sge_current, A_list, delay_N_sum,
                  static_cast<int>(current_times_vec(j)), integeronly,
                  substoch, used_density_input, used_density_index_input,
                  false, sparse_bool, sparse_bool, false, err_check);
                
                // Rcout << "invpre_optim A18          ";
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) { // Could remove later
                  if (!sparse_bool) {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                  } else {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                  }
                } // Could remove later
                // Rcout << "invpre_optim A19          ";
              } else {
                // Rcout << "invpre_optim A20          ";
                double delay_N_sum {0.0};
                
                if (j > 0) {
                  for (int p = 0; p < 2; p++) {
                    int current_variant_index_agg = l;
                    
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
                
                // Rcout << "invpre_optim A21          ";
                arma::vec new_popvec;
                arma::mat new_projmat;
                arma::sp_mat new_projsp;
                
                AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                  running_popvec_mpm, sge_current, A_list, delay_N_sum,
                  static_cast<int>(current_times_vec(j)), integeronly,
                  substoch, used_density_input, used_density_index_input,
                  false, sparse_bool, sparse_bool, false, err_check);
                
                // Rcout << "invpre_optim A22          ";
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) { // Could remove later
                  if (!sparse_bool) {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                  } else {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                  }
                } // Could remove later
              }
              // Rcout << "invpre_optim A23          ";
            } else {
              // Rcout << "invpre_optim A24          ";
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                running_popvec_mpm, sge_current, A_list, 0.0,
                static_cast<int>(current_times_vec(j)), integeronly, substoch,
                used_density_input, used_density_index_input, false,
                sparse_bool, sparse_bool, false, err_check);
              
              // Rcout << "invpre_optim A25          ";
              running_popvec_mpm = new_popvec;
              if (err_check_extreme) { // Could remove later
                if (!sparse_bool) {
                  errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                } else {
                  errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                }
              } // Could remove later
              // Rcout << "invpre_optim A26          ";
            }
          }
          
          // Rcout << "invpre_optim A27          ";
          if (integeronly) running_popvec_mpm = floor(running_popvec_mpm);
          double N_current = arma::sum(running_popvec_mpm);
          N_mpm(m, (j + 1), l) = N_current; // Used to be (k, (j + 1))
          
          running_popvecs(m) = running_popvec_mpm;
          pops_out.col(j + 1) = running_popvec_mpm;
          // Rcout << "invpre_optim A28          ";
        } else {
          // Rcout << "invpre_optim A29          ";
          arma::vec current_zero_vec = as<arma::vec>(zero_stage_vec_list(0));
          pops_out.col(j + 1) = current_zero_vec;
        }
        pop_reps(i) = pops_out;
      } // m loop - var_per_run
      // Rcout << "invpre_optim A30          ";
      if (err_check_extreme) errcheck_mpm_reps_time_vmt(l) = errcheck_mpm_reps_time_vmt_var;
      // Rcout << "invpre_optim A31          ";
      
      { // 0.995 variant run
        int m = 1;
        // Rcout << "invpre_optim A32          ";
        int current_variant_index = l; // Equivalent to index integer k
        
        List sge_list = as<List>(new_stageexpansion_list(current_variant_index));
        DataFrame sge_current = as<DataFrame>(sge_list(m));
        
        List pop_reps = as<List>(all_pops_per_run(m));
        arma::mat pops_out = as<arma::mat>(pop_reps(i));
        
        if (j > (entry_time_vec(m) - 1)) {
          // Rcout << "invpre_optim A33          ";
          List used_times_per_run = as<List>(used_times(0));
          List used_times_current_var = as<List>(used_times_per_run(m));
          IntegerVector current_times_vec = as<IntegerVector>(used_times_current_var(i));
          
          arma::vec running_popvec_mpm;
          if (j == entry_time_vec(m)) {
            running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
            pops_out.col(j) = running_popvec_mpm;
            
          } else {
            running_popvec_mpm = pops_out.col(j);
          }
          
          // Rcout << "invpre_optim A34          ";
          if (!dens_yn_bool) {
            if (!sparse_bool) {
              // Rcout << "invpre_optim A35          ";
              arma::mat current_A = as<arma::mat>(A_list(current_times_vec(j)));
              
              if (A_only) {
                AdaptUtils::Amat_alter(current_A, sge_current); 
              } else {
                arma::mat current_U_unaltered = as<arma::mat>(U_list(current_times_vec(j)));
                arma::mat current_F_unaltered = as<arma::mat>(F_list(current_times_vec(j)));
                AdaptUtils::UFmat_alter(current_A, current_U_unaltered, current_F_unaltered, sge_current);
              }
              
              // Rcout << "invpre_optim A36          ";
              if (err_check_extreme) errcheck_mpm_reps_time_vmt_var(m) = current_A; // Could remove later
              running_popvec_mpm = current_A * running_popvec_mpm; 
            } else {
              // Rcout << "invpre_optim A37          ";
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
              // Rcout << "invpre_optim A38          ";
            }
          } else {
            // Rcout << "invpre_optim A39          ";
            DataFrame used_density_input = density_df;
            DataFrame used_density_index_input = as<DataFrame>(dens_index_df);
            
            IntegerVector ud_delay_vec = as<IntegerVector>(used_density_input["time_delay"]);
            int used_delay = max(ud_delay_vec);
            
            // Rcout << "invpre_optim A40          ";
            if (j >= (used_delay - 1 )) { // Change to allow different delay Ns for different entries
              if (!stages_not_equal) {
                // Rcout << "invpre_optim A41          ";
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
                
                // Rcout << "invpre_optim A42          ";
                running_popvec_mpm = new_popvec;
                if (err_check_extreme) { // Could remove later
                  if (!sparse_bool) {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                  } else {
                    errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                  }
                } // Could remove later
                // Rcout << "invpre_optim A43          ";
              } else {
                // Rcout << "invpre_optim A44          ";
                double delay_N_sum {0.0};
                
                if (j > 0) {
                  for (int p = 0; p < 2; p++) {
                    int current_variant_index_agg = l;
                    
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
                
                // Rcout << "invpre_optim A46          ";
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
              // Rcout << "invpre_optim A47          ";
            } else {
              // Rcout << "invpre_optim A48          ";
              arma::vec new_popvec;
              arma::mat new_projmat;
              arma::sp_mat new_projsp;
              
              AdaptUtils::proj3dens_inv(new_popvec, new_projmat, new_projsp,
                running_popvec_mpm, sge_current, A_list, 0.0,
                static_cast<int>(current_times_vec(j)), integeronly, substoch,
                used_density_input, used_density_index_input, false,
                sparse_bool, sparse_bool, false, err_check);
              
              // Rcout << "invpre_optim A49          ";
              running_popvec_mpm = new_popvec;
              if (err_check_extreme) { // Could remove later
                if (!sparse_bool) {
                  errcheck_mpm_reps_time_vmt_var(m) = new_projmat;
                } else {
                  errcheck_mpm_reps_time_vmt_var(m) = new_projsp;
                }
              } // Could remove later
              // Rcout << "invpre_optim ac          ";
            }
          }
          
          // Rcout << "invpre_optim A50          ";
          if (integeronly) running_popvec_mpm = floor(running_popvec_mpm);
          double N_current = arma::sum(running_popvec_mpm);
          N_mpm(m, (j + 1), l) = N_current; // Used to be (k, (j + 1))
          
          // Rcout << "invpre_optim A51          ";
          running_popvecs(m) = running_popvec_mpm;
          pops_out.col(j + 1) = running_popvec_mpm;
          // Rcout << "invpre_optim A52          ";
        } else {
          // Rcout << "invpre_optim A53          ";
          arma::vec current_zero_vec = as<arma::vec>(zero_stage_vec_list(0));
          pops_out.col(j + 1) = current_zero_vec;
        }
        pop_reps(i) = pops_out;
      } // m loop - var_per_run
      if (err_check_extreme) errcheck_mpm_reps_time_vmt(l) = errcheck_mpm_reps_time_vmt_var;
    } // l loop - var_mat_length
    if (err_check_extreme) errcheck_mpm_reps_time(j) = errcheck_mpm_reps_time_vmt;
    // Rcout << "invpre_optim A54          ";
  } // j loop - time
  if (err_check_extreme) errcheck_mpm_reps(i) = errcheck_mpm_reps_time; // Could remove later
  N_out_pre(i) = N_mpm;
  
  // Rcout << "invpre_optim B" << endl;
}

//' Set-Up Running Invasion Analyses of Existing MPMs
//' 
//' Function \code{invade3_pre_core} is the main function running invasion
//' analyses with existing MPMs supplied by the user.
//' 
//' @name invade3_pre_core
//' 
//' @param Lyapunov The main data frame giving the Lyapunov coefficients
//' estimated, as well as the circumstances resulting in them.
//' @param Lyapunov_optim Main data frame giving Lyapunov coefficients for all
//' trait combinations developed for the ESS optima table. Holds elasticity
//' fitness values.
//' @param ESS_Lyapunov A data frame provided by reference that will hold the
//' ESS optima.
//' @param var_run_mat A matrix giving the order of trait variants in each
//' run.
//' @param N_out The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param N_out_optim The main list of final population sizes from ESS
//' optimization, supplied as a reference and altered by this function.
//' @param comm_out_optim The main list of full projection results for the
//' community resulting from ESS optimization, supplied as a pointer and altered
//' by this function.
//' @param zero_stage_vec_list A list of vectors giving zero stage vectors for
//' each MPM, if entry times are staggered.
//' @param trait_axis A data frame of class \code{adaptAxis} holding the trait
//' data to test.
//' @param new_trait_axis A data frame giving trait axis data post-processing
//' with function \code{ta_reassess()}.
//' @param optim_trait_axis A data frame giving trait axis data processed for
//' ESS optimization.
//' @param optim_trait_axis_995 A data frame giving trait axis data processed
//' for ESS optimization for variants 99.5% the values of the core frame.
//' @param new_stageexpansion_list A list with stage expansions for all trait
//' axis data leading to matrix element changes with each list element
//' corresponding to each respective variant.
//' @param new_stageexpansion_list_optima A list with stage expansions for all
//' variant data used in ESS evaluation.
//' @param errcheck_mpms An optional list of all MPMs post-processing. Only
//' output if \code{err_check = TRUE}.
//' @param errcheck_mpms_optim An optional list of all MPMs used in ESS optima
//' evaluation. Only output if \code{err_check = TRUE}.
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
//' @param main_optim_res An integer giving the number of variants being tested.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param loop_max The maximum number of times to attempt optimization, if
//' convergence does not occur.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param stochastic A Boolean value indicating to perform a temporally
//' stochastic projection.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used in the MPM, given through a \code{lefkoDens} object.
//' @param entry_time_vec_use A Boolean value indicating whether entry times
//' differ among MPMs.
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
//' @param threshold The lower limit for the absolute value of fitness, below
//' which fitness is rounded to 0. Defaults to 0.00000001.
//' @param fitness_table A Boolean value dictating whether to include a data
//' frame giving Lyapunov coefficients for all combinations of variants tested.
//' Necessary for the creation of pairwise invasibility plots (PIPs). Defaults
//' to \code{TRUE}.
//' @param ESS_optima A logical value indicating whether to assess the values of
//' ESS optima for traits that vary among variants in the given \code{adaptAxis}
//' table. Defaults to \code{TRUE}.
//' @param elast_mult A multiplier for traits to assess the elasticity of
//' fitness in trait optimization. Defaults to 0.995.
//' @param zap_min A Boolean value describing whether to round fitness values
//' below the value given in \code{threshold}.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
void invade3_pre_core (DataFrame& Lyapunov, DataFrame& Lyapunov_optim,
  DataFrame& ESS_Lyapunov, const arma::mat& var_run_mat, List& N_out,
  List& comm_out, List& N_out_optim, List& comm_out_optim,
  List& zero_stage_vec_list, DataFrame& trait_axis, DataFrame& new_trait_axis,
  DataFrame& optim_trait_axis, DataFrame& optim_trait_axis_995,
  List& new_stageexpansion_list, List& new_stageexpansion_list_optima,
  List& errcheck_mpms, List& errcheck_mpms_optima, const List chosen_mpm,
  const List tweights_list, const List start_list, const List vrm_list,
  DataFrame stageframe_df, const List allmodels_all, const List allstages_all,
  const DataFrame supplement_df, const CharacterVector chosen_years,
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
  const int substoch, const int format_int, const int preexisting_mpm_size,
  const int firstage_int, const int finalage_int, const int main_optim_res,
  const int opt_res, const double exp_tol, const double theta_tol,
  const int loop_max, const bool integeronly, const bool stages_not_equal,
  const bool stochastic, const bool dens_yn_bool, const bool entry_time_vec_use,
  const bool sparse_bool, const bool historical, const bool pure_leslie,
  const bool A_only, const bool err_check, const bool err_check_extreme,
  const double threshold, const bool fitness_table, const bool ESS_optima,
  double elast_mult, const bool zap_min) {
  
  // Structures for optim
  List N_out_pre_optim (nreps);
  DataFrame ESS_trait_axis;
  IntegerVector ESS_var_traits;
  IntegerVector flipped_traits (20);
  
  int ehrlen_optim {0};
  int style_optim {0};
  int filter_optim {0};
    
  // patches?????
  
  // Rcout << "invade3_pre_core a          ";
  int var_mat_length = static_cast<int>(var_run_mat.n_rows);
  bool opt_res_squared {false};
  int opt_res_true = opt_res;
  
  List A_list = as<List>(chosen_mpm["A"]);
  
  List U_list;
  List F_list;
  
  if (!A_only) {
    U_list = as<List>(chosen_mpm["U"]);
    F_list = as<List>(chosen_mpm["F"]);
  }
  
  // Rcout << "invade3_pre_core b          ";
  // Vectors of axis variants
  IntegerVector axis_variant_vec = as<IntegerVector>(trait_axis["variant"]);
  IntegerVector axis_variants_unique = sort_unique(axis_variant_vec);
  
  DataFrame trait_axis_clone = clone(trait_axis);
  new_trait_axis = AdaptUtils::ta_reassess(stageframe_df, trait_axis_clone, firstage_int,
    historical, false, pure_leslie);
  
  if (ESS_optima) {
    AdaptUtils::optim_ta_setup(new_trait_axis, ESS_trait_axis, optim_trait_axis,
      optim_trait_axis_995, ESS_var_traits, flipped_traits, opt_res,
      elast_mult);
    
    //int var_traits = sum(ESS_var_traits);
    if (opt_res_squared) opt_res_true = opt_res * opt_res; /////
  }
  
  List comm_out_pre_optim (opt_res_true);
  List trait_axis_by_variant_optim (opt_res_true); // Might wish to remove this later
  
  if (ESS_optima) {
    List stageexpansion_by_variant_optim (opt_res_true); // Might wish to remove this later
    // Set up comm_out_pre_optim 
    for (int i = 0; i < opt_res_true; i++) {
      // Rcout << "invade3_pre_core c          ";
      //int chosen_mats_length = static_cast<int>(chosen_mats.length());
      
      List all_pops_per_run (2);
      List all_used_times_per_run (2);
      for (int m = 0; m < 2; m++) {
        // Rcout << "invade3_pre_core d          ";
        //int current_variant_index = var_run_mat(i, m);
        
        List pop_reps (nreps);
        List used_times_reps (nreps);
        for (int j = 0; j < nreps; j++) {
          // Rcout << "invade3_pre_core e          ";
          arma::mat pops_out_pre (stagecounts, (times + 1), fill::zeros);
          pop_reps(j) = pops_out_pre;
        }
        //all_used_times_per_run(m) = used_times_reps;
        all_pops_per_run(m) = pop_reps;
      }
      //used_times(i) = all_used_times_per_run;
      comm_out_pre_optim(i) = all_pops_per_run;
    }
  
    // Rcout << "invade3_pre_core e1          ";
    for (int i = 0; i < opt_res_true; i++) {
      IntegerVector used_i = {i + 1};
      StringVector focused_var = {"variant"};
      // Rcout << "invade3_pre_core e2          ";
      DataFrame current_optim_trait_axis = LefkoUtils::df_subset(optim_trait_axis,
        as<RObject>(used_i), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "invade3_pre_core e2a          ";
      DataFrame current_optim_trait_axis_995 = LefkoUtils::df_subset(optim_trait_axis_995,
        as<RObject>(used_i), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "invade3_pre_core e3          ";
      List tabvo = Rcpp::List::create(_["main"] = current_optim_trait_axis,
        _["e995"] = current_optim_trait_axis_995);
      trait_axis_by_variant_optim(i) = tabvo;
      // Rcout << "invade3_pre_core e4          ";
      
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
      
      ehrlen_optim = ehrlen;
      style_optim = style;
      filter_optim = filter;
      
      // Rcout << "invade3_pre_core e5          ";
      DataFrame stageexpansion = AdaptMats::thenewpizzle(stageframe_df,
        current_optim_trait_axis, firstage_int, finalage_int, ehrlen, style,
        filter);
      DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(stageframe_df,
        current_optim_trait_axis_995, firstage_int, finalage_int, ehrlen, style,
        filter);
      // Rcout << "invade3_pre_core e6          ";
      focused_var = {"mpm_altered"};
      IntegerVector chosen_int = {1};
      DataFrame stageexpansion_reduced = LefkoUtils::df_subset(stageexpansion,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      DataFrame stageexpansion_reduced_995 = LefkoUtils::df_subset(stageexpansion_995,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "invade3_pre_core e7          ";
      List sbvo = Rcpp::List::create(_["main"] = stageexpansion_reduced,
        _["e995"] = stageexpansion_reduced_995);
      stageexpansion_by_variant_optim(i) = sbvo;
      // Rcout << "invade3_pre_core e8          ";
    }
    new_stageexpansion_list_optima = stageexpansion_by_variant_optim;
    // Rcout << "invade3_pre_core e9          ";
  }
  
  // Rcout << "invade3_pre_core f          ";
  // Matrix order set up and creation of zero stage vectors
  IntegerVector chosen_mats; // Eliminated matrix_choice_list
  List zvl (variant_count);
  
  IntegerVector chosen_matrix_vec;
  
  DataFrame chosen_labels = as<DataFrame>(chosen_mpm["labels"]);
  CharacterVector chosen_labels_names = as<CharacterVector>(chosen_labels.attr("names"));
  IntegerVector clm_y2 = index_l3(chosen_labels_names, "year2");
  
  CharacterVector mpm_labels_patch = as<CharacterVector>(chosen_labels["patch"]);
  IntegerVector pvis = index_l3(mpm_labels_patch, patch_vec(0)); // Is patch done properly?
  
  // Rcout << "invade3_pre_core g          ";
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
  
  // Rcout << "invade3_pre_core h          ";
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
  
  // Rcout << "invade3_pre_core i          ";
  // Year order determination
  List comm_out_pre (var_mat_length);
  List used_times (var_mat_length);
  
  // Determine order of matrices for each variant run
  for (int i = 0; i < var_mat_length; i++) {
    // Rcout << "invade3_pre_core k          ";
    int chosen_mats_length = static_cast<int>(chosen_mats.length());
    
    List all_pops_per_run (var_per_run);
    List all_used_times_per_run (var_per_run);
    for (int m = 0; m < var_per_run; m++) {
      // Rcout << "invade3_pre_core l          ";
      //int current_variant_index = var_run_mat(i, m);
      
      List pop_reps (nreps);
      List used_times_reps (nreps);
      for (int j = 0; j < nreps; j++) {
        // Rcout << "invade3_pre_core m          ";
        arma::mat pops_out_pre (stagecounts, (times + 1), fill::zeros);
      
        IntegerVector years_topull;
        
        if (!stochastic) {
          // Rcout << "invade3_pre_core n          ";
          IntegerVector years_topull_pre (times);
          
          int mat_tracker {0};
          for (int k = 0; k < times; k++) {
            if (mat_tracker >= chosen_mats_length) mat_tracker = 0;
            
            years_topull_pre(k) = chosen_mats(mat_tracker);
            mat_tracker++;
          }
          years_topull = years_topull_pre;
          // Rcout << "invade3_pre_core o          ";
        } else {
          // Rcout << "invade3_pre_core p          ";
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
            // Rcout << "invade3_pre_core q          ";
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
  
  // Rcout << "invade3_pre_core r          ";
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
  
  // Rcout << "invade3_pre_core s          ";
  // Main projection
  List N_out_pre (nreps);
  List errcheck_mpm_reps (nreps);
  List errcheck_mpm_reps_optima (nreps);
  
  for (int i = 0; i < nreps; i ++) { // 1st loop - reps i
    // Rcout << "invade3_pre_core t          ";
    invpre_project(var_run_mat, N_out_pre, comm_out_pre, new_stageexpansion_list,
      errcheck_mpm_reps, used_times, zero_stage_vec_list, start_list,
      equivalence_list, A_list, U_list, F_list, density_df, dens_index_df,
      entry_time_vec, var_per_run, times, var_mat_length, format_int, i,
      firstage_int, finalage_int, substoch, exp_tol, theta_tol, err_check,
      err_check_extreme, sparse_bool, A_only, stages_not_equal, integeronly,
      dens_yn_bool);
    
    if (ESS_optima) {
      invpre_optim(N_out_pre_optim, comm_out_pre_optim, new_stageexpansion_list_optima,
        errcheck_mpm_reps_optima, used_times, zero_stage_vec_list, start_list,
        equivalence_list, A_list, U_list, F_list, density_df, dens_index_df,
        entry_time_vec, var_per_run, times, var_mat_length, format_int, i,
        firstage_int, finalage_int, substoch, opt_res_true, opt_res, exp_tol,
        theta_tol, threshold, err_check, err_check_extreme, sparse_bool, A_only,
        stages_not_equal, integeronly, dens_yn_bool);
    }
  } // i loop - reps
  comm_out = comm_out_pre;
  comm_out_optim = comm_out_pre_optim;
  
  N_out = N_out_pre;
  N_out_optim = N_out_pre_optim;
  
  if (err_check_extreme) {
    errcheck_mpms = errcheck_mpm_reps;
    errcheck_mpms_optima = errcheck_mpm_reps_optima;
  }
  
  // Rcout << "invade3_pre_core u     ";
  if (fitness_table) {
    AdaptUtils::Lyapunov_creator (Lyapunov, N_out, entry_time_vec, nreps, var_per_run,
      var_mat_length, times, fitness_times, threshold, 0, false, zap_min);
  }
  
  // Rcout << "invade3_pre_core v     ";
  if (ESS_optima) {
    AdaptUtils::Lyapunov_creator (Lyapunov_optim, N_out_optim, entry_time_vec, nreps, 2,
      opt_res_true, times, fitness_times, threshold, 1, false, zap_min);
    
    ESS_optimizer_pre(ESS_Lyapunov, ESS_trait_axis, Lyapunov_optim,
      optim_trait_axis, ESS_var_traits, flipped_traits, new_stageexpansion_list_optima,
      used_times, zero_stage_vec_list, start_list, equivalence_list, A_list,
      U_list, F_list, density_df, dens_index_df, stageframe_df, entry_time_vec,
      times, fitness_times, format_int, stagecounts, firstage_int, finalage_int,
      substoch, exp_tol, theta_tol, sparse_bool, A_only, stages_not_equal,
      integeronly, dens_yn_bool, threshold, opt_res_true, opt_res, ehrlen_optim,
      style_optim, loop_max, filter_optim, elast_mult, zap_min);
  }
}

//' Core Function-Based Projection Engine
//' 
//' Function \code{invfb_project} runs the projections in function
//' \code{invade3_fb_core}.
//' 
//' @name invfb_project
//' 
//' @param var_run_mat A matrix giving the the variants to be run in each
//' projection, with rows giving the projections and columns giving the
//' variants.
//' @param surv_dev_nta The survival column in the reassessed trait axis.
//' @param obs_dev_nta The observation status column in the reassessed trait
//' axis.
//' @param size_dev_nta The primary size column in the reassessed trait axis.
//' @param sizeb_dev_nta The secondary size column in the reassessed trait axis.
//' @param sizec_dev_nta The tertiary size column in the reassessed trait axis.
//' @param repst_dev_nta The reproductive status column in the reassessed trait
//' axis.
//' @param fec_dev_nta The fecundity column in the reassessed trait axis.
//' @param jsurv_dev_nta The juvenile survival column in the reassessed trait
//' axis.
//' @param jobs_dev_nta The juvenile observation status column in the reassessed
//' trait axis.
//' @param jsize_dev_nta The juvenile primary size column in the reassessed
//' trait axis.
//' @param jsizeb_dev_nta The juvenile secondary size column in the reassessed
//' trait axis.
//' @param jsizec_dev_nta The juvenile tertiary size column in the reassessed
//' trait axis.
//' @param jrepst_dev_nta The juvenile reproductive status column in the
//' reassessed trait axis.
//' @param jmatst_dev_nta The juvenile maturity status column in the reassessed
//' trait axis.
//' @param indcova_nta The individual covariate a column in the reassessed trait
//' axis.
//' @param indcovb_nta The individual covariate b column in the reassessed trait
//' axis.
//' @param indcovc_nta The individual covariate c column in the reassessed trait
//' axis.
//' @param variant_nta The variant column in the reassessed trait axis.
//' @param N_out_pre The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out_pre The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param new_stageexpansion_list A list with stage expansions for all trait
//' axis data leading to matrix element changes with each list element
//' corresponding to each respective variant.
//' @param errcheck_mpm_reps An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param errcheck_mpmout_reps An optional list of all mpm_out matrices from MPM
//' processing. Only output if \code{err_check = "extreme"}.
//' @param mdtl The modified dev terms list.
//' @param used_times A list of year numbers for each time per run.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param vrm_list A list of \code{vrm_input} objects.
//' @param allstages_all The allstages indexing data frame used to produce MPMs.
//' @param dev_terms_list List of deviations for vital rate models.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates.
//' @param stageexpansion_ta_devterms_by_variant A list giving trait axis info
//' by variant, with each variant given a list element.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param density_vr_list Data frame of \code{lefkoDensVR} objects holding
//' density relationships for all 14 vital rate models.
//' @param current_stageframe The main stageframe, including extra stages.
//' @param current_supplement A supplement in \code{lefkoSD} format.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms.
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
//' @param sp_density_counter A matrix counting the use of spatial density terms
//' across variants and times.
//' @param inda_num_terms_previous A vector with numeric individual covariate a
//' terms for time t-1.
//' @param indb_num_terms_previous A vector with numeric individual covariate b
//' terms for time t-1.
//' @param indc_num_terms_previous A vector with numeric individual covariate c
//' terms for time t-1.
//' @param inda_cat_terms_previous A vector with categorical individual
//' covariate a terms for time t-1.
//' @param indb_cat_terms_previous A vector with categorical individual
//' covariate b terms for time t-1.
//' @param indc_cat_terms_previous A vector with categorical individual
//' covariate c terms for time t-1.
//' @param inda_num_terms_counter A vector with numeric individual covariate a
//' terms for time t.
//' @param indb_num_terms_counter A vector with numeric individual covariate b
//' terms for time t.
//' @param indc_num_terms_counter A vector with numeric individual covariate c
//' terms for time t.
//' @param inda_cat_terms_counter A vector with categorical individual covariate
//' a terms for time t.
//' @param indb_cat_terms_counter A vector with categorical individual covariate
//' b terms for time t.
//' @param indc_cat_terms_counter A vector with categorical individual covariate
//' c terms for time t.
//' @param dev_num_counter A vector giving the number of vital rate deviations
//' in each MPM.
//' @param fecmod_vec A numeric vector giving the fecmod values.
//' @param year_vec A vector giving the main years used.
//' @param patch_vec A vector giving the name of each patch used in projection.
//' @param var_per_run An integer giving the number of variants per run.
//' @param times An integer giving the number of occasions to project.
//' @param var_mat_length An integer giving the number of rows in the variant
//' matrix.
//' @param format_int An integer giving the MPM format.
//' @param current_rep The integer giving the current replicate.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param dev_terms_times_int A vector giving the number of occasions over
//' which vital rate y-intercept deviations cycle.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param year_counter An integer for year counts during projection.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
inline void invfb_project (const arma::mat var_run_mat, arma::vec& surv_dev_nta,
  arma::vec& obs_dev_nta, arma::vec& size_dev_nta, arma::vec& sizeb_dev_nta,
  arma::vec& sizec_dev_nta, arma::vec& repst_dev_nta, arma::vec& fec_dev_nta,
  arma::vec& jsurv_dev_nta, arma::vec& jobs_dev_nta, arma::vec& jsize_dev_nta,
  arma::vec& jsizeb_dev_nta, arma::vec& jsizec_dev_nta, arma::vec& jrepst_dev_nta,
  arma::vec& jmatst_dev_nta, arma::vec& indcova_nta, arma::vec& indcovb_nta,
  arma::vec& indcovc_nta, arma::ivec& variant_nta, List& N_out_pre,
  List& comm_out_pre, List& new_stageexpansion_list, List& errcheck_mpm_reps,
  List& errcheck_mpmout_reps, List& mdtl, List& used_times,
  const List allmodels_all, const List vrm_list, const List allstages_all,
  const List dev_terms_list, const List ind_terms_num_list,
  const List ind_terms_cat_list, List& stageexpansion_ta_devterms_by_variant,
  const List sp_density_list, const List start_list, const List equivalence_list,
  const DataFrame density_vr_list, DataFrame& current_stageframe,
  const DataFrame current_supplement, const DataFrame density_df,
  const DataFrame dens_index_df, const IntegerVector entry_time_vec,
  const IntegerVector sp_density_num_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_vr_yn_vec,
  IntegerMatrix& sp_density_counter, IntegerMatrix& inda_num_terms_previous,
  IntegerMatrix& indb_num_terms_previous, IntegerMatrix& indc_num_terms_previous,
  IntegerMatrix& inda_cat_terms_previous, IntegerMatrix& indb_cat_terms_previous,
  IntegerMatrix& indc_cat_terms_previous, IntegerMatrix& inda_num_terms_counter,
  IntegerMatrix& indb_num_terms_counter, IntegerMatrix& indc_num_terms_counter,
  IntegerMatrix& inda_cat_terms_counter, IntegerMatrix& indb_cat_terms_counter,
  IntegerMatrix& indc_cat_terms_counter, IntegerMatrix& dev_num_counter,
  const NumericVector fecmod_vec, const CharacterVector year_vec,
  const CharacterVector patch_vec, const int var_per_run, const int times,
  const int var_mat_length, const int format_int, const int current_rep,
  const int firstage_int, const int finalage_int, const int dev_terms_times_int,
  const int substoch, int& year_counter, const double exp_tol,
  const double theta_tol, const bool err_check, const bool err_check_extreme,
  const bool sparse_bool, const bool A_only, const bool stages_not_equal,
  const bool integeronly, const bool dens_yn_bool) {
  
  // Rcout << "invfb_project A" << endl;
  List running_popvecs; // = clone(start_list);
  List running_popvecs_startonly;
  arma::cube N_mpm (var_per_run, (times + 1), var_mat_length);
  
  List errcheck_mpm_reps_time (times); 
  List errcheck_mpmout_reps_time (times); 
  
  for (int j = 0; j < times; j++) { // 2nd loop - time j
    if (j % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    // Rcout << "invfb_project A1          ";
    List errcheck_mpm_reps_time_vmt (var_mat_length); 
    List errcheck_mpmout_reps_time_vmt (var_mat_length); 
    
    for (int i = 0; i < var_mat_length; i++) { // 3rd loop - permutes i
      List errcheck_mpm_reps_time_vmt_var(var_per_run); 
      List errcheck_mpmout_reps_time_vmt_var(var_per_run); 
      
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
        
        // Rcout << "invfb_project A2          ";
        List pop_reps = as<List>(all_pops_per_run(m));
        arma::mat pops_out = as<arma::mat>(pop_reps(current_rep));
        
        // Rcout << "invfb_project A3          ";
        if (j > (entry_time_vec(m) - 1)) {
          // Rcout << "invfb_project A4          ";
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
          // Rcout << "invfb_project A5          ";
          
          List current_vrm_extract = allmodels_all; // (i)
          List current_vrm_unextract = vrm_list; // (i)
          //int ehrlen_format {1}; // This will need to be dealt with differently later
          
          // Rcout << "invfb_project A6          ";
          DataFrame current_mpm_allstages = allstages_all; // (i)
          
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
          
          // Rcout << "invfb_project A7          ";
          CharacterVector current_mainyears = year_vec;
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
          
          // Rcout << "invfb_project A8          ";
          DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
          DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
          DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
          CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
          CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
          CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
          
          // Rcout << "invfb_project A9          ";
          // Counter resets
          int yearnumber = current_times_vec(j); // year_counter
          // Rcout << "invfb_project A10          ";
          
          // Rcout << "current_mainyears: " << current_mainyears << "               ";
          // Rcout << "yearnumber: " << yearnumber << "               ";
          CharacterVector current_year = as<CharacterVector>(current_mainyears(yearnumber));
          // Rcout << "invfb_project A11          ";
          
          if (inda_num_terms_counter(i, m) >= inda_terms_num_vec(0)) inda_num_terms_counter(i, m) = 0;
          if (indb_num_terms_counter(i, m) >= indb_terms_num_vec(0)) indb_num_terms_counter(i, m) = 0;
          if (indc_num_terms_counter(i, m) >= indc_terms_num_vec(0)) indc_num_terms_counter(i, m) = 0;
          if (inda_cat_terms_counter(i, m) >= inda_terms_cat_vec(0)) inda_cat_terms_counter(i, m) = 0;
          if (indb_cat_terms_counter(i, m) >= indb_terms_cat_vec(0)) indb_cat_terms_counter(i, m) = 0;
          if (indc_cat_terms_counter(i, m) >= indc_terms_cat_vec(0)) indc_cat_terms_counter(i, m) = 0;
          
          List current_ind_terms_num = ind_terms_num_list(0);
          List current_ind_terms_cat = ind_terms_cat_list(0);
          
          // Rcout << "invfb_project A12          ";
          NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
          NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
          NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
          CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
          CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
          CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
          
          // Rcout << "invfb_project A13          ";
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
          // Rcout << "invfb_project A14          ";
          NumericVector dv_terms (17);
          arma::uvec var_corresponding_elems;
          int vce_found {0};
          if (dev_terms_times_int > 0) {
            NumericMatrix used_dv_df = as<NumericMatrix>(dev_terms_list(current_variant_index));
            if (dev_num_counter(i, m) >= dev_terms_times_int) dev_num_counter(i, m) = 0;
            dv_terms = used_dv_df(_, dev_num_counter(i, m));
            
            var_corresponding_elems = find(variant_nta == (i + 1));
            vce_found = static_cast<int>(var_corresponding_elems.n_elem);
            
            // Rcout << "invfb_project A15          ";
            if (vce_found > 0) {
              // Rcout << "invfb_project A16          ";
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
              arma::vec indcova_nta_sub = indcova_nta.elem(var_corresponding_elems);
              arma::vec indcovb_nta_sub = indcovb_nta.elem(var_corresponding_elems);
              arma::vec indcovc_nta_sub = indcovc_nta.elem(var_corresponding_elems);
              
              // Rcout << "invfb_project A17          ";
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
                if(!NumericVector::is_na(indcova_nta_sub(vce_track))) dv_terms(14) =
                    dv_terms(14) + indcova_nta_sub(vce_track);
                if(!NumericVector::is_na(indcovb_nta_sub(vce_track))) dv_terms(15) =
                    dv_terms(15) + indcovb_nta_sub(vce_track);
                if(!NumericVector::is_na(indcovc_nta_sub(vce_track))) dv_terms(16) =
                    dv_terms(16) + indcovc_nta_sub(vce_track);
              }
            }
            dev_num_counter(i, m) = dev_num_counter(i, m) + 1;
          }
          
          // Rcout << "invfb_project A17a            ";
          NumericVector stdbv = as<NumericVector>(stageexpansion_ta_devterms_by_variant(current_variant_index));
          // Rcout << "invfb_project dv_terms: " << dv_terms << endl;
          // Rcout << "stdbv: " << stdbv << endl;
          for (int z = 0; z < 17; z++) {
            dv_terms(z) = dv_terms(z) + stdbv(z);
          }
          // Rcout << "invfb_project A18          ";
          
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
          
          // Rcout << "invfb_project A19          ";
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
          
          // Rcout << "invfb_project A20          ";
          double dens_sp {1.0};
          
          if (sp_density_num_vec(0) > 0) {
            if (sp_density_counter(i, m) >= sp_density_num_vec(0)) sp_density_counter(i, m) = 0;
            
            NumericVector current_sp_density = as<NumericVector>(sp_density_list(0));
            dens_sp = current_sp_density(sp_density_counter(i, m));
            
            sp_density_counter(i, m) = sp_density_counter(i, m) + 1;
          }
          
          List current_mpm;
          if (format_int < 5) {
            // Rcout << "invfb_project A21          ";
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
            
            // Rcout << "invfb_project A22          ";
          } else {
            IntegerVector all_ages = seq(firstage_int, finalage_int);
            //DataFrame current_supplement;
            if (!(current_supplement.length() > 1)) {
              // Rcout << "invfb_project A23          ";
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
                surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
                f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
                r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp,
                fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                exp_tol, theta_tol, sparse_bool);
              // Rcout << "invfb_project A24          ";
              
            } else {
              //current_supplement = as<DataFrame>(supplement_list(0));
              
              // Rcout << "invfb_project A25          ";
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
                surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
                f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
                r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(0),
                finalage_int, true, yearnumber, patchnumber, dvr_bool, dvr_yn,
                dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
                sparse_bool, current_supplement);
              // Rcout << "invfb_project A26          ";
            }
          }
          // Rcout << "invfb_project A27          ";
          
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
            // dens_bool = true
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
          
          // Rcout << "invfb_project A28          ";
          
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
}

//' Core Function-Based Projection Engine for ESS Evaluation
//' 
//' Function \code{invfb_optim} runs the optimization projections in function
//' \code{invade3_fb_core}. These are responsible for the estimation of ESS
//' values.
//' 
//' @name invfb_optim
//' 
//' @param surv_dev_nta The survival column in the reassessed trait axis.
//' @param obs_dev_nta The observation status column in the reassessed trait
//' axis.
//' @param size_dev_nta The primary size column in the reassessed trait axis.
//' @param sizeb_dev_nta The secondary size column in the reassessed trait axis.
//' @param sizec_dev_nta The tertiary size column in the reassessed trait axis.
//' @param repst_dev_nta The reproductive status column in the reassessed trait
//' axis.
//' @param fec_dev_nta The fecundity column in the reassessed trait axis.
//' @param jsurv_dev_nta The juvenile survival column in the reassessed trait
//' axis.
//' @param jobs_dev_nta The juvenile observation status column in the reassessed
//' trait axis.
//' @param jsize_dev_nta The juvenile primary size column in the reassessed
//' trait axis.
//' @param jsizeb_dev_nta The juvenile secondary size column in the reassessed
//' trait axis.
//' @param jsizec_dev_nta The juvenile tertiary size column in the reassessed
//' trait axis.
//' @param jrepst_dev_nta The juvenile reproductive status column in the
//' reassessed trait axis.
//' @param jmatst_dev_nta The juvenile maturity status column in the reassessed
//' trait axis.
//' @param indcova_nta The individual covariate a column in the reassessed
//' trait axis.
//' @param indcovb_nta The individual covariate b column in the reassessed
//' trait axis.
//' @param indcovc_nta The individual covariate c column in the reassessed
//' trait axis.
//' @param variant_nta The variant column in the reassessed 995 trait axis.
//' @param surv_dev_nta_995 The survival column in the reassessed 995 trait
//' axis.
//' @param obs_dev_nta_995 The observation status column in the reassessed 995
//' trait axis.
//' @param size_dev_nta_995 The primary size column in the reassessed 995 trait
//' axis.
//' @param sizeb_dev_nta_995 The secondary size column in the reassessed 995
//' trait axis.
//' @param sizec_dev_nta_995 The tertiary size column in the reassessed 995
//' trait axis.
//' @param repst_dev_nta_995 The reproductive status column in the reassessed
//' 995 trait axis.
//' @param fec_dev_nta_995 The fecundity column in the reassessed 995 trait
//' axis.
//' @param jsurv_dev_nta_995 The juvenile survival column in the reassessed 995
//' trait axis.
//' @param jobs_dev_nta_995 The juvenile observation status column in the
//' reassessed 995 trait axis.
//' @param jsize_dev_nta_995 The juvenile primary size column in the reassessed
//' 995 trait axis.
//' @param jsizeb_dev_nta_995 The juvenile secondary size column in the
//' reassessed 995 trait axis.
//' @param jsizec_dev_nta_995 The juvenile tertiary size column in the
//' reassessed 995 trait axis.
//' @param jrepst_dev_nta_995 The juvenile reproductive status column in the
//' reassessed 995 trait axis.
//' @param jmatst_dev_nta_995 The juvenile maturity status column in the
//' reassessed 995 trait axis.
//' @param indcova_nta_995 The individual covariate a column in the reassessed
//' 995 trait axis.
//' @param indcovb_nta_995 The individual covariate b column in the reassessed
//' 995 trait axis.
//' @param indcovc_nta_995 The individual covariate c column in the reassessed
//' 995 trait axis.
//' @param variant_nta_995 The variant column in the reassessed 995 trait axis.
//' @param N_out_pre The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out_pre The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param new_stageexpansion_list A list with stage expansions for all
//' variant data used in ESS evaluation. This list includes an extra layer of
//' list elements, corresponding to the optim_ta and optim_ta_995 data.
//' @param errcheck_mpm_reps An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param errcheck_mpmout_reps An optional list of all mpm_out matrices from MPM
//' processing. Only output if \code{err_check = "extreme"}.
//' @param mdtl The modified dev terms list.
//' @param used_times A list of year numbers for each time per run.
//' @param allmodels_all A list of extracted vrm inputs for all MPMs.
//' @param vrm_list A list of \code{vrm_input} objects.
//' @param allstages_all The allstages indexing data frame used to produce MPMs.
//' @param dev_terms_list List of deviations for vital rate models.
//' @param ind_terms_num_list List of data frames giving values of numeric
//' individual covariates.
//' @param ind_terms_cat_list List of data frames giving values of factor
//' individual covariates.
//' @param stageexpansion_ta_devterms_by_variant A list giving trait axis info
//' by variant, with each variant given a list element.
//' @param sp_density_list A list of values of spatial density for all MPMs.
//' @param start_list A list of starting information, supplied in \code{lefkoSV}
//' format.
//' @param equivalence_list A list giving the effect of each individual in each
//' stage relative to a reference individual.
//' @param density_vr_list Data frame of \code{lefkoDensVR} objects holding
//' density relationships for all 14 vital rate models.
//' @param current_stageframe The main stageframe, including extra stages.
//' @param current_supplement A supplement in \code{lefkoSD} format.
//' @param density_df A data frame of class \code{lefkoDens}.
//' @param dens_index_df A data frame giving indices for density dependent
//' transitions.
//' @param entry_time_vec An IntegerVector containing the entry time of each
//' mutant, population, or species, as given by each MPM.
//' @param sp_density_num_vec A vector giving the number of spatial density
//' terms.
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
//' @param fecmod_vec A numeric vector giving the fecmod values.
//' @param year_vec A vector giving the main years used.
//' @param patch_vec A vector giving the name of each patch used in projection.
//' @param var_per_run An integer giving the number of variants per run.
//' @param times An integer giving the number of occasions to project.
//' @param var_mat_length An integer giving the number of rows in the variant
//' matrix.
//' @param format_int An integer giving the MPM format.
//' @param current_rep The integer giving the current replicate.
//' @param firstage_int An integer giving the first age in a Leslie or
//' age-by-stage MPM.
//' @param finalage_int  An integer giving the final age in a Leslie or
//' age-by-stage MPM.
//' @param dev_terms_times_int A vector giving the number of occasions over
//' which vital rate y-intercept deviations cycle.
//' @param substoch An integer giving the level of sustochasticity to enforce.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis. Note that the version used in this
//' function is actually equivalent to \code{opt_res_true}.
//' @param opt_res_orig The original value of \code{opt_res}, prior to the
//' determination of the number of variable traits. Equal to \code{opt_res} if
//' the number of variable traits is 1, and to the square root of \code{opt_res}
//' if the number of variable traits is 2.
//' @param year_counter An integer for year counts during projection.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param threshold The lower limit for the absolute value of fitness, below
//' which fitness is rounded to 0. Defaults to 0.00000001.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking.
//' @param err_check_extreme A logical value indicating whether to include an
//' extra list of all matrices projected in the \code{err_check} object.
//' @param sparse_bool A Boolean value indiating whether the MPM is in sparse
//' matrix format.
//' @param A_only A Boolean value indicating whether to export U and F matrices
//' for alteration, or only A matrices.
//' @param stages_not_equal A Boolean value indicating whether equivalence
//' info is supplied suggesting even stages within MPMs are not equal.
//' @param integeronly A Boolean value indicating whether to allow only whole
//' values of individuals or not.
//' @param dens_yn_bool A Boolean value stating whether density dependence is
//' used, given through \code{lefkoDens} objects.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
inline void invfb_optim (arma::vec& surv_dev_nta,
  arma::vec& obs_dev_nta, arma::vec& size_dev_nta, arma::vec& sizeb_dev_nta,
  arma::vec& sizec_dev_nta, arma::vec& repst_dev_nta, arma::vec& fec_dev_nta,
  arma::vec& jsurv_dev_nta, arma::vec& jobs_dev_nta, arma::vec& jsize_dev_nta,
  arma::vec& jsizeb_dev_nta, arma::vec& jsizec_dev_nta,
  arma::vec& jrepst_dev_nta, arma::vec& jmatst_dev_nta, arma::vec& indcova_nta,
  arma::vec& indcovb_nta, arma::vec& indcovc_nta, arma::ivec& variant_nta,
  arma::vec& surv_dev_nta_995, arma::vec& obs_dev_nta_995,
  arma::vec& size_dev_nta_995, arma::vec& sizeb_dev_nta_995,
  arma::vec& sizec_dev_nta_995, arma::vec& repst_dev_nta_995,
  arma::vec& fec_dev_nta_995, arma::vec& jsurv_dev_nta_995,
  arma::vec& jobs_dev_nta_995, arma::vec& jsize_dev_nta_995,
  arma::vec& jsizeb_dev_nta_995, arma::vec& jsizec_dev_nta_995,
  arma::vec& jrepst_dev_nta_995, arma::vec& jmatst_dev_nta_995,
  arma::vec& indcova_nta_995, arma::vec& indcovb_nta_995,
  arma::vec& indcovc_nta_995, arma::ivec& variant_nta_995, List& N_out_pre,
  List& comm_out_pre, List& new_stageexpansion_list, List& errcheck_mpm_reps,
  List& errcheck_mpmout_reps, List& mdtl, List& used_times,
  const List allmodels_all, const List vrm_list, const List allstages_all,
  const List dev_terms_list, const List ind_terms_num_list,
  const List ind_terms_cat_list, List& stageexpansion_ta_devterms_by_variant,
  const List sp_density_list, const List start_list, const List equivalence_list,
  const DataFrame density_vr_list, DataFrame& current_stageframe,
  const DataFrame current_supplement, const DataFrame density_df,
  const DataFrame dens_index_df, const IntegerVector entry_time_vec,
  const IntegerVector sp_density_num_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_vr_yn_vec,
  const NumericVector fecmod_vec, const CharacterVector year_vec,
  const CharacterVector patch_vec, const int var_per_run, const int times,
  const int var_mat_length, const int format_int, const int current_rep,
  const int firstage_int, const int finalage_int, const int dev_terms_times_int,
  const int substoch, const int opt_res, const int opt_res_orig,
  int& year_counter, const double exp_tol, const double theta_tol,
  const double threshold, const bool err_check, const bool err_check_extreme,
  const bool sparse_bool, const bool A_only, const bool stages_not_equal,
  const bool integeronly, const bool dens_yn_bool) {
  
  // Rcout << "invfb_optim A" << endl;
  List running_popvecs; // = clone(start_list);
  List running_popvecs_startonly;
  arma::cube N_mpm (2, (times + 1), opt_res);
  
  List errcheck_mpm_reps_time (times); 
  List errcheck_mpmout_reps_time (times); 
  
  IntegerMatrix inda_num_terms_counter (opt_res, var_per_run);
  IntegerMatrix indb_num_terms_counter (opt_res, var_per_run);
  IntegerMatrix indc_num_terms_counter (opt_res, var_per_run);
  IntegerMatrix inda_cat_terms_counter (opt_res, var_per_run);
  IntegerMatrix indb_cat_terms_counter (opt_res, var_per_run);
  IntegerMatrix indc_cat_terms_counter (opt_res, var_per_run);
  IntegerMatrix inda_num_terms_previous (opt_res, var_per_run);
  IntegerMatrix indb_num_terms_previous (opt_res, var_per_run);
  IntegerMatrix indc_num_terms_previous (opt_res, var_per_run);
  IntegerMatrix inda_cat_terms_previous (opt_res, var_per_run);
  IntegerMatrix indb_cat_terms_previous (opt_res, var_per_run);
  IntegerMatrix indc_cat_terms_previous (opt_res, var_per_run);
  IntegerMatrix dev_num_counter (opt_res, var_per_run);
  IntegerMatrix sp_density_counter (opt_res, var_per_run);
  
  for (int j = 0; j < times; j++) { // 2nd loop - time j
    // Rcout << "invfb_optim A1          ";
    List errcheck_mpm_reps_time_vmt (opt_res); 
    List errcheck_mpmout_reps_time_vmt (opt_res); 
    
    for (int i = 0; i < opt_res; i++) { // 3rd loop - permutes i
      if (i % 10 == 0){
        Rcpp::checkUserInterrupt();
      }
      
      List errcheck_mpm_reps_time_vmt_var(2); 
      List errcheck_mpmout_reps_time_vmt_var(2); 
      
      if (j == 0) {
        List var_popvecs_to_start (2);
        for (int n = 0; n < 2; n++) {
          var_popvecs_to_start(n) = as<arma::vec>(start_list(0));
        }
        running_popvecs = var_popvecs_to_start;
        running_popvecs_startonly = clone(var_popvecs_to_start);
      }
      
      for (int m = 0; m < 2; m++) { // 4th loop - var per run m
        if (j == entry_time_vec(m)) {
          arma::vec running_popvec_mpm = as<arma::vec>(running_popvecs_startonly(m));
          
          double N_current = accu(running_popvec_mpm);
          N_mpm(m, j, i) = N_current; // Used to be (k, j)
        }
      }
      
      List all_pops_per_run = as<List>(comm_out_pre(i));
      List mdtl_1 (var_per_run);
     
      { // Main variant section
        // Rcout << "invfb_optim A2          ";
        int m = 0;
        int current_variant_index = i; // Equivalent to index integer k
        
        List sgec = as<List> (new_stageexpansion_list(current_variant_index));
        DataFrame sge_current = as<DataFrame>(sgec(m));
        
        // Rcout << "invfb_optim A3          ";
        List pop_reps = as<List>(all_pops_per_run(m));
        arma::mat pops_out = as<arma::mat>(pop_reps(current_rep));
        
        // Rcout << "invfb_optim A4          ";
        if (j > (entry_time_vec(m) - 1)) {
          // Rcout << "invfb_optim A5          ";
          List used_times_per_run = as<List>(used_times(0));
          List used_times_current_var = as<List>(used_times_per_run(m));
          IntegerVector current_times_vec = as<IntegerVector>(used_times_current_var(current_rep));
          
          arma::vec running_popvec_vrm;
          if (j == entry_time_vec(m)) {
            running_popvec_vrm = as<arma::vec>(running_popvecs_startonly(m));
            pops_out.col(j) = running_popvec_vrm;
          } else {
            running_popvec_vrm = pops_out.col(j);
          }
          // Rcout << "invfb_optim A6 current_rep: " << current_rep << "      ";
            
          List current_vrm_extract = allmodels_all; // (i)
          List current_vrm_unextract = vrm_list; // (i)
          //int ehrlen_format {1}; // This will need to be dealt with differently later
          
          // Rcout << "invfb_optim A7       ";
          //int mpm_style {1};
          //if (format_int < 3) {
          //  mpm_style = 0;
          //  if (format_int == 2) ehrlen_format = 2;
          //} else if (format_int == 4) {
          //  mpm_style = 2;
          //}
          
          // Rcout << "invfb_optim A8       ";
          DataFrame current_mpm_allstages = allstages_all; // (i)
          
          // Rcout << "invfb_optim A9       ";
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
          
          // Rcout << "invfb_optim A10       ";
          CharacterVector current_mainyears = year_vec;
          //unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
          
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
          
          // Rcout << "invfb_optim A11        ";
          
          // Not sure if we need the next bit
          DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
          DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
          DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
          CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
          CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
          CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
          
          // Rcout << "invfb_optim A12        ";
          // Counter resets
          int yearnumber = current_times_vec(j); // year_counter
          CharacterVector current_year = as<CharacterVector>(current_mainyears(yearnumber));
          // Rcout << "invfb_optim A14        ";
          
          if (inda_num_terms_counter(i, m) >= inda_terms_num_vec(0)) inda_num_terms_counter(i, m) = 0;
          if (indb_num_terms_counter(i, m) >= indb_terms_num_vec(0)) indb_num_terms_counter(i, m) = 0;
          if (indc_num_terms_counter(i, m) >= indc_terms_num_vec(0)) indc_num_terms_counter(i, m) = 0;
          if (inda_cat_terms_counter(i, m) >= inda_terms_cat_vec(0)) inda_cat_terms_counter(i, m) = 0;
          if (indb_cat_terms_counter(i, m) >= indb_terms_cat_vec(0)) indb_cat_terms_counter(i, m) = 0;
          if (indc_cat_terms_counter(i, m) >= indc_terms_cat_vec(0)) indc_cat_terms_counter(i, m) = 0;
          
          List current_ind_terms_num = ind_terms_num_list(0);
          List current_ind_terms_cat = ind_terms_cat_list(0);
          
          // Rcout << "invfb_optim A15        ";
          NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
          NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
          NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
          CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
          CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
          CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
          
          // Rcout << "invfb_optim A16        ";
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
          
          // Rcout << "invfb_optim A17        ";
          // dev_terms and vrm trait axis processing
          NumericVector dv_terms (17);
          arma::uvec var_corresponding_elems;
          int vce_found {0};
          if (dev_terms_times_int > 0) {
            // Rcout << "invfb_optim A18        ";
            NumericMatrix used_dv_df = clone(as<NumericMatrix>(dev_terms_list(current_variant_index)));
            if (dev_num_counter(i, 0) >= dev_terms_times_int) dev_num_counter(i, 0) = 0;
            dv_terms = used_dv_df(_, dev_num_counter(i, 0));
          }
          
          var_corresponding_elems = find(variant_nta == (i + 1));
          vce_found = static_cast<int>(var_corresponding_elems.n_elem);
          
          // Rcout << "invfb_optim A19        ";
          if (vce_found > 0) {
            // Rcout << "invfb_optim A20        ";
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
            arma::vec indcova_nta_sub = indcova_nta.elem(var_corresponding_elems);
            arma::vec indcovb_nta_sub = indcovb_nta.elem(var_corresponding_elems);
            arma::vec indcovc_nta_sub = indcovc_nta.elem(var_corresponding_elems);
            
            // Rcout << "invfb_optim A21        ";
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
              if(!NumericVector::is_na(indcova_nta_sub(vce_track))) dv_terms(14) =
                  dv_terms(14) + indcova_nta_sub(vce_track);
              if(!NumericVector::is_na(indcovb_nta_sub(vce_track))) dv_terms(15) =
                  dv_terms(15) + indcovb_nta_sub(vce_track);
              if(!NumericVector::is_na(indcovc_nta_sub(vce_track))) dv_terms(16) =
                  dv_terms(16) + indcovc_nta_sub(vce_track);
            }
          }
          dev_num_counter(i, 0) = dev_num_counter(i, 0) + 1;
          
          // Rcout << "invfb_optim A22        ";
          List stdbv1 = as<List>(stageexpansion_ta_devterms_by_variant(current_variant_index));
          NumericVector stdbv = as<NumericVector>(stdbv1(0));
          // Rcout << "invfb_optim A23        ";
          
          for (int z = 0; z < 17; z++) {
            dv_terms(z) = dv_terms(z) + stdbv(z);
          }
          // Rcout << "invfb_optim m = 0 dv_terms: " << dv_terms << endl;
          // Rcout << "invfb_optim A24        ";
          
          if (err_check) {
            mdtl(i) = mdtl_1;
          }
          bool dvr_bool {false};
          
          // Rcout << "invfb_optim A25        ";
          LogicalVector dvr_yn = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          IntegerVector dvr_style (14);
          IntegerVector dvr_time_delay (14);
          NumericVector dvr_alpha (14);
          NumericVector dvr_beta (14);
          NumericVector dens_n (14);
          
          // Rcout << "invfb_optim A26        ";
          
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
          
          // Rcout << "invfb_optim A27        ";
          double dens_sp {1.0};
          
          if (sp_density_num_vec(0) > 0) {
            if (sp_density_counter(i, m) >= sp_density_num_vec(0)) sp_density_counter(i, m) = 0;
            
            NumericVector current_sp_density = as<NumericVector>(sp_density_list(0));
            dens_sp = current_sp_density(sp_density_counter(i, m));
            
            sp_density_counter(i, m) = sp_density_counter(i, m) + 1;
          }
          
          List current_mpm;
          if (format_int < 5) {
            // Rcout << "invfb_optim A28        ";
            // Rcout << "dv_terms: " << dv_terms << endl;
            
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
            
          // Rcout << "invfb_optim A29        ";
          } else {
            IntegerVector all_ages = seq(firstage_int, finalage_int);
            //DataFrame current_supplement;
            if (!(current_supplement.length() > 1)) {
              // Rcout << "inv_fb_optim g3        ";
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
                surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
                f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
                r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(0),
                finalage_int, true, yearnumber, patchnumber, dvr_bool, dvr_yn,
                dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
                sparse_bool);
              // Rcout << "invfb_optim A30        ";
              
            } else {
              //current_supplement = as<DataFrame>(supplement_list(0));
              // Rcout << "invfb_optim A31        ";
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages, current_stageframe,
                surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
                f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
                r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp, fecmod_vec(0),
                finalage_int, true, yearnumber, patchnumber, dvr_bool, dvr_yn,
                dvr_style, dvr_alpha, dvr_beta, dens_n, exp_tol, theta_tol,
                sparse_bool, current_supplement);
              // Rcout << "invfb_optim A32        ";
            }
          }
          // Rcout << "invfb_optim A33        ";
          
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
            // dens_bool = true
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
                  for (int l = 0; l < 2; l++) {
                    int current_variant_index_agg = i;
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
          
          // Rcout << "invfb_optim A34        ";
          
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
        
        // Rcout << "invfb_optim A100           ";
        
        pop_reps(current_rep) = pops_out;
        // Rcout << "invfb_optim A101           ";
      } // main variant section
      
      // Rcout << "invfb_optim B1          ";
      { // e995 variant section
        int m = 1;
        int current_variant_index = i; // Equivalent to index integer k
        
        // Rcout << "invfb_optim B2          ";
        List sgec = as<List> (new_stageexpansion_list(current_variant_index));
        DataFrame sge_current = as<DataFrame>(sgec(m));
        
        List pop_reps = as<List>(all_pops_per_run(m));
        arma::mat pops_out = as<arma::mat>(pop_reps(current_rep));
        
        // Rcout << "invfb_optim B3          ";
        if (j > (entry_time_vec(m) - 1)) {
          // Rcout << "invfb_optim B4          ";
          List used_times_per_run = as<List>(used_times(0));
          List used_times_current_var = as<List>(used_times_per_run(m));
          IntegerVector current_times_vec = as<IntegerVector>(used_times_current_var(current_rep));
            
          // Rcout << "invfb_optim B5          ";
          
          arma::vec running_popvec_vrm;
          if (j == entry_time_vec(m)) {
            running_popvec_vrm = as<arma::vec>(running_popvecs_startonly(m));
            pops_out.col(j) = running_popvec_vrm;
          } else {
            running_popvec_vrm = pops_out.col(j);
          }
          // Rcout << "invfb_optim B6 current_rep: " << current_rep << "      ";
          
          List current_vrm_extract = allmodels_all; // (i)
          List current_vrm_unextract = vrm_list; // (i)
          //int ehrlen_format {1}; // This will need to be dealt with differently later
          
          // Rcout << "invfb_optim B7       ";
          //int mpm_style {1};
          //if (format_int < 3) {
          //  mpm_style = 0;
          //  if (format_int == 2) ehrlen_format = 2;
          //} else if (format_int == 4) {
          //  mpm_style = 2;
          //}
          
          // Rcout << "invfb_optim B8       ";
          DataFrame current_mpm_allstages = allstages_all; // (i)
          
          // Rcout << "invfb_optim B9       ";
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
          
          // Rcout << "invfb_optim B10       ";
          CharacterVector current_mainyears = year_vec;
          //unsigned int no_mainyears = static_cast<unsigned int>(current_mainyears.length());
          
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
          
          // Rcout << "invfb_optim B11        ";
          
          // Not sure if we need the next bit
          DataFrame indcova2_frame = as<DataFrame>(current_vrm_unextract["indcova2_frame"]);
          DataFrame indcovb2_frame = as<DataFrame>(current_vrm_unextract["indcovb2_frame"]);
          DataFrame indcovc2_frame = as<DataFrame>(current_vrm_unextract["indcovc2_frame"]);
          CharacterVector current_mainindcova = as<CharacterVector>(indcova2_frame["indcova"]);
          CharacterVector current_mainindcovb = as<CharacterVector>(indcovb2_frame["indcovb"]);
          CharacterVector current_mainindcovc = as<CharacterVector>(indcovc2_frame["indcovc"]);
          
          // Rcout << "invfb_optim B12        ";
          // Counter resets
          int yearnumber = current_times_vec(j); // year_counter
          CharacterVector current_year = as<CharacterVector>(current_mainyears(yearnumber));
          // Rcout << "invfb_optim B14        ";
          
          if (inda_num_terms_counter(i, 0) >= inda_terms_num_vec(0)) inda_num_terms_counter(i, 0) = 0;
          if (indb_num_terms_counter(i, 0) >= indb_terms_num_vec(0)) indb_num_terms_counter(i, 0) = 0;
          if (indc_num_terms_counter(i, 0) >= indc_terms_num_vec(0)) indc_num_terms_counter(i, 0) = 0;
          if (inda_cat_terms_counter(i, 0) >= inda_terms_cat_vec(0)) inda_cat_terms_counter(i, 0) = 0;
          if (indb_cat_terms_counter(i, 0) >= indb_terms_cat_vec(0)) indb_cat_terms_counter(i, 0) = 0;
          if (indc_cat_terms_counter(i, 0) >= indc_terms_cat_vec(0)) indc_cat_terms_counter(i, 0) = 0;
          
          List current_ind_terms_num = ind_terms_num_list(0);
          List current_ind_terms_cat = ind_terms_cat_list(0);
          
          // Rcout << "invfb_optim B15        ";
          NumericVector f_inda_full = as<NumericVector>(current_ind_terms_num(0));
          NumericVector f_indb_full = as<NumericVector>(current_ind_terms_num(1));
          NumericVector f_indc_full = as<NumericVector>(current_ind_terms_num(2));
          CharacterVector r_inda_full = as<CharacterVector>(current_ind_terms_cat(0));
          CharacterVector r_indb_full = as<CharacterVector>(current_ind_terms_cat(1));
          CharacterVector r_indc_full = as<CharacterVector>(current_ind_terms_cat(2));
          
          // Rcout << "invfb_optim B16        ";
          NumericVector f2_inda = {f_inda_full(inda_num_terms_counter(i, 0))}; // i
          NumericVector f1_inda = {f_inda_full(inda_num_terms_previous(i, 0))};
          NumericVector f2_indb = {f_indb_full(indb_num_terms_counter(i, 0))};
          NumericVector f1_indb = {f_indb_full(indb_num_terms_previous(i, 0))};
          NumericVector f2_indc = {f_indc_full(indc_num_terms_counter(i, 0))};
          NumericVector f1_indc = {f_indc_full(indc_num_terms_previous(i, 0))};
          CharacterVector r2_inda = as<CharacterVector>(r_inda_full(inda_cat_terms_counter(i, 0)));
          CharacterVector r1_inda = 
            as<CharacterVector>(r_inda_full(inda_cat_terms_previous(i, 0)));
          CharacterVector r2_indb = as<CharacterVector>
            (r_indb_full(indb_cat_terms_counter(i, 0)));
          CharacterVector r1_indb = as<CharacterVector>
            (r_indb_full(indb_cat_terms_previous(i, 0)));
          CharacterVector r2_indc = as<CharacterVector>
            (r_indc_full(indc_cat_terms_counter(i, 0)));
          CharacterVector r1_indc = 
            as<CharacterVector>(r_indc_full(indc_cat_terms_previous(i, 0)));
          
          // Rcout << "invfb_optim B17        ";
          // dev_terms and vrm trait axis processing
          NumericVector dv_terms (17);
          arma::uvec var_corresponding_elems;
          int vce_found {0};
          if (dev_terms_times_int > 0) {
            // Rcout << "invfb_optim B18        ";
            NumericMatrix used_dv_df = clone(as<NumericMatrix>(dev_terms_list(current_variant_index)));
            if (dev_num_counter(i, 0) >= dev_terms_times_int) dev_num_counter(i, 0) = 0;
            dv_terms = used_dv_df(_, dev_num_counter(i, 0));
          }
          
          var_corresponding_elems = find(variant_nta == (i + 1));
          vce_found = static_cast<int>(var_corresponding_elems.n_elem);
          
          // Rcout << "invfb_optim B19        ";
          if (vce_found > 0) {
            // Rcout << "invfb_optim B20        ";
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
            arma::vec indcova_nta_sub = indcova_nta.elem(var_corresponding_elems);
            arma::vec indcovb_nta_sub = indcovb_nta.elem(var_corresponding_elems);
            arma::vec indcovc_nta_sub = indcovc_nta.elem(var_corresponding_elems);
            
            // Rcout << "invfb_optim B21        ";
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
              if(!NumericVector::is_na(indcova_nta_sub(vce_track))) dv_terms(14) =
                  dv_terms(14) + indcova_nta_sub(vce_track);
              if(!NumericVector::is_na(indcovb_nta_sub(vce_track))) dv_terms(15) =
                  dv_terms(15) + indcovb_nta_sub(vce_track);
              if(!NumericVector::is_na(indcovc_nta_sub(vce_track))) dv_terms(16) =
                  dv_terms(16) + indcovc_nta_sub(vce_track);
            }
          }
          dev_num_counter(i, 0) = dev_num_counter(i, 0) + 1;
          
          // Rcout << "invfb_optim B22        ";
          List stdbv2 = as<List>(stageexpansion_ta_devterms_by_variant(current_variant_index));
          NumericVector stdbv = as<NumericVector>(stdbv2(1));
          // Rcout << "invfb_optim B23        ";
          
          for (int z = 0; z < 17; z++) {
            dv_terms(z) = dv_terms(z) + stdbv(z);
          }
          // Rcout << "invfb_optim B24        ";
          
          if (err_check) {
            mdtl(i) = mdtl_1;
          }
          bool dvr_bool {false};
          
          // Rcout << "invfb_optim B25        ";
          LogicalVector dvr_yn = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          IntegerVector dvr_style (14);
          IntegerVector dvr_time_delay (14);
          NumericVector dvr_alpha (14);
          NumericVector dvr_beta (14);
          NumericVector dens_n (14);
          
          // Rcout << "invfb_optim B26        ";
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
          
          // Rcout << "invfb_optim B27        ";
          
          double dens_sp {1.0};
          
          if (sp_density_num_vec(0) > 0) {
            if (sp_density_counter(i, m) >= sp_density_num_vec(0)) sp_density_counter(i, m) = 0;
            
            NumericVector current_sp_density = as<NumericVector>(sp_density_list(0));
            dens_sp = current_sp_density(sp_density_counter(i, m));
            
            sp_density_counter(i, m) = sp_density_counter(i, m) + 1;
          }
          
          List current_mpm;
          if (format_int < 5) {
            // Rcout << "invfb_optim B28        ";
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
            
            // Rcout << "invfb_optim B29        ";
          } else {
            IntegerVector all_ages = seq(firstage_int, finalage_int);
            if (!(current_supplement.length() > 1)) {
              // Rcout << "invfb_optim B30        ";
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages,
                current_stageframe, surv_proxy, fec_proxy, f2_inda, f1_inda,
                f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                r1_indb, r2_indc, r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp,
                fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                exp_tol, theta_tol, sparse_bool);
              // Rcout << "invfb_optim B31        ";
              
            } else {
              // Rcout << "invfb_optim B32        ";
              NumericVector dvt3 = {dv_terms(14), dv_terms(15), dv_terms(16)};
              current_mpm = AdaptMats::mdabrowskiego(all_ages,
                current_stageframe, surv_proxy, fec_proxy, f2_inda, f1_inda,
                f2_indb, f1_indb, f2_indc, f1_indc, r2_inda, r1_inda, r2_indb,
                r1_indb, r2_indc, r1_indc, dvt3, dv_terms(0), dv_terms(6), dens_sp,
                fecmod_vec(0), finalage_int, true, yearnumber, patchnumber,
                dvr_bool, dvr_yn, dvr_style, dvr_alpha, dvr_beta, dens_n,
                exp_tol, theta_tol, sparse_bool, current_supplement);
              // Rcout << "invfb_optim B33        ";
            }
          }
          // Rcout << "invfb_optim B34        ";
          
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
            // dens_bool = true
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
                  for (int l = 0; l < 2; l++) {
                    int current_variant_index_agg = i;
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
          
          // Rcout << "invfb_optim B33        ";
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
        // Rcout << "invfb_optim B100              ";
        
        pop_reps(current_rep) = pops_out;
        // Rcout << "invfb_optim B101              ";
      }
        
      // Rcout << "invfb_optim C1              ";
      if (err_check_extreme) {
        errcheck_mpmout_reps_time_vmt(i) = errcheck_mpmout_reps_time_vmt_var;
        errcheck_mpm_reps_time_vmt(i) = errcheck_mpm_reps_time_vmt_var;
      }
      // Rcout << "invfb_optim C2              ";
    }
    
    // Rcout << "invfb_optim C3              ";
    if (err_check_extreme) {
      errcheck_mpmout_reps_time(j) = errcheck_mpmout_reps_time_vmt;
      errcheck_mpm_reps_time(j) = errcheck_mpm_reps_time_vmt;
    }
    // Rcout << "invfb_optim C4              ";
    
    year_counter++;
  } // j loop
  // Rcout << "invfb_optim C5              ";
  if (err_check_extreme) {
    errcheck_mpmout_reps(current_rep) = errcheck_mpmout_reps_time;
    errcheck_mpm_reps(current_rep) = errcheck_mpm_reps_time;
  }
  // Rcout << "invfb_optim C6              ";
  N_out_pre(current_rep) = N_mpm;
}

//' Set-Up Function Running Invasion Analyses of Function-based MPMs
//' 
//' Function \code{invade3_fb_core} is the main function running invasion
//' analyses in which matrices must be created at each time step.
//' 
//' @name invade3_fb_core
//' 
//' @param Lyapunov The main data frame giving the Lyapunov coefficients
//' estimated, as well as the circumstances resulting in them.
//' @param Lyapunov_optim Main data frame giving Lyapunov coefficients for all
//' trait combinations developed for the ESS optima table. Holds elasticity
//' fitness values.
//' @param ESS_Lyapunov A data frame provided by reference that will hold the
//' ESS optima.
//' @param var_run_mat A matrix giving the the variants to be run in each
//' projection, with rows giving the projections and columns giving the
//' variants.
//' @param N_out The main list of final population sizes, supplied as a
//' reference and altered by this function.
//' @param comm_out The main list of full projection results for the community,
//' supplied as a pointer and altered by this function.
//' @param N_out_optim The main list of final population sizes from ESS
//' optimization, supplied as a reference and altered by this function.
//' @param comm_out_optim The main list of full projection results for the
//' community resulting from ESS optimization, supplied as a pointer and altered
//' by this function.
//' @param zero_stage_vec_list A list of vectors giving zero stage vectors for
//' each MPM, if entry times are staggered.
//' @param trait_axis A data frame of class \code{adaptAxis} holding the trait
//' data to test.
//' @param new_trait_axis A data frame giving trait axis data post-processing
//' with function \code{ta_reassess()}.
//' @param optim_trait_axis A data frame giving trait axis data processed for
//' ESS optimization.
//' @param optim_trait_axis_995 A data frame giving trait axis data processed
//' for ESS optimization for variants 99.5% the values of the core frame.
//' @param new_stageexpansion_list A list with stage expansions for all trait
//' axis data leading to matrix element changes with each list element
//' corresponding to each respective variant.
//' @param new_stageexpansion_list_optim A list with stage expansions for all
//' variant data used in ESS evaluation.
//' @param modified_dev_terms_list An optional list giving the vital rate
//' y-intercept deviations by variant once data from the \code{trait_axis} data
//' frame has been allocated.
//' @param errcheck_mpms An optional list of all MPMs post-processing. Only
//' output if \code{err_check = "extreme"}.
//' @param errcheck_mpms_optima An optional list of all MPMs used in the first
//' stage of ESS trait evaluation. Only output if \code{err_check = "extreme"}.
//' @param errcheck_mpm_ESS An optional list of all MPMs used in the second
//' stage of ESS trait evaluation. Only output if \code{err_check = "extreme"}.
//' @param errcheck_mpmouts An optional list of all mpm_out matrices from MPM
//' processing. Only output if \code{err_check = "extreme"}.
//' @param errcheck_mpmouts_optima An optional list of all mpm_out matrices from
//' the first stage of ESS trait evaluation. Only output if
//' \code{err_check = "extreme"}.
//' @param errcheck_mpmout_ESS An optional list of all mpm_out matrices from the
//' second stage of ESS trait evaluation. Only output if
//' \code{err_check = "extreme"}.
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
//' @param dev_terms_times_int A vector giving the number of occasions over
//' which vital rate y-intercept deviations cycle.
//' @param main_optim_res An integer giving the number of variants being tested.
//' @param opt_res If evaluating optima, then this integer gives the number
//' of variants to create between each minimum and maximum for each trait found
//' to be variable in the input trait axis.
//' @param exp_tol The maximum tolerated exponent.
//' @param theta_tol The maximum tolerated limit for theta, in non-linear
//' models such as those using the negative binomial.
//' @param loop_max An integer value denoting the number of search cycles
//' allowed per ESS during ESS optimization. Defaults to 150.
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
//' @param threshold The lower limit for the absolute value of fitness, below
//' which fitness is rounded to 0. Defaults to 0.00000001.
//' @param fitness_table A Boolean value dictating whether to include a data
//' frame giving Lyapunov coefficients for all combinations of variants tested.
//' Necessary for the creation of pairwise invasibility plots (PIPs). Defaults
//' to \code{TRUE}.
//' @param ESS_optima A logical value indicating whether to assess the values of
//' ESSoptima for traits that vary among variants in the given \code{adaptAxis}
//' table. Defaults to \code{TRUE}.
//' @param elast_mult A multiplier for traits to assess the elasticity of
//' fitness in trait optimization. Defaults to 0.995.
//' @param zap_min A Boolean value describing whether to round fitness values
//' below the value given in \code{threshold}.
//' 
//' @return The first four arguments are directly manipulated without any
//' values returned.
//' 
//' @keywords internal
//' @noRd
void invade3_fb_core (DataFrame& Lyapunov, DataFrame& Lyapunov_optim,
  DataFrame& ESS_Lyapunov, const arma::mat& var_run_mat, List& N_out,
  List& comm_out, List& N_out_optim, List& comm_out_optim,
  List& zero_stage_vec_list, DataFrame& trait_axis, DataFrame& new_trait_axis,
  DataFrame& optim_trait_axis, DataFrame& optim_trait_axis_995,
  List& new_stageexpansion_list, List& new_stageexpansion_list_optim,
  List& modified_dev_terms_list, List& errcheck_mpms, List& errcheck_mpms_optima,
  List& errcheck_mpm_ESS, List& errcheck_mpmouts, List& errcheck_mpmouts_optima,
  List& errcheck_mpmout_ESS, const List tweights_list, const List start_list,
  const List vrm_list, DataFrame current_stageframe, const List allmodels_all,
  const List allstages_all, const DataFrame current_supplement,
  const CharacterVector year_vec, const List ind_terms_num_list,
  const List ind_terms_cat_list, const List dev_terms_list,
  const DataFrame density_vr_list, const List sp_density_list,
  const DataFrame density_df, const DataFrame dens_index_df,
  const List equivalence_list, const IntegerVector sp_density_num_vec,
  const IntegerVector entry_time_vec, const IntegerVector inda_terms_num_vec,
  const IntegerVector indb_terms_num_vec, const IntegerVector indc_terms_num_vec,
  const IntegerVector inda_terms_cat_vec, const IntegerVector indb_terms_cat_vec,
  const IntegerVector indc_terms_cat_vec, const IntegerVector dens_vr_yn_vec,
  const IntegerVector tweights_type_vec, const NumericVector fecmod_vec,
  const CharacterVector patch_vec, const int variant_count,
  const int var_per_run, const int nreps, const int times,
  const int fitness_times, const int stagecounts, const int substoch,
  const int format_int, const int firstage_int, const int finalage_int,
  const int dev_terms_times_int, const int main_optim_res, const int opt_res,
  const double exp_tol, const double theta_tol, const int loop_max,
  const bool integeronly, const bool stochastic, const bool dens_yn_bool,
  const bool stages_not_equal, const bool sparse_bool, const bool historical,
  const bool pure_leslie, const bool A_only, const bool err_check,
  const bool err_check_extreme, const double threshold,
  const bool fitness_table, const bool ESS_optima, double elast_mult,
  const bool zap_min) {
  
  // Structures for optim
  arma::ivec variant_nta_optim;
  arma::vec surv_dev_nta_optim;
  arma::vec obs_dev_nta_optim;
  arma::vec size_dev_nta_optim;
  arma::vec sizeb_dev_nta_optim;
  arma::vec sizec_dev_nta_optim;
  arma::vec repst_dev_nta_optim;
  arma::vec fec_dev_nta_optim;
  arma::vec jsurv_dev_nta_optim;
  arma::vec jobs_dev_nta_optim;
  arma::vec jsize_dev_nta_optim;
  arma::vec jsizeb_dev_nta_optim;
  arma::vec jsizec_dev_nta_optim;
  arma::vec jrepst_dev_nta_optim;
  arma::vec jmatst_dev_nta_optim;
  arma::vec indcova_nta_optim;
  arma::vec indcovb_nta_optim;
  arma::vec indcovc_nta_optim;
  
  arma::ivec variant_nta_optim_995;
  arma::vec surv_dev_nta_optim_995;
  arma::vec obs_dev_nta_optim_995;
  arma::vec size_dev_nta_optim_995;
  arma::vec sizeb_dev_nta_optim_995;
  arma::vec sizec_dev_nta_optim_995;
  arma::vec repst_dev_nta_optim_995;
  arma::vec fec_dev_nta_optim_995;
  arma::vec jsurv_dev_nta_optim_995;
  arma::vec jobs_dev_nta_optim_995;
  arma::vec jsize_dev_nta_optim_995;
  arma::vec jsizeb_dev_nta_optim_995;
  arma::vec jsizec_dev_nta_optim_995;
  arma::vec jrepst_dev_nta_optim_995;
  arma::vec jmatst_dev_nta_optim_995;
  arma::vec indcova_nta_optim_995;
  arma::vec indcovb_nta_optim_995;
  arma::vec indcovc_nta_optim_995;
  
  List N_out_pre_optim (nreps);
  DataFrame ESS_trait_axis; 
  DataFrame stageframe_df_ESS;
  IntegerVector ESS_var_traits;
  IntegerVector flipped_traits (20);
  
  int opt_res_true = opt_res;
  int ehrlen_optim {0};
  int style_optim {0};
  int filter_optim {0};
  bool opt_res_squared {false};
  
  // Rcout << "invade3_fb_core A" << endl;
  // patch_vec???
  // density_vr_list
  // dens_vr_yn_vec
  
  // Vectors of axis variants
  IntegerVector axis_variant_vec = as<IntegerVector>(trait_axis["variant"]);
  IntegerVector axis_variants_unique = sort_unique(axis_variant_vec);
  
  int var_mat_length = static_cast<int>(var_run_mat.n_rows);
  
  DataFrame trait_axis_clone = clone(trait_axis);
  new_trait_axis = AdaptUtils::ta_reassess(current_stageframe, trait_axis_clone,
    firstage_int, historical, true, pure_leslie);
  // Rcout << "invade3_fb_core B" << endl;
  
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
  arma::vec indcova_nta = as<arma::vec>(new_trait_axis["indcova"]);
  arma::vec indcovb_nta = as<arma::vec>(new_trait_axis["indcovb"]);
  arma::vec indcovc_nta = as<arma::vec>(new_trait_axis["indcovc"]);
  
  // Rcout << "indcova in new trait axis: " << indcova_nta.t() << endl;
  // Rcout << "indcovb in new trait axis: " << indcovb_nta.t() << endl;
  // Rcout << "indcovc in new trait axis: " << indcovc_nta.t() << endl;
  
  // Rcout << "invade3_fb_core C" << endl;
  if (ESS_optima) {
    AdaptUtils::optim_ta_setup(new_trait_axis, ESS_trait_axis, optim_trait_axis,
      optim_trait_axis_995, ESS_var_traits, flipped_traits, opt_res,
      elast_mult);
    
    //int var_traits = sum(ESS_var_traits);
    if (opt_res_squared) opt_res_true = opt_res * opt_res; /////
  }
  
  // Rcout << "invade3_fb_core D" << endl;
  List comm_out_pre_optim (opt_res_true);
  List trait_axis_by_variant_optima (opt_res_true); // Might wish to remove this later
  List stageexpansion_by_variant_optima (opt_res_true); // Might wish to remove this later
  List stageexpansion_ta_devterms_by_variant_optima (opt_res_true);
  
  if (ESS_optima) {
    // Rcout << "invade3_fb_core E" << endl;
    variant_nta_optim = as<arma::ivec>(optim_trait_axis["variant"]);
    surv_dev_nta_optim = as<arma::vec>(optim_trait_axis["surv_dev"]);
    obs_dev_nta_optim = as<arma::vec>(optim_trait_axis["obs_dev"]);
    size_dev_nta_optim = as<arma::vec>(optim_trait_axis["size_dev"]);
    sizeb_dev_nta_optim = as<arma::vec>(optim_trait_axis["sizeb_dev"]);
    sizec_dev_nta_optim = as<arma::vec>(optim_trait_axis["sizec_dev"]);
    repst_dev_nta_optim = as<arma::vec>(optim_trait_axis["repst_dev"]);
    fec_dev_nta_optim = as<arma::vec>(optim_trait_axis["fec_dev"]);
    jsurv_dev_nta_optim = as<arma::vec>(optim_trait_axis["jsurv_dev"]);
    jobs_dev_nta_optim = as<arma::vec>(optim_trait_axis["jobs_dev"]);
    jsize_dev_nta_optim = as<arma::vec>(optim_trait_axis["jsize_dev"]);
    jsizeb_dev_nta_optim = as<arma::vec>(optim_trait_axis["jsizeb_dev"]);
    jsizec_dev_nta_optim = as<arma::vec>(optim_trait_axis["jsizec_dev"]);
    jrepst_dev_nta_optim = as<arma::vec>(optim_trait_axis["jrepst_dev"]);
    jmatst_dev_nta_optim = as<arma::vec>(optim_trait_axis["jmatst_dev"]);
    indcova_nta_optim = as<arma::vec>(optim_trait_axis["indcova"]);
    indcovb_nta_optim = as<arma::vec>(optim_trait_axis["indcovb"]);
    indcovc_nta_optim = as<arma::vec>(optim_trait_axis["indcovc"]);
    // Rcout << "invade3_fb_core F" << endl;
    
    variant_nta_optim_995 = as<arma::ivec>(optim_trait_axis_995["variant"]);
    surv_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["surv_dev"]);
    obs_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["obs_dev"]);
    size_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["size_dev"]);
    sizeb_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["sizeb_dev"]);
    sizec_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["sizec_dev"]);
    repst_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["repst_dev"]);
    fec_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["fec_dev"]);
    jsurv_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jsurv_dev"]);
    jobs_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jobs_dev"]);
    jsize_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jsize_dev"]);
    jsizeb_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jsizeb_dev"]);
    jsizec_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jsizec_dev"]);
    jrepst_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jrepst_dev"]);
    jmatst_dev_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["jmatst_dev"]);
    indcova_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["indcova"]);
    indcovb_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["indcovb"]);
    indcovc_nta_optim_995 = as<arma::vec>(optim_trait_axis_995["indcovc"]);
    // Rcout << "invade3_fb_core G" << endl;
    
    // This next loop sets up the structure for comm_out_pre_optim
    for (int i = 0; i < opt_res_true; i++) {
      //int year_length = static_cast<int>(year_vec.length());
      //IntegerVector year_int_vec = seq(0, (year_length - 1));
      
      List all_pops_per_run (var_per_run);
      List all_used_times_per_run (var_per_run);
      for (int m = 0; m < var_per_run; m++) {
        //int current_variant_index = var_run_mat(i, m);
        
        List pop_reps (nreps);
        List used_times_reps (nreps);
        for (int j = 0; j < nreps; j++) {
          arma::mat pops_out_pre (stagecounts, (times + 1), fill::zeros);
          pop_reps(j) = pops_out_pre;
        }
        //all_used_times_per_run(m) = used_times_reps;
        all_pops_per_run(m) = pop_reps;
      }
      //used_times(i) = all_used_times_per_run;
      comm_out_pre_optim(i) = all_pops_per_run;
    }
    
    // Rcout << "invade3_fb_core H" << endl;
    
    // Stage expansions for ESS optimization
    for (int i = 0; i < opt_res_true; i++) {
      IntegerVector used_i = {i + 1};
      StringVector focused_var = {"variant"};
      // Rcout << "invade3_fb_core I" << endl;
      
      DataFrame current_optim_trait_axis = LefkoUtils::df_subset(optim_trait_axis,
        as<RObject>(used_i), false, true, false, false, true, as<RObject>(focused_var));
      // Rcout << "invade3_fb_core J" << endl;
      
      DataFrame current_optim_trait_axis_995 = LefkoUtils::df_subset(optim_trait_axis_995,
        as<RObject>(used_i), false, true, false, false, true, as<RObject>(focused_var));
      // Rcout << "invade3_fb_core K" << endl;
      
      List tabvo = Rcpp::List::create(_["main"] = current_optim_trait_axis,
        _["e995"] = current_optim_trait_axis_995);
      trait_axis_by_variant_optima(i) = tabvo;
      
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
      
      ehrlen_optim = ehrlen;
      style_optim = style;
      filter_optim = filter;
      
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
      
      // Rcout << "invade3_fb_core L" << endl;
      
      stageframe_df_ESS = new_sf;
      DataFrame stageexpansion = AdaptMats::thenewpizzle(new_sf,
        current_optim_trait_axis, firstage_int, finalage_int, ehrlen, style,
        filter);
      // Rcout << "invade3_fb_core M" << endl;
      
      DataFrame stageexpansion_995 = AdaptMats::thenewpizzle(new_sf,
        current_optim_trait_axis_995, firstage_int, finalage_int, ehrlen, style,
        filter);
      // Rcout << "invade3_fb_core N" << endl;
      
      focused_var = {"mpm_altered"};
      IntegerVector chosen_int = {1};
      DataFrame stageexpansion_red_mpm = LefkoUtils::df_subset(stageexpansion,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "invade3_fb_core O" << endl;
      
      DataFrame stageexpansion_red_mpm_995 = LefkoUtils::df_subset(stageexpansion_995,
        as<RObject>(chosen_int), false, true, false, false, true,
        as<RObject>(focused_var));
      // Rcout << "invade3_fb_core P" << endl;
      
      List sbvo = Rcpp::List::create(_["main"] = stageexpansion_red_mpm,
        _["e995"] = stageexpansion_red_mpm_995);
      stageexpansion_by_variant_optima(i) = sbvo;
      
      NumericVector variant_ta_devterms (17);
      NumericVector variant_ta_devterms_995 (17);
      int current_variant_rows = current_optim_trait_axis.nrows();
      
      // Rcout << "invade3_fb_core Q" << endl;
      for (int j = 0; j < current_variant_rows; j++) {
        NumericVector cv_vt_survdev = as<NumericVector>(current_optim_trait_axis["surv_dev"]);
        NumericVector cv_vt_obsdev = as<NumericVector>(current_optim_trait_axis["obs_dev"]);
        NumericVector cv_vt_sizedev = as<NumericVector>(current_optim_trait_axis["size_dev"]);
        NumericVector cv_vt_sizebdev = as<NumericVector>(current_optim_trait_axis["sizeb_dev"]);
        NumericVector cv_vt_sizecdev = as<NumericVector>(current_optim_trait_axis["sizec_dev"]);
        NumericVector cv_vt_repstdev = as<NumericVector>(current_optim_trait_axis["repst_dev"]);
        NumericVector cv_vt_fecdev = as<NumericVector>(current_optim_trait_axis["fec_dev"]);
        
        NumericVector cv_vt_jsurvdev = as<NumericVector>(current_optim_trait_axis["jsurv_dev"]);
        NumericVector cv_vt_jobsdev = as<NumericVector>(current_optim_trait_axis["jobs_dev"]);
        NumericVector cv_vt_jsizedev = as<NumericVector>(current_optim_trait_axis["jsize_dev"]);
        NumericVector cv_vt_jsizebdev = as<NumericVector>(current_optim_trait_axis["jsizeb_dev"]);
        NumericVector cv_vt_jsizecdev = as<NumericVector>(current_optim_trait_axis["jsizec_dev"]);
        NumericVector cv_vt_jrepstdev = as<NumericVector>(current_optim_trait_axis["jrepst_dev"]);
        NumericVector cv_vt_jmatstdev = as<NumericVector>(current_optim_trait_axis["jmatst_dev"]);
        
        NumericVector cv_vt_indcova = as<NumericVector>(current_optim_trait_axis["indcova"]);
        NumericVector cv_vt_indcovb = as<NumericVector>(current_optim_trait_axis["indcovb"]);
        NumericVector cv_vt_indcovc = as<NumericVector>(current_optim_trait_axis["indcovc"]);
        
        NumericVector cv_vt_survdev_995 = as<NumericVector>(current_optim_trait_axis_995["surv_dev"]);
        NumericVector cv_vt_obsdev_995 = as<NumericVector>(current_optim_trait_axis_995["obs_dev"]);
        NumericVector cv_vt_sizedev_995 = as<NumericVector>(current_optim_trait_axis_995["size_dev"]);
        NumericVector cv_vt_sizebdev_995 = as<NumericVector>(current_optim_trait_axis_995["sizeb_dev"]);
        NumericVector cv_vt_sizecdev_995 = as<NumericVector>(current_optim_trait_axis_995["sizec_dev"]);
        NumericVector cv_vt_repstdev_995 = as<NumericVector>(current_optim_trait_axis_995["repst_dev"]);
        NumericVector cv_vt_fecdev_995 = as<NumericVector>(current_optim_trait_axis_995["fec_dev"]);
        
        NumericVector cv_vt_jsurvdev_995 = as<NumericVector>(current_optim_trait_axis_995["jsurv_dev"]);
        NumericVector cv_vt_jobsdev_995 = as<NumericVector>(current_optim_trait_axis_995["jobs_dev"]);
        NumericVector cv_vt_jsizedev_995 = as<NumericVector>(current_optim_trait_axis_995["jsize_dev"]);
        NumericVector cv_vt_jsizebdev_995 = as<NumericVector>(current_optim_trait_axis_995["jsizeb_dev"]);
        NumericVector cv_vt_jsizecdev_995 = as<NumericVector>(current_optim_trait_axis_995["jsizec_dev"]);
        NumericVector cv_vt_jrepstdev_995 = as<NumericVector>(current_optim_trait_axis_995["jrepst_dev"]);
        NumericVector cv_vt_jmatstdev_995 = as<NumericVector>(current_optim_trait_axis_995["jmatst_dev"]);
        
        NumericVector cv_vt_indcova_995 = as<NumericVector>(current_optim_trait_axis_995["indcova"]);
        NumericVector cv_vt_indcovb_995 = as<NumericVector>(current_optim_trait_axis_995["indcovb"]);
        NumericVector cv_vt_indcovc_995 = as<NumericVector>(current_optim_trait_axis_995["indcovc"]);
        
        if (!NumericVector::is_na(cv_vt_survdev(j))) variant_ta_devterms(0) =
          variant_ta_devterms(0) + cv_vt_survdev(j);
        if (!NumericVector::is_na(cv_vt_obsdev(j))) variant_ta_devterms(1) =
          variant_ta_devterms(1) + cv_vt_obsdev(j);
        if (!NumericVector::is_na(cv_vt_sizedev(j))) variant_ta_devterms(2) =
          variant_ta_devterms(2) + cv_vt_sizedev(j);
        if (!NumericVector::is_na(cv_vt_sizebdev(j))) variant_ta_devterms(3) =
          variant_ta_devterms(3) + cv_vt_sizebdev(j);
        if (!NumericVector::is_na(cv_vt_sizecdev(j))) variant_ta_devterms(4) =
          variant_ta_devterms(4) + cv_vt_sizecdev(j);
        if (!NumericVector::is_na(cv_vt_repstdev(j))) variant_ta_devterms(5) =
          variant_ta_devterms(5) + cv_vt_repstdev(j);
        if (!NumericVector::is_na(cv_vt_fecdev(j))) variant_ta_devterms(6) =
          variant_ta_devterms(6) + cv_vt_fecdev(j);
        
        if (!NumericVector::is_na(cv_vt_jsurvdev(j))) variant_ta_devterms(7) =
          variant_ta_devterms(7) + cv_vt_jsurvdev(j);
        if (!NumericVector::is_na(cv_vt_jobsdev(j))) variant_ta_devterms(8) =
          variant_ta_devterms(8) + cv_vt_jobsdev(j);
        if (!NumericVector::is_na(cv_vt_jsizedev(j))) variant_ta_devterms(9) =
          variant_ta_devterms(9) + cv_vt_jsizedev(j);
        if (!NumericVector::is_na(cv_vt_jsizebdev(j))) variant_ta_devterms(10) =
          variant_ta_devterms(10) + cv_vt_jsizebdev(j);
        if (!NumericVector::is_na(cv_vt_jsizecdev(j))) variant_ta_devterms(11) =
          variant_ta_devterms(11) + cv_vt_jsizecdev(j);
        if (!NumericVector::is_na(cv_vt_jrepstdev(j))) variant_ta_devterms(12) =
          variant_ta_devterms(12) + cv_vt_jrepstdev(j);
        if (!NumericVector::is_na(cv_vt_jmatstdev(j))) variant_ta_devterms(13) =
          variant_ta_devterms(13) + cv_vt_jmatstdev(j);
        if (!NumericVector::is_na(cv_vt_indcova(j))) variant_ta_devterms(14) =
          variant_ta_devterms(14) + cv_vt_indcova(j);
        if (!NumericVector::is_na(cv_vt_indcovb(j))) variant_ta_devterms(15) =
          variant_ta_devterms(15) + cv_vt_indcovb(j);
        if (!NumericVector::is_na(cv_vt_indcovc(j))) variant_ta_devterms(16) =
          variant_ta_devterms(16) + cv_vt_indcovc(j);
        
        if (!NumericVector::is_na(cv_vt_survdev_995(j))) variant_ta_devterms_995(0) =
          variant_ta_devterms_995(0) + cv_vt_survdev_995(j);
        if (!NumericVector::is_na(cv_vt_obsdev_995(j))) variant_ta_devterms_995(1) =
          variant_ta_devterms_995(1) + cv_vt_obsdev_995(j);
        if (!NumericVector::is_na(cv_vt_sizedev_995(j))) variant_ta_devterms_995(2) =
          variant_ta_devterms_995(2) + cv_vt_sizedev_995(j);
        if (!NumericVector::is_na(cv_vt_sizebdev_995(j))) variant_ta_devterms_995(3) =
          variant_ta_devterms_995(3) + cv_vt_sizebdev_995(j);
        if (!NumericVector::is_na(cv_vt_sizecdev_995(j))) variant_ta_devterms_995(4) =
          variant_ta_devterms_995(4) + cv_vt_sizecdev_995(j);
        if (!NumericVector::is_na(cv_vt_repstdev_995(j))) variant_ta_devterms_995(5) =
          variant_ta_devterms_995(5) + cv_vt_repstdev_995(j);
        if (!NumericVector::is_na(cv_vt_fecdev_995(j))) variant_ta_devterms_995(6) =
          variant_ta_devterms_995(6) + cv_vt_fecdev_995(j);
        
        if (!NumericVector::is_na(cv_vt_jsurvdev_995(j))) variant_ta_devterms_995(7) =
          variant_ta_devterms_995(7) + cv_vt_jsurvdev_995(j);
        if (!NumericVector::is_na(cv_vt_jobsdev_995(j))) variant_ta_devterms_995(8) =
          variant_ta_devterms_995(8) + cv_vt_jobsdev_995(j);
        if (!NumericVector::is_na(cv_vt_jsizedev_995(j))) variant_ta_devterms_995(9) =
          variant_ta_devterms_995(9) + cv_vt_jsizedev_995(j);
        if (!NumericVector::is_na(cv_vt_jsizebdev_995(j))) variant_ta_devterms_995(10) =
          variant_ta_devterms_995(10) + cv_vt_jsizebdev_995(j);
        if (!NumericVector::is_na(cv_vt_jsizecdev_995(j))) variant_ta_devterms_995(11) =
          variant_ta_devterms_995(11) + cv_vt_jsizecdev_995(j);
        if (!NumericVector::is_na(cv_vt_jrepstdev_995(j))) variant_ta_devterms_995(12) =
          variant_ta_devterms_995(12) + cv_vt_jrepstdev_995(j);
        if (!NumericVector::is_na(cv_vt_jmatstdev_995(j))) variant_ta_devterms_995(13) =
          variant_ta_devterms_995(13) + cv_vt_jmatstdev_995(j);
        
        if (!NumericVector::is_na(cv_vt_indcova_995(j))) variant_ta_devterms_995(14) =
          variant_ta_devterms_995(14) + cv_vt_indcova_995(j);
        if (!NumericVector::is_na(cv_vt_indcovb_995(j))) variant_ta_devterms_995(15) =
          variant_ta_devterms_995(15) + cv_vt_indcovb_995(j);
        if (!NumericVector::is_na(cv_vt_indcovc_995(j))) variant_ta_devterms_995(16) =
          variant_ta_devterms_995(16) + cv_vt_indcovc_995(j);
      }
      
      List stdbvo = Rcpp::List::create(_["main"] = variant_ta_devterms,
        _["e995"] = variant_ta_devterms_995);
      stageexpansion_ta_devterms_by_variant_optima(i) = stdbvo;
    }
    new_stageexpansion_list_optim = stageexpansion_by_variant_optima;
  }
  
  // Rcout << "invade3_fb_core R" << endl;
  
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
  
  IntegerMatrix inda_num_terms_counter_optim (opt_res_true, 2);
  IntegerMatrix indb_num_terms_counter_optim (opt_res_true, 2);
  IntegerMatrix indc_num_terms_counter_optim (opt_res_true, 2);
  IntegerMatrix inda_cat_terms_counter_optim (opt_res_true, 2);
  IntegerMatrix indb_cat_terms_counter_optim (opt_res_true, 2);
  IntegerMatrix indc_cat_terms_counter_optim (opt_res_true, 2);
  IntegerMatrix inda_num_terms_previous_optim (opt_res_true, 2);
  IntegerMatrix indb_num_terms_previous_optim (opt_res_true, 2);
  IntegerMatrix indc_num_terms_previous_optim (opt_res_true, 2);
  IntegerMatrix inda_cat_terms_previous_optim (opt_res_true, 2);
  IntegerMatrix indb_cat_terms_previous_optim (opt_res_true, 2);
  IntegerMatrix indc_cat_terms_previous_optim (opt_res_true, 2);
  IntegerMatrix dev_num_counter_optim (opt_res_true, 2);
  IntegerMatrix sp_density_counter_optim (opt_res_true, 2);
  
  // Year order determination
  List comm_out_pre (var_mat_length);
  List used_times (var_mat_length);
  
  // Rcout << "invade3_fb_core S" << endl;
  
  // This next loop determines the order of matrices for each variant run
  for (int i = 0; i < var_mat_length; i++) {
    int year_length = static_cast<int>(year_vec.length());
    IntegerVector year_int_vec = seq(0, (year_length - 1));
    
    List all_pops_per_run (var_per_run);
    List all_used_times_per_run (var_per_run);
    for (int m = 0; m < var_per_run; m++) {
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
  
  // Rcout << "invade3_fb_core T" << endl;
  
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
    
    NumericVector variant_ta_devterms (17);
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
      
      NumericVector cv_vt_indcova = as<NumericVector>(current_trait_axis["indcova"]);
      NumericVector cv_vt_indcovb = as<NumericVector>(current_trait_axis["indcovb"]);
      NumericVector cv_vt_indcovc = as<NumericVector>(current_trait_axis["indcovc"]);
      
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
      
      if (!NumericVector::is_na(cv_vt_indcova(j))) variant_ta_devterms(14) = variant_ta_devterms(14) + cv_vt_indcova(j);
      if (!NumericVector::is_na(cv_vt_indcovb(j))) variant_ta_devterms(15) = variant_ta_devterms(15) + cv_vt_indcovb(j);
      if (!NumericVector::is_na(cv_vt_indcovc(j))) variant_ta_devterms(16) = variant_ta_devterms(16) + cv_vt_indcovc(j);
    }
    stageexpansion_ta_devterms_by_variant(i) = variant_ta_devterms;
  }
  new_stageexpansion_list = stageexpansion_by_variant;
  
  // Rcout << "invade3_fb_core U" << endl;
  
  // Main projection
  List N_out_pre (nreps);
  List mdtl (var_mat_length);
  List mdtl_optim (opt_res_true);
  List errcheck_mpm_reps (nreps);
  List errcheck_mpm_reps_optima (nreps);
  List errcheck_mpmout_reps (nreps);
  List errcheck_mpmout_reps_optima (nreps);
  
  // Rcout << "indcova_nta: " << indcova_nta.t() << endl;
  // Rcout << "indcovb_nta: " << indcovb_nta.t() << endl;
  // Rcout << "indcovc_nta: " << indcovc_nta.t() << endl;
  
  for (int current_rep = 0; current_rep < nreps; current_rep++) { // 1st loop - reps current_rep
    invfb_project(var_run_mat, surv_dev_nta, obs_dev_nta, size_dev_nta,
      sizeb_dev_nta, sizec_dev_nta, repst_dev_nta, fec_dev_nta, jsurv_dev_nta,
      jobs_dev_nta, jsize_dev_nta, jsizeb_dev_nta, jsizec_dev_nta,
      jrepst_dev_nta, jmatst_dev_nta, indcova_nta, indcovb_nta, indcovc_nta,
      variant_nta, N_out_pre, comm_out_pre, new_stageexpansion_list,
      errcheck_mpm_reps, errcheck_mpmout_reps, mdtl, used_times, allmodels_all,
      vrm_list, allstages_all, dev_terms_list, ind_terms_num_list, ind_terms_cat_list,
      stageexpansion_ta_devterms_by_variant, sp_density_list, start_list,
      equivalence_list, density_vr_list, current_stageframe, current_supplement,
      density_df, dens_index_df, entry_time_vec, sp_density_num_vec,
      inda_terms_num_vec, indb_terms_num_vec, indc_terms_num_vec, inda_terms_cat_vec,
      indb_terms_cat_vec, indc_terms_cat_vec, dens_vr_yn_vec, sp_density_counter,
      inda_num_terms_previous, indb_num_terms_previous, indc_num_terms_previous,
      inda_cat_terms_previous, indb_cat_terms_previous, indc_cat_terms_previous,
      inda_num_terms_counter, indb_num_terms_counter, indc_num_terms_counter,
      inda_cat_terms_counter, indb_cat_terms_counter, indc_cat_terms_counter,
      dev_num_counter, fecmod_vec, year_vec, patch_vec, var_per_run, times,
      var_mat_length, format_int, current_rep, firstage_int, finalage_int,
      dev_terms_times_int, substoch, year_counter, exp_tol, theta_tol,
      err_check, err_check_extreme, sparse_bool, A_only, stages_not_equal,
      integeronly, dens_yn_bool);
    
    // Rcout << "obs_dev_nta_optim: " << obs_dev_nta_optim.t() << endl;
    // Rcout << "obs_dev_nta_optim_995: " << obs_dev_nta_optim_995.t() << endl;
    
    if (ESS_optima) {
      invfb_optim(surv_dev_nta_optim, obs_dev_nta_optim, size_dev_nta_optim,
        sizeb_dev_nta_optim, sizec_dev_nta_optim, repst_dev_nta_optim,
        fec_dev_nta_optim, jsurv_dev_nta_optim, jobs_dev_nta_optim,
        jsize_dev_nta_optim, jsizeb_dev_nta_optim, jsizec_dev_nta_optim,
        jrepst_dev_nta_optim, jmatst_dev_nta_optim, indcova_nta_optim,
        indcovb_nta_optim, indcovc_nta_optim, variant_nta_optim,
        surv_dev_nta_optim_995, obs_dev_nta_optim_995, size_dev_nta_optim_995,
        sizeb_dev_nta_optim_995, sizec_dev_nta_optim_995,
        repst_dev_nta_optim_995, fec_dev_nta_optim_995, jsurv_dev_nta_optim_995,
        jobs_dev_nta_optim_995, jsize_dev_nta_optim_995,
        jsizeb_dev_nta_optim_995, jsizec_dev_nta_optim_995,
        jrepst_dev_nta_optim_995, jmatst_dev_nta_optim_995,
        indcova_nta_optim_995, indcovb_nta_optim_995, indcovc_nta_optim_995,
        variant_nta_optim_995, N_out_pre_optim, comm_out_pre_optim,
        new_stageexpansion_list_optim, errcheck_mpm_reps_optima,
        errcheck_mpmout_reps_optima, mdtl_optim, used_times, allmodels_all,
        vrm_list, allstages_all, dev_terms_list, ind_terms_num_list,
        ind_terms_cat_list, stageexpansion_ta_devterms_by_variant_optima,
        sp_density_list, start_list, equivalence_list, density_vr_list,
        current_stageframe, current_supplement, density_df, dens_index_df,
        entry_time_vec, sp_density_num_vec, inda_terms_num_vec,
        indb_terms_num_vec, indc_terms_num_vec, inda_terms_cat_vec,
        indb_terms_cat_vec, indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec,
        year_vec, patch_vec, var_per_run, times, var_mat_length, format_int,
        current_rep, firstage_int, finalage_int, dev_terms_times_int, substoch,
        opt_res_true, opt_res, year_counter, exp_tol, theta_tol, threshold,
        err_check, err_check_extreme, sparse_bool, A_only, stages_not_equal,
        integeronly, dens_yn_bool);
    }
  } // current_rep loop
  comm_out = comm_out_pre;
  comm_out_optim = comm_out_pre_optim;
  
  N_out = N_out_pre;
  N_out_optim = N_out_pre_optim;
  
  if (err_check) {
    errcheck_mpmouts = errcheck_mpmout_reps;
    errcheck_mpmouts_optima = errcheck_mpmout_reps_optima;
    errcheck_mpms = errcheck_mpm_reps;
    errcheck_mpms_optima = errcheck_mpm_reps_optima;
    modified_dev_terms_list = mdtl;
  }
  
  // Rcout << "invade3_fb_core V        ";
  
  if (fitness_table) {
    AdaptUtils::Lyapunov_creator (Lyapunov, N_out, entry_time_vec, nreps, var_per_run,
      var_mat_length, times, fitness_times, threshold, 0, false, zap_min);
  }
  
  // Rcout << "invade3_fb_core W right before call to ESS_optimizer_fb       " << endl;
  
  if (ESS_optima) {
    AdaptUtils::Lyapunov_creator (Lyapunov_optim, N_out_optim, entry_time_vec, nreps, 2,
      opt_res_true, times, fitness_times, threshold, 1, false, zap_min);
    
    ESS_optimizer_fb(ESS_Lyapunov, ESS_trait_axis, Lyapunov_optim,
      optim_trait_axis, ESS_var_traits, flipped_traits, surv_dev_nta_optim,
      obs_dev_nta_optim, size_dev_nta_optim, sizeb_dev_nta_optim,
      sizec_dev_nta_optim, repst_dev_nta_optim, fec_dev_nta_optim,
      jsurv_dev_nta_optim, jobs_dev_nta_optim, jsize_dev_nta_optim,
      jsizeb_dev_nta_optim, jsizec_dev_nta_optim, jrepst_dev_nta_optim,
      jmatst_dev_nta_optim, indcova_nta_optim, indcovb_nta_optim,
      indcovc_nta_optim, variant_nta_optim, new_stageexpansion_list_optim,
      used_times, errcheck_mpm_ESS, errcheck_mpmout_ESS, allmodels_all, vrm_list,
      allstages_all, dev_terms_list, ind_terms_num_list, ind_terms_cat_list,
      stageexpansion_ta_devterms_by_variant_optima, sp_density_list, start_list,
      equivalence_list, density_vr_list, current_stageframe, current_supplement,
      density_df, dens_index_df, stageframe_df_ESS, entry_time_vec,
      sp_density_num_vec, inda_terms_num_vec, indb_terms_num_vec,
      indc_terms_num_vec, inda_terms_cat_vec, indb_terms_cat_vec,
      indc_terms_cat_vec, dens_vr_yn_vec, fecmod_vec, year_vec, patch_vec, times,
      fitness_times, format_int, stagecounts, firstage_int, finalage_int,
      dev_terms_times_int, substoch, exp_tol, theta_tol, sparse_bool, A_only,
      stages_not_equal, integeronly, dens_yn_bool, threshold, opt_res_true,
      opt_res, ehrlen_optim, style_optim, loop_max, filter_optim,
      err_check_extreme, elast_mult, zap_min);
  }
  
  // Rcout << "invade3_fb_core X right after call to ESS_optimizer_fb        " << endl;
}

//' Run Pairwise and Multiple Invasibility Analysis
//' 
//' Function \code{invade3} runs pairwise and multiple invasisibility analyses,
//' and provides output to analyze all aspects of the dynamics of generic
//' variants. The output includes invasion fitness for all variants included, as
//' well as resident fitness in all simulations run. Forms the core data plotted
//' by the \code{\link{plot.adaptInv}()} function to plot pairwise
//' invasibility plots (PIPs) and elasticity plots. Also estimates trait optima,
//' if any occur.
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
//' competition terms, such as are used in multiple species competition models.
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
//' @param tweights An optional numeric vector or matrix denoting the
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
//' argument \code{var_per_run}. Note that if two variants are to be run at a
//' time, as in a pairwise invasion analysis, then the length of the vector
//' should be equal to 2.
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
//' @param dev_terms An optional data frame including 14 or 17 columns and up to
//' \code{times} rows showing the values of the deviation terms to be added to
//' each linear vital rate. The column order should be: 1: survival,
//' 2: observation, 3: primary size, 4: secondary size, 5: tertiary size,
//' 6: reproduction, 7: fecundity, 8: juvenile survival, 9: juvenile
//' observation, 10: juvenile primary size, 11: juvenile secondary size,
//' 12: juvenile tertiary size, 13: juvenile reproduction, and 14: juvenile
//' maturity transition. Optionally, the user may add three columns designating
//' additive deviations to: 15: individual covariate a, 16: individual covariate
//' b, and 17: individual covariate c. Unused terms must be set to \code{0} (use
//' of \code{NA} will produce errors). Single or small numbers of values per
//' vital rate model are also allowed, and if the number of rows is less than
//' \code{times}, then the terms will be cycled.
//' @param fb_sparse A logical value indicating whether function-based MPMs
//' should be produced in sparse matrix format. Defaults to \code{FALSE}.
//' @param firstage An optional integer vector used for function-based Leslie
//' and age-by-stage MPMs giving the starting ages in such MPMs. Use only if at
//' least one MPM is both function-based and has age structure. Typically,
//' the starting age in such MPMs should be set to \code{0} if post-breeding and
//' \code{1} if pre-breeding. All other MPMs should be set to \code{0}. Do not
//' use if no MPM has age structure. Defaults to \code{1} in Leslie and
//' age-by-stage MPMs.
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
//' @param stochastic A logical value indicating whether the projection will be
//' run as a temporally stochastic projection. Defaults to \code{FALSE}.
//' @param A_only A logical value indicating whether to alter survival and
//' fecundity matrix elements separately prior to creating the overall \code{A}
//' matrix, or whether to alter elements directly on \code{A} matrices. Defaults
//' to \code{TRUE}, and should be kept to that setting unless some matrix
//' elements to be altered are sums of survival and fecundity transitions.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion down to the next lower
//' integer. Defaults to \code{FALSE}.
//' @param fitness_table A logical value dictating whether to include a data
//' frame giving Lyapunov coefficients for all combinations of variants tested.
//' Necessary for the creation of pairwise invasibility plots (PIPs). Defaults
//' to \code{TRUE}.
//' @param trait_optima A logical value indicating whether to assess the optimal
//' values of traits, generally as kinds of evolutionary stage equilibrium (ESS)
//' points. Trait optimization is conducted via elasticity analysis of traits
//' that are variable within the \code{trait_axis} table. Defaults to
//' \code{FALSE}.
//' @param zap_min A logical value indicating whether to treat traits and
//' fitness as 0 when their absolute values are less than the value given in
//' argument \code{threshold}.
//' @param converged_only A logical value indicating whether to show predicted
//' trait optima only in cases where the elasticity analysis has converged.
//' Defaults to \code{FALSE}.
//' @param err_check A logical value indicating whether to include an extra list
//' of output objects for error checking. Can also be set to the text value
//' \code{"extreme"}, in which case all \code{err_check} output plus a multiple
//' level list with each MPM used in each time step will be output.
//' @param var_per_run The number of variants to run in each simulation.
//' Defaults to \code{2}, resulting in pairwise invasibility analysis. See
//' \code{Notes} for details.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent and density
//' independent simulations. Defaults to \code{0}, which does not enforce
//' substochasticity. Alternatively, \code{1} forces all survival-transition
//' elements to range from 0.0 to 1.0, and forces fecundity to be non-negative;
//' and \code{2} forces all column rows in the survival-transition matrices to
//' total no more than 1.0, in addition to the actions outlined for option
//' \code{1}. Both settings \code{1} and \code{2} change negative fecundity
//' elements to \code{0.0}, while setting \code{0} does not alter fecundity.
//' @param elast_mult A multiplier for traits to assess the elasticity of
//' fitness in trait optimization. Defaults to \code{0.995}.
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
//' @param threshold The lower limit for the absolute value of fitness, below
//' which fitness is rounded to 0. Defaults to 0.00000001.
//' @param loop_max An integer value denoting the number of search cycles
//' allowed per ESS during ESS optimization. Defaults to 150.
//' 
//' @return A list of class \code{adaptInv}, with the following elements:
//' \item{fitness}{A data frame giving the resident and invasion fitness
//' estimated for each variant, per replicate. In typical cases, the user will
//' be interested in the invasion fitness, which corresponds to the Lyapunov
//' coefficient of the invader in invasibility analysis.}
//' \item{variants_out}{A two-level list with the top level list having number of
//' elements equal to the number of variants, and the lower level corresponding
//' to the number of replicates. Each element of the lower level list is a
//' matrix showing the number of individuals in each stage (row) at each time
//' (column).}
//' \item{N_out}{A list with the number of elements equal to the number of
//' replicates. Each element within this list is data frame showing the number
//' of individuals of each species or genotype alive at each time. The number of
//' rows are equal to the number of MPMs used, and the columns correspond to the
//' time steps.}
//' \item{trait_axis}{An edited version of the input trait axis.}
//' \item{stageframe_list}{A list in which each element is the stageframe for
//' each MPM used.}
//' \item{hstages_list}{A list giving the used \code{hstages} data frames, which
//' identify the correct stage pairing for each row / column in each historical
//' MPM utilized.}
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
//' The argument \code{var_per_run} establishes the style of simulation to run.
//' Entering \code{var_per_run = 1} runs each variant singly. Entering
//' \code{var_per_run = 2} runs pairwise invasibility analysis, trying each pair
//' permutation of variants. Greater values will lead to multiple invasibility
//' analysis with different permutations of groups. For example,
//' \code{var_per_run = 3} runs each permutation of groups of three. The integer
//' set must be positive, and must not be larger than the number of variants.
//' 
//' When \code{optima = TRUE}, ESS values for traits that vary in the input
//' \code{adaptAxis} data frame are evaluated. The methodology is that
//' originally developed in Benton and Grant (1999, Evolution 53:677-688), as
//' communicated in Roff (2010, Modeling evolution: an introduction to numerical
//' methods, Oxford University Press). In essence, function \code{invade3}
//' determines which traits vary among all traits noted in the input trait axis.
//' A new trait axis is then created with values of variable traits multiplied
//' by \code{0.995}, although this can be reset with argument \code{elast_mult}.
//' This new trait axis is composed entirely of invaders that will be paired
//' against each respective row in the original trait axis. These two trait axis
//' frames are then used to conduct pairwise invasibility elasticity analyses,
//' particularly noting where fitness values and trends invert. Elasticity plots
//' show the results of these elasticity simulations, using R's line segment
//' default in the \code{plot.default()} function. Results are then used to
//' determine golden sections for further evaluation, within which the bounded
//' Nelder-Mead algorithm is used to determine the optimum value. Note that this
//' optimization approach typically works properly with only one variable trait.
//' Also note that if a trait spans both negative and positive values, then
//' elasticity calculation is altered such that invader values for negative
//' values are calculated as the sum of the old trait value and the difference
//' between that value and the value when multiplied by the elasticity value.
//' This is done to ensure that all elasticity values are consistently lower
//' than the original values.
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
  Nullable<RObject> stochastic = R_NilValue, Nullable<RObject> A_only = R_NilValue,
  Nullable<RObject> integeronly = R_NilValue, Nullable<RObject> fitness_table = R_NilValue,
  Nullable<RObject> trait_optima = R_NilValue, Nullable<RObject> zap_min = R_NilValue,
  Nullable<RObject> converged_only = R_NilValue, Nullable<RObject> err_check = R_NilValue,
  
  int var_per_run = 2,
  
  int substoch = 0, double elast_mult = 0.995,
  int nreps = 1, int times = 10000, int fitness_times = 100,
  double exp_tol = 700.0, double theta_tol = 100000000.0,
  double threshold = 0.00000001, int loop_max = 150) {
  
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
  int opt_res {100};
  int main_optim_res {0};
  
  bool preexisting {false}; // Are preexisting MPMs being used?
  bool funcbased {false}; // Will function-based MPMs be created?
  bool entry_time_vec_use {false}; // Are any elements in entry_time greater than 0?
  bool stages_not_equal {false}; // Are equivalence vectors supplied separating even stages?
  bool pure_leslie {false}; // Is core MPM Leslie?
  bool dens_yn_bool {false}; // density input for each MPM (0 = no, 1 = yes)
  bool prebreeding_bool {true};
  bool stochastic_bool {false};
  bool integeronly_bool {false};
  bool A_only_bool {true};
  bool fitness_table_bool {true};
  bool trait_optima_bool {false};
  
  //bool inda_char {false};
  //bool indb_char {false};
  //bool indc_char {false};
  bool historical {false};
  bool sparse_bool {false};
  bool zap_min_bool {true};
  bool converged_only_bool {false};
  bool err_check_bool {false};
  bool err_check_extreme {false};
  
  // Boolean variable processing
  LefkoInputs::RObj_TF_input_check("err_check", "extreme", err_check_bool,
    err_check_extreme, true, true, err_check);
  
  if(stochastic.isNotNull()) {
    stochastic_bool = LefkoInputs::yesno_to_logic(as<RObject>(stochastic), "stochastic");
  }
  
  if(integeronly.isNotNull()) {
    integeronly_bool = LefkoInputs::yesno_to_logic(as<RObject>(integeronly), "integeronly");
  }
  
  if(A_only.isNotNull()) {
    A_only_bool = LefkoInputs::yesno_to_logic(as<RObject>(A_only), "A_only");
  }
  
  if(fitness_table.isNotNull()) {
    fitness_table_bool = LefkoInputs::yesno_to_logic(as<RObject>(fitness_table), "fitness_table");
  }
  
  if(trait_optima.isNotNull()) {
    trait_optima_bool = LefkoInputs::yesno_to_logic(as<RObject>(trait_optima), "trait_optima");
  }
  
  if(zap_min.isNotNull()) {
    zap_min_bool = LefkoInputs::yesno_to_logic(as<RObject>(zap_min), "zap_min");
  }
  
  if(converged_only.isNotNull()) {
    converged_only_bool = LefkoInputs::yesno_to_logic(as<RObject>(converged_only), "converged_only");
  }
  
  // Age-by-stage and Leslie MPM settings
  int firstage_int;
  int finalage_int;
  IntegerVector cont_vec;
  
  // Leslie MPM only
  IntegerVector fecage_min_vec;
  IntegerVector fecage_max_vec;
  
  // Main data frames
  DataFrame final_stageframe;
  DataFrame stageframe_df; // List ending in _fb are only used in function-based cases
  DataFrame supplement_df;
  DataFrame supplement_list_fb;
  DataFrame density_df;
  DataFrame dens_index_df;  // Holds element index vectors for density_frames in density_df
  DataFrame density_vr_list;
  DataFrame trait_axis; // DataFrame used to house trait axis
  DataFrame new_trait_axis; // Data frame housing re-assessed trait_axis
  DataFrame optim_trait_axis; // Data frame housing re-assessed trait axis for optimization
  DataFrame optim_trait_axis_995; // Data frame housing re-assessed trait axis for optimization
  DataFrame Lyapunov; // DataFrame showing fitness of all genotypes under all circumstances
  DataFrame Lyapunov_optim; // DataFrame showing fitness of all genotypes in ESS optima table
  DataFrame ESS_Lyapunov; // DataFrame showing optimized (ESS) values of traits
  DataFrame labels;  // Data frame to provide order of MPMs
  
  // Main lists
  List comm_out; // List of mats of pop vecs (top lvl: mpm, lower lvl: reps, mats: stages x times)
  List comm_out_optim; // List of optim mats of pop vecs (top lvl: mpm, lower lvl: reps, mats: stages x times)
  List N_out;  // List of pop size mats (top level: reps, mats: mpm by times)
  List N_out_optim;  // List of optim pop size mats (top level: reps, mats: mpm by times)
  List mpm_list;
  List A_list;
  List vrm_list;
  List labels_list;  // List to hold data for data frame labels
  List repmatrix_list;
  List equivalence_list;
  List hstages_list;
  List agestages_list;
  List start_list;
  List tweights_list;
  List sp_density_list;
  List ind_terms_num_list;
  List ind_terms_cat_list;
  List dev_terms_list;
  List modified_dev_terms_list; // Only used in fb processing
  List allstages_all; // Used in fbMPM processing
  List allmodels_all; // Used in fbMPM processirg
  List stageexpansion_list; // Classic stage expansion data frame for each variant
  List stageexpansion_list_optim; // Classic stage expansion data frame for each variant in ESS evaluation
  List errcheck_mpm_list;
  List errcheck_mpm_list_optima;
  List errcheck_mpm_list_ESS;
  List errcheck_mpmout_list;
  List errcheck_mpmout_list_optima;
  List errcheck_mpmout_list_ESS;
  
  CharacterVector patch_vec; // choice of patch in MPM
  CharacterVector year_vec; // choice of years in MPM
  NumericVector equivalence_vec; // equivalence vector if !stages_not_equal
  NumericVector fecmod_vec; // fecundity multipliers for multiple offspring stages
  
  IntegerVector matrowcounts; // # rows in each MPM
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
  
  IntegerVector axis_variants_unique; // Vector giving the variant numbers
  arma::mat var_run_mat; // Holds indices for variants to run per run in invasion analysis
  List zero_stage_vec_list; // Holds zero stage vectors, if entry times are staggered
  
  // Rcout << "invade3 a ";
  
  if (fitness_times > times) {
    Rf_warningcall(R_NilValue,
      "Argument fitness_times is too large. Reseting to value of argument times.");
    
    fitness_times = times;
  }
  
  // Rcout << "invade3 b ";
  
  if (axis.isNotNull()) {
    if (is<DataFrame>(axis)) {
      trait_axis = as<DataFrame>(axis);
      CharacterVector ta_vars = trait_axis.names();
      int ta_vars_num = static_cast<int>(ta_vars.length());
      
      if (ta_vars_num != 33) {
        throw Rcpp::exception("Argument axis is not recognized.", false);
      }
      
      IntegerVector axis_variants = as<IntegerVector>(trait_axis["variant"]);
      axis_variants_unique = sort_unique(axis_variants);
      variant_count = static_cast<unsigned int>(axis_variants_unique.length());
      
      if (!(variant_count >= var_per_run)) {
        throw Rcpp::exception("Argument var_per_run may not be greater than the number of variants.", false);
      }
      if (var_per_run < 1) {
        throw Rcpp::exception("Argument var_per_run must be a positive integer.", false);
      }
      //double exp_permutes = pow(variant_count, var_per_run);
      
      var_run_mat = AdaptUtils::exp_grid_single(variant_count, var_per_run); // Rows = runs, Cols = variants
      
    } else {
      AdaptUtils::pop_error2("axis", "a data frame of class adaptAxis", "", 1);
    }
  } else {
    throw Rcpp::exception("Argument axis is required to run invasion analyses.", false);
  }
  
  //DataFrame final_Lyapunov = clone(Lyapunov); // Remove later
  AdaptUtils::Lyapunov_df_maker(Lyapunov, var_run_mat, var_per_run, axis_variants_unique);
  
  List cleaned_input = cleanup3_inv(mpm, vrm, stageframe, supplement, format,
    firstage, finalage, fecage_min, fecage_max, cont, fecmod, starts, patches,
    years, tweights, density, entry_time, density_vr, sp_density, ind_terms,
    dev_terms, fb_sparse, equivalence, prebreeding, exp_tol, theta_tol,
    substoch, variant_count, var_per_run);
  
  // Rcout << "invade3 c ";
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
  
  IntegerVector variant_ta = trait_axis["variant"];
  IntegerVector variant_ta_unique = unique(variant_ta);
  main_optim_res = static_cast<int>(variant_ta_unique.length());
  opt_res = static_cast<int>(variant_ta_unique.length());
  
  if (total_mpms != 1) {
    throw Rcpp::exception("Function invade3 only allows a single lefkoMat or vrm_input object.",
      false);
  }
  
  // Rcout << "invade3 d ";
  // Projection runs
  if (preexisting) {
    // Rcout << "invade3 e ";
    invade3_pre_core (Lyapunov, Lyapunov_optim, ESS_Lyapunov, var_run_mat,
      N_out, comm_out, N_out_optim, comm_out_optim, zero_stage_vec_list,
      trait_axis, new_trait_axis, optim_trait_axis, optim_trait_axis_995,
      stageexpansion_list, stageexpansion_list_optim, errcheck_mpm_list,
      errcheck_mpm_list_optima, mpm_list, tweights_list, start_list, vrm_list,
      stageframe_df, allmodels_all, allstages_all, supplement_df, year_vec,
      sp_density_list, density_df, dens_index_df, equivalence_list,
      sp_density_num_vec, entry_time_vec, inda_terms_num_vec,
      indb_terms_num_vec, indc_terms_num_vec, inda_terms_cat_vec,
      indb_terms_cat_vec, indc_terms_cat_vec, dens_vr_yn_vec, tweights_type_vec,
      fecmod_vec, patch_vec, variant_count, var_per_run, nreps, times,
      fitness_times, stagecounts, substoch, format_int, preexisting_mpm_size,
      firstage_int, finalage_int, main_optim_res, opt_res, exp_tol, theta_tol,
      loop_max, integeronly_bool, stages_not_equal, stochastic_bool, dens_yn_bool,
      entry_time_vec_use, sparse_bool, historical, pure_leslie, A_only_bool,
      err_check_bool, err_check_extreme, threshold, fitness_table_bool,
      trait_optima_bool, elast_mult, zap_min_bool);
  } else if (funcbased) {
    // Rcout << "invade3 f ";
    invade3_fb_core (Lyapunov, Lyapunov_optim, ESS_Lyapunov, var_run_mat, N_out,
      comm_out, N_out_optim, comm_out_optim, zero_stage_vec_list, trait_axis,
      new_trait_axis, optim_trait_axis, optim_trait_axis_995,
      stageexpansion_list, stageexpansion_list_optim, modified_dev_terms_list,
      errcheck_mpm_list, errcheck_mpm_list_optima, errcheck_mpm_list_ESS,
      errcheck_mpmout_list, errcheck_mpmout_list_optima,
      errcheck_mpmout_list_ESS, tweights_list, start_list, vrm_list,
      final_stageframe, allmodels_all, allstages_all, supplement_df, year_vec,
      ind_terms_num_list, ind_terms_cat_list, dev_terms_list, density_vr_list,
      sp_density_list, density_df, dens_index_df, equivalence_list,
      sp_density_num_vec, entry_time_vec, inda_terms_num_vec,
      indb_terms_num_vec, indc_terms_num_vec, inda_terms_cat_vec,
      indb_terms_cat_vec, indc_terms_cat_vec, dens_vr_yn_vec, tweights_type_vec,
      fecmod_vec, patch_vec, variant_count, var_per_run, nreps, times,
      fitness_times, stagecounts, substoch, format_int, firstage_int,
      finalage_int, dev_terms_times_int, main_optim_res, opt_res, exp_tol,
      theta_tol, loop_max, integeronly_bool, stochastic_bool, dens_yn_bool,
      stages_not_equal, sparse_bool, historical, pure_leslie, A_only_bool,
      err_check_bool, err_check_extreme, threshold, fitness_table_bool,
      trait_optima_bool, elast_mult, zap_min_bool);
  } // funcbased processing
  
  // Rcout << "invade3 g ";
  int ESS_Lyap_nrows = static_cast<int>(ESS_Lyapunov.nrows());
  
  // Zap tiny values
  if (zap_min_bool && ESS_Lyap_nrows > 0) {
    for (int i = 11; i < 14; i++) { 
      NumericVector ESS_Lyapunov_current_var = as<NumericVector>(ESS_Lyapunov(i));
      for (int j = 0; j < ESS_Lyap_nrows; j++) { 
        if (!NumericVector::is_na(ESS_Lyapunov_current_var(j))) {
          if (abs(ESS_Lyapunov_current_var(j)) < threshold) ESS_Lyapunov_current_var(j) = 0.;
        }
      }
    }
    
    for (int i = 16; i < 30; i++) { 
      NumericVector ESS_Lyapunov_current_var = as<NumericVector>(ESS_Lyapunov(i));
      for (int j = 0; j < ESS_Lyap_nrows; j++) { 
        if (!NumericVector::is_na(ESS_Lyapunov_current_var(j))) {
          if (abs(ESS_Lyapunov_current_var(j)) < threshold) ESS_Lyapunov_current_var(j) = 0.;
        }
      }
    }
  }
  
  // Restrict optima presented to converged fitness values
  if (converged_only_bool) {
    StringVector focused_var = {"converged"};
    IntegerVector chosen_value = {1};
    
    if (ESS_Lyap_nrows > 0) ESS_Lyapunov = LefkoUtils::df_subset(ESS_Lyapunov,
      as<RObject>(chosen_value), false, true, true, false, true,
      as<RObject>(focused_var));
  }
  
  List optim_list = Rcpp::List::create(_["ESS_values"] = ESS_Lyapunov,
    _["fitness"] = Lyapunov_optim, _["comm_elas_out"] = comm_out_optim,
    _["N_elas_out"] = N_out_optim);
  
  int out_dim = 8;
  if (trait_optima_bool) out_dim++;
  if (err_check_bool) out_dim++;
  List output (out_dim);
  
  output(0) = Lyapunov; // Needed in final output
  output(1) = comm_out; // Needed in final output
  output(2) = N_out; // Needed in final output
  output(3) = new_trait_axis;
  output(4) = stageframe_df; // Needed in final output
  output(5) = hstages_list;
  output(6) = agestages_list;
  output(7) = labels;
  if (trait_optima_bool) output(8) = optim_list;
  
  if (err_check_bool) {
    if (err_check_extreme) {
      List optim (7);
      
      optim(0) = optim_trait_axis;
      optim(1) = optim_trait_axis_995;
      optim(2) = stageexpansion_list_optim;
      optim(3) = errcheck_mpm_list_optima;
      optim(4) = errcheck_mpmout_list_optima;
      optim(5) = errcheck_mpm_list_ESS;
      optim(6) = errcheck_mpmout_list_ESS;
      
      CharacterVector optim_errcheck_names = {"optim_trait_axis",
        "optim_trait_axis_995", "stage_expansion_list_optim",
        "errcheck_mpm_list_optima", "errcheck_mpmout_list_optima",
        "errcheck_mpm_list_ESS", "errcheck_mpmout_list_ESS"};
      optim.attr("names") = optim_errcheck_names;
      
      List output_errcheck (15);
      
      output_errcheck(0) = allstages_all;
      output_errcheck(1) = allmodels_all;
      output_errcheck(2) = equivalence_list;
      output_errcheck(3) = density_df;
      output_errcheck(4) = dens_index_df;
      output_errcheck(5) = density_vr_list;
      output_errcheck(6) = stageexpansion_list;
      output_errcheck(7) = dev_terms_list;
      output_errcheck(8) = modified_dev_terms_list;
      output_errcheck(9) = errcheck_mpm_list;
      output_errcheck(10) = errcheck_mpmout_list;
      output_errcheck(11) = var_run_mat;
      output_errcheck(12) = start_list;
      output_errcheck(13) = final_stageframe;
      output_errcheck(14) = optim;
      
      CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
        "equivalence_list", "density_df", "dens_index_df", "density_vr_list",
        "stageexpansion_list", "dev_terms_list", "modified_dev_terms_list",
        "modified_mpms", "fb_mpm_out_matrices", "var_run_mat", "start_list",
        "final_stageframe", "optim"};
      output_errcheck.attr("names") = output_errcheck_names;
      
      output(out_dim - 1) = output_errcheck;
      
      CharacterVector output_main_names;
      if (trait_optima_bool) {
        output_main_names = {"fitness", "variants_out", "N_out", "trait_axis",
          "stageframe", "hstages", "agestages", "labels", "optim", "err_check"};
      } else {
        output_main_names = {"fitness", "variants_out", "N_out", "trait_axis",
          "stageframe", "hstages", "agestages", "labels", "err_check"};
      }
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
      output_errcheck(11) = start_list;
      
      CharacterVector output_errcheck_names = {"allstages_all", "allmodels_all",
        "equivalence_list", "density_df", "dens_index_df", "density_vr_list",
        "trait_axis_reassessed", "stageexpansion_list", "dev_terms_list",
        "modified_dev_terms_list", "var_run_mat", "start_list"};
      output_errcheck.attr("names") = output_errcheck_names;
      
      output(out_dim - 1) = output_errcheck;
      
      CharacterVector output_main_names;
      if (trait_optima_bool) {
        output_main_names = {"fitness", "variants_out", "N_out", "trait_axis",
          "stageframe", "hstages", "agestages", "labels", "optim", "err_check"};
      } else {
        output_main_names = {"fitness", "variants_out", "N_out", "trait_axis",
          "stageframe", "hstages", "agestages", "labels", "err_check"};
      }
      output.attr("names") = output_main_names;
    }
  } else {
    CharacterVector output_main_names;
    if (trait_optima_bool) {
      output_main_names = {"fitness", "variants_out", "N_out", "trait_axis",
        "stageframe", "hstages", "agestages", "labels", "optim"};
    } else {
      output_main_names = {"fitness", "variants_out", "N_out", "trait_axis",
        "stageframe", "hstages", "agestages", "labels"};
    }
    output.attr("names") = output_main_names;
  }
  output.attr("class") = "adaptInv";
  
  return output;
}

