#ifndef ADAPTUTILS_main_utils_H
#define ADAPTUTILS_main_utils_H

#include <RcppArmadillo.h>
#define BOOST_DISABLE_ASSERTS

#include <LefkoUtils.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;





// Index of functions
// 1. void Amat_alter Alter MPM A Matrices with Trait Axis Inputs, Standard Matrix Format
// 2. void UFmat_alter Alter MPM U and F Matrices with Trait Axis Inputs, Standard Matrix Format
// 3. void sp_Amat_alter Alter MPM A Matrices with Trait Axis Inputs, Sparse Matrix Format
// 4. void sp_UFmat_alter Alter MPM U and F Matrices with Trait Axis Inputs, Sparse Matrix Format
// 5. arma::mat exp_grid_single  Create Expanded Matrix Giving Permutations with Replacement
// 6. void fastLm_sl  Fast Linear Regression, Slopes Only
// 7. void proj3dens_ad  Project Forward By One Time Step
// 8. void proj3dens_inv  Project Forward By One Time Step in Invastion Run
// 9. void Lyapunov_df_maker  Create Data Frame to Hold Fitness Output from Function invade3()
// 10. DataFrame ta_reassess  Expand Trait Axis Table Given User Input
// 11. void pop_error2  Standardized Error Messages
// 12. DataFrame df_rbind  Bind Two Data Frames By Row
// 13. List df_indices  Subset A Data Frame By Row Index
// 14. void Lyapunov_creator  Creates Final Table of Lyapunov Estimates
// 15. void optim_ta_setup  Create Trait Axis Reassessed for Trait Optimization


namespace AdaptUtils {
  
  //' Alter MPM A Matrices with Trait Axis Inputs, Standard Matrix Format
  //' 
  //' This function alters A matrices in lefkoMat objects with information
  //' supplied in \code{adaptAxis} objects.
  //' 
  //' @name Amat_alter
  //' 
  //' @param Amat The A matrix to modify, in arma::mat format.
  //' @param stageexpansion The input trait_axis data frame, post processing by
  //' function \code{thenewpizzle}.
  //' 
  //' @return This function alters object \code{Amat} by reference, and so does
  //' not return any objects.
  //' 
  //' @keywords internal
  //' @noRd
  inline void Amat_alter (arma::mat& Amat, DataFrame& stageexpansion) {
    
    NumericVector tagiven_t = as<NumericVector>(stageexpansion["tagiven_t"]);
    NumericVector tagiven_f = as<NumericVector>(stageexpansion["tagiven_f"]);
    NumericVector taoffset_t = as<NumericVector>(stageexpansion["taoffset_t"]);
    NumericVector taoffset_f = as<NumericVector>(stageexpansion["taoffset_f"]);
    NumericVector tasurvmult = as<NumericVector>(stageexpansion["tasurvmult"]);
    NumericVector tafecmult = as<NumericVector>(stageexpansion["tafecmult"]);
    IntegerVector aliveequal = as<IntegerVector>(stageexpansion["aliveandequal"]);
    IntegerVector aliveequal_proxy = as<IntegerVector>(stageexpansion["aliveandequal_proxy"]);
    IntegerVector taeststage3 = as<IntegerVector>(stageexpansion["eststage3"]);
    IntegerVector taconvtype = as<IntegerVector>(stageexpansion["taconvtype"]);
    IntegerVector index321 = as<IntegerVector>(stageexpansion["index321"]);
    IntegerVector mpm_altered = as<IntegerVector>(stageexpansion["mpm_altered"]);
    int num_alterations = static_cast<int>(stageexpansion.nrows());
    
    for (int i = 0; i < num_alterations; i++) {
      if (mpm_altered(i) > 0) {
        if(aliveequal(i) > -1) {
          unsigned int properindex = static_cast<unsigned int>(aliveequal(i));
          
          if (!NumericVector::is_na(tagiven_t(i))) {
            if (tagiven_t(i) > -1.) {
              Amat(properindex) = tagiven_t(i);
            }
          }
          if (!NumericVector::is_na(tagiven_f(i))) {
            if (tagiven_f(i) > -1.) {
              Amat(properindex) = tagiven_f(i);
            }
          }
          
          if (!NumericVector::is_na(taoffset_t(i))) {
            if (taoffset_t(i) != 0.) {
              Amat(properindex) = Amat(properindex) + taoffset_t(i);
            }
          }
          if (!NumericVector::is_na(taoffset_f(i))) {
            if (taoffset_f(i) != 0.) {
              Amat(properindex) = Amat(properindex) + taoffset_f(i);
            }
          }
          
          if (!NumericVector::is_na(tasurvmult(i))) {
            if (tasurvmult(i) > -1.) {
              Amat(properindex) = Amat(properindex) * tasurvmult(i);
            }
          }
          if (!NumericVector::is_na(tafecmult(i))) {
            if (tafecmult(i) > -1.) {
              Amat(properindex) = Amat(properindex) * tafecmult(i);
            }
          }
          
          if (!IntegerVector::is_na(taeststage3(i)) && aliveequal_proxy(i) > -1) {
            unsigned int proxyindex = static_cast<unsigned int>(aliveequal_proxy(i));
            if (taeststage3(i) > 0) {
              Amat(properindex) = Amat(proxyindex);
            }
          }
        }
      }
    }
  }

  //' Alter MPM U and F Matrices with Trait Axis Inputs, Standard Matrix Format
  //' 
  //' This function alters U and F matrices in lefkoMat objects with information
  //' supplied in \code{adaptAxis} objects.
  //' 
  //' @name UFmat_alter
  //' 
  //' @param Amat The A matrix to modify, in arma::mat format.
  //' @param Umat The U matrix to modify, in arma::mat format.
  //' @param Fmat The F matrix to modify, in arma::mat format.
  //' @param stageexpansion The input trait_axis data frame, post processing by
  //' function \code{thenewpizzle}.
  //' 
  //' @return This function alters objects \code{Amat}, \code{Umat}, and
  //' \code{Fmat} by reference, and so does not return any objects.
  //' 
  //' @keywords internal
  //' @noRd
  inline void UFmat_alter (arma::mat& Amat, arma::mat& Umat, arma::mat& Fmat,
    DataFrame& stageexpansion) {
    
    NumericVector tagiven_t = as<NumericVector>(stageexpansion["tagiven_t"]);
    NumericVector tagiven_f = as<NumericVector>(stageexpansion["tagiven_f"]);
    NumericVector taoffset_t = as<NumericVector>(stageexpansion["taoffset_t"]);
    NumericVector taoffset_f = as<NumericVector>(stageexpansion["taoffset_f"]);
    NumericVector tasurvmult = as<NumericVector>(stageexpansion["tasurvmult"]);
    NumericVector tafecmult = as<NumericVector>(stageexpansion["tafecmult"]);
    IntegerVector aliveequal = as<IntegerVector>(stageexpansion["aliveandequal"]);
    IntegerVector aliveequal_proxy = as<IntegerVector>(stageexpansion["aliveandequal_proxy"]);
    IntegerVector taeststage3 = as<IntegerVector>(stageexpansion["eststage3"]);
    IntegerVector taconvtype = as<IntegerVector>(stageexpansion["taconvtype"]);
    IntegerVector index321 = as<IntegerVector>(stageexpansion["index321"]);
    IntegerVector mpm_altered = as<IntegerVector>(stageexpansion["mpm_altered"]);
    int num_alterations = static_cast<int>(stageexpansion.nrows());
    
    for (int i = 0; i < num_alterations; i++) {
      if (mpm_altered(i) > 0) {
        if(aliveequal(i) > -1) {
          unsigned int properindex = static_cast<unsigned int>(aliveequal(i));
          
          if (!NumericVector::is_na(tagiven_t(i))) {
            if (tagiven_t(i) > -1.) {
              Umat(properindex) = tagiven_t(i);
            }
          }
          if (!NumericVector::is_na(tagiven_f(i))) {
            if (tagiven_f(i) > -1.) {
              Fmat(properindex) = tagiven_f(i);
            }
          }
          
          if (!NumericVector::is_na(taoffset_t(i))) {
            if (taoffset_t(i) != 0.) {
              Umat(properindex) = Umat(properindex) + taoffset_t(i);
            }
          }
          if (!NumericVector::is_na(taoffset_f(i))) {
            if (taoffset_f(i) != 0.) {
              Fmat(properindex) = Fmat(properindex) + taoffset_f(i);
            }
          }
          
          if (!NumericVector::is_na(tasurvmult(i))) {
            if (tasurvmult(i) > -1.) {
              Umat(properindex) = Umat(properindex) * tasurvmult(i);
            }
          }
          if (!NumericVector::is_na(tafecmult(i))) {
            if (tafecmult(i) > -1.) {
              Fmat(properindex) = Fmat(properindex) * tafecmult(i);
            }
          }
          
          if (!IntegerVector::is_na(taeststage3(i)) && aliveequal_proxy(i) > -1) {
            unsigned int proxyindex = static_cast<unsigned int>(aliveequal_proxy(i));
            if (taeststage3(i) > 0) {
              if (taconvtype(i) == 1) {
                Umat(properindex) = Umat(proxyindex);
              } else {
                Fmat(properindex) = Fmat(proxyindex);
              }
            }
          }
        }
      }
    }
    arma::mat A_final = Umat + Fmat;
    Amat = A_final;
  }

  //' Alter MPM A Matrices with Trait Axis Inputs, Sparse Matrix Format
  //' 
  //' This function alters A matrices in lefkoMat objects with information
  //' supplied in \code{adaptAxis} objects.
  //' 
  //' @name sp_Amat_alter
  //' 
  //' @param Amat The A matrix to modify, in arma::sp_mat format.
  //' @param stageexpansion The input trait_axis data frame, post processing by
  //' function \code{thenewpizzle}.
  //' 
  //' @return This function alters object \code{Amat} by reference, and so does
  //' not return any objects.
  //' 
  //' @keywords internal
  //' @noRd
  inline void sp_Amat_alter (arma::sp_mat& Amat, DataFrame& stageexpansion) {
    
    NumericVector tagiven_t = as<NumericVector>(stageexpansion["tagiven_t"]);
    NumericVector tagiven_f = as<NumericVector>(stageexpansion["tagiven_f"]);
    NumericVector taoffset_t = as<NumericVector>(stageexpansion["taoffset_t"]);
    NumericVector taoffset_f = as<NumericVector>(stageexpansion["taoffset_f"]);
    NumericVector tasurvmult = as<NumericVector>(stageexpansion["tasurvmult"]);
    NumericVector tafecmult = as<NumericVector>(stageexpansion["tafecmult"]);
    IntegerVector aliveequal = as<IntegerVector>(stageexpansion["aliveandequal"]);
    IntegerVector aliveequal_proxy = as<IntegerVector>(stageexpansion["aliveandequal_proxy"]);
    IntegerVector taeststage3 = as<IntegerVector>(stageexpansion["eststage3"]);
    IntegerVector taconvtype = as<IntegerVector>(stageexpansion["taconvtype"]);
    IntegerVector index321 = as<IntegerVector>(stageexpansion["index321"]);
    IntegerVector mpm_altered = as<IntegerVector>(stageexpansion["mpm_altered"]);
    int num_alterations = static_cast<int>(stageexpansion.nrows());
    
    for (int i = 0; i < num_alterations; i++) {
      if (mpm_altered(i) > 0) {
        if(aliveequal(i) > -1) {
          unsigned int properindex = static_cast<unsigned int>(aliveequal(i));
          
          if (!NumericVector::is_na(tagiven_t(i))) {
            if (tagiven_t(i) > -1.) {
              Amat(properindex) = tagiven_t(i);
            }
          }
          if (!NumericVector::is_na(tagiven_f(i))) {
            if (tagiven_f(i) > -1.) {
              Amat(properindex) = tagiven_f(i);
            }
          }
          
          if (!NumericVector::is_na(taoffset_t(i))) {
            if (taoffset_t(i) != 0.) {
              Amat(properindex) = Amat(properindex) + taoffset_t(i);
            }
          }
          if (!NumericVector::is_na(taoffset_f(i))) {
            if (taoffset_f(i) != 0.) {
              Amat(properindex) = Amat(properindex) + taoffset_f(i);
            }
          }
          
          if (!NumericVector::is_na(tasurvmult(i))) {
            if (tasurvmult(i) > -1.) {
              Amat(properindex) = Amat(properindex) * tasurvmult(i);
            }
          }
          if (!NumericVector::is_na(tafecmult(i))) {
            if (tafecmult(i) > -1.) {
              Amat(properindex) = Amat(properindex) * tafecmult(i);
            }
          }
          
          if (!IntegerVector::is_na(taeststage3(i)) && aliveequal_proxy(i) > -1) {
            unsigned int proxyindex = static_cast<unsigned int>(aliveequal_proxy(i));
            if (taeststage3(i) > 0) {
              Amat(properindex) = Amat(proxyindex);
            }
          }
        }
      }
    }
  }
  
  //' Alter MPM U and F Matrices with Trait Axis Inputs, Sparse Matrix Format
  //' 
  //' This function alters U and F matrices in lefkoMat objects with information
  //' supplied in \code{adaptAxis} objects.
  //' 
  //' @name sp_UFmat_alter
  //' 
  //' @param Amat The A matrix to modify, in arma::sp_mat format.
  //' @param Umat The U matrix to modify, in arma::sp_mat format.
  //' @param Fmat The F matrix to modify, in arma::sp_mat format.
  //' @param stageexpansion The input trait_axis data frame, post processing by
  //' function \code{thenewpizzle}.
  //' 
  //' @return This function alters objects \code{Amat}, \code{Umat}, and
  //' \code{Fmat} by reference, and so does not return any objects.
  //' 
  //' @keywords internal
  //' @noRd
  inline void sp_UFmat_alter (arma::sp_mat& Amat, arma::sp_mat& Umat, arma::sp_mat& Fmat,
    DataFrame& stageexpansion) {
    
    NumericVector tagiven_t = as<NumericVector>(stageexpansion["tagiven_t"]);
    NumericVector tagiven_f = as<NumericVector>(stageexpansion["tagiven_f"]);
    NumericVector taoffset_t = as<NumericVector>(stageexpansion["taoffset_t"]);
    NumericVector taoffset_f = as<NumericVector>(stageexpansion["taoffset_f"]);
    NumericVector tasurvmult = as<NumericVector>(stageexpansion["tasurvmult"]);
    NumericVector tafecmult = as<NumericVector>(stageexpansion["tafecmult"]);
    IntegerVector aliveequal = as<IntegerVector>(stageexpansion["aliveandequal"]);
    IntegerVector aliveequal_proxy = as<IntegerVector>(stageexpansion["aliveandequal_proxy"]);
    IntegerVector taeststage3 = as<IntegerVector>(stageexpansion["eststage3"]);
    IntegerVector taconvtype = as<IntegerVector>(stageexpansion["taconvtype"]);
    IntegerVector index321 = as<IntegerVector>(stageexpansion["index321"]);
    IntegerVector mpm_altered = as<IntegerVector>(stageexpansion["mpm_altered"]);
    int num_alterations = static_cast<int>(stageexpansion.nrows());
    
    for (int i = 0; i < num_alterations; i++) {
      if (mpm_altered(i) > 0) {
        if(aliveequal(i) > -1) {
          unsigned int properindex = static_cast<unsigned int>(aliveequal(i));
          
          if (!NumericVector::is_na(tagiven_t(i))) {
            if (tagiven_t(i) > -1.) {
              Umat(properindex) = tagiven_t(i);
            }
          }
          if (!NumericVector::is_na(tagiven_f(i))) {
            if (tagiven_f(i) > -1.) {
              Fmat(properindex) = tagiven_f(i);
            }
          }
          
          if (!NumericVector::is_na(taoffset_t(i))) {
            if (taoffset_t(i) != 0.) {
              Umat(properindex) = Umat(properindex) + taoffset_t(i);
            }
          }
          if (!NumericVector::is_na(taoffset_f(i))) {
            if (taoffset_f(i) != 0.) {
              Fmat(properindex) = Fmat(properindex) + taoffset_f(i);
            }
          }
          
          if (!NumericVector::is_na(tasurvmult(i))) {
            if (tasurvmult(i) > -1.) {
              Umat(properindex) = Umat(properindex) * tasurvmult(i);
            }
          }
          if (!NumericVector::is_na(tafecmult(i))) {
            if (tafecmult(i) > -1.) {
              Fmat(properindex) = Fmat(properindex) * tafecmult(i);
            }
          }
          
          if (!IntegerVector::is_na(taeststage3(i)) && aliveequal_proxy(i) > -1) {
            unsigned int proxyindex = static_cast<unsigned int>(aliveequal_proxy(i));
            if (taeststage3(i) > 0) {
              if (taconvtype(i) == 1) {
                Umat(properindex) = Umat(proxyindex);
              } else {
                Fmat(properindex) = Fmat(proxyindex);
              }
            }
          }
        }
      }
    }
    Amat = Umat + Fmat;
  }
  
  //' Create Expanded Matrix Giving Permutations with Replacement
  //' 
  //' @name exp_grid_single
  //' 
  //' Function \code{exp_grid_single()} creates a matrix in which a series of
  //' indices populate the columns. A vector is created of sequential integer
  //' indices starting from 0 and going to the length given in the first
  //' argument. All permutations with replacement of the number of these indices
  //' given in the second argument are then created, laid out by row in the
  //' resulting matrix.
  //' 
  //' @param num_indices An integer giving the numer of indices. Results in an
  //' integer vector of indices, starting with 0.
  //' @param num_cols An integer giving the number of columns to prepagate the
  //' values in the vector generated.
  //' 
  //' @return An arma::mat matrix with \code{num_cols} columns and number of rows
  //' equal to \code{num_indices} to the power of \code{num_cols}.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::mat exp_grid_single (int num_indices, int num_cols) {
    IntegerVector int_limits (num_cols);
    IntegerVector int_minlimits (num_cols);
    double exp_permutes = pow(num_indices, num_cols);
    
    arma::mat exp_mat = arma::mat(static_cast<int>(exp_permutes), num_cols); // Rows = runs, Cols = variants
    
    int var_counter {0};
    int base_num {0};
    
    for (int i = 0; i < num_cols; i++) {
      int_minlimits(i) = static_cast<int>(pow(num_indices, (i)));
      int_limits(i) = static_cast<int>(pow(num_indices, (i+1)));
      
      for (int j = 0; j < static_cast<int>(exp_permutes); j++) {
        if (var_counter >= int_limits(i) || j == 0) var_counter = 0;
        base_num = floor(var_counter / int_minlimits(i));
        if (base_num >= int_limits(i)) {
          var_counter = 0;
          base_num = floor(var_counter / int_minlimits(i));
        }
        exp_mat(j, i) = base_num;
        var_counter++;
      }
    }
    return exp_mat;
  }

  //' Fast Linear Regression, Slopes Only
  //' 
  //' This function performs a simple linear / multiple regression quickly and
  //' accurately using RcppArmadillo. Based on code by Dirk Eddelbuettel,
  //' derived from 
  //' https://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
  //' 
  //' @name fastLm_sl
  //' 
  //' @param slopes A column vector to hold the slope coefficients.
  //' @param y A column vector holding the values of the response variable, in
  //' order.
  //' @param x A matrix holding the values of the independent variables, with the
  //' first column corresponding to the y-intercept.
  //' 
  //' @return A vector named holding the values of the slope coefficients in the
  //' order of the columns in the \code{x} matrix.
  //' 
  //' @keywords internal
  //' @noRd
  inline void fastLm_sl(arma::colvec& slopes, const arma::vec & y,
    const arma::mat & X) {
    
    arma::colvec coef = arma::solve(X, y);
    slopes = coef;
  }
  
  //' Project Forward By One Time Step
  //' 
  //' This function projects the community forward by one time step.
  //' 
  //' @name proj3dens_ad
  //' 
  //' @param proj_vec The output projected vector.
  //' @param prophesized_mat The projection matrix modified for density
  //' dependence, in standard matrix format.
  //' @param prophesized_sp The projection matrix modified for density
  //' dependence, in sparse matrix format.
  //' @param start_vec The starting population vector for the projection.
  //' @param core_list A list of full projection matrices, corresponding to the 
  //' \code{A} list within a \code{lefkoMat} object.
  //' @param the_foreseen The C++ integer index of the appropriate A matrix.
  //' @param delay_N The total number of individuals across all populations, after
  //' accounting for equivalency of individuals in different populations. Should
  //' be a single double-precision floating point value.
  //' @param growthonly A logical value stating whether to output only a matrix
  //' showing the change in population size from one year to the next for use in
  //' stochastic population growth rate estimation (TRUE), or a larger matrix also
  //' containing the w and v projections for stochastic perturbation analysis,
  //' stage distribution estimation, and reproductive value estimation.
  //' @param integeronly A logical value indicating whether to round all projected
  //' numbers of individuals to the nearest integer.
  //' @param substoch An integer value indicating whether to force survival-
  //' transition matrices to be substochastic in density dependent simulations.
  //' Defaults to \code{0}, which does not force substochasticity. Alternatively,
  //' \code{1} forces all survival-transition elements to range from 0.0 to 1.0
  //' and fecundity to be non-negative, and \code{2} forces all column rows to
  //' total no more than 1.0.
  //' @param dens_input The original \code{lefkoDens} data frame supplied through
  //' the \code{\link{density_input}()} function.
  //' @param dens_index A list giving the indices of elements in object
  //' \code{dens_input}.
  //' @param sparse_auto A logical value indicating whether to determine whether
  //' to use sparse matrix encoding automatically.
  //' @param sparse A logical value indicating whether to use sparse matrix
  //' encoding if \code{sparse_auto = FALSE}.
  //' @param sparse_input A logical value indicating whether matrices in the input
  //' MPM are in sparse format (class \code{dgCMatrix}). If so, then all
  //' projection will be handled in sparse format. Defaults to \code{FALSE}.
  //' @param allow_warnings A logical value indicating whether the function should
  //' send warnings if estimated values fall outside of the realm of possibility.
  //' @param err_check A Boolean value indicating whether to store projection
  //' matrices themselves in memory.
  //' 
  //' @return A one-column matrix in which showing the number of individuals in
  //' each stage across time.
  //' 
  //' @section Notes:
  //' There is no option to standardize population vectors here, because density
  //' dependence requires the full population size to be tracked.
  //' 
  //' @keywords internal
  //' @noRd
  inline void proj3dens_ad(arma::vec& proj_vec, arma::mat& prophesized_mat,
    arma::sp_mat& prophesized_sp, const arma::vec& start_vec,
    const List& core_list, const double delay_N, const int the_foreseen,
    bool integeronly, int substoch, const Rcpp::DataFrame& dens_input,
    const Rcpp::List& dens_index, bool sparse_auto, bool sparse,
    bool sparse_input = false, bool allow_warnings = false,
    bool err_check = false) {
    
    int sparse_switch {0};
    int time_delay {1};
    double pop_size {0};
    bool warn_trigger_neg = false;
    bool warn_trigger_1 = false;
    
    int nostages = static_cast<int>(start_vec.n_elem);
    arma::vec theseventhson;
    arma::mat theprophecy;
    arma::sp_mat theprophecy_sp;
    
    // Density dependence
    arma::uvec dyn_index321 = as<arma::uvec>(dens_index["index321"]);
    arma::uvec dyn_index_col = as<arma::uvec>(dens_index[1]);
    arma::uvec dyn_style = as<arma::uvec>(dens_input["style"]);
    arma::vec dyn_alpha = as<arma::vec>(dens_input["alpha"]);
    arma::vec dyn_beta = as<arma::vec>(dens_input["beta"]);
    arma::uvec dyn_delay = as<arma::uvec>(dens_input["time_delay"]);
    arma::uvec dyn_type = as<arma::uvec>(dens_input["type"]);
    int n_dyn_elems = static_cast<int>(dyn_index321.n_elem);
    
    arma::vec popproj(nostages, fill::zeros); // Population vector
    theseventhson = start_vec;
    
    // Check if matrix is large and sparse
    if (sparse_auto && !sparse_input) {
      int test_elems = static_cast<int>(as<arma::mat>(core_list(0)).n_elem);
      arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
      int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
      double sparse_check = static_cast<double>(all_nonzeros) /
        static_cast<double>(test_elems);
      if (sparse_check <= 0.5 && theseventhson.n_elem > 50) {
        sparse_switch = 1;
      } else sparse_switch = 0;
    } else sparse_switch = sparse;
    
    double changing_element {0.0};
    
    if (sparse_switch == 0 && !sparse_input) {
      // Dense matrix projection
      theprophecy = as<arma::mat>(core_list[(the_foreseen)]);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        pop_size = delay_N;
        
        if (dyn_style(j) == 1) { // Ricker
          changing_element = theprophecy(dyn_index321(j)) * 
            dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
          
        } else if (dyn_style(j) == 2) { // Beverton-Holt
          changing_element = theprophecy(dyn_index321(j)) * 
            dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
          
        } else if (dyn_style(j) == 3) { // Usher function
          changing_element = theprophecy(dyn_index321(j)) * 
            (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
          
        } else if (dyn_style(j) == 4) { // Logistic function
          double used_popsize = pop_size;
          if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
            used_popsize = dyn_alpha(j);
          }
          changing_element = theprophecy(dyn_index321(j)) * 
            (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
        }
        
        if (substoch == 1 && dyn_type(j) == 1) {
          if (changing_element > 1.0) {
            changing_element = 1.0;
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch == 2 && dyn_type(j) == 1) {
          arma::vec given_col = theprophecy.col(dyn_index_col(j));
          arma::uvec gc_negs = find(given_col < 0.0);
          arma::uvec gc_pos = find(given_col > 0.0);
          
          double barnyard_antics = sum(given_col(gc_pos)) -
            theprophecy(dyn_index321(j)) + changing_element;
          
          if (barnyard_antics > 1.0 && changing_element > 0.0) {
            double proposed_element = changing_element - barnyard_antics *
              (changing_element / barnyard_antics);
            
            if (proposed_element >= 0.0) {
              changing_element = proposed_element;
            } else {
              changing_element = 0.0;
            }
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch > 0 && dyn_type(j) == 2) {
          if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        }
        theprophecy(dyn_index321(j)) = changing_element;
        
        if (allow_warnings) {
          if (dyn_type(j) == 1 && theprophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
            warn_trigger_1 = true;
            Rf_warningcall(R_NilValue,
              "Some probabilities with value > 1.0 produced during density adjustment.");
            
          } else if (theprophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
            warn_trigger_neg = true;
            Rf_warningcall(R_NilValue,
              "Some matrix elements with value < 0.0 produced during density adjustment.");
          }
        }
      }
      theseventhson = theprophecy * theseventhson;
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      if (err_check) prophesized_mat = theprophecy;
      
      proj_vec = theseventhson;
      
    } else {
      // Sparse matrix projection
      arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
      int matlist_length = static_cast<int>(core_list.size());
      Rcpp::List sparse_list(matlist_length);
      
      if (!sparse_input) {
        arma::mat first_mat = core_list(0);
        arma::sp_mat new_sparse = arma::sp_mat(first_mat);
        sparse_list(0) = new_sparse;
        
        if(matlist_length > 1) {
          for (int i = 1; i < matlist_length; i++) {
            first_mat = as<arma::mat>(core_list(i));
            new_sparse = arma::sp_mat(first_mat);
            sparse_list(i) = new_sparse;
          }
        }
      } else sparse_list = core_list;
      
      arma::sp_mat sparse_prophecy;
      sparse_prophecy = as<arma::sp_mat>(sparse_list[the_foreseen]);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        pop_size = delay_N;
        
        if (dyn_style(j) == 1) { // Ricker
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
          
        } else if (dyn_style(j) == 2) { // Beverton-Holt
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
          
        } else if (dyn_style(j) == 3) { // Usher function
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
          
        } else if (dyn_style(j) == 4) { // Logistic function
          double used_popsize = pop_size;
          if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
            used_popsize = dyn_alpha(j);
          }
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
        }
        
        if (substoch == 1 && dyn_type(j) == 1) {
          if (changing_element > 1.0) {
            changing_element = 1.0;
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch == 2 && dyn_type(j) == 1) {
          arma::vec given_col = arma::vec(sparse_prophecy.col(dyn_index_col(j)));
          arma::uvec gc_negs = find(given_col < 0.0);
          arma::uvec gc_pos = find(given_col > 0.0);
          
          double barnyard_antics = sum(given_col(gc_pos)) -
            sparse_prophecy(dyn_index321(j)) + changing_element;
          
          if (barnyard_antics > 1.0 && changing_element > 0.0) {
            double proposed_element = changing_element - barnyard_antics *
              (changing_element / barnyard_antics);
            
            if (proposed_element >= 0.0) {
              changing_element = proposed_element;
            } else {
              changing_element = 0.0;
            }
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch > 0 && dyn_type(j) == 2) {
          if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        }
        sparse_prophecy(dyn_index321(j)) = changing_element;
        
        if (allow_warnings) {
          if (dyn_type(j) == 1 && sparse_prophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
            warn_trigger_1 = true;
            Rf_warningcall(R_NilValue,
              "Some probabilities with value > 1.0 produced during density adjustment.");
              
          } else if (sparse_prophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
            warn_trigger_neg = true;
            Rf_warningcall(R_NilValue,
              "Some matrix elements with value < 0.0 produced during density adjustment.");
          }
        }
      }
      
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      if (err_check) prophesized_sp = sparse_prophecy;
      
      popproj = arma::vec(arma::mat(sparse_seventhson));
      proj_vec = popproj;
    }
  }
  
  //' Project Forward By One Time Step in Invastion Run
  //' 
  //' This function projects the community forward by one time step.
  //' 
  //' @name proj3dens_inv
  //' 
  //' @param proj_vec The output projected vector.
  //' @param prophesized_mat The projection matrix modified for density
  //' dependence, in standard matrix format.
  //' @param prophesized_sp The projection matrix modified for density
  //' dependence, in sparse matrix format.
  //' @param start_vec The starting population vector for the projection.
  //' @param core_list A list of full projection matrices, corresponding to the 
  //' \code{A} list within a \code{lefkoMat} object.
  //' @param the_foreseen The C++ integer index of the appropriate A matrix.
  //' @param delay_N The total number of individuals across all populations, after
  //' accounting for equivalency of individuals in different populations. Should
  //' be a single double-precision floating point value.
  //' @param growthonly A logical value stating whether to output only a matrix
  //' showing the change in population size from one year to the next for use in
  //' stochastic population growth rate estimation (TRUE), or a larger matrix also
  //' containing the w and v projections for stochastic perturbation analysis,
  //' stage distribution estimation, and reproductive value estimation.
  //' @param integeronly A logical value indicating whether to round all projected
  //' numbers of individuals to the nearest integer.
  //' @param substoch An integer value indicating whether to force survival-
  //' transition matrices to be substochastic in density dependent simulations.
  //' Defaults to \code{0}, which does not force substochasticity. Alternatively,
  //' \code{1} forces all survival-transition elements to range from 0.0 to 1.0
  //' and fecundity to be non-negative, and \code{2} forces all column rows to
  //' total no more than 1.0.
  //' @param dens_input The original \code{lefkoDens} data frame supplied through
  //' the \code{\link{density_input}()} function.
  //' @param dens_index A list giving the indices of elements in object
  //' \code{dens_input}.
  //' @param sparse_auto A logical value indicating whether to determine whether
  //' to use sparse matrix encoding automatically.
  //' @param sparse A logical value indicating whether to use sparse matrix
  //' encoding if \code{sparse_auto = FALSE}.
  //' @param sparse_input A logical value indicating whether matrices in the input
  //' MPM are in sparse format (class \code{dgCMatrix}). If so, then all
  //' projection will be handled in sparse format. Defaults to \code{FALSE}.
  //' @param allow_warnings A logical value indicating whether the function should
  //' send warnings if estimated values fall outside of the realm of possibility.
  //' @param err_check A Boolean value indicating whether to store projection
  //' matrices themselves in memory.
  //' 
  //' @return A one-column matrix in which showing the number of individuals in
  //' each stage across time.
  //' 
  //' @section Notes:
  //' There is no option to standardize population vectors here, because density
  //' dependence requires the full population size to be tracked.
  //' 
  //' @keywords internal
  //' @noRd
  inline void proj3dens_inv(arma::vec& proj_vec, arma::mat& prophesized_mat,
    arma::sp_mat& prophesized_sp, const arma::vec& start_vec,
    DataFrame& sge_current, const List& core_list, const double delay_N,
    const int the_foreseen,
    bool integeronly, int substoch, const Rcpp::DataFrame& dens_input,
    const Rcpp::List& dens_index, bool sparse_auto, bool sparse,
    bool sparse_input = false, bool allow_warnings = false,
    bool err_check = false) {
    
    int sparse_switch {0};
    int time_delay {1};
    double pop_size {0};
    bool warn_trigger_neg = false;
    bool warn_trigger_1 = false;
    
    int nostages = static_cast<int>(start_vec.n_elem);
    arma::vec theseventhson;
    arma::mat theprophecy;
    arma::sp_mat theprophecy_sp;
    
    // Density dependence
    arma::uvec dyn_index321 = as<arma::uvec>(dens_index["index321"]);
    arma::uvec dyn_index_col = as<arma::uvec>(dens_index[1]);
    arma::uvec dyn_style = as<arma::uvec>(dens_input["style"]);
    arma::vec dyn_alpha = as<arma::vec>(dens_input["alpha"]);
    arma::vec dyn_beta = as<arma::vec>(dens_input["beta"]);
    arma::uvec dyn_delay = as<arma::uvec>(dens_input["time_delay"]);
    arma::uvec dyn_type = as<arma::uvec>(dens_input["type"]);
    int n_dyn_elems = static_cast<int>(dyn_index321.n_elem);
    
    arma::vec popproj(nostages, fill::zeros); // Population vector
    theseventhson = start_vec;
    
    // Check if matrix is large and sparse
    if (sparse_auto && !sparse_input) {
      int test_elems = static_cast<int>(as<arma::mat>(core_list(0)).n_elem);
      arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
      int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
      double sparse_check = static_cast<double>(all_nonzeros) /
        static_cast<double>(test_elems);
      if (sparse_check <= 0.5 && theseventhson.n_elem > 50) {
        sparse_switch = 1;
      } else sparse_switch = 0;
    } else sparse_switch = sparse;
    
    double changing_element {0.0};
    
    if (sparse_switch == 0 && !sparse_input) {
      // Dense matrix projection
      theprophecy = as<arma::mat>(core_list[(the_foreseen)]);
      Amat_alter(theprophecy, sge_current);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        pop_size = delay_N;
        
        if (dyn_style(j) == 1) { // Ricker
          changing_element = theprophecy(dyn_index321(j)) * 
            dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
          
        } else if (dyn_style(j) == 2) { // Beverton-Holt
          changing_element = theprophecy(dyn_index321(j)) * 
            dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
          
        } else if (dyn_style(j) == 3) { // Usher function
          changing_element = theprophecy(dyn_index321(j)) * 
            (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
          
        } else if (dyn_style(j) == 4) { // Logistic function
          double used_popsize = pop_size;
          if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
            used_popsize = dyn_alpha(j);
          }
          changing_element = theprophecy(dyn_index321(j)) * 
            (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
        }
        
        if (substoch == 1 && dyn_type(j) == 1) {
          if (changing_element > 1.0) {
            changing_element = 1.0;
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch == 2 && dyn_type(j) == 1) {
          arma::vec given_col = theprophecy.col(dyn_index_col(j));
          arma::uvec gc_negs = find(given_col < 0.0);
          arma::uvec gc_pos = find(given_col > 0.0);
          
          double barnyard_antics = sum(given_col(gc_pos)) -
            theprophecy(dyn_index321(j)) + changing_element;
          
          if (barnyard_antics > 1.0 && changing_element > 0.0) {
            double proposed_element = changing_element - barnyard_antics *
              (changing_element / barnyard_antics);
            
            if (proposed_element >= 0.0) {
              changing_element = proposed_element;
            } else {
              changing_element = 0.0;
            }
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch > 0 && dyn_type(j) == 2) {
          if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        }
        theprophecy(dyn_index321(j)) = changing_element;
        
        if (allow_warnings) {
          if (dyn_type(j) == 1 && theprophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
            warn_trigger_1 = true;
            Rf_warningcall(R_NilValue,
              "Some probabilities with value > 1.0 produced during density adjustment.");
            
          } else if (theprophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
            warn_trigger_neg = true;
            Rf_warningcall(R_NilValue,
              "Some matrix elements with value < 0.0 produced during density adjustment.");
          }
        }
      }
      theseventhson = theprophecy * theseventhson;
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      if (err_check) prophesized_mat = theprophecy;
      
      proj_vec = theseventhson;
      
    } else {
      // Sparse matrix projection
      arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
      int matlist_length = static_cast<int>(core_list.size());
      Rcpp::List sparse_list(matlist_length);
      
      if (!sparse_input) {
        arma::mat first_mat = core_list(0);
        arma::sp_mat new_sparse = arma::sp_mat(first_mat);
        sparse_list(0) = new_sparse;
        
        if(matlist_length > 1) {
          for (int i = 1; i < matlist_length; i++) {
            first_mat = as<arma::mat>(core_list(i));
            new_sparse = arma::sp_mat(first_mat);
            sparse_list(i) = new_sparse;
          }
        }
      } else sparse_list = core_list;
      
      arma::sp_mat sparse_prophecy;
      sparse_prophecy = as<arma::sp_mat>(sparse_list[the_foreseen]);
      sp_Amat_alter(sparse_prophecy, sge_current);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        pop_size = delay_N;
        
        if (dyn_style(j) == 1) { // Ricker
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
          
        } else if (dyn_style(j) == 2) { // Beverton-Holt
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
          
        } else if (dyn_style(j) == 3) { // Usher function
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
          
        } else if (dyn_style(j) == 4) { // Logistic function
          double used_popsize = pop_size;
          if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
            used_popsize = dyn_alpha(j);
          }
          changing_element = sparse_prophecy(dyn_index321(j)) * 
            (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
        }
        
        if (substoch == 1 && dyn_type(j) == 1) {
          if (changing_element > 1.0) {
            changing_element = 1.0;
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch == 2 && dyn_type(j) == 1) {
          arma::vec given_col = arma::vec(sparse_prophecy.col(dyn_index_col(j)));
          arma::uvec gc_negs = find(given_col < 0.0);
          arma::uvec gc_pos = find(given_col > 0.0);
          
          double barnyard_antics = sum(given_col(gc_pos)) -
            sparse_prophecy(dyn_index321(j)) + changing_element;
          
          if (barnyard_antics > 1.0 && changing_element > 0.0) {
            double proposed_element = changing_element - barnyard_antics *
              (changing_element / barnyard_antics);
            
            if (proposed_element >= 0.0) {
              changing_element = proposed_element;
            } else {
              changing_element = 0.0;
            }
          } else if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        } else if (substoch > 0 && dyn_type(j) == 2) {
          if (changing_element < 0.0) {
            changing_element = 0.0;
          }
        }
        sparse_prophecy(dyn_index321(j)) = changing_element;
        
        if (allow_warnings) {
          if (dyn_type(j) == 1 && sparse_prophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
            warn_trigger_1 = true;
            Rf_warningcall(R_NilValue,
              "Some probabilities with value > 1.0 produced during density adjustment.");
              
          } else if (sparse_prophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
            warn_trigger_neg = true;
            Rf_warningcall(R_NilValue,
              "Some matrix elements with value < 0.0 produced during density adjustment.");
          }
        }
      }
      
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      if (err_check) prophesized_sp = sparse_prophecy;
      
      popproj = arma::vec(arma::mat(sparse_seventhson));
      proj_vec = popproj;
    }
  }
  
  //' Create Data Frame to Hold Fitness Output from Function invade3()
  //' 
  //' This function performs a simple linear / multiple regression quickly and
  //' accurately using RcppArmadillo. Based on code by Dirk Eddelbuettel,
  //' derived from 
  //' https://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
  //' 
  //' @name Lyapunov_df_maker
  //' 
  //' @param Lyapunov An empty data frame to be modified.
  //' @param var_run_mat A matrix holding the index order of variants to run,
  //' with columns specifying variants and rows specifying combinations.
  //' @param var_per_run An integer specifying the number of variants to run in
  //' each simulation.
  //' @param axis_variants_unique An integer vector giving the exact integers
  //' that serve as names of variants.
  //' 
  //' @return A data frame with an initial row number vector, followed by
  //' \code{var_per_run} integer vectors specifying the variants run in each
  //' simulation, followed by \code{var_per_run} empty numeric vectors to hold
  //' estimated Lyapunov coefficients.
  //' 
  //' @keywords internal
  //' @noRd
  inline void Lyapunov_df_maker(DataFrame& Lyapunov, const arma::mat& var_run_mat,
    const int var_per_run, const IntegerVector& axis_variants_unique) {
    
    int var_length = static_cast<int>(var_run_mat.n_rows);
    IntegerVector Lyap_rows = seq(1, var_length);
    CharacterVector df_var_names ((2 * var_per_run) + 1);
    
    df_var_names(0) = "simulation_num";
    String var_name_base = "variant";
    String var_fit_base = "fitness_variant";
    
    List output ((2 * var_per_run) + 1);
    output(0) = Lyap_rows;
    
    for (int i = 0; i < var_per_run; i++) {
      String current_call = var_name_base;
      current_call += (i+1);
      df_var_names(1 + i) = current_call;
      
      String next_call = var_fit_base;
      next_call += (i+1);
      df_var_names(1 + var_per_run + i) = next_call;
    }
    
    for (int i = 0; i < var_per_run; i++) {
      IntegerVector Lyap_index_current (var_length);
      NumericVector Lyap_fitness_current (var_length);
      
      for (int j = 0; j < var_length; j++) {
        Lyap_index_current(j) = axis_variants_unique(var_run_mat(j, i));
      }
      
      output(1 + i) = Lyap_index_current;
      output(1 + var_per_run + i) = Lyap_fitness_current;
    }
    
    output.attr("names") = df_var_names;
    output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, var_length);
    StringVector needed_class {"data.frame"};
    output.attr("class") = needed_class;
    
    Lyapunov = output;
  }
  
  //' Expand Trait Axis Table Given User Input
  //' 
  //' The function takes a supplemental table as input and produces an edited
  //' and expanded version for calculation.
  //' 
  //' @name ta_reassess
  //' 
  //' @param stageframe The stageframe used for assessment.
  //' @param trait_axis The input trait_axis data frame.
  //' @param first_age_int The first age in the MPM, if using Leslie MPMs.
  //' @param historical A Boolean value indicating whether the MPM should be
  //' historical (\code{TRUE}) or not (\code{FALSE}).
  //' @param functionbased A Boolean value indicating whether the MPM is
  //' function-based.
  //' @param pure_leslie A Boolean value indicating that the MPM is a Leslie
  //' MPM.
  //' 
  //' @return This function returns a new data frame that acts as the expanded
  //' supplemental table without shorthand codes. It uses only stage
  //' designations from the stageframe used.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::DataFrame ta_reassess (DataFrame& stageframe,
    DataFrame& trait_axis, const int first_age_int, const bool historical,
    const bool functionbased, const bool pure_leslie) {
    
    IntegerVector variant_ta;
    StringVector stage3_ta;
    StringVector stage2_ta;
    StringVector stage1_ta;
    IntegerVector age3_ta;
    IntegerVector age2_ta;
    StringVector eststage3_ta;
    StringVector eststage2_ta;
    StringVector eststage1_ta;
    IntegerVector estage3_ta;
    IntegerVector estage2_ta;
    NumericVector givenrate_ta;
    NumericVector offset_ta;
    NumericVector multiplier_ta;
    IntegerVector convtype_ta;
    IntegerVector convtype_t12_ta;
    
    NumericVector surv_dev_ta;
    NumericVector obs_dev_ta;
    NumericVector size_dev_ta;
    NumericVector sizeb_dev_ta;
    NumericVector sizec_dev_ta;
    NumericVector repst_dev_ta;
    NumericVector fec_dev_ta;
    NumericVector jsurv_dev_ta;
    NumericVector jobs_dev_ta;
    NumericVector jsize_dev_ta;
    NumericVector jsizeb_dev_ta;
    NumericVector jsizec_dev_ta;
    NumericVector jrepst_dev_ta;
    NumericVector jmatst_dev_ta;
    
    StringVector pop_ta;
    StringVector patch_ta;
    StringVector year2_ta;
    
    variant_ta = as<IntegerVector>(trait_axis["variant"]);
    stage3_ta = as<StringVector>(trait_axis["stage3"]);
    stage2_ta = as<StringVector>(trait_axis["stage2"]);
    stage1_ta = as<StringVector>(trait_axis["stage1"]);
    age3_ta = as<IntegerVector>(trait_axis["age3"]);
    age2_ta = as<IntegerVector>(trait_axis["age2"]);
    eststage3_ta = as<StringVector>(trait_axis["eststage3"]);
    eststage2_ta = as<StringVector>(trait_axis["eststage2"]);
    eststage1_ta = as<StringVector>(trait_axis["eststage1"]);
    estage3_ta = as<IntegerVector>(trait_axis["estage3"]);
    estage2_ta = as<IntegerVector>(trait_axis["estage2"]);
    givenrate_ta = as<NumericVector>(trait_axis["givenrate"]);
    offset_ta = as<NumericVector>(trait_axis["offset"]);
    multiplier_ta = as<NumericVector>(trait_axis["multiplier"]);
    convtype_ta = as<IntegerVector>(trait_axis["convtype"]);
    convtype_t12_ta = as<IntegerVector>(trait_axis["convtype_t12"]);
    
    surv_dev_ta = as<NumericVector>(trait_axis["surv_dev"]);
    obs_dev_ta = as<NumericVector>(trait_axis["obs_dev"]);
    size_dev_ta = as<NumericVector>(trait_axis["size_dev"]);
    sizeb_dev_ta = as<NumericVector>(trait_axis["sizeb_dev"]);
    sizec_dev_ta = as<NumericVector>(trait_axis["sizec_dev"]);
    repst_dev_ta = as<NumericVector>(trait_axis["repst_dev"]);
    fec_dev_ta = as<NumericVector>(trait_axis["fec_dev"]);
    jsurv_dev_ta = as<NumericVector>(trait_axis["jsurv_dev"]);
    jobs_dev_ta = as<NumericVector>(trait_axis["jobs_dev"]);
    jsize_dev_ta = as<NumericVector>(trait_axis["jsize_dev"]);
    jsizeb_dev_ta = as<NumericVector>(trait_axis["jsizeb_dev"]);
    jsizec_dev_ta = as<NumericVector>(trait_axis["jsizec_dev"]);
    jrepst_dev_ta = as<NumericVector>(trait_axis["jrepst_dev"]);
    jmatst_dev_ta = as<NumericVector>(trait_axis["jmatst_dev"]);
    int ta_rows = static_cast<int>(stage3_ta.length());
    
    // Prep based on stages in stageframe
    StringVector stagevec = as<StringVector>(stageframe["stage"]);
    
    arma::ivec groupvec = as<arma::ivec>(stageframe["group"]);
    int stageframe_length {static_cast<int>(stagevec.length())};
    IntegerVector stage_id = seq(1, stageframe_length);
    
    if (pure_leslie) {
      IntegerVector sf_min_age = as<IntegerVector>(stageframe["min_age"]);
      int base_age = first_age_int;
      
      for (int i = 0; i < stageframe_length; i++) {
        if (IntegerVector::is_na(sf_min_age(i)) || NumericVector::is_na(sf_min_age(i))) {
          sf_min_age(i) = first_age_int + i;
        }
        
        if (StringVector::is_na(stagevec(i))) {
          if (!IntegerVector::is_na(sf_min_age(i)) && !NumericVector::is_na(sf_min_age(i))) {
            base_age = static_cast<int>(sf_min_age(i));
            String base_age_text = "Age";
            base_age_text += base_age;
            
            stagevec(i) = base_age_text;
          }
        }
      }
    }
    
    // Identify all groups
    arma::ivec all_groups = unique(groupvec);
    int no_groups {static_cast<int>(all_groups.n_elem)};
    StringVector group_text(no_groups);
    
    for (int i = 0; i < no_groups; i++) {
      group_text(i) = "group";
      group_text(i) += std::to_string(all_groups(i));
    }
    
    StringVector unique_stages = unique(stagevec);
    StringVector extra_terms = {"rep", "nrep", "immat", "mat", "prop", "npr", "all", "obs", "nobs"};
    
    int no_newstages {static_cast<int>(unique_stages.length())};
    int no_extraterms {static_cast<int>(extra_terms.length())};
    
    StringVector all_possible_stage_terms(no_newstages + no_extraterms + no_groups);
    for (int i = 0; i < no_newstages; i++) {
      all_possible_stage_terms(i) = unique_stages(i);
    }
    for (int i = 0; i < no_extraterms; i++) {
      all_possible_stage_terms(i + no_newstages) = extra_terms(i);
    }
    for (int i = 0; i < no_groups; i++) {
      all_possible_stage_terms(i + no_newstages + no_extraterms) = group_text(i);
    }
    
    if (trait_axis.containsElementNamed("year2")) {
      year2_ta = as<StringVector>(trait_axis["year2"]);
    } else {
      StringVector year2_ta_ (stage2_ta.length(), NA_STRING);
      year2_ta = year2_ta_;
    }
    
    arma::uvec vrm_altered (ta_rows, fill::zeros); // vrm alteration (binary)  by trait axis row
    arma::uvec mpm_altered (ta_rows, fill::zeros); // mpm alteration (binary)  by trait axis row
    
    for (int i = 0; i < ta_rows; i++) {
      if (!StringVector::is_na(eststage3_ta(i))) {
        mpm_altered(i) = 1;
      }
      if (!IntegerVector::is_na(estage2_ta(i))) {
        mpm_altered(i) = 1;
      }
      if (!NumericVector::is_na(givenrate_ta(i))) {
        mpm_altered(i) = 1;
      }
      if (!NumericVector::is_na(offset_ta(i))) {
        if (offset_ta(i) != 0.) mpm_altered(i) = 1;
      }
      if (!NumericVector::is_na(multiplier_ta(i))) {
        if (multiplier_ta(i) != 1.) mpm_altered(i) = 1;
      }
      
      if (!NumericVector::is_na(surv_dev_ta(i))) {
        if (surv_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(obs_dev_ta(i))) {
        if (obs_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(size_dev_ta(i))) {
        if (size_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(sizeb_dev_ta(i))) {
        if (sizeb_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(sizec_dev_ta(i))) {
        if (sizec_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(repst_dev_ta(i))) {
        if (repst_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(fec_dev_ta(i))) {
        if (fec_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      
      if (!NumericVector::is_na(jsurv_dev_ta(i))) {
        if (jsurv_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(jobs_dev_ta(i))) {
        if (jobs_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(jsize_dev_ta(i))) {
        if (jsize_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(jsizeb_dev_ta(i))) {
        if (jsizeb_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(jsizec_dev_ta(i))) {
        if (jsizec_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(jrepst_dev_ta(i))) {
        if (jrepst_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
      if (!NumericVector::is_na(jmatst_dev_ta(i))) {
        if (jmatst_dev_ta(i) != 0.) vrm_altered(i) = 1;
      }
    }
    
    // Check entries in trait_axis table
    for (int i = 0; i < static_cast<int>(ta_rows); i++) {
      int s3ta_count {0};
      int s2ta_count {0};
      int s1ta_count {0};
      
      bool ests3_used {false};
      bool ests2_used {false};
      bool ests1_used {false};
      
      bool s3_used {false};
      bool s2_used {false};
      bool s1_used {false};
      
      for (int j = 0; j < static_cast<int>(all_possible_stage_terms.length()); j++) {
        if (stage3_ta(i) == all_possible_stage_terms(j)) s3ta_count++;
        if (stage2_ta(i) == all_possible_stage_terms(j)) s2ta_count++;
        
        if (!StringVector::is_na(stage3_ta(i))) {
          s3_used = true;
        }
        if (!StringVector::is_na(stage2_ta(i))) {
          s2_used = true;
        }
        
        if (!StringVector::is_na(eststage3_ta(i))) {
          ests3_used = true;
        }
        if (!StringVector::is_na(eststage2_ta(i))) {
          ests2_used = true;
        }
        
        if (historical) {
          if (stage1_ta(i) == all_possible_stage_terms(j)) s1ta_count++;
          if (stage1_ta(i) == "NotAlive") s1ta_count++;
          
          if (!StringVector::is_na(stage1_ta(i))) {
            s1_used = true;
          }
          
          if (!StringVector::is_na(eststage1_ta(i))) {
            ests1_used = true;
          }
        } 
      }
      
      if (pure_leslie) {
        if (!IntegerVector::is_na(age2_ta(i)) && !NumericVector::is_na(age2_ta(i))) {
          if (StringVector::is_na(stage3_ta(i))) {
            String base_age_text = "Age";
            
            if (!IntegerVector::is_na(age3_ta(i)) && !NumericVector::is_na(age3_ta(i))) {
              base_age_text += static_cast<int>(age3_ta(i));
            } else {
              base_age_text += (static_cast<int>(age2_ta(i)) + 1);
            }
            
            bool found_stage3 {false};
            for (int j = 0; j < stageframe_length; j++) {
              if (base_age_text == stagevec(j)) found_stage3 = true;
            }
            
            if (found_stage3) {
              stage3_ta(i) = base_age_text;
              s3ta_count++;
              s3_used = true;
            }
          }
          
          if (StringVector::is_na(stage2_ta(i))) {
            String base_age_text = "Age";
            
            if (!IntegerVector::is_na(age2_ta(i)) && !NumericVector::is_na(age2_ta(i))) {
              base_age_text += static_cast<int>(age2_ta(i));
            } else {
              throw Rcpp::exception("Insufficient age information provided in trait_axis.", false);
            }
            
            bool found_stage2 {false};
            for (int j = 0; j < stageframe_length; j++) {
              if (base_age_text == stagevec(j)) found_stage2 = true;
            }
            
            if (found_stage2) {
              stage2_ta(i) = base_age_text;
              s2ta_count++;
              s2_used = true;
            }
          }
        }
        
        if (!IntegerVector::is_na(estage2_ta(i)) && !NumericVector::is_na(estage2_ta(i))) {
          if (StringVector::is_na(eststage3_ta(i))) {
            String base_age_text = "Age";
            
            if (!IntegerVector::is_na(estage3_ta(i)) && !NumericVector::is_na(estage3_ta(i))) {
              base_age_text += static_cast<int>(estage3_ta(i));
            } else {
              base_age_text += (static_cast<int>(estage2_ta(i)) + 1);
            }
            
            bool found_eststage3 {false};
            for (int j = 0; j < stageframe_length; j++) {
              if (base_age_text == stagevec(j)) found_eststage3 = true;
            }
            
            if (found_eststage3) {
              eststage3_ta(i) = base_age_text;
              s3ta_count++;
              ests3_used = true;
            }
          }
          
          if (StringVector::is_na(eststage2_ta(i))) {
            String base_age_text = "Age";
            
            if (!IntegerVector::is_na(estage2_ta(i)) && !NumericVector::is_na(estage2_ta(i))) {
              base_age_text += static_cast<int>(estage2_ta(i));
            } else {
              throw Rcpp::exception("Insufficient proxy age information provided in trait_axis.", false);
            }
            
            bool found_eststage2 {false};
            for (int j = 0; j < stageframe_length; j++) {
              if (base_age_text == stagevec(j)) found_eststage2 = true;
            }
            
            if (found_eststage2) {
              eststage2_ta(i) = base_age_text;
              s2ta_count++;
              ests2_used = true;
            }
          }
        }
      }
      
      if (s3ta_count == 0 && (s3_used || !functionbased) && mpm_altered(i) > 0) {
        String eat_my_shorts = "Stage names in trait_axis variable ";
        String eat_my_shorts1 = "stage3 must match stageframe.";
        eat_my_shorts += eat_my_shorts1;
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
      if (s2ta_count == 0 && (s2_used || !functionbased) && mpm_altered(i) > 0) {
        String eat_my_shorts = "Stage names in trait_axis variable ";
        String eat_my_shorts1 = "stage2 must match stageframe.";
        eat_my_shorts += eat_my_shorts1;
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
      if (ests3_used && mpm_altered(i) > 0) {
        if (s3ta_count == 0) {
          String eat_my_shorts = "Stage names in trait_axis variable ";
          String eat_my_shorts1 = "eststage3 must match stageframe.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      if (ests2_used && mpm_altered(i) > 0) {
        if (s2ta_count == 0) {
          String eat_my_shorts = "Stage names in trait_axis variable ";
          String eat_my_shorts1 = "eststage2 must match stageframe.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      if (historical && mpm_altered(i) > 0) {
        if (s1ta_count == 0 && (s1_used || !functionbased)) {
          String eat_my_shorts = "Stage names in trait_axis variable ";
          String eat_my_shorts1 = "stage1 must match stageframe.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        if (ests1_used && mpm_altered(i) > 0) {
          if (s1ta_count == 0) {
            String eat_my_shorts = "Stage names in trait_axis variable ";
            String eat_my_shorts1 = "eststage1 must match stageframe.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
      }
    }
    
    IntegerVector s1_calls (ta_rows, 1);
    IntegerVector s2_calls (ta_rows, 1);
    IntegerVector s3_calls (ta_rows, 1);
    IntegerVector ests1_calls (ta_rows, 1);
    IntegerVector ests2_calls (ta_rows, 1);
    IntegerVector ests3_calls (ta_rows, 1);
    IntegerVector s3_planned (ta_rows, 1);
    IntegerVector s2_planned (ta_rows, 1);
    IntegerVector s1_planned (ta_rows, 1);
    
    IntegerVector found_est_calls (ta_rows);
    
    IntegerVector s123_calls (ta_rows, 1);
    
    // Create indices for edited trait_axis table
    arma::uvec alive;
    if (stageframe.containsElementNamed("alive")) {
      arma::uvec alive_temp = as<arma::uvec>(stageframe["alive"]);
      alive = alive_temp;
    } else {
      arma::uvec alive_temp (stageframe_length, fill::ones);
      alive = alive_temp;
    }
    
    arma::uvec repvec = as<arma::uvec>(stageframe["repstatus"]);
    arma::uvec obsvec = as<arma::uvec>(stageframe["obsstatus"]);
    arma::uvec propvec = as<arma::uvec>(stageframe["propstatus"]);
    arma::uvec immvec = as<arma::uvec>(stageframe["immstatus"]);
    arma::uvec matvec = as<arma::uvec>(stageframe["matstatus"]);
    arma::uvec indvec = as<arma::uvec>(stageframe["indataset"]);
    
    arma::uvec newprop_stages = find(propvec);
    arma::uvec newprop0_stages = find(propvec == 0);
    arma::uvec newimm_stages = find(immvec);
    arma::uvec newalive_stages = find(alive);
    arma::uvec newmat_stages1 = find(matvec);
    arma::uvec newmat_stages = intersect(newalive_stages, newmat_stages1);
    arma::uvec newrep_stages = find(repvec);
    arma::uvec newrep0_stages = find(repvec == 0);
    arma::uvec newmat_rep0_stages = intersect(newmat_stages, newrep0_stages);
    arma::uvec newobs_stages = find(obsvec);
    arma::uvec newobs0_stages = find(obsvec == 0);
    arma::uvec all_stages = find(alive);
    
    int np_s = static_cast<int>(newprop_stages.n_elem);
    int np0_s = static_cast<int>(newprop0_stages.n_elem);
    int ni_s = static_cast<int>(newimm_stages.n_elem);
    int nm_s = static_cast<int>(newmat_stages.n_elem);
    int nr_s = static_cast<int>(newrep_stages.n_elem);
    int nmr0_s = static_cast<int>(newmat_rep0_stages.n_elem);
    int no_s = static_cast<int>(newobs_stages.n_elem);
    int no0_s = static_cast<int>(newobs0_stages.n_elem);
    int a_s = static_cast<int>(all_stages.n_elem);
    
    // Build expanded trait_axis table
    for (int i = 0; i < ta_rows; i++) {
      s3_calls(i) = LefkoMats::supp_decision1(as<std::string>(stage3_ta(i)),
        np_s, np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups,
        groupvec, group_text);
      
      ests3_calls(i) = LefkoMats::supp_decision1(as<std::string>(eststage3_ta(i)),
        np_s, np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups,
        groupvec, group_text);
      
      s2_calls(i) = LefkoMats::supp_decision1(as<std::string>(stage2_ta(i)),
        np_s, np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups,
        groupvec, group_text);
      
      ests2_calls(i) = LefkoMats::supp_decision1(as<std::string>(eststage2_ta(i)),
        np_s, np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups,
        groupvec, group_text);
      
      s1_calls(i) = LefkoMats::supp_decision1(as<std::string>(stage1_ta(i)),
        np_s, np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups,
        groupvec, group_text);
      
      ests1_calls(i) = LefkoMats::supp_decision1(as<std::string>(eststage1_ta(i)),
        np_s, np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups,
        groupvec, group_text);
      
      String eat_my_shorts_gse = "If stage group shorthand is used to designate ";
      eat_my_shorts_gse += "both a transition and a proxy, then the ";
      eat_my_shorts_gse += "shorthand group must be the same in both cases.";
      
      if (!StringVector::is_na(eststage3_ta(i))) {
        if (eststage3_ta(i) != stage3_ta(i)) {
          if (s3_calls(i) == 1 && ests3_calls(i) > 1) {
            s3_planned(i) = ests3_calls(i);
            found_est_calls(i) = 1;
          } else if (s3_calls(i) > 1 && ests3_calls(i) > 1) {
            throw Rcpp::exception(eat_my_shorts_gse.get_cstring(), false);
          }
        } else {
          s3_planned(i) = s3_calls(i);
        }
      } else {
        s3_planned(i) = s3_calls(i);
      }
      
      if (!StringVector::is_na(eststage2_ta(i))) {
        if (eststage2_ta(i) != stage2_ta(i)) {
          if (s2_calls(i) == 1 && ests2_calls(i) > 1) {
            s2_planned(i) = ests2_calls(i);
          } else if (s2_calls(i) > 1 && ests2_calls(i) > 1) {
            throw Rcpp::exception(eat_my_shorts_gse.get_cstring(), false);
          }
        } else {
          s2_planned(i) = s2_calls(i);
        }
      } else {
        s2_planned(i) = s2_calls(i);
      }
      
      if (!StringVector::is_na(eststage1_ta(i))) {
        if (historical && eststage1_ta(i) != stage1_ta(i)) {
          if (s1_calls(i) == 1 && ests1_calls(i) > 1) {
            s1_planned(i) = ests1_calls(i);
          } else if (s1_calls(i) > 1 && ests1_calls(i) > 1) {
            throw Rcpp::exception(eat_my_shorts_gse.get_cstring(), false);
          }
        } else if (historical) {
          s1_planned(i) = s1_calls(i);
        } else if (!historical) {
          s1_planned(i) = 1;
        }
      } else {
        s1_planned(i) = s1_calls(i);
      }
      
      s123_calls(i) = s3_planned(i) * s2_planned(i) * s1_planned(i);
    }
    
    NumericVector basepoints(ta_rows, 0.0);
    for (int i = 0; i < (ta_rows - 1); i++) {
      basepoints(i+1) = basepoints(i) + s123_calls(i);
    }
    
    // Basic format-checking of vectors - should include check of eststages
    IntegerVector vital_rate_modifiers_by_variant (ta_rows);
    IntegerVector mat_elem_modifiers_by_variant (ta_rows);
    
    for (int i = 0; i < ta_rows; i++) {
      NumericVector mat_elem_dev = {givenrate_ta(i), offset_ta(i),
        multiplier_ta(i)};
      arma::vec mat_elem_dev_arma = as<arma::vec>(mat_elem_dev);
      arma::uvec mat_elem_dev_nonzeros = find(mat_elem_dev_arma);
      int mat_elem_found_nonzeros = static_cast<int>(mat_elem_dev_nonzeros.n_elem);
      mat_elem_modifiers_by_variant(i) = mat_elem_found_nonzeros;
      
      NumericVector current_dev = {surv_dev_ta(i), obs_dev_ta(i), size_dev_ta(i),
        sizeb_dev_ta(i), sizec_dev_ta(i), repst_dev_ta(i), fec_dev_ta(i),
        jsurv_dev_ta(i), jobs_dev_ta(i), jsize_dev_ta(i), jsizeb_dev_ta(i),
        jsizec_dev_ta(i), jrepst_dev_ta(i), jmatst_dev_ta(i)};
      arma::vec current_dev_arma = as<arma::vec>(current_dev);
      arma::uvec current_dev_nonzeros = find(current_dev_arma);
      int current_found_nonzeros = static_cast<int>(current_dev_nonzeros.n_elem);
      vital_rate_modifiers_by_variant(i) = current_found_nonzeros;
    }
    
    {
      if (!functionbased) {
        IntegerVector check_calls_mpm = mat_elem_modifiers_by_variant + found_est_calls;
        arma::ivec check_calls_mpm_arma = as<arma::ivec>(check_calls_mpm);
        arma::uvec cca_foundzeros_mpm = find(check_calls_mpm_arma == 0);
        int ccafz_length_mpm = cca_foundzeros_mpm.n_elem;
        
        if (ccafz_length_mpm > 0) {
          Rf_warningcall(R_NilValue,
            "Some variants have no alterations to matrices in argument trait_axis");
        }
      }
      
      IntegerVector check_calls = mat_elem_modifiers_by_variant + 
        vital_rate_modifiers_by_variant + found_est_calls;
      arma::ivec check_calls_arma = as<arma::ivec>(check_calls);
      arma::uvec cca_foundzeros = find(check_calls_arma == 0);
      int ccafz_length = cca_foundzeros.n_elem;
      
      if (ccafz_length > 0 && functionbased) {
        Rf_warningcall(R_NilValue,
          "Some variants have no alterations to matrices and vital rate models in argument trait_axis");
      }
    }
    
    // New trait_axis set-up
    int newta_rows = sum(s123_calls);
    
    IntegerVector variant_newta (newta_rows);
    StringVector stage3_newta (newta_rows);
    StringVector stage2_newta (newta_rows);
    StringVector stage1_newta (newta_rows);
    IntegerVector age3_newta (newta_rows);
    IntegerVector age2_newta (newta_rows);
    StringVector eststage3_newta (newta_rows);
    StringVector eststage2_newta (newta_rows);
    StringVector eststage1_newta (newta_rows);
    IntegerVector estage3_newta (newta_rows);
    IntegerVector estage2_newta (newta_rows);
    NumericVector givenrate_newta (newta_rows);
    NumericVector offset_newta (newta_rows);
    NumericVector multiplier_newta (newta_rows);
    IntegerVector convtype_newta (newta_rows);
    IntegerVector convtype_t12_newta (newta_rows);
    NumericVector surv_dev_newta (newta_rows);
    NumericVector obs_dev_newta (newta_rows);
    NumericVector size_dev_newta (newta_rows);
    NumericVector sizeb_dev_newta (newta_rows);
    NumericVector sizec_dev_newta (newta_rows);
    NumericVector repst_dev_newta (newta_rows);
    NumericVector fec_dev_newta (newta_rows);
    NumericVector jsurv_dev_newta (newta_rows);
    NumericVector jobs_dev_newta (newta_rows);
    NumericVector jsize_dev_newta (newta_rows);
    NumericVector jsizeb_dev_newta (newta_rows);
    NumericVector jsizec_dev_newta (newta_rows);
    NumericVector jrepst_dev_newta (newta_rows);
    NumericVector jmatst_dev_newta (newta_rows);
    NumericVector mpm_altered_newta (newta_rows);
    NumericVector vrm_altered_newta (newta_rows);
    
    StringVector year2_newta (newta_rows);
    
    int overall_counter {0};
    int group_check {0};
    
    int group_baseline3 {0};
    int group_baseline2 {0};
    int group_baseline1 {0};
    int group_baselinee3 {0};
    int group_baselinee2 {0};
    int group_baselinee1 {0};
    
    int group_ratchet3 {0};
    int group_ratchet2 {0};
    int group_ratchet1 {0};
    int group_ratchete3 {0};
    int group_ratchete2 {0};
    int group_ratchete1 {0};
    
    int prevl3 {0};
    int prevl2 {0};
    int prevl1 {0};
    int prevle3 {0};
    int prevle2 {0};
    int prevle1 {0};
    
    for (int i = 0; i < ta_rows; i++) {
      overall_counter = 0;
      for (int j = 0; j < s1_planned(i); j++) {
        for (int k = 0; k < s2_planned(i); k++) {
          for (int l = 0; l < s3_planned(i); l++) {
            stage3_newta(basepoints(i) + overall_counter) =
              LefkoMats::supp_decision2(as<std::string>(stage3_ta(i)),
                newprop_stages, newprop0_stages, newimm_stages, newmat_stages,
                newrep_stages, newmat_rep0_stages, newobs_stages,
                newobs0_stages, all_stages, no_groups, groupvec, group_text,
                stagevec, l, group_check, group_ratchet3, group_baseline3,
                prevl3);
            
            variant_newta(basepoints(i) + overall_counter) = variant_ta(i);
            givenrate_newta(basepoints(i) + overall_counter) = givenrate_ta(i);
            offset_newta(basepoints(i) + overall_counter) = offset_ta(i);
            multiplier_newta(basepoints(i) + overall_counter) = multiplier_ta(i);
            convtype_newta(basepoints(i) + overall_counter) = convtype_ta(i);
            convtype_t12_newta(basepoints(i) + overall_counter) = convtype_t12_ta(i);
            surv_dev_newta(basepoints(i) + overall_counter) = surv_dev_ta(i);
            obs_dev_newta(basepoints(i) + overall_counter) = obs_dev_ta(i);
            size_dev_newta(basepoints(i) + overall_counter) = size_dev_ta(i);
            sizeb_dev_newta(basepoints(i) + overall_counter) = sizeb_dev_ta(i);
            sizec_dev_newta(basepoints(i) + overall_counter) = sizec_dev_ta(i);
            repst_dev_newta(basepoints(i) + overall_counter) = repst_dev_ta(i);
            fec_dev_newta(basepoints(i) + overall_counter) = fec_dev_ta(i);
            jsurv_dev_newta(basepoints(i) + overall_counter) = jsurv_dev_ta(i);
            jobs_dev_newta(basepoints(i) + overall_counter) = jobs_dev_ta(i);
            jsize_dev_newta(basepoints(i) + overall_counter) = jsize_dev_ta(i);
            jsizeb_dev_newta(basepoints(i) + overall_counter) = jsizeb_dev_ta(i);
            jsizec_dev_newta(basepoints(i) + overall_counter) = jsizec_dev_ta(i);
            jrepst_dev_newta(basepoints(i) + overall_counter) = jrepst_dev_ta(i);
            jmatst_dev_newta(basepoints(i) + overall_counter) = jmatst_dev_ta(i);
            mpm_altered_newta(basepoints(i) + overall_counter) = mpm_altered(i);
            vrm_altered_newta(basepoints(i) + overall_counter) = vrm_altered(i);
            
            eststage3_newta(basepoints(i) + overall_counter) =
              LefkoMats::supp_decision2(as<std::string>(eststage3_ta(i)),
                newprop_stages, newprop0_stages, newimm_stages, newmat_stages,
                newrep_stages, newmat_rep0_stages, newobs_stages,
                newobs0_stages, all_stages, no_groups, groupvec, group_text,
                stagevec, l, group_check, group_ratchete3, group_baselinee3,
                prevle3);
            
            stage2_newta(basepoints(i) + overall_counter) =
              LefkoMats::supp_decision2(as<std::string>(stage2_ta(i)),
                newprop_stages, newprop0_stages, newimm_stages, newmat_stages,
                newrep_stages, newmat_rep0_stages, newobs_stages,
                newobs0_stages, all_stages, no_groups, groupvec, group_text,
                stagevec, k, group_check, group_ratchet2, group_baseline2,
                prevl2);
            
            eststage2_newta(basepoints(i) + overall_counter) =
              LefkoMats::supp_decision2(as<std::string>(eststage2_ta(i)),
                newprop_stages, newprop0_stages, newimm_stages, newmat_stages,
                newrep_stages, newmat_rep0_stages, newobs_stages,
                newobs0_stages, all_stages, no_groups, groupvec, group_text,
                stagevec, k, group_check, group_ratchete2, group_baselinee2,
                prevle2);
            
            stage1_newta(basepoints(i) + overall_counter) =
              LefkoMats::supp_decision2(as<std::string>(stage1_ta(i)),
                newprop_stages, newprop0_stages, newimm_stages, newmat_stages,
                newrep_stages, newmat_rep0_stages, newobs_stages,
                newobs0_stages, all_stages, no_groups, groupvec, group_text,
                stagevec, j, group_check, group_ratchet1, group_baseline1,
                prevl1);
            
            eststage1_newta(basepoints(i) + overall_counter) =
              LefkoMats::supp_decision2(as<std::string>(eststage1_ta(i)),
                newprop_stages, newprop0_stages, newimm_stages, newmat_stages,
                newrep_stages, newmat_rep0_stages, newobs_stages,
                newobs0_stages, all_stages, no_groups, groupvec, group_text,
                stagevec, j, group_check, group_ratchete1, group_baselinee1,
                prevle1);
            
            age3_newta(basepoints(i) + overall_counter) = age3_ta(i);
            age2_newta(basepoints(i) + overall_counter) = age2_ta(i);
            estage3_newta(basepoints(i) + overall_counter) = estage3_ta(i);
            estage2_newta(basepoints(i) + overall_counter) = estage2_ta(i);
            
            year2_newta(basepoints(i) + overall_counter) = year2_ta(i);
            
            overall_counter++;
          }
        }
      }
    }
    
    Rcpp::List newtraitaxis(33);
    
    newtraitaxis(0) = variant_newta;
    newtraitaxis(1) = stage3_newta;
    newtraitaxis(2) = stage2_newta;
    newtraitaxis(3) = stage1_newta;
    newtraitaxis(4) = age3_newta;
    newtraitaxis(5) = age2_newta;
    newtraitaxis(6) = eststage3_newta;
    newtraitaxis(7) = eststage2_newta;
    newtraitaxis(8) = eststage1_newta;
    newtraitaxis(9) = estage3_newta;
    newtraitaxis(10) = estage2_newta;
    newtraitaxis(11) = givenrate_newta;
    newtraitaxis(12) = offset_newta;
    newtraitaxis(13) = multiplier_newta;
    newtraitaxis(14) = convtype_newta;
    newtraitaxis(15) = convtype_t12_newta;
    newtraitaxis(16) = surv_dev_newta;
    newtraitaxis(17) = obs_dev_newta;
    newtraitaxis(18) = size_dev_newta;
    newtraitaxis(19) = sizeb_dev_newta;
    newtraitaxis(20) = sizec_dev_newta;
    newtraitaxis(21) = repst_dev_newta;
    newtraitaxis(22) = fec_dev_newta;
    newtraitaxis(23) = jsurv_dev_newta;
    newtraitaxis(24) = jobs_dev_newta;
    newtraitaxis(25) = jsize_dev_newta;
    newtraitaxis(26) = jsizeb_dev_newta;
    newtraitaxis(27) = jsizec_dev_newta;
    newtraitaxis(28) = jrepst_dev_newta;
    newtraitaxis(29) = jmatst_dev_newta;
    newtraitaxis(30) = year2_newta;
    newtraitaxis(31) = mpm_altered_newta;
    newtraitaxis(32) = vrm_altered_newta;
    
    CharacterVector ta_namevec = {"variant", "stage3", "stage2", "stage1",
      "age3", "age2", "eststage3", "eststage2", "eststage1", "estage3",
      "estage2", "givenrate", "offset", "multiplier", "convtype",
      "convtype_t12", "surv_dev", "obs_dev", "size_dev", "sizeb_dev",
      "sizec_dev", "repst_dev", "fec_dev", "jsurv_dev", "jobs_dev", "jsize_dev",
      "jsizeb_dev", "jsizec_dev", "jrepst_dev", "jmatst_dev", "year2",
      "mpm_altered", "vrm_altered"};
    CharacterVector ta_newclasses = {"data.frame", "adaptAxis"};
    newtraitaxis.attr("names") = ta_namevec;
    newtraitaxis.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newta_rows);
    newtraitaxis.attr("class") = ta_newclasses;
    
    return newtraitaxis;
  }

  //' Standardized Error Messages
  //' 
  //' Function \code{pop_error()} produces standardized error messages.
  //' 
  //' @name pop_error
  //' 
  //' @param input1 First string input.
  //' @param input2 Second string input.
  //' @param input3 Third string input.
  //' @param type Designates output message type.
  //' 
  //' @return Stops R and produces an error message.
  //' 
  //' @section Notes:
  //' Pop errors 3 and 6 were merged with 1.
  //' 
  //' @keywords internal
  //' @noRd
  inline void pop_error2 (String input1, String input2, String input3,
    int type = 1) {
    
    String eat_my_shorts;
    if (type == 1) { // Very useful
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " should be entered as ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 2) {
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be the same length as the number of ";
      eat_my_shorts += input2;
      eat_my_shorts += " entered in argument ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    } else if (type == 24) {
      eat_my_shorts = "Do not use arguments ";
      eat_my_shorts += input1;
      eat_my_shorts += " and ";
      eat_my_shorts += input2;
      eat_my_shorts += " if ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    } else if (type == 25) {
      eat_my_shorts = input1;
      eat_my_shorts += " are not allowed in argument ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 26) {
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " is required to ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 27) {
      eat_my_shorts = "Arguments ";
      eat_my_shorts += input1;
      eat_my_shorts += " and ";
      eat_my_shorts += input2;
      eat_my_shorts += " must be ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    } else if (type == 28) {
      eat_my_shorts = "Arguments ";
      eat_my_shorts += input1;
      eat_my_shorts += " can only be used in ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 29) {
      eat_my_shorts = "Vector ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be the same length as ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 30) {
      eat_my_shorts = "Elements in argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " may not be negative.";
      
    }
    
    throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
    
    return;
  }
  
  //' Bind Two Data Frames By Row
  //' 
  //' This function takes two data frames, which must be composed of the same
  //' variables in the same order. Although the variables may have different
  //' names, the variables must be of the same type. The names of the variables in
  //' the new merged data frame will match those of the first data frame.
  //' Developed with the help of Microsoft Gemini AI.
  //' 
  //' @name df_rbind
  //' 
  //' @param df1 The first data frame, which will form the top of the new data
  //' frame and will be used as a reference for variable names.
  //' @param df2 The second data frame, which will be attached below data frame
  //' \code{df1}.
  //' 
  //' @return A new data frame composed of the merged data frames.
  //' 
  //' @keywords internal
  //' @noRd
  inline DataFrame df_rbind(DataFrame df1, DataFrame df2) {
    
    int nrows1 = df1.nrows();
    int nrows2 = df2.nrows();
    int ncols = df1.size();
  
    if (ncols != df2.size()) {
      stop("Data frames must have the same number of columns.");
    }
    
    CharacterVector df_names = df1.names();
  
    List out_df (ncols);
  
    for (int i = 0; i < ncols; ++i) {
      SEXP col1 = df1[i];
      SEXP col2 = df2[i];
  
      if (TYPEOF(col1) != TYPEOF(col2)) {
        stop("Columns must have the same type");
      }
  
      switch(TYPEOF(col1)){
        case INTSXP: {
          IntegerVector combinedCol(nrows1 + nrows2);
          IntegerVector col1Vector = as<IntegerVector>(col1);
          IntegerVector col2Vector = as<IntegerVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          bool current_int_class1 = col1Vector.hasAttribute("levels");
          if (current_int_class1) {
            CharacterVector col1lvls = col1Vector.attr("levels");
            CharacterVector col2lvls = col2Vector.attr("levels");
            
            CharacterVector col_mergedLevels = LefkoUtils::concat_str(col1lvls, col2lvls);
            CharacterVector unique_levels = unique(col_mergedLevels);
            
            combinedCol.attr("levels") = unique_levels;
          }
          out_df(i) = combinedCol;
          break;
        }
  
        case LGLSXP: {
          LogicalVector combinedCol(nrows1 + nrows2);
          LogicalVector col1Vector = as<LogicalVector>(col1);
          LogicalVector col2Vector = as<LogicalVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          out_df(i) = combinedCol;
          break;
        }
  
        case REALSXP: {
          NumericVector combinedCol(nrows1 + nrows2);
          NumericVector col1Vector = as<NumericVector>(col1);
          NumericVector col2Vector = as<NumericVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          out_df(i) = combinedCol;
          break;
        }
  
        case STRSXP: {
          CharacterVector combinedCol(nrows1 + nrows2);
          CharacterVector col1Vector = as<CharacterVector>(col1);
          CharacterVector col2Vector = as<CharacterVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          out_df(i) = combinedCol;
          break;
        }
        default:
          stop("Unsupported column type.");
      }
    }
    
    out_df.attr("names") = df_names;
    out_df.attr("class") = "data.frame";
    StringVector row_names(nrows1 + nrows2);
    for (int i = 0; i < (nrows1 + nrows2); i++) {
      row_names(i) = std::to_string(i+1);
    }
    out_df.attr("row.names") = row_names;
    
    return out_df;
  }
  
  //' Subset A Data Frame By Row Index
  //' 
  //' This function takes a data frame and subsets it according to row indices.
  //' The default is to use C++ indexing.
  //' 
  //' @name df_indices
  //' 
  //' @param x The data frame to subset.
  //' @param indices Integer vector giving the rows to keep.
  //' 
  //' @return A new data frame subset from the old.
  //' 
  //' @keywords internal
  //' @noRd
  inline List df_indices(const DataFrame& x, IntegerVector indices) {
    
    StringVector var_names = x.attr("names");
    StringVector df_class = x.attr("class");
    int no_vars = x.length();
    int no_rows = x.nrows();
    
    arma::uvec useable_indices;
    double chosen_level_d {0.0};
    String chosen_level_s;
    bool string_used {false};
    
    int max_asked_for = max(indices);
    if (max_asked_for >= no_rows) {
      throw Rcpp::exception("Row index too high for subsetting data frame chosen.", false);
    }
    
    int new_rows = static_cast<int>(indices.length());
    List new_df (no_vars);
    
    for (int i = 0; i < no_vars; i++) {
      if (is<NumericVector>(x[i])) {
        NumericVector old_var_i = as<NumericVector>(x[i]);
        NumericVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(indices(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<IntegerVector>(x[i])) {
        IntegerVector old_var_i = as<IntegerVector>(x[i]);
        IntegerVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(indices(j));
        }
        
        if (old_var_i.hasAttribute("levels")) {
          StringVector int_class = old_var_i.attr("class");
          if (stringcompare_simple(String(int_class(0)), "fact", false)) {
            new_var_i.attr("levels") = old_var_i.attr("levels");
            new_var_i.attr("class") = "factor";
          }
        }
        
        new_df(i) = new_var_i;
        
      } else if (is<LogicalVector>(x[i])) {
        LogicalVector old_var_i = as<LogicalVector>(x[i]);
        LogicalVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(indices(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<StringVector>(x[i])) {
        StringVector old_var_i = as<StringVector>(x[i]);
        StringVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(indices(j));
        }
        new_df(i) = new_var_i;
        
      } else {
        throw Rcpp::exception("Variable found of unrecognized type.", false);
      }
    }
    
    new_df.attr("names") = var_names;
    new_df.attr("class") = df_class;
    
    StringVector row_names(new_rows);
    for (int i = 0; i < new_rows; i++) {
      row_names(i) = std::to_string(i+1);
    }
    new_df.attr("row.names") = row_names;
    
    return new_df;
  }
  
  //' Creates Final Table of Lyapunov Estimates
  //' 
  //' @name Lyapunov_creator
  //' 
  //' @param in_df A pointer to a data frame that will be formed by this function.
  //' @param N_out The main list of population / variant size matrices.
  //' @param entry_time_vec An integer vector giving the entry time of each
  //' variant.
  //' @param nreps The number of replicates performed.
  //' @param var_per_run The number of variants per simulation.
  //' @param var_mat_length The total number of simulations.
  //' @param times The length of time to run simulations.
  //' @param fitness_times The time from which fitness should be estimated.
  //' @param threshold The threshold value to use to determine whether a fitness
  //' value is essentially zero.
  //' @param style An integer coding for whether the output should be for the
  //' full pairwise / multiple invasibility analyses (\code{0}), or for the ESS
  //' optimization (\code{1}). Defaults to \code{0}.
  //' @param optim_run A Boolean variable denoting whether the function is being
  //' used to optimize ESS values in the final search (not the initial opt_res
  //' search). If \code{true}, then assumes \code{N_out} is a list of matrices,
  //' rather than a list of cubes. Defaults to \code{false}.
  //' @param round A Boolean variable denoting whether to round fitness values
  //' below \code{threshold} to zero.
  //' 
  //' @return A data frame is formed and modified via pointer, so no actual
  //' output results from this function.
  //' 
  //' @keywords internal
  //' @noRd
  inline void Lyapunov_creator (DataFrame& in_df, List& N_out,
    IntegerVector entry_time_vec, const int nreps, const int var_per_run,
    const int var_mat_length, const int times, const int fitness_times,
    const double threshold, int style = 0, bool optim_run = false,
    bool round = false) {
    
    IntegerVector Lyap_rows = seq(1, nreps * var_mat_length);
    CharacterVector Lyap_df_names ((3 * var_per_run) + 2);
    Lyap_df_names(0) = "simulation_num";
    Lyap_df_names(1) = "rep";
    
    String Lyap_df_names_base1 = "variant";
    String Lyap_df_names_base2 = "entrytime";
    String Lyap_df_names_base3 = "fitness_variant";
    
    List output ((3 * var_per_run) + 2);
    output(0) = Lyap_rows;
    
    if (style == 1) {
      IntegerVector int_runs = seq(1, var_mat_length);
      IntegerVector int_runs2 = seq(1, var_mat_length);
      NumericVector fitness (var_mat_length);
      
      DataFrame coarse_agency = Rcpp::DataFrame::create(_["simulation_run"] = int_runs,
        _["variant"] = int_runs2, _["fitness_variant"] = fitness);
      in_df = coarse_agency;
    }
    
    int Lyap_counter {0};
    for (int m = 0; m < var_per_run; m++) {
      Lyap_counter = 0;
      
      IntegerVector Lyap_var_orig;
      if (style == 0 ) {
        Lyap_var_orig = as<IntegerVector>(in_df(1+m));
      } else {
        Lyap_var_orig = as<IntegerVector>(in_df(1));
      }
      IntegerVector Lyap_rep (nreps * var_mat_length);
      IntegerVector Lyap_var (nreps * var_mat_length);
      IntegerVector Lyap_etime (nreps * var_mat_length);
      NumericVector Lyap_fitness (nreps * var_mat_length);
      
      for (int i = 0; i < nreps; i++) {
        arma::cube current_N_out_cube;
        if (!optim_run) current_N_out_cube = as<arma::cube>(N_out(i));
        
        for (int l = 0; l < var_mat_length; l++) {
          int used_fitness_times = fitness_times;
          if (times - (entry_time_vec(m)) < fitness_times) {
            used_fitness_times = times - (entry_time_vec(m));
            Rf_warningcall(R_NilValue, "Truncating fitness_times due to late entry time.");
          }
          
          arma::mat current_N_out;
          if (!optim_run) {
            current_N_out = current_N_out_cube.slice(l);
          } else {
            current_N_out = as<arma::mat>(N_out(i));
          }
          int time_length = static_cast<int>(current_N_out.n_cols);
          int start_time = time_length - used_fitness_times;
          
          arma::vec running_fitness(used_fitness_times, fill::zeros);
          arma::vec generations(used_fitness_times, fill::zeros);
          arma::vec intercept_ones(used_fitness_times, fill::ones);
          
          for (int k = 0; k < used_fitness_times; k++) {
            running_fitness(k) = current_N_out(m, k+start_time);
            generations(k) = static_cast<double>(k);
          }
          
          arma::mat xmat = join_rows(intercept_ones, generations);
          arma::vec Lyap_regr;
          
          AdaptUtils::fastLm_sl(Lyap_regr, running_fitness, xmat);
          
          double Lyapunov_estimate = Lyap_regr(1);
          if (round && (abs(Lyapunov_estimate) <= threshold)) Lyapunov_estimate = 0.;
          
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
      if (style == 1 && m == 1) new_col_name += "_e995";
      Lyap_df_names(2 + m) = new_col_name;
      
      String next_col_name = Lyap_df_names_base2;
      next_col_name += (m + 1);
      if (style == 1 && m == 1) next_col_name += "_e995";
      Lyap_df_names(2 + var_per_run + m) = next_col_name;
      
      String last_col_name = Lyap_df_names_base3;
      last_col_name += (m + 1);
      if (style == 1 && m == 1) last_col_name += "_e995";
      Lyap_df_names(2 + (2 * var_per_run) + m) = last_col_name;
    }
    
    output.attr("names") = Lyap_df_names;
    output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      nreps * var_mat_length);
    StringVector needed_class {"data.frame"};
    output.attr("class") = needed_class;
    
    in_df = output;
  }
  
  //' Create Trait Axis Reassessed for Trait Optimization
  //' 
  //' Given some core choices by the user, this function will create a new trait
  //' axis for ESS trait value esitmation.
  //' 
  //' @name optim_ta_setup
  //' 
  //' @param ta_reassessed A data frame giving trait axis data post-processing
  //' with function \code{ta_reassess()}.
  //' @param ESS_ta An initially empty data frame provided as a reference, which
  //' is modified by this function to yield a trait axis for ESS optimization.
  //' @param optim_ta A data frame giving trait axis data processed for ESS
  //' optimization. This is one of two data frames, and gives the core info.
  //' @param optim_ta_995 A data frame giving trait axis data processed for ESS
  //' optimization. This is one of two data frames, and gives info for 99.5% the
  //' values of variable traits in \code{optim_ta}.
  //' @param variable_traits An integer vector modifed by this function by
  //' reference, indicating the actual traits that vary. The element order is:
  //' 1, givenrate; 2, offset; 3, multiplier; 4, surv_dev; 5, obs_dev;
  //' 6, size_dev; 7, sizeb_dev; 8, sizec_dev; 9, repst_dev; 10, fec_dev;
  //' 11, jsurv_dev; 12, jobs_dev; 13, jsize_dev; 14, jsizeb_dev; 15, jsizec_dev;
  //' 16, jrepst_dev; and 17, jmatst_dev.
  //' @param opt_res Relic integer value denoting the number of rows in the
  //' trait_axis used in trait optimization.
  //' @param elast_multiplier A double precision floating point number to
  //' multiply traits by to assess elasticity. Defaults to 0.995. 
  //' 
  //' @return Creates a new data frame with all variants to be tested, in order,
  //' and places that at the \code{optim_ta} reference. Also creates a list with
  //' the first element corresponding to the first and last elements for
  //' optimization, and the second element a vector giving the variable traits.
  //' 
  //' @keywords internal
  //' @noRd
  inline void optim_ta_setup (DataFrame& ta_reassessed, DataFrame& ESS_ta,
    DataFrame& optim_ta, DataFrame& optim_ta_995, IntegerVector& variable_traits,
    int opt_res, double elast_multiplier = 0.995) {
    
    DataFrame ESS_optim_ta;
    DataFrame cloned_ta_reassess;
    IntegerVector ESS_variable_traits (17, 0);
    
    // Vectors for optim_ta
    IntegerVector optim_variant;
    
    StringVector optim_stage3;
    StringVector optim_stage2;
    StringVector optim_stage1;
    IntegerVector optim_age3;
    IntegerVector optim_age2;
    StringVector optim_eststage3;
    StringVector optim_eststage2;
    StringVector optim_eststage1;
    IntegerVector optim_estage3;
    IntegerVector optim_estage2;
    NumericVector optim_givenrate;
    NumericVector optim_offset;
    NumericVector optim_multiplier;
    IntegerVector optim_convtype;
    IntegerVector optim_convtype_t12;
    StringVector optim_year2;
    IntegerVector optim_mpm_altered;
    IntegerVector optim_vrm_altered;
    
    NumericVector optim_surv_dev;
    NumericVector optim_obs_dev;
    NumericVector optim_size_dev;
    NumericVector optim_sizeb_dev;
    NumericVector optim_sizec_dev;
    NumericVector optim_repst_dev;
    NumericVector optim_fec_dev;
    
    NumericVector optim_jsurv_dev;
    NumericVector optim_jobs_dev;
    NumericVector optim_jsize_dev;
    NumericVector optim_jsizeb_dev;
    NumericVector optim_jsizec_dev;
    NumericVector optim_jrepst_dev;
    NumericVector optim_jmatst_dev;
    
    // Vectors for optim_ta_995
    IntegerVector optim_variant_995;
    
    StringVector optim_stage3_995;
    StringVector optim_stage2_995;
    StringVector optim_stage1_995;
    IntegerVector optim_age3_995;
    IntegerVector optim_age2_995;
    StringVector optim_eststage3_995;
    StringVector optim_eststage2_995;
    StringVector optim_eststage1_995;
    IntegerVector optim_estage3_995;
    IntegerVector optim_estage2_995;
    NumericVector optim_givenrate_995;
    NumericVector optim_offset_995;
    NumericVector optim_multiplier_995;
    IntegerVector optim_convtype_995;
    IntegerVector optim_convtype_t12_995;
    StringVector optim_year2_995;
    IntegerVector optim_mpm_altered_995;
    IntegerVector optim_vrm_altered_995;
    
    NumericVector optim_surv_dev_995;
    NumericVector optim_obs_dev_995;
    NumericVector optim_size_dev_995;
    NumericVector optim_sizeb_dev_995;
    NumericVector optim_sizec_dev_995;
    NumericVector optim_repst_dev_995;
    NumericVector optim_fec_dev_995;
    
    NumericVector optim_jsurv_dev_995;
    NumericVector optim_jobs_dev_995;
    NumericVector optim_jsize_dev_995;
    NumericVector optim_jsizeb_dev_995;
    NumericVector optim_jsizec_dev_995;
    NumericVector optim_jrepst_dev_995;
    NumericVector optim_jmatst_dev_995;
    
    // Vectors for ESS_trait_axis
    IntegerVector ESS_variant = {1, 2};
    
    StringVector ESS_stage3 (2, NA_STRING);
    StringVector ESS_stage2 (2, NA_STRING);
    StringVector ESS_stage1 (2, NA_STRING);
    IntegerVector ESS_age3 (2, NA_INTEGER);
    IntegerVector ESS_age2 (2, NA_INTEGER);
    StringVector ESS_eststage3 (2, NA_STRING);
    StringVector ESS_eststage2 (2, NA_STRING);
    StringVector ESS_eststage1 (2, NA_STRING);
    IntegerVector ESS_estage3 (2, NA_INTEGER);
    IntegerVector ESS_estage2 (2, NA_INTEGER);
    IntegerVector ESS_convtype (2, NA_INTEGER);
    IntegerVector ESS_convtype_t12 (2, NA_INTEGER);
    StringVector ESS_year2 (2, NA_STRING);
    IntegerVector ESS_mpm_altered (2, NA_INTEGER);
    IntegerVector ESS_vrm_altered (2, NA_INTEGER);
    
    NumericVector ESS_givenrate (2, NA_REAL);
    NumericVector ESS_offset (2, NA_REAL);
    NumericVector ESS_multiplier (2, NA_REAL);
    
    NumericVector ESS_surv_dev (2, NA_REAL);
    NumericVector ESS_obs_dev (2, NA_REAL);
    NumericVector ESS_size_dev (2, NA_REAL);
    NumericVector ESS_sizeb_dev (2, NA_REAL);
    NumericVector ESS_sizec_dev (2, NA_REAL);
    NumericVector ESS_repst_dev (2, NA_REAL);
    NumericVector ESS_fec_dev (2, NA_REAL);
    
    NumericVector ESS_jsurv_dev (2, NA_REAL);
    NumericVector ESS_jobs_dev (2, NA_REAL);
    NumericVector ESS_jsize_dev (2, NA_REAL);
    NumericVector ESS_jsizeb_dev (2, NA_REAL);
    NumericVector ESS_jsizec_dev (2, NA_REAL);
    NumericVector ESS_jrepst_dev (2, NA_REAL);
    NumericVector ESS_jmatst_dev (2, NA_REAL);
    
    // Vectors supplied in ta_reassessed
    IntegerVector variant = as<IntegerVector>(ta_reassessed["variant"]);
    
    NumericVector givenrate = as<NumericVector>(ta_reassessed["givenrate"]);
    NumericVector offset = as<NumericVector>(ta_reassessed["offset"]);
    NumericVector multiplier = as<NumericVector>(ta_reassessed["multiplier"]);
    
    NumericVector surv_dev = as<NumericVector>(ta_reassessed["surv_dev"]);
    NumericVector obs_dev = as<NumericVector>(ta_reassessed["obs_dev"]);
    NumericVector size_dev = as<NumericVector>(ta_reassessed["size_dev"]);
    NumericVector sizeb_dev = as<NumericVector>(ta_reassessed["sizeb_dev"]);
    NumericVector sizec_dev = as<NumericVector>(ta_reassessed["sizec_dev"]);
    NumericVector repst_dev = as<NumericVector>(ta_reassessed["repst_dev"]);
    NumericVector fec_dev = as<NumericVector>(ta_reassessed["fec_dev"]);
    NumericVector jsurv_dev = as<NumericVector>(ta_reassessed["jsurv_dev"]);
    NumericVector jobs_dev = as<NumericVector>(ta_reassessed["jobs_dev"]);
    NumericVector jsize_dev = as<NumericVector>(ta_reassessed["jsize_dev"]);
    NumericVector jsizeb_dev = as<NumericVector>(ta_reassessed["jsizeb_dev"]);
    NumericVector jsizec_dev = as<NumericVector>(ta_reassessed["jsizec_dev"]);
    NumericVector jrepst_dev = as<NumericVector>(ta_reassessed["jrepst_dev"]);
    NumericVector jmatst_dev = as<NumericVector>(ta_reassessed["jmatst_dev"]);
    
    StringVector stage3 = as<StringVector>(ta_reassessed["stage3"]);
    StringVector stage2 = as<StringVector>(ta_reassessed["stage2"]);
    StringVector stage1 = as<StringVector>(ta_reassessed["stage1"]);
    IntegerVector age3 = as<IntegerVector>(ta_reassessed["age3"]);
    IntegerVector age2 = as<IntegerVector>(ta_reassessed["age2"]);
    StringVector eststage3 = as<StringVector>(ta_reassessed["eststage3"]);
    StringVector eststage2 = as<StringVector>(ta_reassessed["eststage2"]);
    StringVector eststage1 = as<StringVector>(ta_reassessed["eststage1"]);
    IntegerVector estage3 = as<IntegerVector>(ta_reassessed["estage3"]);
    IntegerVector estage2 = as<IntegerVector>(ta_reassessed["estage2"]);
    IntegerVector convtype = as<IntegerVector>(ta_reassessed["convtype"]);
    IntegerVector convtype_t12 = as<IntegerVector>(ta_reassessed["convtype_t12"]);
    StringVector year2 = as<StringVector>(ta_reassessed["year2"]);
    IntegerVector mpm_altered = as<IntegerVector>(ta_reassessed["mpm_altered"]);
    IntegerVector vrm_altered = as<IntegerVector>(ta_reassessed["vrm_altered"]);
    
    NumericVector givenrate_noNA = na_omit(givenrate);
    NumericVector offset_noNA = na_omit(offset);
    NumericVector multiplier_noNA = na_omit(multiplier);
    NumericVector surv_dev_noNA = na_omit(surv_dev);
    NumericVector obs_dev_noNA = na_omit(obs_dev);
    NumericVector size_dev_noNA = na_omit(size_dev);
    NumericVector sizeb_dev_noNA = na_omit(sizeb_dev);
    NumericVector sizec_dev_noNA = na_omit(sizec_dev);
    NumericVector repst_dev_noNA = na_omit(repst_dev);
    NumericVector fec_dev_noNA = na_omit(fec_dev);
    NumericVector jsurv_dev_noNA = na_omit(jsurv_dev);
    NumericVector jobs_dev_noNA = na_omit(jobs_dev);
    NumericVector jsize_dev_noNA = na_omit(jsize_dev);
    NumericVector jsizeb_dev_noNA = na_omit(jsizeb_dev);
    NumericVector jsizec_dev_noNA = na_omit(jsizec_dev);
    NumericVector jrepst_dev_noNA = na_omit(jrepst_dev);
    NumericVector jmatst_dev_noNA = na_omit(jmatst_dev);
    
    bool givenrate_NA {false};
    bool multiplier_NA {false};
    bool offset_NA {false};
    bool surv_dev_NA {false};
    bool obs_dev_NA {false};
    bool size_dev_NA {false};
    bool sizeb_dev_NA {false};
    bool sizec_dev_NA {false};
    bool repst_dev_NA {false};
    bool fec_dev_NA {false};
    bool jsurv_dev_NA {false};
    bool jobs_dev_NA {false};
    bool jsize_dev_NA {false};
    bool jsizeb_dev_NA {false};
    bool jsizec_dev_NA {false};
    bool jrepst_dev_NA {false};
    bool jmatst_dev_NA {false};
    
    int var1_length {0};
    
    for (int i = 0; i < static_cast<int>(variant.length()); i++) {
      if (variant(i) == 1) {
        var1_length++;
      }
    }
    IntegerVector var1_bounds (var1_length);
    int var1_bounds_counter {0};
    for (int i = 0; i < static_cast<int>(variant.length()); i++) {
      if (variant(i) == 1) {
        var1_bounds(var1_bounds_counter) = i;
        var1_bounds_counter++;
      }
    }
    
    if (var1_length == 0) {
      throw Rcpp::exception("Cannot locate variant 1. Cannot determine variant dimensions.", false);
    }
    
    double givenrate_min {0.};
    double givenrate_max {0.};
    double offset_min {0.};
    double offset_max {0.};
    double multiplier_min {0.};
    double multiplier_max {0.};
    
    double surv_dev_min {0.};
    double surv_dev_max {0.};
    double obs_dev_min {0.};
    double obs_dev_max {0.};
    double size_dev_min {0.};
    double size_dev_max {0.};
    double sizeb_dev_min {0.};
    double sizeb_dev_max {0.};
    double sizec_dev_min {0.};
    double sizec_dev_max {0.};
    double repst_dev_min {0.};
    double repst_dev_max {0.};
    double fec_dev_min {0.};
    double fec_dev_max {0.};
    
    double jsurv_dev_min {0.};
    double jsurv_dev_max {0.};
    double jobs_dev_min {0.};
    double jobs_dev_max {0.};
    double jsize_dev_min {0.};
    double jsize_dev_max {0.};
    double jsizeb_dev_min {0.};
    double jsizeb_dev_max {0.};
    double jsizec_dev_min {0.};
    double jsizec_dev_max {0.};
    double jrepst_dev_min {0.};
    double jrepst_dev_max {0.};
    double jmatst_dev_min {0.};
    double jmatst_dev_max {0.};
    
    bool givenrate_reversed {false};
    bool offset_reversed {false};
    bool multiplier_reversed {false};
    
    bool surv_dev_reversed {false};
    bool obs_dev_reversed {false};
    bool size_dev_reversed {false};
    bool sizeb_dev_reversed {false};
    bool sizec_dev_reversed {false};
    bool repst_dev_reversed {false};
    bool fec_dev_reversed {false};
    
    bool jsurv_dev_reversed {false};
    bool jobs_dev_reversed {false};
    bool jsize_dev_reversed {false};
    bool jsizeb_dev_reversed {false};
    bool jsizec_dev_reversed {false};
    bool jrepst_dev_reversed {false};
    bool jmatst_dev_reversed {false};
    
    if (static_cast<int>(givenrate_noNA.length()) > 0) {
      givenrate_min = min(givenrate_noNA);
      givenrate_max = max(givenrate_noNA);
      
      unsigned int givenrate_min_pos = which_min(givenrate_noNA);
      unsigned int givenrate_max_pos = which_max(givenrate_noNA);
      
      if (givenrate_min_pos < givenrate_max_pos) {
        ESS_givenrate(0) = givenrate_min;
        ESS_givenrate(1) = givenrate_max;
      } else {
        ESS_givenrate(0) = givenrate_max;
        ESS_givenrate(1) = givenrate_min;
        
        givenrate_reversed = true;
      }
    } else {
      givenrate_NA = true;
    }
    
    if (static_cast<int>(offset_noNA.length()) > 0) {
      offset_min = min(offset_noNA);
      offset_max = max(offset_noNA);
      
      unsigned int offset_min_pos = which_min(offset_noNA);
      unsigned int offset_max_pos = which_max(offset_noNA);
      
      if (offset_min_pos < offset_max_pos) {
        ESS_offset(0) = offset_min;
        ESS_offset(1) = offset_max;
      } else {
        ESS_offset(0) = offset_max;
        ESS_offset(1) = offset_min;
        
        offset_reversed = true;
      }
    } else {
      offset_NA = true;
    }
    
    if (static_cast<int>(multiplier_noNA.length()) > 0) {
      multiplier_min = min(multiplier_noNA);
      multiplier_max = max(multiplier_noNA);
      
      unsigned int multiplier_min_pos = which_min(multiplier_noNA);
      unsigned int multiplier_max_pos = which_max(multiplier_noNA);
      
      if (multiplier_min_pos < multiplier_max_pos) {
        ESS_multiplier(0) = multiplier_min;
        ESS_multiplier(1) = multiplier_max;
      } else {
        ESS_multiplier(0) = multiplier_max;
        ESS_multiplier(1) = multiplier_min;
        
        multiplier_reversed = true;
      }
    } else {
      multiplier_NA = true;
    }
    
    if (static_cast<int>(surv_dev_noNA.length()) > 0) {
      surv_dev_min = min(surv_dev_noNA);
      surv_dev_max = max(surv_dev_noNA);
      
      unsigned int surv_dev_min_pos = which_min(surv_dev_noNA);
      unsigned int surv_dev_max_pos = which_max(surv_dev_noNA);
      
      if (surv_dev_min_pos < surv_dev_max_pos) {
        ESS_surv_dev(0) = surv_dev_min;
        ESS_surv_dev(1) = surv_dev_max;
      } else {
        ESS_surv_dev(0) = surv_dev_max;
        ESS_surv_dev(1) = surv_dev_min;
        
        surv_dev_reversed = true;
      }
    } else {
      surv_dev_NA = true;
    }
    
    if (static_cast<int>(obs_dev_noNA.length()) > 0) {
      obs_dev_min = min(obs_dev_noNA);
      obs_dev_max = max(obs_dev_noNA);
      
      unsigned int obs_dev_min_pos = which_min(obs_dev_noNA);
      unsigned int obs_dev_max_pos = which_max(obs_dev_noNA);
      
      if (obs_dev_min_pos < obs_dev_max_pos) {
        ESS_obs_dev(0) = obs_dev_min;
        ESS_obs_dev(1) = obs_dev_max;
      } else {
        ESS_obs_dev(0) = obs_dev_max;
        ESS_obs_dev(1) = obs_dev_min;
        
        obs_dev_reversed = true;
      }
    } else {
      obs_dev_NA = true;
    }
    
    
    if (static_cast<int>(size_dev_noNA.length()) > 0) {
      size_dev_min = min(size_dev_noNA);
      size_dev_max = max(size_dev_noNA);
      
      unsigned int size_dev_min_pos = which_min(size_dev_noNA);
      unsigned int size_dev_max_pos = which_max(size_dev_noNA);
      
      if (size_dev_min_pos < size_dev_max_pos) {
        ESS_size_dev(0) = size_dev_min;
        ESS_size_dev(1) = size_dev_max;
      } else {
        ESS_size_dev(0) = size_dev_max;
        ESS_size_dev(1) = size_dev_min;
        
        size_dev_reversed = true;
      }
    } else {
      size_dev_NA = true;
    }
    
    if (static_cast<int>(sizeb_dev_noNA.length()) > 0) {
      sizeb_dev_min = min(sizeb_dev_noNA);
      sizeb_dev_max = max(sizeb_dev_noNA);
      
      unsigned int sizeb_dev_min_pos = which_min(sizeb_dev_noNA);
      unsigned int sizeb_dev_max_pos = which_max(sizeb_dev_noNA);
      
      if (sizeb_dev_min_pos < sizeb_dev_max_pos) {
        ESS_sizeb_dev(0) = sizeb_dev_min;
        ESS_sizeb_dev(1) = sizeb_dev_max;
      } else {
        ESS_sizeb_dev(0) = sizeb_dev_max;
        ESS_sizeb_dev(1) = sizeb_dev_min;
        
        sizeb_dev_reversed = true;
      }
    } else {
      sizeb_dev_NA = true;
    }
    
    if (static_cast<int>(sizec_dev_noNA.length()) > 0) {
      sizec_dev_min = min(sizec_dev_noNA);
      sizec_dev_max = max(sizec_dev_noNA);
      
      unsigned int sizec_dev_min_pos = which_min(sizec_dev_noNA);
      unsigned int sizec_dev_max_pos = which_max(sizec_dev_noNA);
      
      if (sizec_dev_min_pos < sizec_dev_max_pos) {
        ESS_sizec_dev(0) = sizec_dev_min;
        ESS_sizec_dev(1) = sizec_dev_max;
      } else {
        ESS_sizec_dev(0) = sizec_dev_max;
        ESS_sizec_dev(1) = sizec_dev_min;
        
        sizec_dev_reversed = true;
      }
    } else {
      sizec_dev_NA = true;
    }
    
    if (static_cast<int>(repst_dev_noNA.length()) > 0) {
      repst_dev_min = min(repst_dev_noNA);
      repst_dev_max = max(repst_dev_noNA);
      
      unsigned int repst_dev_min_pos = which_min(repst_dev_noNA);
      unsigned int repst_dev_max_pos = which_max(repst_dev_noNA);
      
      if (repst_dev_min_pos < repst_dev_max_pos) {
        ESS_repst_dev(0) = repst_dev_min;
        ESS_repst_dev(1) = repst_dev_max;
      } else {
        ESS_repst_dev(0) = repst_dev_max;
        ESS_repst_dev(1) = repst_dev_min;
        
        repst_dev_reversed = true;
      }
    } else {
      repst_dev_NA = true;
    }
    
    if (static_cast<int>(fec_dev_noNA.length()) > 0) {
      fec_dev_min = min(fec_dev_noNA);
      fec_dev_max = max(fec_dev_noNA);
      
      unsigned int fec_dev_min_pos = which_min(fec_dev_noNA);
      unsigned int fec_dev_max_pos = which_max(fec_dev_noNA);
      
      if (fec_dev_min_pos < fec_dev_max_pos) {
        ESS_fec_dev(0) = fec_dev_min;
        ESS_fec_dev(1) = fec_dev_max;
      } else {
        ESS_fec_dev(0) = fec_dev_max;
        ESS_fec_dev(1) = fec_dev_min;
        
        fec_dev_reversed = true;
      }
    } else {
      fec_dev_NA = true;
    }
    
    if (static_cast<int>(jsurv_dev_noNA.length()) > 0) {
      jsurv_dev_min = min(jsurv_dev_noNA);
      jsurv_dev_max = max(jsurv_dev_noNA);
      
      unsigned int jsurv_dev_min_pos = which_min(jsurv_dev_noNA);
      unsigned int jsurv_dev_max_pos = which_max(jsurv_dev_noNA);
      
      if (jsurv_dev_min_pos < jsurv_dev_max_pos) {
        ESS_jsurv_dev(0) = jsurv_dev_min;
        ESS_jsurv_dev(1) = jsurv_dev_max;
      } else {
        ESS_jsurv_dev(0) = jsurv_dev_max;
        ESS_jsurv_dev(1) = jsurv_dev_min;
        
        jsurv_dev_reversed = true;
      }
    } else {
      jsurv_dev_NA = true;
    }
    
    if (static_cast<int>(jobs_dev_noNA.length()) > 0) {
      jobs_dev_min = min(jobs_dev_noNA);
      jobs_dev_max = max(jobs_dev_noNA);
      
      unsigned int jobs_dev_min_pos = which_min(jobs_dev_noNA);
      unsigned int jobs_dev_max_pos = which_max(jobs_dev_noNA);
      
      if (jobs_dev_min_pos < jobs_dev_max_pos) {
        ESS_jobs_dev(0) = jobs_dev_min;
        ESS_jobs_dev(1) = jobs_dev_max;
      } else {
        ESS_jobs_dev(0) = jobs_dev_max;
        ESS_jobs_dev(1) = jobs_dev_min;
        
        jobs_dev_reversed = true;
      }
    } else {
      jobs_dev_NA = true;
    }
    
    if (static_cast<int>(jsize_dev_noNA.length()) > 0) {
      jsize_dev_min = min(jsize_dev_noNA);
      jsize_dev_max = max(jsize_dev_noNA);
      
      unsigned int jsize_dev_min_pos = which_min(jsize_dev_noNA);
      unsigned int jsize_dev_max_pos = which_max(jsize_dev_noNA);
      
      if (jsize_dev_min_pos < jsize_dev_max_pos) {
        ESS_jsize_dev(0) = jsize_dev_min;
        ESS_jsize_dev(1) = jsize_dev_max;
      } else {
        ESS_jsize_dev(0) = jsize_dev_max;
        ESS_jsize_dev(1) = jsize_dev_min;
        
        jsize_dev_reversed = true;
      }
    } else {
      jsize_dev_NA = true;
    }
    
    if (static_cast<int>(jsizeb_dev_noNA.length()) > 0) {
      jsizeb_dev_min = min(jsizeb_dev_noNA);
      jsizeb_dev_max = max(jsizeb_dev_noNA);
      
      unsigned int jsizeb_dev_min_pos = which_min(jsizeb_dev_noNA);
      unsigned int jsizeb_dev_max_pos = which_max(jsizeb_dev_noNA);
      
      if (jsizeb_dev_min_pos < jsizeb_dev_max_pos) {
        ESS_jsizeb_dev(0) = jsizeb_dev_min;
        ESS_jsizeb_dev(1) = jsizeb_dev_max;
      } else {
        ESS_jsizeb_dev(0) = jsizeb_dev_max;
        ESS_jsizeb_dev(1) = jsizeb_dev_min;
        
        jsizeb_dev_reversed = true;
      }
    } else {
      jsizeb_dev_NA = true;
    }
    
    if (static_cast<int>(jsizec_dev_noNA.length()) > 0) {
      jsizec_dev_min = min(jsizec_dev_noNA);
      jsizec_dev_max = max(jsizec_dev_noNA);
      
      unsigned int jsizec_dev_min_pos = which_min(jsizec_dev_noNA);
      unsigned int jsizec_dev_max_pos = which_max(jsizec_dev_noNA);
      
      if (jsizec_dev_min_pos < jsizec_dev_max_pos) {
        ESS_jsizec_dev(0) = jsizec_dev_min;
        ESS_jsizec_dev(1) = jsizec_dev_max;
      } else {
        ESS_jsizec_dev(0) = jsizec_dev_max;
        ESS_jsizec_dev(1) = jsizec_dev_min;
        
        jsizec_dev_reversed = true;
      }
    } else {
      jsizec_dev_NA = true;
    }
    
    if (static_cast<int>(jrepst_dev_noNA.length()) > 0) {
      jrepst_dev_min = min(jrepst_dev_noNA);
      jrepst_dev_max = max(jrepst_dev_noNA);
      
      unsigned int jrepst_dev_min_pos = which_min(jrepst_dev_noNA);
      unsigned int jrepst_dev_max_pos = which_max(jrepst_dev_noNA);
      
      if (jrepst_dev_min_pos < jrepst_dev_max_pos) {
        ESS_jrepst_dev(0) = jrepst_dev_min;
        ESS_jrepst_dev(1) = jrepst_dev_max;
      } else {
        ESS_jrepst_dev(0) = jrepst_dev_max;
        ESS_jrepst_dev(1) = jrepst_dev_min;
        
        jrepst_dev_reversed = true;
      }
    } else {
      jrepst_dev_NA = true;
    }
    
    if (static_cast<int>(jmatst_dev_noNA.length()) > 0) {
      jmatst_dev_min = min(jmatst_dev_noNA);
      jmatst_dev_max = max(jmatst_dev_noNA);
      
      unsigned int jmatst_dev_min_pos = which_min(jmatst_dev_noNA);
      unsigned int jmatst_dev_max_pos = which_max(jmatst_dev_noNA);
      
      if (jmatst_dev_min_pos < jmatst_dev_max_pos) {
        ESS_jmatst_dev(0) = jmatst_dev_min;
        ESS_jmatst_dev(1) = jmatst_dev_max;
      } else {
        ESS_jmatst_dev(0) = jmatst_dev_max;
        ESS_jmatst_dev(1) = jmatst_dev_min;
        
        jmatst_dev_reversed = true;
      }
    } else {
      jmatst_dev_NA = true;
    }
    
    
    int found_variables_all {0};
    int found_variables {0};
    bool first_var_found {false};
    
    arma::ivec mpm_altered_arma = as<arma::ivec>(mpm_altered);
    arma::uvec found_mpm_altered = find(mpm_altered_arma);
    int input_ta_length = static_cast<int>(stage3.length());
    
    StringVector core_stage_index (input_ta_length);
    if (static_cast<int>(found_mpm_altered.n_elem) > 0) {
      //Rcout << "using mpm_altered if-else route A" << endl;
      for (int i = 0; i < input_ta_length; i++) {
        String a3 = String(stage3(i));
        String a2 = String(stage2(i));
        String a1 = String(stage1(i));
        String ag2;
        if (IntegerVector::is_na(age2(i))) {
          ag2 = "0 ";
        } else {
          std::string ag2_conv = std::to_string(static_cast<int>(age2(i)));
          ag2 = String(ag2_conv);
        }
        if (a3 == "NA_STRING") a3 = "0 ";
        if (a2 == "NA_STRING") a2 = "0 ";
        if (a1 == "NA_STRING") a1 = "0 ";
        
        String new_string_code = a3;
        new_string_code += " ";
        new_string_code += a2;
        new_string_code += " ";
        new_string_code += a1;
        new_string_code += " ";
        new_string_code += ag2;
        
        core_stage_index(i) = new_string_code;
      }
    } else {
      //Rcout << "using mpm_altered if-else route B" << endl;
      for (int i = 0; i < input_ta_length; i++) {
        core_stage_index(i) = "1";
      }
    }
    
    //Rcout << "core_stage_index: " << core_stage_index << endl;
    
    StringVector unique_core_stages = unique(core_stage_index);
    int found_core_stage_indices = unique_core_stages.length();
    
    arma::mat ta_rows_per_change_per_variant (found_core_stage_indices, 17, fill::zeros); // rows = var_rows, cols = vital rates
    IntegerVector place_holder_for_mat (found_core_stage_indices, 0);
    
    StringVector new_stage3_found_variables (found_core_stage_indices);
    StringVector new_stage2_found_variables (found_core_stage_indices);
    StringVector new_stage1_found_variables (found_core_stage_indices);
    IntegerVector new_age3_found_variables (found_core_stage_indices);
    IntegerVector new_age2_found_variables (found_core_stage_indices);
    
    StringVector new_eststage3_found_variables (found_core_stage_indices);
    StringVector new_eststage2_found_variables (found_core_stage_indices);
    StringVector new_eststage1_found_variables (found_core_stage_indices);
    IntegerVector new_estage3_found_variables (found_core_stage_indices);
    IntegerVector new_estage2_found_variables (found_core_stage_indices);
    
    IntegerVector new_convtype_found_variables (found_core_stage_indices);
    IntegerVector new_convtype_t12_found_variables (found_core_stage_indices);
    StringVector new_year2_found_variables (found_core_stage_indices);
    IntegerVector new_mpm_altered_found_variables (found_core_stage_indices);
    IntegerVector new_vrm_altered_found_variables (found_core_stage_indices);
    
    IntegerVector found_core_stage_index_vector (input_ta_length);
    for (int i = 0; i < input_ta_length; i++) {
      for (int j = 0; j < found_core_stage_indices; j++) {
        if (core_stage_index(i) == unique_core_stages(j)) {
          new_stage3_found_variables(j) = stage3(i);
          new_stage2_found_variables(j) = stage2(i);
          new_stage1_found_variables(j) = stage1(i);
          new_age3_found_variables(j) = age3(i);
          new_age2_found_variables(j) = age2(i);
          
          new_eststage3_found_variables(j) = eststage3(i);
          new_eststage2_found_variables(j) = eststage2(i);
          new_eststage1_found_variables(j) = eststage1(i);
          new_estage3_found_variables(j) = estage3(i);
          new_estage2_found_variables(j) = estage2(i);
          
          new_convtype_found_variables(j) = convtype(i);
          new_convtype_t12_found_variables(j) = convtype_t12(i);
          new_year2_found_variables(j) = year2(i);
          new_mpm_altered_found_variables(j) = mpm_altered(i);
          new_vrm_altered_found_variables(j) = vrm_altered(i);
          
          found_core_stage_index_vector(i) = j + 1;
          place_holder_for_mat(j) = place_holder_for_mat(j) + 1;
        }
      }
    }
    
    // mpm-altered only
    if (givenrate_min != givenrate_max) {
      ESS_variable_traits(0) = found_core_stage_indices;
      found_variables_all = found_variables_all + found_core_stage_indices;
      found_variables++;
      
      for (int i = 0; i < found_core_stage_indices; i++) {
        ta_rows_per_change_per_variant(i, 0) = place_holder_for_mat(i);
      }
    }
    if (offset_min != offset_max) {
      ESS_variable_traits(1) = found_core_stage_indices;
      found_variables_all = found_variables_all + found_core_stage_indices;
      found_variables++;
      
      for (int i = 0; i < found_core_stage_indices; i++) {
        ta_rows_per_change_per_variant(i, 1) = place_holder_for_mat(i);
      }
    }
    if (multiplier_min != multiplier_max) {
      ESS_variable_traits(2) = found_core_stage_indices;
      found_variables_all = found_variables_all + found_core_stage_indices;
      found_variables++;
      
      for (int i = 0; i < found_core_stage_indices; i++) {
        ta_rows_per_change_per_variant(i, 2) = place_holder_for_mat(i);
      }
    }
    
    // vrm-altered only
    if (surv_dev_min != surv_dev_max) {
      ESS_variable_traits(3) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 3) = 1;
    }
    if (obs_dev_min != obs_dev_max) {
      ESS_variable_traits(4) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 4) = 1;
    }
    if (size_dev_min != size_dev_max) {
      ESS_variable_traits(5) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 5) = 1;
    }
    if (sizeb_dev_min != sizeb_dev_max) {
      ESS_variable_traits(6) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 6) = 1;
    }
    if (sizec_dev_min != sizec_dev_max) {
      ESS_variable_traits(7) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 7) = 1;
    }
    if (repst_dev_min != repst_dev_max) {
      ESS_variable_traits(8) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 8) = 1;
    }
    if (fec_dev_min != fec_dev_max) {
      ESS_variable_traits(9) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 9) = 1;
    }
    if (jsurv_dev_min != jsurv_dev_max) {
      ESS_variable_traits(10) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 10) = 1;
    }
    if (jobs_dev_min != jobs_dev_max) {
      ESS_variable_traits(11) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 11) = 1;
    }
    if (jsize_dev_min != jsize_dev_max) {
      ESS_variable_traits(12) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 12) = 1;
    }
    if (jsizeb_dev_min != jsizeb_dev_max) {
      ESS_variable_traits(13) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 13) = 1;
    }
    if (jsizec_dev_min != jsizec_dev_max) {
      ESS_variable_traits(14) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 14) = 1;
    }
    if (jrepst_dev_min != jrepst_dev_max) {
      ESS_variable_traits(15) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 15) = 1;
    }
    if (jmatst_dev_min != jmatst_dev_max) {
      ESS_variable_traits(16) = 1;
      found_variables_all++;
      found_variables++;
      
      ta_rows_per_change_per_variant(0, 16) = 1;
    }
    
    int data_frame_length = opt_res * var1_length;
    
    int ta_re_vars = ta_reassessed.length();
    int ta_re_nrows = ta_reassessed.nrows();
    int ESS_ta_vars = ESS_ta.length();
    int ESS_ta_nrows = ESS_ta.nrows();
    int optim_ta_vars = optim_ta.length();
    int optim_ta_nrows = optim_ta.nrows();
    int optim_ta_995_vars = optim_ta_995.length();
    int optim_ta_995_nrows = optim_ta_995.nrows();
    
    cloned_ta_reassess = clone(ta_reassessed);
    
    optim_variant = as<IntegerVector>(cloned_ta_reassess["variant"]);
    optim_stage3 = as<CharacterVector>(cloned_ta_reassess["stage3"]);
    optim_stage2 = as<CharacterVector>(cloned_ta_reassess["stage2"]);
    optim_stage1 = as<CharacterVector>(cloned_ta_reassess["stage1"]);
    optim_age3 = as<IntegerVector>(cloned_ta_reassess["age3"]);
    optim_age2 = as<IntegerVector>(cloned_ta_reassess["age2"]);
    optim_eststage3 = as<CharacterVector>(cloned_ta_reassess["eststage3"]);
    optim_eststage2 = as<CharacterVector>(cloned_ta_reassess["eststage2"]);
    optim_eststage1 = as<CharacterVector>(cloned_ta_reassess["eststage1"]);
    optim_estage3 = as<IntegerVector>(cloned_ta_reassess["estage3"]);
    optim_estage2 = as<IntegerVector>(cloned_ta_reassess["estage2"]);
    optim_convtype = as<IntegerVector>(cloned_ta_reassess["convtype"]);
    optim_convtype_t12 = as<IntegerVector>(cloned_ta_reassess["convtype_t12"]);
    optim_year2 = as<CharacterVector>(cloned_ta_reassess["year2"]);
    optim_mpm_altered = as<IntegerVector>(cloned_ta_reassess["mpm_altered"]);
    optim_vrm_altered = as<IntegerVector>(cloned_ta_reassess["vrm_altered"]);
    
    opt_res = static_cast<int>(optim_stage3.length());
    
    for (int i = 0; i < 2; i++) {
      ESS_stage3(i) = stage3(0);
      ESS_stage2(i) = stage2(0);
      ESS_stage1(i) = stage1(0);
      ESS_age3(i) = age3(0);
      ESS_age2(i) = age2(0);
      ESS_eststage3(i) = eststage3(0);
      ESS_eststage2(i) = eststage2(0);
      ESS_eststage1(i) = eststage1(0);
      ESS_estage3(i) = estage3(0);
      ESS_estage2(i) = estage2(0);
      ESS_convtype(i) = convtype(0);
      ESS_convtype_t12(i) = convtype_t12(0);
      ESS_year2(i) = year2(0);
      ESS_mpm_altered(i) = mpm_altered(0);
      ESS_vrm_altered(i) = vrm_altered(0);
    }
    
    // Establish core values in variable vs. constant traits
    optim_givenrate = as<NumericVector>(cloned_ta_reassess["givenrate"]);
    optim_offset = as<NumericVector>(cloned_ta_reassess["offset"]);
    optim_multiplier = as<NumericVector>(cloned_ta_reassess["multiplier"]);
    optim_surv_dev = as<NumericVector>(cloned_ta_reassess["surv_dev"]);
    optim_obs_dev = as<NumericVector>(cloned_ta_reassess["obs_dev"]);
    optim_size_dev = as<NumericVector>(cloned_ta_reassess["size_dev"]);
    optim_sizeb_dev = as<NumericVector>(cloned_ta_reassess["sizeb_dev"]);
    optim_sizec_dev = as<NumericVector>(cloned_ta_reassess["sizec_dev"]);
    optim_repst_dev = as<NumericVector>(cloned_ta_reassess["repst_dev"]);
    optim_fec_dev = as<NumericVector>(cloned_ta_reassess["fec_dev"]);
    optim_jsurv_dev = as<NumericVector>(cloned_ta_reassess["jsurv_dev"]);
    optim_jobs_dev = as<NumericVector>(cloned_ta_reassess["jobs_dev"]);
    optim_jsize_dev = as<NumericVector>(cloned_ta_reassess["jsize_dev"]);
    optim_jsizeb_dev = as<NumericVector>(cloned_ta_reassess["jsizeb_dev"]);
    optim_jsizec_dev = as<NumericVector>(cloned_ta_reassess["jsizec_dev"]);
    optim_jrepst_dev = as<NumericVector>(cloned_ta_reassess["jrepst_dev"]);
    optim_jmatst_dev = as<NumericVector>(cloned_ta_reassess["jmatst_dev"]);
    
    optim_givenrate_995 = clone(optim_givenrate);
    optim_offset_995 = clone(optim_offset);
    optim_multiplier_995 = clone(optim_multiplier);
    optim_surv_dev_995 = clone(optim_surv_dev);
    optim_obs_dev_995 = clone(optim_obs_dev);
    optim_size_dev_995 = clone(optim_size_dev);
    optim_sizeb_dev_995 = clone(optim_sizeb_dev);
    optim_sizec_dev_995 = clone(optim_sizec_dev);
    optim_repst_dev_995 = clone(optim_repst_dev);
    optim_fec_dev_995 = clone(optim_fec_dev);
    optim_jsurv_dev_995 = clone(optim_jsurv_dev);
    optim_jobs_dev_995 = clone(optim_jobs_dev);
    optim_jsize_dev_995 = clone(optim_jsize_dev);
    optim_jsizeb_dev_995 = clone(optim_jsizeb_dev);
    optim_jsizec_dev_995 = clone(optim_jsizec_dev);
    optim_jrepst_dev_995 = clone(optim_jrepst_dev);
    optim_jmatst_dev_995 = clone(optim_jmatst_dev);
    
    
    if (ESS_variable_traits(0) > 0) optim_givenrate_995 = optim_givenrate_995 * elast_multiplier;
    if (ESS_variable_traits(1) > 0) optim_offset_995 = optim_offset_995 * elast_multiplier;
    if (ESS_variable_traits(2) > 0) optim_multiplier_995 = optim_multiplier_995 * elast_multiplier;
    if (ESS_variable_traits(3) > 0) optim_surv_dev_995 = optim_surv_dev_995 * elast_multiplier;
    if (ESS_variable_traits(4) > 0) optim_obs_dev_995 = optim_obs_dev_995 * elast_multiplier;
    if (ESS_variable_traits(5) > 0) optim_size_dev_995 = optim_size_dev_995 * elast_multiplier;
    if (ESS_variable_traits(6) > 0) optim_sizeb_dev_995 = optim_sizeb_dev_995 * elast_multiplier;
    if (ESS_variable_traits(7) > 0) optim_sizec_dev_995 = optim_sizec_dev_995 * elast_multiplier;
    if (ESS_variable_traits(8) > 0) optim_repst_dev_995 = optim_repst_dev_995 * elast_multiplier;
    if (ESS_variable_traits(9) > 0) optim_fec_dev_995 = optim_fec_dev_995 * elast_multiplier;
    if (ESS_variable_traits(10) > 0) optim_jsurv_dev_995 = optim_jsurv_dev_995 * elast_multiplier;
    if (ESS_variable_traits(11) > 0) optim_jobs_dev_995 = optim_jobs_dev_995 * elast_multiplier;
    if (ESS_variable_traits(12) > 0) optim_jsize_dev_995 = optim_jsize_dev_995 * elast_multiplier;
    if (ESS_variable_traits(13) > 0) optim_jsizeb_dev_995 = optim_jsizeb_dev_995 * elast_multiplier;
    if (ESS_variable_traits(14) > 0) optim_jsizec_dev_995 = optim_jsizec_dev_995 * elast_multiplier;
    if (ESS_variable_traits(15) > 0) optim_jrepst_dev_995 = optim_jrepst_dev_995 * elast_multiplier;
    if (ESS_variable_traits(16) > 0) optim_jmatst_dev_995 = optim_jmatst_dev_995 * elast_multiplier;
    
    optim_variant_995 = clone(optim_variant);
    
    optim_stage3_995 = clone(optim_stage3);
    optim_stage2_995 = clone(optim_stage2);
    optim_stage1_995 = clone(optim_stage1);
    optim_age3_995 = clone(optim_age3);
    optim_age2_995 = clone(optim_age2);
    optim_eststage3_995 = clone(optim_eststage3);
    optim_eststage2_995 = clone(optim_eststage2);
    optim_eststage1_995 = clone(optim_eststage1);
    optim_estage3_995 = clone(optim_estage3);
    optim_estage2_995 = clone(optim_estage2);
    optim_convtype_995 = clone(optim_convtype);
    optim_convtype_t12_995 = clone(optim_convtype_t12);
    optim_year2_995 = clone(optim_year2);
    optim_mpm_altered_995 = clone(optim_mpm_altered);
    optim_vrm_altered_995 = clone(optim_vrm_altered);
    
    List output (33);
    
    output(0) = optim_variant;
    output(1) = optim_stage3;
    output(2) = optim_stage2;
    output(3) = optim_stage1;
    output(4) = optim_age3;
    output(5) = optim_age2;
    output(6) = optim_eststage3;
    output(7) = optim_eststage2;
    output(8) = optim_eststage1;
    output(9) = optim_estage3;
    output(10) = optim_estage2;
    output(11) = optim_givenrate;
    output(12) = optim_offset;
    output(13) = optim_multiplier;
    output(14) = optim_convtype;
    output(15) = optim_convtype_t12;
    output(16) = optim_surv_dev;
    output(17) = optim_obs_dev;
    output(18) = optim_size_dev;
    output(19) = optim_sizeb_dev;
    output(20) = optim_sizec_dev;
    output(21) = optim_repst_dev;
    output(22) = optim_fec_dev;
    output(23) = optim_jsurv_dev;
    output(24) = optim_jobs_dev;
    output(25) = optim_jsize_dev;
    output(26) = optim_jsizeb_dev;
    output(27) = optim_jsizec_dev;
    output(28) = optim_jrepst_dev;
    output(29) = optim_jmatst_dev;
    output(30) = optim_year2;
    output(31) = optim_mpm_altered;
    output(32) = optim_vrm_altered;
    
    CharacterVector ta_names = as<CharacterVector>(ta_reassessed.attr("names"));
    output.attr("names") = clone(ta_names);
    output.attr("class") = "data.frame";
    output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, data_frame_length);
    
    optim_ta = output;
    
    List output_995 (33);
    
    output_995(0) = optim_variant_995;
    output_995(1) = optim_stage3_995;
    output_995(2) = optim_stage2_995;
    output_995(3) = optim_stage1_995;
    output_995(4) = optim_age3_995;
    output_995(5) = optim_age2_995;
    output_995(6) = optim_eststage3_995;
    output_995(7) = optim_eststage2_995;
    output_995(8) = optim_eststage1_995;
    output_995(9) = optim_estage3_995;
    output_995(10) = optim_estage2_995;
    output_995(11) = optim_givenrate_995;
    output_995(12) = optim_offset_995;
    output_995(13) = optim_multiplier_995;
    output_995(14) = optim_convtype_995;
    output_995(15) = optim_convtype_t12_995;
    output_995(16) = optim_surv_dev_995;
    output_995(17) = optim_obs_dev_995;
    output_995(18) = optim_size_dev_995;
    output_995(19) = optim_sizeb_dev_995;
    output_995(20) = optim_sizec_dev_995;
    output_995(21) = optim_repst_dev_995;
    output_995(22) = optim_fec_dev_995;
    output_995(23) = optim_jsurv_dev_995;
    output_995(24) = optim_jobs_dev_995;
    output_995(25) = optim_jsize_dev_995;
    output_995(26) = optim_jsizeb_dev_995;
    output_995(27) = optim_jsizec_dev_995;
    output_995(28) = optim_jrepst_dev_995;
    output_995(29) = optim_jmatst_dev_995;
    output_995(30) = optim_year2_995;
    output_995(31) = optim_mpm_altered_995;
    output_995(32) = optim_vrm_altered_995;
    
    output_995.attr("names") = clone(ta_names);
    output_995.attr("class") = "data.frame";
    output_995.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, data_frame_length);
    
    optim_ta_995 = output_995;
    
    List output_ESS (33);
    
    output_ESS(0) = ESS_variant;
    output_ESS(1) = ESS_stage3;
    output_ESS(2) = ESS_stage2;
    output_ESS(3) = ESS_stage1;
    output_ESS(4) = ESS_age3;
    output_ESS(5) = ESS_age2;
    output_ESS(6) = ESS_eststage3;
    output_ESS(7) = ESS_eststage2;
    output_ESS(8) = ESS_eststage1;
    output_ESS(9) = ESS_estage3;
    output_ESS(10) = ESS_estage2;
    output_ESS(11) = ESS_givenrate;
    output_ESS(12) = ESS_offset;
    output_ESS(13) = ESS_multiplier;
    output_ESS(14) = ESS_convtype;
    output_ESS(15) = ESS_convtype_t12;
    output_ESS(16) = ESS_surv_dev;
    output_ESS(17) = ESS_obs_dev;
    output_ESS(18) = ESS_size_dev;
    output_ESS(19) = ESS_sizeb_dev;
    output_ESS(20) = ESS_sizec_dev;
    output_ESS(21) = ESS_repst_dev;
    output_ESS(22) = ESS_fec_dev;
    output_ESS(23) = ESS_jsurv_dev;
    output_ESS(24) = ESS_jobs_dev;
    output_ESS(25) = ESS_jsize_dev;
    output_ESS(26) = ESS_jsizeb_dev;
    output_ESS(27) = ESS_jsizec_dev;
    output_ESS(28) = ESS_jrepst_dev;
    output_ESS(29) = ESS_jmatst_dev;
    output_ESS(30) = ESS_year2;
    output_ESS(31) = ESS_mpm_altered;
    output_ESS(32) = ESS_vrm_altered;
    
    output_ESS.attr("names") = clone(ta_names);
    output_ESS.attr("class") = "data.frame";
    output_ESS.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, 2);
    
    ESS_ta = output_ESS;
    
    variable_traits = ESS_variable_traits;
  }

}

#endif
