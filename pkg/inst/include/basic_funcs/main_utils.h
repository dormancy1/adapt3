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
// 11. void pop_error2


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
    
    //Rcout << "ta_reassess a          ";
    
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
    
    //Rcout << "ta_reassess b          ";
    
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
    
    //Rcout << "ta_reassess c          ";
    
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
    
    //Rcout << "stagevec: " << stagevec << endl;
    //Rcout << "ta_reassess d          ";
    
    // Identify all groups
    arma::ivec all_groups = unique(groupvec);
    int no_groups {static_cast<int>(all_groups.n_elem)};
    StringVector group_text(no_groups);
    
    for (int i = 0; i < no_groups; i++) {
      group_text(i) = "group";
      group_text(i) += std::to_string(all_groups(i));
    }
    
    //Rcout << "ta_reassess e          ";
    
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
    
    //Rcout << "ta_reassess f          ";
    
    if (trait_axis.containsElementNamed("year2")) {
      year2_ta = as<StringVector>(trait_axis["year2"]);
    } else {
      StringVector year2_ta_ (stage2_ta.length(), NA_STRING);
      year2_ta = year2_ta_;
    }
    
    arma::uvec vrm_altered (ta_rows, fill::zeros); // vrm alteration (binary)  by trait axis row
    arma::uvec mpm_altered (ta_rows, fill::zeros); // mpm alteration (binary)  by trait axis row
    
    //Rcout << "ta_reassess g          ";
    
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
    
    //Rcout << "ta_reassess h          ";
    
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
      
      //Rcout << "ta_reassess i          ";
      
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
      
      //Rcout << "ta_reassess j          ";
      
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
    
    //Rcout << "ta_reassess k          ";
    
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
    
    //Rcout << "ta_reassess l          ";
    
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
    
    //Rcout << "ta_reassess m          ";
    
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
      //Rcout << "ta_reassess n  i: " << i << "          ";
      
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
      
      //Rcout << "ta_reassess o          ";
      
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
    
    //Rcout << "ta_reassess p          ";
    
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
    
    //Rcout << "ta_reassess q          ";
    
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
    
    //Rcout << "ta_reassess r          ";
    
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
    
    //Rcout << "ta_reassess s          ";
    
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
    
    //Rcout << "ta_reassess t          ";
    
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
}

#endif
