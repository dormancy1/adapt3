#ifndef ADAPTUTILS_input_utils_H
#define ADAPTUTILS_input_utils_H

#include <RcppArmadillo.h>
#define BOOST_DISABLE_ASSERTS

#include <LefkoUtils.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using namespace arma;





// Index of functions
// 1. void character_vectorizer  Create Standardized CharacterVectors Based on Non-Standard Input
// 2. void integer_vectorizer  Create Standardized IntegerVectors Based on Non-Standard Input
// 3. void numeric_vectorizer  Create Standardized NumericVectors Based on Non-Standard Input
// 4. int list_vector_length  Find Length of Longest Vector in a List


namespace AdaptInputs {
  
  //' Create Standardized CharacterVectors Based on Non-Standard Input
  //' 
  //' @name character_vectorizer
  //' 
  //' @param output The output reference, passed by reference.
  //' @param input The input vector, treated as an \code{RObject}.
  //' @param argument_name The name of the argument used as \code{input}, given as
  //' a String.
  //' @param overall_length An integer giving the desired length of the vector.
  //' @param change_value The string to set \code{NA} values to if
  //' \code{NAasOther = TRUE}.
  //' @param NAasOther A Boolean value indicating whether to treat \code{NA}
  //' values as the value specified in \code{change_value}.
  //' 
  //' @return This function modifies an input vector by reference, given as
  //' argument \code{output}. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void character_vectorizer (CharacterVector& output,
    Nullable<RObject> input, String argument_name, int overall_length,
    String change_value, bool NAasOther = false) {
    
    CharacterVector input_;
    int input_length {0};
    
    if (input.isNotNull()) {
      if (is<CharacterVector>(input)) {
        input_ = as<CharacterVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (CharacterVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (input_length != overall_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as the longest input vector.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else {
        if (overall_length != 0) {
          if (NAasOther) {
            CharacterVector input_temp (overall_length, change_value);
            input_ = input_temp;
          } else {
            CharacterVector input_temp (overall_length, NA_STRING);
            input_ = input_temp;
          }
        }
      }
    } else {
      if (overall_length != 0) {
        if (NAasOther) {
          CharacterVector input_temp (overall_length, change_value);
          input_ = input_temp;
        } else {
          CharacterVector input_temp (overall_length, NA_STRING);
          input_ = input_temp;
        }
      }
    }
    output = input_;
  }
  
  //' Create Standardized IntegerVectors Based on Non-Standard Input
  //' 
  //' @name integer_vectorizer
  //' 
  //' @param output The output reference, passed by reference.
  //' @param input The input vector, treated as an \code{RObject}.
  //' @param argument_name The name of the argument used as \code{input}, given as
  //' a String.
  //' @param overall_length An integer giving the desired length of the vector.
  //' @param min_limit The smallest integer to allow, if using limits.
  //' @param max_limit The largest integer to allow, if using limits.
  //' @param use_limits A Boolean variable indicating whether to limit allowable
  //' values.
  //' @param NAasOther A Boolean value indicating whether to treat \code{NA}
  //' values as the value specified in \code{change_value}.
  //' @param change_value The integer to set \code{NA} values to if
  //' \code{NAasOther = TRUE}.
  //' 
  //' @return This function modifies an input vector by reference, given as
  //' argument \code{output}. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void integer_vectorizer (IntegerVector& output, Nullable<RObject> input,
    String argument_name, int overall_length, int min_limit, int max_limit,
    bool use_limits = false, bool NAasOther = false, int change_value = 0) {
    
    IntegerVector input_;
    int input_length {0};
    
    if (input.isNotNull()) {
      if (is<IntegerVector>(input)) {
        input_ = as<IntegerVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (IntegerVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (use_limits) {
          for (int i = 0; i < input_length; i++) { 
            if (!IntegerVector::is_na(input_(i))) {
              if (input_(i) < min_limit || input_(i) > max_limit) {
                String eat_my_shorts = "Argument ";
                eat_my_shorts += argument_name;
                eat_my_shorts += " must be an integer between ";
                eat_my_shorts += min_limit;
                eat_my_shorts += " and ";
                eat_my_shorts += max_limit;
                eat_my_shorts += ".";
                
                throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
              }
            }
          }
        }
        
        if (input_length != overall_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as the longest input vector.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<NumericVector>(input)) {
        input_ = as<IntegerVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (NumericVector::is_na(input_(i)) || IntegerVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (use_limits) {
          for (int i = 0; i < input_length; i++) { 
            if (!IntegerVector::is_na(input_(i)) && !NumericVector::is_na(input_(i))) {
              if (input_(i) < min_limit || input_(i) > max_limit) {
                String eat_my_shorts = "Argument ";
                eat_my_shorts += argument_name;
                eat_my_shorts += " must be an integer between ";
                eat_my_shorts += min_limit;
                eat_my_shorts += " and ";
                eat_my_shorts += max_limit;
                eat_my_shorts += ".";
                
                throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
              }
            }
          }
        }
        
        if (input_length != overall_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as the longest input vector.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<LogicalVector>(input)) {
          if (overall_length != 0) {
            if (NAasOther) {
              IntegerVector input_temp (overall_length, change_value);
              input_ = input_temp;
            } else {
              IntegerVector input_temp (overall_length, NA_INTEGER);
              input_ = input_temp;
            }
          }
      } else {
        String eat_my_shorts = "Please enter argument ";
        eat_my_shorts += argument_name;
        eat_my_shorts += " in integer format.";
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
    } else {
      if (overall_length != 0) {
        if (NAasOther) {
          IntegerVector input_temp (overall_length, change_value);
          input_ = input_temp;
        } else {
          IntegerVector input_temp (overall_length, NA_INTEGER);
          input_ = input_temp;
        }
      }
    }
    output = input_;
  }
  
  //' Create Standardized NumericVectors Based on Non-Standard Input
  //' 
  //' @name numeric_vectorizer
  //' 
  //' @param output The output reference, passed by reference.
  //' @param input The input vector, treated as an \code{RObject}.
  //' @param argument_name The name of the argument used as \code{input}, given as
  //' a String.
  //' @param overall_length An integer giving the desired length of the vector.
  //' @param NAasOther A Boolean value indicating whether to treat \code{NA}
  //' values as value given in \code{change_value}.
  //' @param change_value The numeric value to change \code{NA}s to, if
  //' \code{NAasOther = TRUE}.
  //' 
  //' @return This function modifies an input vector by reference, given as
  //' argument \code{output}. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void numeric_vectorizer (NumericVector& output, Nullable<RObject> input,
    String argument_name, int overall_length, bool NAasOther = false,
    double change_value = 0) {
    
    NumericVector input_;
    int input_length {0};
    
    if (input.isNotNull()) {
      if (is<NumericVector>(input)) {
        input_ = as<NumericVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (NumericVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (input_length != overall_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as the longest input vector.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<LogicalVector>(input)) {
        if (overall_length != 0) {
          if (NAasOther) {
            NumericVector input_temp (overall_length, change_value);
            input_ = input_temp;
          } else {
            NumericVector input_temp (overall_length, NA_REAL);
            input_ = input_temp;
          }
        }
      } else {
        String eat_my_shorts = "Please enter argument ";
        eat_my_shorts += argument_name;
        eat_my_shorts += " in numeric format.";
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
    } else {
      if (overall_length != 0) {
        if (NAasOther) {
          NumericVector input_temp (overall_length, change_value);
          input_ = input_temp;
        } else {
          NumericVector input_temp (overall_length, NA_REAL);
          input_ = input_temp;
        }
      }
    }
    output = input_;
  }
  
  //' Find Length of Longest Vector in a List
  //' 
  //' Function \code{list_vector_length()} finds the length of the longest
  //' vector in a list of vectors.
  //' 
  //' @name list_vector_length
  //' 
  //' @param all_entered_vectors A List of R style vectors, of any type.
  //' 
  //' @return An integer corresponding to the length of the longest vector in
  //' the entered list.
  //' 
  //' @keywords internal
  //' @noRd
  inline int list_vector_length (Rcpp::List& all_entered_vectors) {
    int list_vec_vars = static_cast<int>(all_entered_vectors.length()); 
    
    int found_max_length {0};
    for (int i = 0; i < list_vec_vars; i++) {
      if (is<NumericVector>(all_entered_vectors(i))) {
        NumericVector temp_vec = as<NumericVector>(all_entered_vectors(i));
        int temp_vec_length = static_cast<int>(temp_vec.length());
        
        if (temp_vec_length > found_max_length) found_max_length = temp_vec_length;
      } else if (is<IntegerVector>(all_entered_vectors(i))) {
        IntegerVector temp_vec = as<IntegerVector>(all_entered_vectors(i));
        int temp_vec_length = static_cast<int>(temp_vec.length());
        
        if (temp_vec_length > found_max_length) found_max_length = temp_vec_length;
      } else if (is<CharacterVector>(all_entered_vectors(i))) {
        CharacterVector temp_vec = as<CharacterVector>(all_entered_vectors(i));
        int temp_vec_length = static_cast<int>(temp_vec.length());
        
        if (temp_vec_length > found_max_length) found_max_length = temp_vec_length;
      } else if (is<LogicalVector>(all_entered_vectors(i))) {
        LogicalVector temp_vec = as<LogicalVector>(all_entered_vectors(i));
        int temp_vec_length = static_cast<int>(temp_vec.length());
        
        if (temp_vec_length > found_max_length) found_max_length = temp_vec_length;
      }
    }
    
    return found_max_length;
  }

}

#endif
