#' Summarize adaptProj Objects
#' 
#' Function \code{summary.adaptProj()} summarizes \code{adaptProj} objects.
#' 
#' @name summary.adaptProj
#' 
#' @param object An \code{adaptProj} object.
#' @param threshold A threshold population size to be searched for in
#' projections. Defaults to 1.
#' @param inf_alive A logical value indicating whether to treat infinitely
#' large population size as indicating that the population is still extant.
#' If \code{FALSE}, then the population is considered extinct. Defaults to
#' \code{TRUE}.
#' @param milepost A numeric vector indicating at which points in the projection
#' to assess detailed results. Can be input as integer values, in which case
#' each number must be between 1 and the total number of occasions projected in
#' each projection, or decimals between 0 and 1, which would then be translated
#' into the corresponding projection steps of the total. Defaults to
#' \code{c(0, 0.25, 0.50, 0.75, 1.00)}.
#' @param ext_time A logical value indicating whether to output extinction times
#' per population-patch. Defaults to \code{FALSE}.
#' @param finalN_mean A logical value indicating whether to take the arithmetic
#' mean of the final population sizes for each MPM across all replicates.
#' Defaults to \code{FALSE}.
#' @param finalN_used An integer value indicating the number of final population
#' sizes in the arithmetic mean noted in argument \code{finalN_mean}. Defaults
#' to \code{100}, unless the projections are for fewer time steps, in which case
#' defaults to \code{10}.
#' @param ... Other parameters currently not utilized.
#' 
#' @return Apart from a statement of the results, this function outputs a list
#' with the following elements:
#' \item{milepost_sums}{A data frame showing the number of replicates at each
#' of the milepost times that is above the threshold population/patch size.}
#' \item{extinction_times}{A dataframe showing the numbers of replicates going
#' extinct (\code{ext_reps}) and mean extinction time (\code{ext_time}) per
#' population-patch. If \code{ext_time = FALSE}, then only outputs \code{NA}.}
#' \item{final_N}{The final population size of each MPM across all replicates.
#' Only given if \code{finalN_mean = TRUE}.}
#' 
#' @section Notes:
#' The \code{inf_alive} and \code{ext_time} options both assess whether
#' replicates have reached a value of \code{NaN} or \code{Inf}. If
#' \code{inf_alive = TRUE} or \code{ext_time = TRUE} and one of these values is
#' found, then the replicate is counted in the \code{milepost_sums} object if
#' the last numeric value in the replicate is above the \code{threshold} value,
#' and is counted as extant and not extinct if the last numeric value in the
#' replicate is above the extinction threshold of a single individual.
#' 
#' Extinction time is calculated on the basis of whether the replicate ever
#' falls below a single individual. A replicate with a positive population size
#' below 0.0 that manages to rise above 1.0 individual is still considered to
#' have gone extinct the first time it crossed below 1.0.
#' 
#' If the input \code{lefkoProj} object is a mixture of two or more other
#' \code{lefkoProj} objects, then mileposts will be given relative to the
#' maximum number of time steps noted.
#' 
#' @examples
#' library(lefko3)
#' data(cypdata)
#' 
#' data(cypa_data)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cycaraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#'   
#' cyparaw_v1 <- verticalize3(data = cypa_data, noyears = 18, firstyear = 1994,
#'   individcol = "plant_id", blocksize = 2, sizeacol = "Inf.94",
#'   sizebcol = "Veg.94", repstracol = "Inf.94", fecacol = "Inf.94",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' cyp_supp_list1 <- list(cypsupp2r, cypsupp2r)
#' 
#' cycamatrix2r <- rlefko2(data = cycaraw_v1, stageframe = cypframe_raw, 
#'   year = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' cypamatrix2r <- rlefko2(data = cyparaw_v1, stageframe = cypframe_raw, 
#'   year = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' cyp_mpm_list <- list(cycamatrix2r, cypamatrix2r)
#' 
#' cyca2_start <- start_input(cycamatrix2r, stage2 = c("SD", "P1", "P2"),
#'   value = c(500, 100, 200))
#' cypa2_start <- start_input(cypamatrix2r, stage2 = c("SD", "P1", "P2"),
#'   value = c(5000, 1000, 2000))
#' cyp_start_list <- list(cyca2_start, cypa2_start)
#' 
#' cyp2_dv <- density_input(cypamatrix2r, stage3 = c("SD", "P1"),
#'   stage2 = c("rep", "rep"), style = c(1, 1), alpha = c(0.5, 1.2),
#'   beta = c(1.0, 2.0), type = c(2, 1))
#' cyp_dv_list <- list(cyp2_dv, cyp2_dv)
#' 
#' cyp_comm_proj <- project3(mpms = cyp_mpm_list, starts = cyp_start_list,
#'   density = cyp_dv_list, times = 10)
#'   
#' summary(cyp_comm_proj)
#' 
#' 
#' @export
summary.adaptProj <- function(object, threshold = 1, inf_alive = TRUE,
  milepost = c(0, 0.25, 0.50, 0.75, 1.00), ext_time = FALSE, finalN_mean = FALSE,
  finalN_used = 100, ...) {
  
  num_reps <- num_times <- 0
  appended <- FALSE
  max_times <- max_reps <- 1L
  ave_times <- ave_reps <- 1.0
  step_text <- "steps"
  pop_text <- "populations"
  rep_text <- "replicates"
  
  num_reps <- length(object$N_out)
  num_pops <- dim(object$N_out[[1]])[1]
  num_times <- dim(object$N_out[[1]])[2]
  
  mean_throw_error <- FALSE
  
  mean_N_vec <- rep(NA, num_pops)
  
  if (any(milepost < 0)) {
    stop("Option milepost may not take negative values.", call. = FALSE)
  }
  if (any(milepost > num_times)) {
    stop("Option milepost may not take values higher than the number of actual
      number of projected occasions.", call. = FALSE)
  }
  
  if (inf_alive | ext_time) {
    for (i in c(1:num_reps)) {
      for (j in c(1:num_pops)) {
        for (k in c(1:num_times)) {
          if ((is.nan(object$N_out[[i]][j, k]) | is.infinite(object$N_out[[i]][j, k])) & k > 1) {
            object$N_out[[i]][j, k] <- object$N_out[[i]][j, (k - 1)]   #. max_found
          }
        }
      }
    }
  }
  
  if (ext_time) {
    the_numbers <- apply(as.matrix(c(1:num_reps)), 1, function(X) {
        freemasonry <- apply(as.matrix(c(1:num_pops)), 1, function(Y) {
            ext_points <- which(object$N_out[[X]][Y,] < 1)
            if (length(ext_points) > 0) return (min(ext_points)) else return (NA)
          }
        )
        ext_varmints <- length(which(!is.na(freemasonry) & !is.nan(freemasonry)))
        if (ext_varmints > 0) {
          ext_time <- mean(freemasonry, na.rm = TRUE)
        } else {
          ext_time <- NA
        }
        return (c(ext_varmints, ext_time))
      }
    )
    
    the_numbers <- t(the_numbers)
    the_numbers <- as.data.frame(the_numbers)
    colnames(the_numbers) <- c("ext_reps", "ext_time")
    
  } else {
    the_numbers <- NA
  }
  
  for (i in c(1:num_pops)) {
    if (any(milepost > num_times)) {
      stop("Entered milepost values are outside the allowable range.", call. = FALSE)
    }
  }
  
  if (finalN_mean) {
    if (finalN_used > num_times) {
      if (num_times > 10) {
        finalN_used <- 10
        warning("Mean N will be calculated over the final 10 time steps.", call. = FALSE)
      } else {
        warning("Unable to generate mean given the length of projection and the set value of finalN_used.", call. = FALSE)
        mean_throw_error <- TRUE
      }
    }
    
    if (!mean_throw_error) {
      useable_N_mat <- Reduce("+", object$N_out)
      useable_N_mat <- useable_N_mat / num_reps
      
      mean_N_vec <- apply(useable_N_mat[, c((num_times - finalN_used + 1):num_times)], 1, mean)
      names(mean_N_vec) <- c(1:num_pops)
    }
  }
  
  if (num_pops == 1) pop_text <- "population"
  if (num_reps == 1) rep_text <- "replicate"
  if (num_times == 1) step_text <- "step"
  
  writeLines(paste0("\nThe input adaptProj object covers ", num_pops, " ",
      pop_text, ", ", num_times, " projected ", step_text, " and ", num_reps,
      " projected ", rep_text, "."), con = stdout())
  
  writeLines(paste0("The number of replicates with population size above the threshold size of ",
    threshold, " is as in"), con = stdout())
  writeLines(paste0("the following matrix, with populations given by row and milepost times given by column: \n"),
    con = stdout())
  
  used_milepost <- milepost
  
  if (all(milepost >= 0) & all(milepost <= 1)) {
    used_milepost <- floor(used_milepost * num_times)
  } else if (any(milepost == 0)) {
    used_milepost <- used_milepost
  }
  if (any(used_milepost == 0)) {
    used_milepost[which(used_milepost == 0)] <- 1
  }
  used_milepost <- unique(used_milepost)
  used_milepost <- sort(used_milepost)
  
  milepost_sums <- matrix(0, nrow = num_pops, ncol = length(used_milepost))
  
  for (i in c(1:num_reps)) {
    for (j in c(1:num_pops)) {
      for (k in c(1:length(used_milepost))) {
        if (object$N_out[[i]][j,used_milepost[k]] >= threshold) milepost_sums[j, k] <- milepost_sums[j, k] + 1
      }
    }
  }
  colnames(milepost_sums) <- used_milepost
  
  if (finalN_mean) {
    output <- list(milepost_sums = milepost_sums, extinction_times = the_numbers,
      final_N = mean_N_vec)
  } else {
    output <- list(milepost_sums = milepost_sums, extinction_times = the_numbers)
  }
  
  return (output)
}

#' Summarize adaptProjBatch Objects
#' 
#' Function \code{summary.adaptProjBatch()} summarizes \code{adaptProjBatch}
#' objects.
#' 
#' @name summary.adaptProjBatch
#' 
#' @param object An \code{adaptProjBatch} object.
#' @param finalN_mean A logical value indicating whether to take the arithmetic
#' mean of the final population sizes for each MPM in each projection. Defaults
#' to \code{FALSE}, in which case only the final population sizes are reported.
#' @param finalN_used An integer value indicating the number of final population
#' sizes in the arithmetic mean noted in argument \code{finalN_mean}. Defaults
#' to \code{100}, unless the projections are for fewer time steps, in which case
#' defaults to \code{10}.
#' @param threshold A threshold population size to be searched for in
#' projections. Defaults to 1.
#' @param inf_alive A logical value indicating whether to treat infinitely
#' large population size as indicating that the population is still extant.
#' If \code{FALSE}, then the population is considered extinct. Defaults to
#' \code{TRUE}.
#' @param ext_time A logical value indicating whether to output extinction times
#' per population-patch. Defaults to \code{FALSE}.
#' @param print_output A logical value indicating whether to print the output
#' data frame to the screen. Defaults to \code{FALSE}.
#' @param ... Other parameters currently not utilized.
#' 
#' @return Apart from a statement of the results, this function outputs a data
#' frame with the following elements:
#' \item{projection}{The identity of the current projection in the original
#' \code{adaptProjBatch} object.}
#' \item{target_mpm}{The identity of the MPM targeted for alteration in the
#' batch projection.}
#' \item{stage3}{Stage at occasion \emph{t}+1 in the transition replaced.}
#' \item{stage2}{Stage at occasion \emph{t} in the transition replaced.}
#' \item{stage1}{Stage at occasion \emph{t}-1 in the transition replaced.}
#' \item{age2}{Age at occasion \emph{t} in the transition replaced.}
#' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and
#' \code{stage1}.}
#' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and
#' \code{stage1}.}
#' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and
#' \code{stage1}.}
#' \item{estage2}{Age at occasion \emph{t} in the transition to replace
#' the transition designated by \code{age2}.}
#' \item{givenrate}{A constant to be used as the value of the transition.}
#' \item{offset}{A constant value to be added to the transition or proxy
#' transition.}
#' \item{multiplier}{A multiplier for proxy transitions or for fecundity.}
#' \item{convtype}{Designates whether the transition from occasion \emph{t} to
#' occasion \emph{t}+1 is a survival transition probability (1), a fecundity
#' rate (2), or a fecundity multiplier (3).}
#' \item{convtype_t12}{Designates whether the transition from occasion
#' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), or
#' a fecundity rate (2).}
#' \item{rep}{The identity of the replicate being summarized, within the
#' current projection.}
#' \item{mpm}{The identity of the MPM for which the population summary
#' corresponding to the row in question is being given.}
#' \item{final_N}{The final population size, meaning the population size given
#' for the current MPM in the current replicate in the current projection, in
#' the final time recorded.}
#' \item{extinct_by}{The first time by which the population size goes below the
#' extinction threshold, or hits 0.}
#' \item{final_N_mean}{The mean population size during the final
#' \code{finalN_used} times for the current MPM in the current replicate in the
#' current projection.}
#' 
#' \item{extinction_times}{A dataframe showing the numbers of replicates going
#' extinct (\code{ext_reps}) and mean extinction time (\code{ext_time}) per
#' population-patch. If \code{ext_time = FALSE}, then only outputs \code{NA}.}
#' 
#' @section Notes:
#' The \code{inf_alive} and \code{ext_time} options both assess whether
#' replicates have reached a value of \code{NaN} or \code{Inf}. If
#' \code{inf_alive = TRUE} or \code{ext_time = TRUE} and one of these values is
#' found, then the replicate is counted in the \code{milepost_sums} object if
#' the last numeric value in the replicate is above the \code{threshold} value,
#' and is counted as extant and not extinct if the last numeric value in the
#' replicate is above the extinction threshold of a single individual.
#' 
#' Extinction time is calculated on the basis of whether the replicate ever
#' falls below a single individual. A replicate with a positive population size
#' below 0.0 that manages to rise above 1.0 individual is still considered to
#' have gone extinct the first time it crossed below 1.0.
#' 
#' If the input \code{lefkoProj} object is a mixture of two or more other
#' \code{lefkoProj} objects, then mileposts will be given relative to the
#' maximum number of time steps noted.
#' 
#' @examples
#' library(lefko3)
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' sizevector <- c(0, 0, 3.0, 15)
#' stagevector <- c("P1", "D", "Sm", "Lg")
#' repvector <- c(0, 0, 1, 1)
#' obsvector <- c(0, 0, 1, 1)
#' matvector <- c(0, 1, 1, 1)
#' immvector <- c(1, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1)
#' binvec <- c(0, 0.5, 2.5, 9.5)
#' 
#' cypframe_small_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypraw_v2 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_small_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypraw_v3 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   NAas0 = TRUE, NRasRep = TRUE)
#' 
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 1500, 500),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' cypmean <- lmean(cypmatrix2r)
#' 
#' cypsupp2r_small <- supplemental(stage3 = c("D", "Sm", "Lg", "P1"),
#'   stage2 = c("P1", "P1", "P1", "rep"), eststage3 = c(NA, "Sm", "Lg", NA),
#'   eststage2 = c(NA, "D", "D", NA), givenrate = c(0.05, NA, NA, NA),
#'   offset = c(NA, NA, -0.1, NA), multiplier = c(NA, NA, NA, 0.5),
#'   type =c(1, 1, 1, 3), stageframe = cypframe_small_raw, historical = FALSE)
#' cypmatrix2r_small <- rlefko2(data = cypraw_v2, stageframe = cypframe_small_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r_small,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' cypmean_small <- lmean(cypmatrix2r_small)
#' 
#' cypmatrixL_small <- rleslie(data = cypraw_v3, start_age = 1, last_age = 4,
#'   continue = TRUE, fecage_min = 3, year = "all", pop = NA, patch = "all",
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cyp_mpms1 <- list(cypmatrix2r, cypmatrix2r_small, cypmatrixL_small)
#' 
#' c2d_4 <- density_input(cypmean, stage3 = c("P1", "P1"), stage2= c("SD", "rep"),
#'   style = 1, time_delay = 1, alpha = 1, beta = 0.0005, type = c(2, 2))
#' c2d_4a <- density_input(cypmean_small, stage3 = c("P1", "P1"), stage2= c("P1", "rep"),
#'   style = 1, time_delay = 1, alpha = 1, beta = 0.0005, type = c(2, 2))
#' cypL_dv <- density_input(cypmatrixL_small, stage3 = c("Age1"), stage2 = c("rep"),
#'   style = c(1), alpha = c(0.5), beta = c(1.0), type = c(2))
#' cyp_density <- list(c2d_4, c2d_4a, cypL_dv)
#' 
#' cyp_start1 <- start_input(cypmatrix2r, stage2 = c("SD", "P1", "D"),
#'   value = c(100, 200, 4))
#' cyp_start2 <- start_input(cypmatrix2r_small, stage2 = c("P1", "D"),
#'   value = c(10, 2000))
#' cypL_start_1 <- start_input(cypmatrixL_small, stage2 = c("Age1"),
#'   value = c(200))
#' cyp_start <- list(cyp_start1, cyp_start2, cypL_start_1)
#' 
#' new_supplement_cyp2_small <- sup_skeleton(2)
#' new_supplement_cyp2_small$stage3 <- c("D", "Sm")
#' new_supplement_cyp2_small$stage2 <- c("Lg", "Lg")
#' new_supplement_cyp2_small$convtype <- c(1, 1)
#' used_supplements <- list(new_supplement_cyp2_small,
#'   new_supplement_cyp2_small, NULL)
#' 
#' aaa1_prj_batch2 <- batch_project3(used_mpms = "all", all_elems = FALSE,
#'   mpms =  cyp_mpms1, entry_time = c(0, 5, 8), times = 15, nreps = 3,
#'   supplement = used_supplements, integeronly = TRUE, density = cyp_density)
#'   
#' summary(aaa1_prj_batch2)
#' 
#' @export
summary.adaptProjBatch <- function(object, finalN_mean = FALSE,
  finalN_used = 100, threshold = 1, inf_alive = TRUE, ext_time = FALSE,
  print_output = FALSE, ...) {
  
  if (!("proj_out" %in% names(object))) {
    stop("Argument object does not appear to be an adaptProjBatch object.")
  }
  if (!("ref" %in% names(object))) {
    stop("Argument object does not appear to be an adaptProjBatch object.")
  }
  if (!("control" %in% names(object))) {
    stop("Argument object does not appear to be an adaptProjBatch object.")
  }
  
  total_time_steps <- dim(object$proj_out[[1]][[1]]$N_out[[1]])[2] - 1
  final_time <- total_time_steps + 1
  
  start_time <- final_time - finalN_used
  if (start_time < 1 & finalN_mean) {
    if (final_time > 10) {
      start_time <- final_time - 10
      warning("Mean N will be calculated over the final 10 time steps.", call. = FALSE)
    } else {
      finalN_mean <- FALSE
      warning("Unable to generate mean given the length of projection and the set value of finalN_used.", call. = FALSE)
    }
  }
  
  total_lengths <- apply(as.matrix(seq_along(object$ref)), 1, function(X) {
    return (length(object$ref[[X]]))
  })
  
  used_mpm_vec <- object$control
  used_mpm_vec <- used_mpm_vec[which(total_lengths > 0)]
  all_mpm_vec <- object$proj_out[[1]][[1]]$labels$mpm
  
  targeted_mpms_num <- length(used_mpm_vec)
  all_mpms_num <- length(all_mpm_vec)
  all_reps <- length(object$proj_out[[1]][[1]]$N_out)
  all_rep_vec <- seq(from = 1, to = all_reps)
  
  writeLines(paste0("This adaptProjBatch has ", sum(total_lengths),
    " projections covering ", all_mpms_num, " mpms, of which ", targeted_mpms_num,
    " have been targeted for alterations."))
  writeLines(paste0("Assuming the first projection is representative, each projection is ",
    total_time_steps, " time steps in length.\n\n"))
  
  new_list <- lapply(used_mpm_vec, function(i) {
    list_of_indices_df <- if (i == 1) {
      seq_along(object$ref[[i]])
    } else {
      current_start <- 0
      for (j in c(1:(i-1))) {
        current_start = current_start + length(object$ref[[j]])
      }
      seq((current_start + 1), (current_start + length(object$ref[[i]])))
    }
    
    vector_of_core_indices <- seq_along(object$ref[[i]])
    
    new_list_sub <- lapply(vector_of_core_indices, function(X) {
      new_df <- cbind.data.frame(list_of_indices_df[X], object$ref[[i]][[X]])
      colnames(new_df)[1] <- "projection"
      colnames(new_df)[2] <- "target_mpm"
      
      new_df_list <- vector(mode = "list", length = all_mpms_num)
      for (m in c(1:all_mpms_num)) new_df_list[[m]] <- new_df
      
      new_df_list_sub <- lapply(all_rep_vec, function(Y) {
        new_df_list_sub_sub <- lapply(all_mpm_vec, function(Z){
          new_df_list_df <- cbind.data.frame(new_df_list[[Z]], Y)
          new_df_list_df <- cbind.data.frame(new_df_list_df, Z)
          
          new_finalN <- object$proj_out[[i]][[X]]$N_out[[Y]][Z, final_time]
          new_df_list_df <- cbind.data.frame(new_df_list_df, new_finalN)
          
          new_all_N <- object$proj_out[[i]][[X]]$N_out[[Y]][Z, ]
          first_entry <- min(which(new_all_N > 0))
          new_all_N_subset <- new_all_N[first_entry:length(new_all_N)]
          extinct_times <- which(new_all_N_subset == 0 | new_all_N_subset < threshold)
          
          extinct_by <- NA
          if (length(extinct_times > 0)) extinct_by <- min(extinct_times) + (first_entry - 1)
          
          new_df_list_df <- cbind.data.frame(new_df_list_df, extinct_by)
          
          colnames(new_df_list_df)[16] <- "rep"
          colnames(new_df_list_df)[17] <- "mpm"
          colnames(new_df_list_df)[18] <- "final_N"
          colnames(new_df_list_df)[19] <- "extinct_by"
          
          if (finalN_mean) {
            finalN_vec <- object$proj_out[[i]][[X]]$N_out[[Y]][Z, c(start_time:final_time)]
            finalN_mean_value <- mean(finalN_vec, na.rm = TRUE)
            
            new_df_list_df <- cbind.data.frame(new_df_list_df, finalN_mean_value)
            colnames(new_df_list_df)[20] <- "final_N_mean"
          }
          return (new_df_list_df)
        })
        
        return (do.call(rbind.data.frame, new_df_list_sub_sub))
      })
      
      return (do.call(rbind.data.frame, new_df_list_sub))
    })
    return (do.call(rbind.data.frame, new_list_sub))
  })
  
  out_df <- do.call(rbind.data.frame, new_list)
  
  if (print_output) print(out_df, digits = 3)
  
  return (out_df)
}

#' Summarize adaptInv Objects
#' 
#' Function \code{summary.adaptInv()} summarizes \code{adaptInv} objects.
#' 
#' @name summary.adaptInv
#' 
#' @param object An \code{adaptInv} object.
#' @param ... Other parameters currently not utilized.
#' 
#' @return This function only produces text summarizing the numbers of variants,
#' time steps, replicates, ESS optima, etc.
#' 
#' @examples
#' library(lefko3)
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.40, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, 1000, 1000),
#'   type =c(1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' cypmean <- lmean(cypmatrix2r)
#' 
#' cyp_start <- start_input(cypmean, stage2 = c("SD", "P1", "D"),
#'   value = c(1000, 200, 4))
#' 
#' c2d_4 <- density_input(cypmean, stage3 = c("P1", "P1"), stage2= c("SD", "rep"),
#'   style = 2, time_delay = 1, alpha = 0.005, beta = 0.000005, type = c(2, 2))
#' 
#' # A simple projection allows us to find a combination of density dependence
#' # and running time that produces a stable quasi-equilibrium
#' cyp_proj <- projection3(cypmean, times = 250, start_frame = cyp_start,
#'   density = c2d_4, integeronly = TRUE)
#' plot(cyp_proj)
#' 
#' cyp_ta <- trait_axis(stageframe = cypframe_raw,
#'   stage3 = rep("P1", 15),
#'   stage2 = rep("rep", 15),
#'   multiplier = seq(from = 0.1, to = 10.0, length.out = 15),
#'   type = rep(2, 15))
#' 
#' cyp_inv <- invade3(axis = cyp_ta, mpm = cypmean, density = c2d_4, times = 350,
#'   starts = cyp_start, entry_time = c(0, 250), fitness_times = 30,
#'   var_per_run = 2)
#' summary(cyp_inv)
#' 
#' @export
summary.adaptInv <- function(object, ...) {
  
  num_reps <- num_times <- 0
  appended <- FALSE
  max_times <- max_reps <- 1L
  ave_times <- ave_reps <- 1.0
  step_text <- "steps"
  run_variant_text <- "variants per run"
  rep_text <- "replicates"
  time_text <- "steps"
  run_text <- "runs"
  variant_text <- "variants"
  
  num_reps <- length(object$N_out)
  num_run_variants <- dim(object$N_out[[1]])[1]
  num_times <- dim(object$N_out[[1]])[2]
  num_runs <- dim(object$N_out[[1]])[3]
  
  all_fitness_vars <- names(object$fitness)
  fitness_variant1 <- object$fitness$variant1
  total_variants <- length(unique(fitness_variant1))
  found_entrytime_vars <- grep("entry", all_fitness_vars)
  found_fitness_vars <- grep("fitness", all_fitness_vars)
  
  if (total_variants == 1) variant_text <- "variant"
  if (num_run_variants == 1) run_variant_text <- "variant per run"
  if (num_runs == 1) pop_text <- "run"
  if (num_reps == 1) rep_text <- "replicate"
  if (num_times == 1) time_text <- "step"
  
  writeLines(paste0("\nThe input adaptInv object covers ", total_variants, " ",
      variant_text, ", ", num_times, " projected ", time_text, ", ",
      num_run_variants, " ", run_variant_text, ", and ", max_reps,
      " projected ", rep_text, "."), con = stdout())
  
  if (is.element("optim", names(object))) {
    
    if (length(object$optim) > 0) {
      if (is.element("ESS_values", names(object$optim))) {
        found_ESS_frame <- object$optim$ESS_values
        
        if (length(found_ESS_frame) > 1) {
          found_ESS_values <- nrow(found_ESS_frame)
          
          optimum_text <- "optima"
          if (found_ESS_values == 1) optimum_text <- "optimum"
          writeLines(paste0("It includes optimization data suggesting ", found_ESS_values,
            " purported ESS ", optimum_text, "."), con = stdout())
        } else {
          writeLines("It includes optimization data but suggests no purported ESS optima.",
            con = stdout())
        }
      }
    } else {
      writeLines("No trait optima found.")
    }
  }
  writeLines("", con = stdout())
}
