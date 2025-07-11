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
#' @param ... Other parameters currently not utilized.
#' 
#' @return Apart from a statement of the results, this function outputs a list
#' with the following elements:
#' \item{milepost_sums}{A data frame showing the number of replicates at each
#' of the milepost times that is above the threshold population/patch size.}
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
#'   individcol = "plant_id", blocksize = 3, sizeacol = "Inf2.94",
#'   sizebcol = "Inf.94", sizeccol = "Veg.94", repstracol = "Inf.94",
#'   repstrbcol = "Inf2.94", fecacol = "Inf.94", stageassign = cypframe_raw,
#'   stagesize = "sizeadded", NAas0 = TRUE, NRasRep = TRUE)
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
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypamatrix2r <- rlefko2(data = cyparaw_v1, stageframe = cypframe_raw, 
#'   year = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
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
  milepost = c(0, 0.25, 0.50, 0.75, 1.00), ext_time = FALSE, ...) {
  
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
  
  #used_milepost <- milepost * num_times
  
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
    
    #if (dim(object$labels)[1] > 1) {
    #  row_labels <- apply(object$labels, 1, function(X) {
    #    paste(X[1], X[2])
    #  })
    #  rownames(the_numbers) <- row_labels
    #}
  } else {
    the_numbers <- NA
  }
  
  for (i in c(1:num_pops)) {
    if (any(milepost > num_times)) {
      stop("Entered milepost values are outside the allowable range.", call. = FALSE)
    }
  }
  
  if (num_pops == 1) pop_text <- "population"
  if (num_reps == 1) rep_text <- "replicate"
  if (num_times == 1) step_text <- "step"
  
  writeLines(paste0("\nThe input adaptProj object covers ", num_pops, " ",
      pop_text, ", ", num_times, " projected ", step_text, " and ", max_reps,
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
  
  output <- list(milepost_sums = milepost_sums, extinction_times = the_numbers)
  
  return (output)
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
  found_entrytime_vars <- grep("entry", all_fitness_vars)
  found_fitness_vars <- grep("fitness", all_fitness_vars)
  
  total_variants <- length(found_entrytime_vars)
  
  if (total_variants == 1) variant_text <- "variant"
  if (num_run_variants == 1) run_variant_text <- "variant per run"
  if (num_runs == 1) pop_text <- "run"
  if (num_reps == 1) rep_text <- "replicate"
  if (num_times == 1) time_text <- "step"
  
  writeLines(paste0("\nThe input adaptInv object covers ", total_variants, " ",
      variant_text, ", ", num_times, " projected ", time_text, ", ",
      num_run_variants, " ", run_variant_text, ", and ", max_reps,
      " projected ", rep_text, "."), con = stdout())
  
  if (exists("optim", object)) {
    if (exists("ESS_values", object$optim)) {
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
  }
  writeLines("", con = stdout())
}
