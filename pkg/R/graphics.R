#' Create Plot of Community Projection
#' 
#' Function \code{plot.adaptProj} plots community projections.
#' 
#' @name plot.adaptProj
#' 
#' @param x An \code{adaptProj} object.
#' @param repl The replicate to plot. Defaults to \code{1}, in which case the
#' first replicate is plotted.
#' @param auto_ylim A logical value indicating whether the maximum of the y axis
#' should be determined automatically. Defaults to \code{TRUE}, but reverts to
#' \code{FALSE} if any setting for \code{ylim} is given.
#' @param auto_col A logical value indicating whether to shift the color of
#' lines associated with each patch automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{col} is given.
#' @param auto_lty A logical value indicating whether to shift the line type
#' associated with each replicate automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{lty} is given.
#' @param auto_title A logical value indicating whether to add a title to each
#' plot. The plot is composed of the concatenated population and patch names.
#' Defaults to \code{FALSE}.
#' @param ... Other parameters used by functions \code{plot.default()} and
#' \code{lines()}.
#' 
#' @return A plot of the results of a \code{\link{project3}()} run.
#' 
#' @section Notes:
#' Output plots are currently limited to time series of population size.
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
#' cyparaw_v1 <- verticalize3(data = cypa_data, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
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
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypamatrix2r <- rlefko2(data = cyparaw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r_alt,
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
#' plot(cyp_comm_proj, lwd = 2, bty = "n")
#' 
#' @export
plot.adaptProj <- function(x, repl = 1, auto_ylim = TRUE, auto_col = TRUE,
  auto_lty = TRUE, auto_title = FALSE, ...) {
  
  used_N_mat <- NULL
  appended <- FALSE
  num_pops <- 0
  
  further_args <- list(...)
  
  if (length(further_args) == 0) further_args <- list()
  
  if (!is.element("type", names(further_args))) {
    further_args$type <- "l"
  }
  if (is.element("col", names(further_args))) {
    auto_col <- FALSE
  }
  if (is.element("lty", names(further_args))) {
    auto_lty <- FALSE
  }
  if (is.element("ylim", names(further_args))) {
    auto_ylim <- FALSE
  }
  if (is.element("main", names(further_args))) {
    auto_title <- FALSE
  }
  basal_args <- further_args
  
  if (!is.element("ylab", names(further_args))) {
    further_args$ylab <- "Population size"
  }
  if (!is.element("xlab", names(further_args))) {
    further_args$xlab <- "Time"
  }

  if (is.character(repl)) {
    stop("Argument repl not understood.", call. = FALSE)
  }
  
  used_col <- 1
  
  if (auto_title) {
    used_string <- paste("pop", x$labels[i, 1], "patch", x$labels[i, 2])
    further_args$main <- used_string
  }
  
  used_N_mat <- x$N_out[[repl]]
  
  used_col <- 1
  used_lty <- 1
  
  if (auto_ylim) {
    further_args$ylim <- c(0, max(used_N_mat, na.rm = TRUE))
  }
  
  if (auto_col) {
    further_args$col <- used_col
    basal_args$col <- used_col
  }
  if (auto_lty) {
    further_args$lty <- used_lty
  }
  
  num_pops <- dim(used_N_mat)[1]
  c_xy <- xy.coords(x = c(1:length(used_N_mat[1,])), y = used_N_mat[1,])
  further_args$x <- c_xy
  
  do.call("plot.default", further_args)
  
  if (num_pops > 1) {
    for (j in c(2:num_pops)) {
      if (auto_lty) {
        used_lty <- used_lty + 1;
        basal_args$lty <- used_lty
      }
      if (auto_col) {
        used_col <- used_col + 1;
        if (used_col > length(palette())) used_col <- 1;
        basal_args$col <- used_col
      }
      used_col <- used_col + 1;
      
      basal_args$x <- c(1:length(used_N_mat[j,]))
      basal_args$y <- used_N_mat[j,]
      
      do.call("lines", basal_args)
    }
  }
}

#' Create Contour Plot of Pairwise Invasibility Analysis Results
#' 
#' Function \code{plot.adaptInv} plots pairwise invasibility contour plots. This
#' function is based on code derived from Roff's Modeling Evolution: An
#' Introduction to Numerical Methods (2010, Oxford University Press).
#' 
#' @name plot.adaptInv
#' 
#' @param x An \code{adaptInv} object, created with function
#' \code{\link{invade3}()}.
#' @param xlab The x axis label forthe contour plot. Defaults to
#' \code{Resident}.
#' @param ylab The y axis label forthe contour plot. Defaults to
#' \code{Invader}.
#' @param res_variant The number of the variant representing the resident
#' subpopulation.
#' @param inv_variant The number of the variant representing the mutant
#' subpopulation.
#' @param repl The replicate number to plot, in the \code{fitness} data frame
#' within the \code{adaptInv} object entered in argument \code{x}.
#' @param pip A logical value indicating whether to produce a pairwise
#' invasibility plot. If \code{FALSE}, then will produce a diagnostic population
#' size plot. Defaults to \code{TRUE}.
#' @param elast A logical value indicating whether to produce an elasticity
#' plot. Such plots can only be produced when trait optimization is performed
#' during invasibility analysis. Defaults to \code{FALSE}.
#' @param run An integer giving the run to plot if \code{pip = FALSE}.
#' @param filled A logical value indicating whether to produce a filled contour
#' plot, or a standard contour plot. Defaults to \code{TRUE}, but reverts if
#' invader fitness is consistently positive, or consistently negative, relative
#' to the resident.
#' @param plot.title A title for the plot.
#' @param plot.axes 
#' @param axes A logical value indicating whether to include axis lines.
#' Defaults to \code{TRUE}.
#' @param frame.plot A logical value indicating whether to frame the plot.
#' @param auto_ylim A logical value indicating whether the maximum of the y axis
#' should be determined automatically. Defaults to \code{TRUE}, but reverts to
#' \code{FALSE} if any setting for \code{ylim} is given. Used only if
#' \code{pip = FALSE}.
#' @param auto_col A logical value indicating whether to shift the color of
#' lines associated with each patch automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{col} is given. Used only if
#' \code{pip = FALSE}.
#' @param auto_lty A logical value indicating whether to shift the line type
#' associated with each replicate automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{lty} is given. Used only if
#' \code{pip = FALSE}.
#' @param auto_title A logical value indicating whether to add a title to each
#' plot. The plot is composed of the concatenated population and patch names.
#' Defaults to \code{FALSE}. Used only if \code{pip = FALSE}.
#' @param ... Other parameters used by functions \code{plot.default()}.
#' 
#' @return A contour plot showing the overall fitness dynamics of the invader
#' variant, assuming a pairwise invasibility analysis.
#' 
#' @section Notes:
#' By default, function \code{plot.adaptInv} produces a filled contour plot in
#' which grey regions show where the invader has positive fitness relative to
#' the resident, and white regions show where the invader has negative fitness
#' relative to the resident. Fitness here refers to the Lyapunov coefficient,
#' calculated over the final \code{fitness_times} in the original call to
#' function \code{\link{invade3}()}.
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
#' plot(cyp_inv)
#' 
#' @export
plot.adaptInv <- function(x, xlab = "Resident", ylab = "Invader",
  res_variant = 1, inv_variant = 2, repl = 1, pip = TRUE, elast = FALSE,
  run = 1, filled = TRUE, plot.title, plot.axes, axes = TRUE, frame.plot = TRUE,
  auto_ylim = TRUE, auto_col = TRUE, auto_lty = TRUE, auto_title = FALSE, ...) {
  
  asp <- NA
  las <- 1
  xaxs <- yaxs <- "i"
  xlab <- ylab <- NULL
  chosen_colours <- c("white", "darkgrey")
  
  if (!axes) frame.plot <- FALSE
  
  if (length(pip) > 1) pip <- pip[1]
  if (length(elast) > 1) elast <- elast[1]
  
  further_args <- list(...)
  if (length(further_args) == 0) further_args <- list()
  found_terms <- names(further_args)
  
  if (is.character(repl)) {
    stop("Argument repl not understood.", call. = FALSE)
  }
  
  if (!pip & !elast) {
    if (!is.element("type", found_terms)) {
      further_args$type <- "l"
    }
    if (is.element("col", found_terms)) {
      auto_col <- FALSE
    }
    if (is.element("lty", found_terms)) {
      auto_lty <- FALSE
    }
    if (is.element("ylim", found_terms)) {
      auto_ylim <- FALSE
    }
    if (is.element("main", found_terms)) {
      auto_title <- FALSE
    }
    basal_args <- further_args
    
    if (!is.element("ylab", names(further_args))) {
      further_args$ylab <- "Population size"
    }
    if (!is.element("xlab", names(further_args))) {
      further_args$xlab <- "Time"
    }
    
    used_col <- 1
    
    if (auto_title) {
      used_string <- paste("pop", x$labels[i, 1], "patch", x$labels[i, 2])
      further_args$main <- used_string
    }
    
    N_length <- length(x$N_out)
    if (repl < 1 | repl > N_length) {
      stop("Invalid replciate chosen.", call. = FALSE)
    }
    
    N_dims <- dim(x$N_out[[repl]])
    if (run < 1 | run > N_dims[3]) {
      stop("Invalid run chosen.", call. = FALSE)
    }
    
    used_N_mat <- x$N_out[[repl]][,,run]
    
    used_col <- 1
    used_lty <- 1
    
    if (auto_ylim) {
      further_args$ylim <- c(0, max(used_N_mat, na.rm = TRUE))
    }
    
    if (auto_col) {
      further_args$col <- used_col
      basal_args$col <- used_col
    }
    if (auto_lty) {
      further_args$lty <- used_lty
    }
    
    num_pops <- dim(used_N_mat)[1]
    
    start_1 <- used_N_mat[1,c(min(which(used_N_mat[1,] > 0.)): max(which(used_N_mat[1,] > 0.)))]
    c_xy <- xy.coords(x = c(min(which(used_N_mat[1,] > 0.)): max(which(used_N_mat[1,] > 0.))), y = start_1)
    further_args$x <- c_xy
    
    do.call("plot.default", further_args)
    
    if (num_pops > 1) {
      for (j in c(2:num_pops)) {
        if (auto_lty) {
          used_lty <- used_lty + 1;
          basal_args$lty <- used_lty
        }
        if (auto_col) {
          used_col <- used_col + 1;
          if (used_col > length(palette())) used_col <- 1;
          basal_args$col <- used_col
        }
        used_col <- used_col + 1;
        
        start_j <- used_N_mat[j,c(min(which(used_N_mat[j,] > 0.)): max(which(used_N_mat[j,] > 0.)))]
        #basal_args$x <- xy.coords(x = c(min(which(used_N_mat[j,] > 0.)):max(which(used_N_mat[j,] > 0.))), y = start_j)
        
        basal_args$x <- c(min(which(used_N_mat[j,] > 0.)): max(which(used_N_mat[j,] > 0.)))
        basal_args$y <- start_j
        
        do.call("lines", basal_args)
      }
    }
  }
  
  if (pip) {
    if (is.element("las", found_terms)) las <- further_args$las
    if (is.element("xaxs", found_terms)) xaxs <- further_args$xaxs
    if (is.element("yaxs", found_terms)) yaxs <- further_args$yaxs
    if (is.element("asp", found_terms)) asp <- further_args$asp
    
    if (is.element("xlab", found_terms)) {
      xlab <- further_args$xlab
    } else xlab <- "Resident Value"
    if (is.element("ylab", found_terms)) {
      ylab <- further_args$ylab
    } else ylab <- "Invader Value"
    
    if (is.element("col", found_terms)) {
      found_colours <- further_args$col
      if (length(found_colours) < 2) stop("Argument col needs 2 color choices.", call. = FALSE)
      chosen_colours <- found_colours[c(1,2)]
    }
    
    if (!is.element("fitness", names(x))) {
      stop("Argument x does not appear to be an adaptInv object.", call. = FALSE)
    }
    
    if (res_variant == inv_variant) {
      stop("Arguments res_variant and inv_variant cannot be equal.", call. = FALSE)
    }
    fit_var_name_res <- paste0("variant", res_variant)
    fit_var_name_fitres <- paste0("fitness_variant", res_variant)
    fit_var_name_inv <- paste0("variant", inv_variant)
    fit_var_name_fitinv <- paste0("fitness_variant", inv_variant)
    
    if (!is.element(fit_var_name_inv, names(x$fitness)) | !is.element(fit_var_name_fitinv, names(x$fitness))) {
      stop("Pairwise invasibility analysis requires data on two variants.", call. = FALSE)
    }
    if (!is.element(fit_var_name_res, names(x$fitness)) | !is.element(fit_var_name_fitres, names(x$fitness))) {
      stop("Pairwise invasibility analysis requires data on two variants.", call. = FALSE)
    }
    
    correct_indices <- which(x$fitness$rep == repl)
    if (length(correct_indices) == 0) {
      stop("Replicate entered in argument repl could not be found.", call. = FALSE)
    }
    
    res_column <- which(names(x$fitness) == fit_var_name_res)[1]
    resfit_column <- which(names(x$fitness) == fit_var_name_fitres)[1]
    inv_column <- which(names(x$fitness) == fit_var_name_inv)[1]
    invfit_column <- which(names(x$fitness) == fit_var_name_fitinv)[1]
    
    resident_index <- x$fitness[correct_indices, res_column]
    invader_index <- x$fitness[correct_indices, inv_column]
    invader_fitness <- x$fitness[correct_indices, invfit_column]
    
    unique_variants <- unique(resident_index)
    num_variants <- length(unique_variants)
    uv_range <- range(unique_variants, finite = TRUE)
    
    ifmat <- matrix(invader_fitness, ncol = num_variants)
    ifmat_simplified <- ifmat
    ifmat_simplified[which(ifmat < 0)] <- -1
    ifmat_simplified[which(ifmat > 0)] <- 1
    found_levels <- unique(as.vector(ifmat_simplified))
    xylim <- range(found_levels, finite = TRUE)
    
    if (filled & length(found_levels) == 1) {
      filled <- FALSE
      message("Invader fitness is consistently either negative or positive, relative to the resident.")
      message("Cannot produce a filled contour plot.")
    }
    
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    
    if (!filled) {
      contour(unique_variants, unique_variants, ifmat, col = "black")
    } else {
      w <- (3 + mar.orig[2L]) * par("csi") * 2.54
      layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
      par(las = las)
      mar <- mar.orig
      par(mar=mar)
      plot.new()
      plot.window(xlim = c(0, 1), ylim = uv_range, xaxs = "i", yaxs = "i")
      
      plot.new()
      plot.window(xlim = uv_range, ylim = uv_range, "", xaxs = xaxs,
        yaxs = yaxs, asp = asp, xlab = xlab, ylab = ylab)
      
      .filled.contour(unique_variants, unique_variants, ifmat_simplified,
        levels = c(-1, 0, 1), col = chosen_colours)
    }
    
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(unique_variants, side = 1)
        Axis(unique_variants, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) box()
    if (missing(plot.title)) title(...)
    else plot.title
    invisible()
  }
  
  if (elast) {
    if (!exists("optim", x)) {
      stop("No optimization data included in object.", call. = FALSE)
    }
    
    used_x <- x$optim$fitness
    
    if (!exists("fitness_variant2_e995", used_x)) {
      stop("No optimization fitness data included in object.", call. = FALSE)
    }
    
    if (!is.element("type", found_terms)) {
      further_args$type <- "l"
    }
    if (is.element("col", found_terms)) {
      auto_col <- FALSE
    }
    if (is.element("lty", found_terms)) {
      auto_lty <- FALSE
    }
    if (is.element("ylim", found_terms)) {
      auto_ylim <- FALSE
    }
    if (is.element("main", found_terms)) {
      auto_title <- FALSE
    }
    basal_args <- further_args
    
    if (!is.element("ylab", names(further_args))) {
      further_args$ylab <- "Invader fitness"
    }
    if (!is.element("xlab", names(further_args))) {
      further_args$xlab <- "Trait value"
    }
    
    used_col <- 1
    
    if (auto_title) {
      used_string <- "Trait optimization elasticity"
      further_args$main <- used_string
    }
    
    used_fitness <- used_x$fitness_variant2_e995
    
    used_col <- 1
    used_lty <- 1
    
    if (auto_ylim) {
      further_args$ylim <- c(min(used_fitness, na.rm = TRUE), max(used_fitness, na.rm = TRUE))
    }
    
    if (auto_col) {
      further_args$col <- used_col
      basal_args$col <- used_col
    }
    if (auto_lty) {
      further_args$lty <- used_lty
    }
    
    c_xy <- xy.coords(x = c(1:length(used_fitness)), y = used_fitness)
    further_args$x <- c_xy
    
    do.call("plot.default", further_args)
    
    basal_args$x <- c(1:length(used_fitness))
    basal_args$y <- rep(0, length(used_fitness))
    basal_args$lty = 2
    
    do.call("lines", basal_args)
  }
}
