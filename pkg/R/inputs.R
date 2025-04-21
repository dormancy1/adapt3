#' Create an Equivalence Vector for Each Population
#' 
#' Function \code{equiv_input()} creates a data frame summarizing the degree to
#' which an individual in each stage of a life history is equivalent to a
#' standard individual.
#' 
#' @name equiv_input
#' 
#' @param mpm The lefkoMat object to be used in projection. Can be an example
#' MPM if function-based projection is planned.
#' @param stage2 A vector showing the name or number of a stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Abbreviations for groups of stages are also usable (see Notes).
#' This input is required for all stage-based and age-by-stage MPMs. Defaults to
#' \code{NA}.
#' @param stage1 A vector showing the name or number of a stage in occasion
#' \emph{t}-1 that should be set to a positive number of individuals in the
#' start vector. Abbreviations for groups of stages are also usable (see Notes).
#' This is only used for historical MPMs, since the rows of hMPMs correspond to
#' stage-pairs in times \emph{t} and \emph{t}-1 together. Only required for
#' historical MPMs, and will result in errors if otherwise used.
#' @param age2 A vector showing the age of each respective stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Only used for Leslie and age-by-stage MPMs. Defaults to \code{NA}.
#' @param value A vector showing the values, in order, of the number of
#' individuals set for the stage or stage-pair in question. Defaults to
#' \code{1}.
#' 
#' @return A list of class \code{adaptEq}, with four objects, which can be
#' used as input in function \code{\link{project3}()}. The last three include
#' the \code{ahstages}, \code{hstages}, and \code{agestages} objects from the
#' \code{lefkoMat} object supplied in \code{mpm}. The first element in the list
#' is a data frame with the following variables:
#' 
#' \item{stage2}{Stage at occasion \emph{t}.}
#' \item{stage_id_2}{The stage number associated with \code{stage2}.}
#' \item{stage1}{Stage at occasion \emph{t}-1, if historical. Otherwise NA.}
#' \item{stage_id_1}{The stage number associated with \code{stage1}.}
#' \item{age2}{The age of individuals in \code{stage2} and, if applicable,
#' \code{stage1}. Only used in age-by-stage MPMs.}
#' \item{row_num}{A number indicating the respective starting vector element.}
#' \item{value}{Number of individuals in corresponding stage or stage-pair.}
#' 
#' @section Notes:
#' Entries in \code{stage2}, and \code{stage1} can include abbreviations for
#' groups of stages. Use \code{rep} if all reproductive stages are to be used,
#' \code{nrep} if all mature but non-reproductive stages are to be used,
#' \code{mat} if all mature stages are to be used, \code{immat} if all immature
#' stages are to be used, \code{prop} if all propagule stages are to be used,
#' \code{npr} if all non-propagule stages are to be used, \code{obs} if all
#' observable stages are to be used, \code{nobs} if all unobservable stages are
#' to be used, and leave empty or use \code{all} if all stages in stageframe are
#' to be used.
#' 
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' 
#' e3m_sv <- start_input(ehrlen3mean, stage2 = "Sd", stage1 = "Sd", value = 1000)
#' 
#' lathproj <- projection3(ehrlen3, nreps = 5, times = 100, stochastic = TRUE,
#'   start_frame = e3m_sv)
#' 
#' @export
equiv_input <- function(mpm, stage2 = NA, stage1 = NA, age2 = NA, value = 1.0) {
  
  mpmrows <- stage2_id <- stage1_id <- start_vec <- full_length <- NULL
  
  if (all(!is(mpm, "lefkoMat"))) {
    stop("A regular lefkoMat object is required as input.", call. = FALSE)
  }
  
  if (all(is.na(stage2)) & all(is.na(age2))) {
    stop("Options stage2 and age2 cannot both be set to NA.", call. = FALSE)
  }
  if (all(is.null(stage2)) & all(is.null(age2))) {
    stop("Options stage2 and age2 cannot both be empty.", call. = FALSE)
  }
  
  if (!is.element("stage", names(mpm$ahstages))) {
    stop("Stageframe appears to be modified. Please make sure that a stage
      column exists holding stage names.", call. = FALSE)
  }
  
  if (all(is.na(mpm$hstages)) | all(is.null(mpm$hstages))) {
    historical <- FALSE
  } else {
    historical <- TRUE
  }
  if (all(is.na(mpm$agestages)) | all(is.null(mpm$agestages))) {
    agebystage <- FALSE
  } else {
    agebystage <- TRUE
  }
  
  if (historical & all(is.na(stage1))) {
    stop("Historical projection analysis requires that stage in time t-1 be
      designated for all stage pairs.", call. = FALSE)
  } else if (!historical & !all(is.na(stage1))) {
    stop("Ahistorical projection analysis cannot include designated stages in
      time t-1.", call. = FALSE)
  }
  
  full_length <- max(length(stage2), length(stage1), length(age2), length(value))
  
  if (length(value) == 1 & full_length > 1) {
    value <- rep(value, full_length)
  }
  
  if (length(stage2) == 1 & full_length > 1) {
    stage2 <- rep(stage2, full_length)
  }
  
  if (length(stage1) == 1 & full_length > 1) {
    stage1 <- rep(stage1, full_length)
  }
  
  if (length(age2) == 1 & full_length > 1) {
    age2 <- rep(age2, full_length)
  }
  
  if ((all(is.na(stage2)) | all(is.null(stage2))) & (all(is.na(age2)) | all(is.null(age2)))) {
    stop("Either stage2 or age2 must be provided.", call. = FALSE)
  }
  
  if (all(is.character(stage2))) {
    unknown_stage2 <- which(!is.element(tolower(stage2), c(tolower(mpm$ahstages$stage),
        c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "obs", "nobs"))))
    if (length(unknown_stage2) > 0) {
      stop(paste0("Unknown stage designations used in stage2: ",
        stage2[unknown_stage2]), call. = FALSE)
    }
    
    reassessed <- apply(as.matrix(c(1:length(stage2))), 1, function(X) {
      if (!is.na(stage2[X])) {
        if (is.element(stage2[X], mpm$ahstages$stage)) {
          shrubbery.small <- cbind.data.frame(stage2 = stage2[X], stage1 = stage1[X],
            age2 = age2[X], value = value[X], stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (is.element(stage2[X], as.character(mpm$ahstages$stage_id))) {
          shrubbery.small <- cbind.data.frame(stage2 = as.numeric(stage2[X]),
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "rep") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "all") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage,
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nrep") {
          shrubbery.small <- cbind.data.frame(
            stage2 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                which(mpm$ahstages$matstatus == 1))],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) =="mat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "immat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "prop") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "npr") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "obs") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nobs") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        }
      }
    })
    shrubbery <- do.call(rbind.data.frame, reassessed)
    
  } else if (all(is.numeric(stage2)) & !any(is.na(stage2))) {
    stage2_id <- stage2
    
    if (any(stage2_id > max(mpm$ahstages$stage_id)) | any(stage2_id < min(mpm$ahstages$stage_id))) {
      stop("Unknown stage2 codes used.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(stage2_id), 1, function(X) {
      return(mpm$ahstages$stage[X])
    })
    
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
    
  } else if (all(is.numeric(age2)) & !any(is.na(age2))) {
    if (!(all(is.na(stage2)))) {
      stop("Leslie MPMs should not be entered with the stage2 option set to
        values other than NA.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(age2), 1, function(X) {
      cross_ref <- which(mpm$ahstages$stage_id == X)
      return(mpm$ahstages$stage[cross_ref])
    })
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
    
  } else {
    stop("Input stage2 codes do not conform to accepted inputs.", call. = FALSE)
  }
  
  if (historical) {
    if (all(is.character(shrubbery$stage1)) & !all(is.na(shrubbery$stage1))) {
      
      unknown_stage1 <- which(!is.element(tolower(stage1), c(tolower(mpm$ahstages$stage),
          c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "obs", "nobs", "almostborn"))))
      if (length(unknown_stage1) > 0) {
      stop(paste0("Unknown stage designations used in stage1: ",
        stage1[unknown_stage1]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$stage1[X])) {
          if (is.element(shrubbery$stage1[X], mpm$ahstages$stage)) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (is.element(stage1[X], as.character(mpm$ahstages$stage_id))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = as.numeric(shrubbery$stage1[X]), age2 = age2[X],
              value = value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "rep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage, age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nrep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                  which(mpm$ahstages$matstatus == 1))],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) =="mat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "immat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "prop") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "npr") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "obs") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nobs") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      })
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
      
    } else if (all(is.numeric(shrubbery$stage1)) & !any(is.na(shrubbery$stage1))) {
      stage1_id <- shrubbery$stage1
      
      if (any(stage1_id > max(mpm$ahstages$stage_id)) | any(stage1_id < min(mpm$ahstages$stage_id))) {
        stop("Unknown stage1 codes used.", call. = FALSE)
      }
      
      stage1 <- apply(as.matrix(stage1_id), 1, function(X) {
        return(mpm$ahstages$stage[X])
      })
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (agebystage) {
    if (all(is.character(shrubbery$age2)) & !all(is.na(shrubbery$age2))) {
      
      unknown_age2 <- which(!is.element(tolower(age2), c(tolower(mpm$agestages$age),
          "all")))
      if (length(unknown_age2) > 0) {
        stop(paste0("Unknown age designations used in age2: ",
          stage1[unknown_age2]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$age2[X])) {
          common_ages <- unique(mpm$agestages$age[which(mpm$agestages$stage == shrubbery$stage2[X])])
          
          if (is.element(shrubbery$age2[X], as.character(common_ages))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = as.numeric(shrubbery$age2[X]),
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$age2[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = common_ages,
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      }) 
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
    } else if (all(is.numeric(shrubbery$age2)) & !any(is.na(shrubbery$age2))) {
      if (any(shrubbery$age2 > max(mpm$agestages$age)) | any(shrubbery$age2 < min(mpm$agestages$age))) {
        stop("Unknown age2 values used.", call. = FALSE)
      }
      
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = shrubbery$stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (!all(is.numeric(value))) {
    stop("Object value must be composed only of valid numbers.", call. = FALSE)
  }
  
  shrubbery$stage2_id <- apply(as.matrix(shrubbery$stage2), 1, function(X) {
    return(mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)])
  })
  shrubbery$stage1_id <- apply(as.matrix(shrubbery$stage1), 1, function(X) {
    possible_option <- mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)]
    if (length(possible_option) > 0) return(possible_option) else return(NA)
  })
  
  full_length <- dim(shrubbery)[1]
  
  if (!historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$ahstages)[1]) {
      stop("This ahistorical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- shrubbery$stage2_id
    
  } else if (agebystage & !historical) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$agestages)[1]) {
      stop("This age-by-stage mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    if (any(is.na(shrubbery$age2)) | any(!is.numeric(shrubbery$age2))) {
      stop("Option age2 must include only numbers for age-by-stage MPMs.",
        call. = FALSE)
    }
    
    if (any(shrubbery$age2 < min(mpm$agestages$age)) | any(shrubbery$age2 > max(mpm$agestages$age))) {
      stop("Option age2 can only take ages shown in element agestages within the input MPM.",
        call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$agestages$stage_id == shrubbery$stage2_id[X])
      vec1 <- which(mpm$agestages$age == shrubbery$age2[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else if (historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$hstages)[1]) {
      stop("This historical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$hstages$stage_id_2 == shrubbery$stage2_id[X])
      vec1 <- which(mpm$hstages$stage_id_1 == shrubbery$stage1_id[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else {
    stop("Format of mpm not recognized.", call. = FALSE)
  }
  
  output_tab <- cbind.data.frame(shrubbery$stage2, shrubbery$stage2_id,
    shrubbery$stage1, shrubbery$stage1_id, shrubbery$age2, start_vec,
    shrubbery$value, stringsAsFactors = FALSE)
  
  names(output_tab) <- c("stage2", "stage_id_2", "stage1", "stage_id_1", "age2",
    "row_num", "value")
  
  if (!historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2")])
  } else if (historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "stage1", "stage_id_1")])
  } else if (!historical & agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "age2")])
  }
  if (dim(out_check)[1] < dim(output_tab)[1]) {
    warning("Some stages, stage-pairs, or age-stages appear to be listed
      multiple times. This may cause errors in analysis.", call. = FALSE)
  }

  class(output_tab) <- append(class(output_tab), "adaptEq")
  
  return(output_tab)
}
