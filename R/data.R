#' Demographic Dataset of \emph{Cypripedium parviflorum} Population, in Horizontal
#' Format
#' 
#' A dataset containing the states and fates of a population of the terrestrial
#' orchid species \emph{Cypripedium parviflorum} (small yellow lady's slipper
#' orchids) in Illinois, USA. The dataset resulted from annual monitoring from
#' 1994 and 2011.
#' 
#' @docType data
#' 
#' @usage data(cypa_data)
#' 
#' @format A data frame with 1119 individuals and 37 variables. Each row 
#' corresponds to an unique individual, and each variable from \code{Inf.94} 
#' on refers to the state of the individual in a particular year.
#' 
#' \describe{
#'   \item{plant_id}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{Inf.94}{Number of inflorescences in 1994.}
#'   \item{Veg.94}{Number of stems without inflorescences in 1994.}
#'   \item{Inf.95}{Number of inflorescences in 1995.}
#'   \item{Veg.95}{Number of stems without inflorescences in 1995.}
#'   \item{Inf.96}{Number of inflorescences in 1996.}
#'   \item{Veg.96}{Number of stems without inflorescences in 1996.}
#'   \item{Inf.97}{Number of inflorescences in 1997.}
#'   \item{Veg.97}{Number of stems without inflorescences in 1997.}
#'   \item{Inf.98}{Number of inflorescences in 1998.}
#'   \item{Veg.98}{Number of stems without inflorescences in 1998.}
#'   \item{Inf.99}{Number of inflorescences in 1999.}
#'   \item{Veg.99}{Number of stems without inflorescences in 1999.}
#'   \item{Inf.00}{Number of inflorescences in 2000.}
#'   \item{Veg.00}{Number of stems without inflorescences in 2000.}
#'   \item{Inf.01}{Number of inflorescences in 2001.}
#'   \item{Veg.01}{Number of stems without inflorescences in 2001.}
#'   \item{Inf.02}{Number of inflorescences in 2002.}
#'   \item{Veg.02}{Number of stems without inflorescences in 2002.}
#'   \item{Inf.03}{Number of inflorescences in 2003.}
#'   \item{Veg.03}{Number of stems without inflorescences in 2003.}
#'   \item{Inf.04}{Number of inflorescences in 2004.}
#'   \item{Veg.04}{Number of stems without inflorescences in 2004.}
#'   \item{Inf.05}{Number of inflorescences in 2005.}
#'   \item{Veg.05}{Number of stems without inflorescences in 2005.}
#'   \item{Inf.06}{Number of inflorescences in 2006.}
#'   \item{Veg.06}{Number of stems without inflorescences in 2006.}
#'   \item{Inf.07}{Number of inflorescences in 2007.}
#'   \item{Veg.07}{Number of stems without inflorescences in 2007.}
#'   \item{Inf.08}{Number of inflorescences in 2008.}
#'   \item{Veg.08}{Number of stems without inflorescences in 2008.}
#'   \item{Inf.09}{Number of inflorescences in 2009.}
#'   \item{Veg.09}{Number of stems without inflorescences in 2009.}
#'   \item{Inf.10}{Number of inflorescences in 2010.}
#'   \item{Veg.10}{Number of stems without inflorescences in 2010.}
#'   \item{Inf.11}{Number of inflorescences in 2011.}
#'   \item{Veg.11}{Number of stems without inflorescences in 2011.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' library(lefko3)
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
#' cypraw_v1 <- verticalize3(data = cypa_data, noyears = 18, firstyear = 1994,
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
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", indivcol = "individ")
#'                        
#' lambda3(cypmatrix2r)
"cypa_data"
