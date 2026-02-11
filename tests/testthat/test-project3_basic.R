test_that("project3() basic example works", {
  library(lefko3)
  data(cypdata)
  
  data(cypa_data)
  
  sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
  stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
    "XLg")
  repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
  immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
  binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
  
  cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
    repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
    propstatus = propvector, immstatus = immvector, indataset = indataset,
    binhalfwidth = binvec)
  
  cycaraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
    patchidcol = "patch", individcol = "plantid", blocksize = 4,
    sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
    repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
    stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
    NRasRep = TRUE)
    
  cyparaw_v1 <- verticalize3(data = cypa_data, noyears = 18, firstyear = 1994,
    individcol = "plant_id", blocksize = 2, sizeacol = "Inf.94",
    sizebcol = "Veg.94", repstracol = "Inf.94", fecacol = "Inf.94",
    stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
    NRasRep = TRUE)
  
  cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
      "XSm", "Sm", "SD", "P1"),
    stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
      "rep"),
    eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
    eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
    givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
    multiplier = c(1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5),
    type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
    stageframe = cypframe_raw, historical = FALSE)
  cyp_supp_list1 <- list(cypsupp2r, cypsupp2r)
  
  cycamatrix2r <- rlefko2(data = cycaraw_v1, stageframe = cypframe_raw, 
    year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
    size = c("size3added", "size2added"), supplement = cypsupp2r,
    yearcol = "year2", patchcol = "patchid", indivcol = "individ")
  
  cypamatrix2r <- rlefko2(data = cyparaw_v1, stageframe = cypframe_raw, 
    year = "all", stages = c("stage3", "stage2", "stage1"),
    size = c("size3added", "size2added"), supplement = cypsupp2r,
    yearcol = "year2", patchcol = "patchid", indivcol = "individ")
  
  cyp_mpm_list <- list(cycamatrix2r, cypamatrix2r)
  
  cyca2_start <- start_input(cycamatrix2r, stage2 = c("SD", "P1", "P2"),
    value = c(500, 100, 200))
  cypa2_start <- start_input(cypamatrix2r, stage2 = c("SD", "P1", "P2"),
    value = c(5000, 1000, 2000))
  cyp_start_list <- list(cyca2_start, cypa2_start)
  
  cyp2_dv <- density_input(cypamatrix2r, stage3 = c("SD", "P1"),
    stage2 = c("rep", "rep"), style = c(1, 1), alpha = c(0.5, 1.2),
    beta = c(1.0, 2.0), type = c(2, 1))
  cyp_dv_list <- list(cyp2_dv, cyp2_dv)
  
  cyp_comm_proj <- project3(mpms = cyp_mpm_list, starts = cyp_start_list,
    density = cyp_dv_list, times = 10)
  
  expect_true(class(cyp_comm_proj) == "adaptProj")
  expect_true(cyp_comm_proj$agg_density[1,11] > 1)
})
