test_that("project3() and batch_project3() advanced examples work", {
  # Trials with raw and function-based MPMs in both ahistorical and Leslie formats to test matrix community model projection
  library(lefko3)
  data(cypdata)
  
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
  
  sizevector <- c(0, 0, 3.0, 15)
  stagevector <- c("P1", "D", "Sm", "Lg")
  repvector <- c(0, 0, 1, 1)
  obsvector <- c(0, 0, 1, 1)
  matvector <- c(0, 1, 1, 1)
  immvector <- c(1, 0, 0, 0)
  indataset <- c(0, 1, 1, 1)
  binvec <- c(0, 0.5, 2.5, 9.5)
  
  cypframe_small_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
    repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
    immstatus = immvector, indataset = indataset, binhalfwidth = binvec)
  
  cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
    patchidcol = "patch", individcol = "plantid", blocksize = 4,
    sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
    repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
    stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
    NRasRep = TRUE)
  
  cypraw_v2 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
    patchidcol = "patch", individcol = "plantid", blocksize = 4,
    sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
    repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
    stageassign = cypframe_small_raw, stagesize = "sizeadded", NAas0 = TRUE,
    NRasRep = TRUE)
  
  cypraw_v3 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
    patchidcol = "patch", individcol = "plantid", blocksize = 4,
    sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
    repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
    NAas0 = TRUE, NRasRep = TRUE)
  
  # Here we use supplemental() to provide overwrite and reproductive info
  suppressWarnings(cypsupp2r <- supplemental(
    stage3 = c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "SD", "P1"),
    stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep"),
    eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
    eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
    givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
    multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 1500, 500),
    type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
    stageframe = cypframe_raw, historical = FALSE))
  
  cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
    year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
    size = c("size3added", "size2added"), supplement = cypsupp2r,
    yearcol = "year2", patchcol = "patchid", indivcol = "individ")
  cypmean <- lmean(cypmatrix2r)
  
  suppressWarnings(cypsupp2r_small <- supplemental(stage3 = c("D", "Sm", "Lg", "P1"),
    stage2 = c("P1", "P1", "P1", "rep"), eststage3 = c(NA, "Sm", "Lg", NA),
    eststage2 = c(NA, "D", "D", NA), givenrate = c(0.05, NA, NA, NA),
    offset = c(NA, NA, -0.1, NA), multiplier = c(NA, NA, NA, 0.5),
    type =c(1, 1, 1, 3), stageframe = cypframe_small_raw, historical = FALSE))
  cypmatrix2r_small <- rlefko2(data = cypraw_v2, stageframe = cypframe_small_raw, 
    year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
    size = c("size3added", "size2added"), supplement = cypsupp2r_small,
    yearcol = "year2", patchcol = "patchid", indivcol = "individ")
  cypmean_small <- lmean(cypmatrix2r_small)
  
  cypmatrixL_small <- rleslie(data = cypraw_v3, start_age = 1, last_age = 4,
    continue = TRUE, fecage_min = 3, year = "all", pop = NA, patch = "all",
    yearcol = "year2", patchcol = "patchid", indivcol = "individ")
  
  cypmodels2 <- modelsearch(cypraw_v1, historical = FALSE, approach = "mixed", 
    vitalrates = c("surv", "obs", "size", "repst", "fec"),
    sizedist = "negbin", size.trunc = TRUE, fecdist = "poisson", fec.zero = TRUE,
    suite = "main", size = c("size3added", "size2added"), quiet = "partial")
  cypmodels2_patch <- modelsearch(cypraw_v1, historical = FALSE, approach = "mixed", 
    vitalrates = c("surv", "obs", "size", "repst", "fec"),
    sizedist = "negbin", size.trunc = TRUE, fecdist = "poisson", fec.zero = TRUE,
    suite = "main", size = c("size3added", "size2added"), patch = "patchid",
    quiet = "partial")
  
  cypmatrix2f <- flefko2(stageframe = cypframe_raw, supplement = cypsupp2r, 
    modelsuite = cypmodels2, data = cypraw_v1, err_check = TRUE)
  cypmatrix2f_A <- flefko2(stageframe = cypframe_raw, supplement = cypsupp2r, 
    modelsuite = cypmodels2_patch, data = cypraw_v1, patch = "A", err_check = TRUE)
  cypmatrix2f_patch <- flefko2(stageframe = cypframe_raw, supplement = cypsupp2r, 
    modelsuite = cypmodels2_patch, data = cypraw_v1)
  cypmatrix2f_small <- flefko2(stageframe = cypframe_small_raw, supplement = cypsupp2r_small, 
    modelsuite = cypmodels2, data = cypraw_v1, err_check = TRUE)
  cypmatrix2f_A <- flefko2(stageframe = cypframe_small_raw, supplement = cypsupp2r_small, 
    modelsuite = cypmodels2_patch, data = cypraw_v1, patch = "A", err_check = TRUE)
  cypmatrix2f_patch <- flefko2(stageframe = cypframe_small_raw, supplement = cypsupp2r_small, 
    modelsuite = cypmodels2_patch, data = cypraw_v1)
  cypmodelsL <- modelsearch(cypraw_v3, approach = "mixed", suite = "age",
    vitalrates = c("surv", "fec"), age = "obsage", fecdist = "poisson",
    fec.zero = TRUE, quiet = "partial")
  
  cyp_vrm1 <- miniMod(cypmodels2, hfv_data = cypraw_v1)
  cyp_vrm2p <- miniMod(cypmodels2_patch, hfv_data = cypraw_v2)
  cypL_vrm <- miniMod(cypmodelsL, hfv_data = cypraw_v3)
  cyp_vrms1 <- list(cyp_vrm1, cyp_vrm2p, cypL_vrm)
  
  cypmatrix2fvrm <- flefko2(stageframe = cypframe_raw, supplement = cypsupp2r, 
    modelsuite = cyp_vrm1, data = cypraw_v1, err_check = TRUE)
  cypmatrix2f_patchvrm <- flefko2(stageframe = cypframe_small_raw, supplement = cypsupp2r_small, 
    modelsuite = cyp_vrm2p, data = cypraw_v2)
  cypmatrixLf <- fleslie(last_age = 4, fecage_min = 3, continue = TRUE,
    modelsuite = cypmodelsL, data = cypraw_v3)
  
  cyp_mpms1 <- list(cypmatrix2r, cypmatrix2r_small, cypmatrixL_small, cypmatrixLf)
  
  trial_supplement_cyp_age_small <- sup_skeleton(2)
  trial_supplement_cyp_age_small$stage3 <- c("Age2", "")
  trial_supplement_cyp_age_small$stage2 <- c("Age1", NA)
  trial_supplement_cyp_age_small$age2 <- c(1, 3)
  trial_supplement_cyp_age_small$offset <- c(0.003, 0.006)
  trial_supplement_cyp_age_small$convtype <- c(1, 1)
  cypmatrixLf_trial <- fleslie(last_age = 4, fecage_min = 3, continue = TRUE,
    modelsuite = cypmodelsL, data = cypraw_v3, supplement = trial_supplement_cyp_age_small)
  
  c2d_4 <- density_input(cypmean, stage3 = c("P1", "P1"), stage2= c("SD", "rep"),
    style = 1, time_delay = 1, alpha = 1, beta = 0.0005, type = c(2, 2))
  c2d_4a <- density_input(cypmean_small, stage3 = c("P1", "P1"), stage2= c("P1", "rep"),
    style = 1, time_delay = 1, alpha = 1, beta = 0.0005, type = c(2, 2))
  cypL_dv_1 <- density_input(cypmatrixLf, stage3 = c("Age1"), stage2 = c("rep"),
    style = c(1), alpha = c(0.5), beta = c(1.0), type = c(2))
  cypL_dv_2 <- density_input(cypmatrixLf, age2 = 1,
    alpha = c(0.5), beta = c(1.0), type = c(2))
  
  cypL_dvr <- density_vr(c(T, F, F, F, F, F, T, F, F, F, F, F, F, F),
    style = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    alpha = c(0.5, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0),
    beta = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
  
  cyp_start1 <- start_input(cypmatrix2r, stage2 = c("SD", "P1", "D"),
    value = c(100, 200, 4))
  cyp_start2 <- start_input(cypmatrix2r, stage2 = c("SD", "P1", "D"),
    value = c(10, 2000, 40))
  cypL_start_1 <- start_input(cypmatrixLf, stage2 = c("Age1"), value = c(200))
  cypL_start_2 <- start_input(cypmatrixLf, age2 = c(1), value = c(200))
  
  cyp_start <- list(cyp_start1, cyp_start2, cypL_start_1, cypL_start_2)
  
  cyp_dv <- density_input(cypmatrix2r, stage3 = c("P1", "P2"),
    stage2 = c("rep", "P1"), style = c(1, 1), alpha = c(0.5, 1.2),
    beta = c(1.0, 0), type = c(2, 1))
  
  cyp_dv_small <- density_input(cypmatrix2r_small, stage3 = c("P1", "Sm"),
    stage2 = c("rep", "P1"), style = c(1, 1), alpha = c(0.5, 1.2),
    beta = c(1.0, 0), type = c(2, 1))
  cyp_dv_age <- density_input(cypmatrixL_small, stage3 = c("Age1", "Age2"),
    stage2 = c("rep", "Age1"), style = c(1, 1), alpha = c(0.5, 1.2),
    beta = c(1.0, 0), type = c(2, 1))
  
  cyp_dvr <- density_vr(c(T, F, F, F, F, F, T, F, F, F, F, F, F, F),
    style = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    alpha = c(0.5, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0),
    beta = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
  
  cyp1_eq <- stage_weight(cypmatrix2r,
    stage2 = c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg"),
    value = c(0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1, 1))
  cyp2_eq <- stage_weight(cypmatrix2r_small, stage2 = c("P1", "D", "Sm", "Lg"),
    value = c(0, 1, 1, 1))
  cyp3_eq <- stage_weight(cypmatrixL_small, age2 = c(1, 2, 3, 4),
    value = c(0, 1, 1, 1))
  cyp4_eq <- stage_weight(cypmatrixLf, age2 = c(1, 2, 3, 4),
    value = c(0, 1, 1, 1))
  
  cyp_mpm_eq <- list(cyp1_eq, cyp2_eq, cyp3_eq, cyp4_eq)
  cyp_vrm_eq <- list(cyp1_eq, cyp1_eq, cyp3_eq)
  
  # Simple raw ahistorical
  aaa1_prj <- project3(mpms =  cyp_mpms1, entry_time = c(0, 5, 8, 10), times = 15,
    integeronly = TRUE, nreps = 3)
  expect_true(dim(aaa1_prj$agg_density)[1] == 3)
  expect_true(dim(aaa1_prj$agg_density)[2] == 16)
  
  # Raw ahistorical with stage weights
  aaa1a_prj <- project3(mpms =  cyp_mpms1, entry_time = c(0, 5, 8, 10), times = 15,
    equivalence = cyp_mpm_eq, integeronly = TRUE, err_check = TRUE)
  expect_true(dim(aaa1a_prj$agg_density)[1] == 1)
  expect_true(dim(aaa1a_prj$agg_density)[2] == 16)
  expect_true(aaa1a_prj$agg_density[1,16] > 1)
  
  # Raw ahistorical with subsets of annual matrices used
  aaa1b_prj <- project3(mpms =  cyp_mpms1, entry_time = c(0, 5, 8, 10), times = 15,
    years = list(c("2005", "2006"), c("2006", "2007", "2008"), c("2005", "2006"), c("2005", "2008")),
    integeronly = TRUE)
  expect_true(dim(aaa1b_prj$agg_density)[1] == 1)
  expect_true(dim(aaa1b_prj$agg_density)[2] == 16)
  
  cyp_tweights <- c(0.01, 0.2, 0.1, 0.2, 0.49)
  cyp_tw_list1 <- list(cyp_tweights, cyp_tweights, cyp_tweights)
  cyp_tw_list2 <- list(cyp_tweights, cyp_tweights, cyp_tweights, cyp_tweights)
  
  # Raw stochastic ahistorical with altered environmental transitions
  aaa1d_prj <- project3(mpms =  cyp_mpms1, entry_time = c(0, 5, 8, 10), times = 15,
    integeronly = TRUE, stochastic = TRUE, tweights = cyp_tw_list2, nreps = 3)
  expect_true(dim(aaa1d_prj$agg_density)[1] == 3)
  expect_true(dim(aaa1d_prj$agg_density)[2] == 16)
  
  # Raw ahistorical with matrix element density dependence
  aaa1e_prj <- project3(mpms =  cyp_mpms1, entry_time = c(0, 5, 8, 10), times = 15,
    density = list(c2d_4, c2d_4a, cypL_dv_1, cypL_dv_2), integeronly = TRUE)
  expect_true(dim(aaa1e_prj$agg_density)[1] == 1)
  expect_true(dim(aaa1e_prj$agg_density)[2] == 16)
  expect_true(aaa1e_prj$agg_density[1,16] > 0)
  
  # Supplements used in post-processing
  cypsupp2r_post <- supplemental(stage3 = c("P2"), stage2 = c("P1"),
    offset = c(0.10), type =c(1), stageframe = cypframe_raw, historical = FALSE)
  cypsuppLr_post <- supplemental(stage3 = c("Age3"), stage2 = c("Age2"),
    offset = c(0.10), type =c(1), stageframe = cypmatrixLf$ahstages, historical = FALSE)
  
  # Test of post-processing supplements
  aaa1f_prj <- project3(mpms =  cyp_mpms1, entry_time = c(0, 5, 8, 10),
    times = 15, density = list(c2d_4, c2d_4a, cypL_dv_1, cypL_dv_2),
    supplements = list(cypsupp2r_post, NULL, NULL, cypsuppLr_post),
    integeronly = TRUE, err_check = "extreme")
  expect_true(dim(aaa1f_prj$agg_density)[1] == 1)
  expect_true(dim(aaa1f_prj$agg_density)[2] == 16)
  expect_true(aaa1f_prj$agg_density[1,16] > 1)
  
  # Batch projection
  new_supplement_cyp2_small <- sup_skeleton(2)
  new_supplement_cyp2_small$stage3 <- c("D", "Sm")
  new_supplement_cyp2_small$stage2 <- c("D", "D")
  new_supplement_cyp2_small$convtype <- c(1, 1)
  
  new_supplement_cyp_age_small <- sup_skeleton(2)
  new_supplement_cyp_age_small$stage3 <- c("Age2", "")
  new_supplement_cyp_age_small$stage2 <- c("Age1", NA)
  new_supplement_cyp_age_small$age2 <- c(1, 3)
  new_supplement_cyp_age_small$convtype <- c(1, 1)
  
  suppressWarnings(cypsupp2r_alt <- supplemental(
    stage3 = c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "SD", "P1", "Lg"),
    stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep", "SD"),
    eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA, NA),
    eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA, NA),
    givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA, NA),
    multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 1500, 500, NA),
    type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 1),
    stageframe = cypframe_raw, historical = FALSE))
  
  aaa1_prj_batch4 <- batch_project3(used_mpms = "all", all_elems = FALSE,
    mpms =  list(cypmatrix2r_small, cypmatrixL_small), times = 15,
    supplements = list(new_supplement_cyp2_small, new_supplement_cyp_age_small),
    integeronly = TRUE, nreps = 1)
  aaa14 <- summary(aaa1_prj_batch4, finalN_mean = TRUE, finalN_used = 10)
  expect_true(length(aaa14$final_N) == 40)
  expect_true(mean(aaa14$extinct_by, na.rm = TRUE) > 1)
  
  aaa1_prj_batch5 <- batch_project3(used_mpms = "all", all_elems = FALSE,
    mpms =  list(cypmatrix2r_small, cypmatrixL_small), times = 15,
    supplements = list(new_supplement_cyp2_small, new_supplement_cyp_age_small),
    density = list(cyp_dv_small, cyp_dv_age), integeronly = TRUE, nreps = 1)
  aaa15 <- summary(aaa1_prj_batch5, finalN_mean = TRUE, finalN_used = 10)
  expect_true(length(aaa15$final_N) == 40)
  expect_true(mean(aaa15$extinct_by, na.rm = TRUE) > 1)
})

