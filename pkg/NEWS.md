# adapt3 2.0.0 (2026-02-09)

## NEW FEATURES

* Trait axes now include individual covariate values that can be used to create
  variants.
  
* Density dependence through the additive limit and absolute limit functions has
  been added to projection and invasion functions.
  
* Equivalence frames and density frames can now be added per annual matrix.
  
* Community projections with projection3() can now handle post-processing of
  pre-existing MPMs with argument `supplemental`.
  
* Function `batch_project3()` now runs batches of community projections using
  vectors of alterations input as matrix element offsets, given rate changes,
  and multipliers.
  
* Function `summary.adaptProjBatch()` now summarizes batch community
  projections.
  
* Function `project3()` now includes aggregate density in the output.

* Function `summary.adaptProj()` now also outputs mean final population sizes
  for each MPM across replicates.
  
* Function `plot.adaptInv()` now processes the first variable trait in the trait
  axis as the source of axis labels by default to PIPs, and to the x axis in
  elasticity plots.
  
## USER-VISIBLE CHANGES

* The output list structure for function project3() has changed, with item
  `comm_out` renamed to `structure`, and with its own structure changing such
  that the top-level list refers to replicates, and the lower-level list refers
  to input MPMs.

## BUGS FIXED

* Projection via vrm_input with deviance terms now no longer crashes.

* Function project3() no longer always incorporates the current MPMs start
  vector as the vector at time 0, but instead incorporates it only as the
  correct start time.
  
* Corrected issue resulting in Leslie-format function-based MPMs yielding
  matrices full of NAs if the format vector is missing.
  
* Fixed fatal error in project3() occurring when the number of MPMs entered is
  smaller than the times projected using unequal stage weights across MPMs.
  
* Fixed bug that led function `summary.adaptProj()` to acknowledge only a single
  replicate.

# adapt3 1.0.2 (2025-08-27)

## NEW FEATURES

* ESS trait value optimization now handled via Golden section Nelder-Mead
  search algorithm.
  
* Standard error bootstrapping example added to invasibility analysis vignette.

## USER-VISIBLE CHANGES

* Corrections to documentation.

* Corrected default axis labels in function plot.adaptInv().

## BUGS FIXED

* Corrected optimization trait axis calculations in cases where traits span both
  positive and negative values.

* Corrected ESS trait value calculations in cases where traits span both
  positive and negative values.

* Function plot.adaptInv() now correctly labels axes with user-supplied values.

* Function summary.adaptInv() now gives the correct number of variants tested.

# adapt3 1.0.1 (2025-07-11)

## USER-VISIBLE CHANGES

* Corrections to documentation.

# adapt3 1.0.0 (2025-07-10)

## NEW FEATURES

* New community projection vignette.

## USER-VISIBLE CHANGES

* Corrected help files for certain functions.

# adapt3 0.9.9 (2025-07-08)

## USER-VISIBLE CHANGES

* Invasibility analysis vignette fully created.

## BUGS FIXED

* Function invade3() now correctly subsets converged trait optima rows.

# adapt 0.9.8 (2025-07-07)

## NEW FEATURES

* Functions plot.adaptInv() and summary.adaptInv() created to plot and
  summarize invasibility analyses, respectively.

# adapt3 0.9.7 (2025-07-05)

## NEW FEATURES

* Functions plot.adaptProj() and summary.adaptProj() created to plot and
  summarize community projections, respectively.

* ESS optimization routine enabled for function invade3().

* New invasibility analysis vignette added.

* Threshold values for 0s created for function invade3().

# adapt3 0.9.6 (2025-04-12)

## NEW FEATURES

* Function invade3() added to run invasibility analyses.

* Function plot.adaptInv() created to plot the results of pairwise invasibility
  analyses.

* Functions ta_skeletion() and trait_axis() created to develop quick trait axis
  data frames.

# adapt3 0.9.5 (2025-03-09)

## NEW FEATURES

* Full project3() function developed. First full pre-invasion version.

## USER-VISIBLE CHANGES

* Updated error messages.

# adapt3 0.9.3 (2025-02-19)

## USER-VISIBLE CHANGES

* Underlying functions adapted to newest LefkoUtils.

# adapt3 0.9.2 (2024-05-12)

## NEW FEATURES

* Function project3() split up to allow creation of main invasion functions.

# adapt3 0.9.1 (2024-05-02)

## NEW FEATURES

* Cypripedium parviflorum dataset added.

# adapt3 0.9 (2024-04-28)

## NEW FEATURES

* First development version. Core function addyn_proj() function and
  accessory function equiv_input() developed.


