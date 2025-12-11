# adapt3 1.1.0 (2025-XX-XX)

## NEW FEATURES

* Trait axes now include individual covariate values that can be used to create
  variants.
  
* Density dependence through the additive limit function has been added to
  projection and invasion functions.
  
* Equivalence frames and density frames can now be added per annual matrix.

## BUGS FIXED

* Projection via vrm_input with deviance terms now no longer crashes.

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


