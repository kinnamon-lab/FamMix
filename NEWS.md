# FamModel 0.2.0

## Major Changes

* Renamed package to FamModel.
* Added new classes for working with model results:
  * `FamModelFit` abstract class.
  * `FamLMMFit` subclass specialized for linear mixed model results.
  * `Contrast` class for results of general contrasts for a `FamModelFit`
      object.
* Modified linear mixed modeling framework to allow:
  * Variance component heterogeneity across families from
    different populations.
  * Return of various types of residuals.
* Added utility functions for printing coefficient matrices.
* Made various other improvements for robustness.
* Updated and expanded documentation with additional methodological details and
  references.

# FamMix 0.1.3

## Minor Changes

* Modified `FamData$lmm()` method by adding code to transform likelihood
  maximization problem by scaling parameters by square root of their Hessian
  diagonal elements at the initial estimates, which should improve conditioning
  of the Hessian and therefore numerical convergence.

# FamMix 0.1.2

## Minor Changes

* Changed `FamData` constructors to use character string column names and updated
  documentation.

# FamMix 0.1.1

## Bug Fixes

* Fixed bug in `FamData$new()` using `phi` signature.

# FamMix 0.1.0

* Initial release.
