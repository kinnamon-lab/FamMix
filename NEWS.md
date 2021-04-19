# FamModel 1.0.0

## R code changes and improvements

### `Contrast` class
* Changed to store estimate of $L\theta - m$ rather than $L\theta$.
* Changed chi-square df calculation to use rank of $L$ rather than number of
  rows.
* Improved output from `print` method.

### `FamData` class
* Added checks to constructor to make sure that individual order and identifiers
  are identical in `data` and `phi` members.
* Added pedigree check to constructor.
* Converted all integer storage type numeric columns in input data set to double
  storage type to avoid later comparison errors with `identical()`.
* Added identification of families with consanguinity and `get_consang` method.
* Modified `get_data` method to return copy of `data` member to avoid
  inadvertent modification by reference.
* Added `print` method.
* Modified `make_model_mats_lmm` to track included observations from data member
  by `fmid` and `id` and include these identifiers as names or row names to
  increase robustness.
* Modified `lmm` method to use only non-probands in `lm()` fit used to obtain
  initial parameter estimates for optimization.

### `FamLMMFit` subclass
* Disallowed construction outside of internal `FamModel` methods and functions.
* Made sure that analytical, not numeric, Hessian was being used in covariance
  matrix calculation.
* Added data mapping checks to constructor to protect against bugs due to future
  API changes and ensure accurate diagnostic output.
* Added `get_model_res` method with additional diagnostic features.
* Added scaled gradient convergence criterion and negative Hessian diagnostics
  to `print` and `get_h2_a_lrts` methods.
* Modified `get_h2_a_lrts` to cache heritability LRT results after first call
  and use cached values on subsequent calls unless the user explicitly
  overrides.

### Other
* Moved common LMM optimization code to utility function and switched
  optimization to use `optim()` without relative function convergence criterion.
* Switched from variance to SD paramaterization for linear mixed model to
  improve scaling properties of optimization problem.
* Added rigorous argument checks to `make_coef_mat()` and `print_ests()`.
* Added title argument to `print_ests()`.
* Replaced `all.equal()` with `identical()` for all object comparisons where
  floating point imprecision should not be an issue.

## C++ template infrastructure changes
* Fixed bug in `sing_asc_lmm` model that would yield misspecified trait variance
  for inbred individuals.
* Moved diagnostic calculations, which are not necessary for optimization,
  to `FamLMMFit$get_model_res()` method to speed up iterations.
* Transitioned to shared library with same name as package to avoid issues with
  newer `pkgload` versions during development.
* Refactored to allow easy incorporation of multiple models into a single shared
  library as per https://github.com/mlysy/TMBtools.

## Documentation changes
* Added `vignette("linear_mixed_models")`.
  * Moved methodological details from individual class and method documentation
    to vignette.
  * Discussed conditions required for valid inferences when true values fall on
    boundaries of the parameter space.
  * Discussed properties of estimated residuals and diagnostics relevant to
    assessing goodness of fit.

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

* Changed `FamData` constructors to use character string column names and
  updated documentation.

# FamMix 0.1.1

## Bug Fixes

* Fixed bug in `FamData$new()` using `phi` signature.

# FamMix 0.1.0

* Initial release.
