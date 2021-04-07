
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FamModel

<!-- badges: start -->

[![R-CMD-check](https://github.com/kinnamon-lab/FamModel/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/kinnamon-lab/FamModel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This package implements regression models for family-based genetic
studies, including mixed models. While such models can be fit with other
existing software, this package tries to provide a unified solution that
leverages the flexibility of the R environment. It also leverages modern
automatic differentiation techniques to fit likelihood-based models
quickly.

## Installation

You can install the latest released version of FamModel from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kinnamon-lab/FamModel", ref = "[tag]")
```

where `[tag]` is the tag of the most recent release.

## Example

The following code fits a univariate version of the linear mixed model
presented in [Cowan et
al. (2018)](https://doi.org/10.1161/CIRCGEN.117.002038) with an
equivalent parameterization using the `lmna_nonseg` example data
provided with the package. Note how the package permits flexible use of
`formula` constructs and provides an appropriate likelihood ratio test
for a null narrow-sense heritability on the boundary of the parameter
space.

``` r
lmna_data <- copy(lmna_nonseg)[,
  `:=`(
    female = as.integer(sex == 2),
    # Use N rather than N-1 for SD divisor in standardization (like Mendel 16.0)
    age_echo_std = (age_echo_yrs - mean(age_echo_yrs, na.rm = TRUE)) /
      (
        sd(age_echo_yrs, na.rm = TRUE) *
          sqrt((sum(!is.na(age_echo_yrs)) - 1) / sum(!is.na(age_echo_yrs)))
      )
  )
]
lmna_fd <- FamData$new(
  lmna_data,
  family_id = "family_ID",
  indiv_id = "individual_ID",
  proband = "proband",
  sex = "sex",
  maternal_id = "maternal_ID",
  paternal_id = "paternal_ID",
  mzgrp = "mzpair",
  dzgrp = "dzpair"
)
lmna_lvef_model <- lmna_fd$lmm(
  lvef ~ female + age_echo_std + I(n_lmna_vars > 0) + I(n_oth_vars > 0)
)
lmna_lvef_model$print()
```


    ===LINEAR MIXED MODEL RESULTS===
    DATA: lmna_data
    MEAN MODEL: lvef ~ female + age_echo_std + I(n_lmna_vars > 0) + I(n_oth_vars > 
        0)
    VARIANCE PARAMETER GROUPS: ~1

    FAMILIES USED: 5
    SUBJECTS USED: 36
    PROBANDS: 5
    FAMILY SIZE DISTRIBUTION:
     2 4 7 19
             
     1 2 1  1
    CONVERGENCE ACHIEVED AT -2 LL = 251.2582
    EVALUATIONS:
    function gradient 
          22       22 
    MESSAGE: CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL
    MAX ABSOLUTE ELEMENT OF LL GRADIENT (g) AT SOLUTION: 2.111349e-06
    NEGATIVE LL HESSIAN (-H) CHARACTERISTICS AT SOLUTION:
       SMALLEST EIGENVALUE: 2.071621e-02
       RECIPROCAL CONDITION NUMBER: 2.204742e-03
    SCALED LL GRADIENT (-g' * H^-1 * g) CRITERION AT SOLUTION: 5.399669e-12

    VARIANCE PARAMETERS

    Parameter Estimates
    -------------------
          Estimate       SE
    h2_a   0.35479  0.59385
    sigma 14.23753  2.05032

    Likelihood Ratio Tests
    ----------------------
           Ho      Max |g| -g' * H^-1 * g Min lambda(-H) 1 / kappa(-H)  LR X^2 Pr(> X^2)
     h2_a = 0 1.788213e-08   1.015742e-15   2.038201e-02  6.474284e-02 0.48601   0.24286

    MEAN MODEL

    Parameter Estimates
    -------------------
                           Estimate       SE  95% LCL  95% UCL Z value  Pr(>|Z|)    
    (Intercept)             50.6429   5.9038  39.0717  62.2142  8.5780 < 2.2e-16 ***
    female                   3.4080   5.2904  -6.9609  13.7770  0.6442  0.519446    
    age_echo_std            -5.6955   3.1419 -11.8535   0.4626 -1.8127  0.069873 .  
    I(n_lmna_vars > 0)TRUE   3.3815   5.7393  -7.8674  14.6303  0.5892  0.555747    
    I(n_oth_vars > 0)TRUE  -13.3117   5.1573 -23.4199  -3.2035 -2.5811  0.009848 ** 
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    NOTE: Wald tests and CIs are displayed in the above output

## Acknowledgements

Development of this software was supported by the National Heart, Lung,
and Blood Institute and National Human Genome Research Institute of the
National Institutes of Health under award numbers R01HL128857 and
R01HL149423. The content is solely the responsibility of the author and
does not necessarily represent the official views of the National
Institutes of Health.
