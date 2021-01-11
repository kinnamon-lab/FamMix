
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FamModel

<!-- badges: start -->

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
presented in [Cowan et al.
(2018)](https://doi.org/10.1161/CIRCGEN.117.002038) with an equivalent
parameterization using the `lmna_nonseg` example data provided with the
package. Note how the package permits flexible use of `formula`
constructs and provides an appropriate likelihood ratio test for a null
narrow-sense heritability on the boundary of the parameter space.

``` r
lmna_data <- copy(lmna_nonseg)[,
  `:=`(
    female = as.integer(sex == 2),
    age_echo_std = (age_echo_yrs - mean(age_echo_yrs, na.rm = TRUE)) /
      sd(age_echo_yrs, na.rm = TRUE)
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
  lvef ~ age_echo_std + female + I(n_lmna_vars > 0) + I(n_oth_vars > 0)
)
lmna_lvef_model$print()
```

``` 

===LINEAR MIXED MODEL RESULTS===
DATA: lmna_data
MEAN MODEL: lvef ~ age_echo_std + female + I(n_lmna_vars > 0) + I(n_oth_vars > 
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
      17       17 
MESSAGE: CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
MAX ABSOLUTE ELEMENT OF LL GRADIENT (g) AT SOLUTION: 1.796447e-05
NEGATIVE LL HESSIAN (-H) CHARACTERISTICS AT SOLUTION:
   SMALLEST EIGENVALUE: 2.071124e-02
   RECIPROCAL CONDITION NUMBER: 2.204337e-03
SCALED LL GRADIENT (-g' * H^-1 * g) CRITERION AT SOLUTION: 2.967911e-09

VARIANCE PARAMETERS

Parameter Estimates
-------------------
      Estimate       SE
h2_a   0.35480  0.59386
sigma 14.23758  2.05036

Likelihood Ratio Tests
----------------------
       Ho      Max |g| -g' * H^-1 * g Min lambda(-H) 1 / kappa(-H)  LR X^2 DF Pr(> X^2)
 h2_a = 0 9.973149e-06   1.922893e-09   2.038053e-02  6.473839e-02 0.48601  1   0.24286

NOTE: P-values are calculated from a 50:50 mixture of chi-square(0)
      and chi-square(1) per Self and Liang (1987)

MEAN MODEL

Parameter Estimates
-------------------
                        Estimate        SE   95% LCL   95% UCL Z value  Pr(>|Z|)    
(Intercept)             50.64281   5.90390  39.07138  62.21425  8.5779 < 2.2e-16 ***
age_echo_std            -5.75977   3.17743 -11.98741   0.46788 -1.8127  0.069876 .  
female                   3.40821   5.29036  -6.96070  13.77712  0.6442  0.519426    
I(n_lmna_vars > 0)TRUE   3.38168   5.73935  -7.86724  14.63059  0.5892  0.555721    
I(n_oth_vars > 0)TRUE  -13.31161   5.15733 -23.41980  -3.20342 -2.5811  0.009849 ** 
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

NOTE: Wald tests and CIs are displayed in the above output
```

## Acknowledgements

Development of this software was supported by the National Heart, Lung,
and Blood Institute and National Human Genome Research Institute of the
National Institutes of Health under award numbers R01HL128857 and
R01HL149423. The content is solely the responsibility of the author and
does not necessarily represent the official views of the National
Institutes of Health.
