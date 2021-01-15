# Can reproduce univariate model results from Mendel 16.0

    
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
          17       17 
    MESSAGE: CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
    MAX ABSOLUTE ELEMENT OF LL GRADIENT (g) AT SOLUTION: 1.796447e-05
    NEGATIVE LL HESSIAN (-H) CHARACTERISTICS AT SOLUTION:
       SMALLEST EIGENVALUE: 2.0716e-02
       RECIPROCAL CONDITION NUMBER: 2.20471e-03
    SCALED LL GRADIENT (-g' * H^-1 * g) CRITERION AT SOLUTION: 2.967911e-09
    
    VARIANCE PARAMETERS
    
    Parameter Estimates
    -------------------
          Estimate       SE
    h2_a   0.35480  0.59386
    sigma 14.23758  2.05036
    
    Likelihood Ratio Tests
    ----------------------
           Ho      Max |g| -g' * H^-1 * g Min lambda(-H) 1 / kappa(-H)  LR X^2 DF
     h2_a = 0 9.973149e-06   1.922893e-09   2.038198e-02  6.474299e-02 0.48601  1
     Pr(> X^2)
       0.24286
    
    NOTE: P-values are calculated from a 50:50 mixture of chi-square(0)
          and chi-square(1) per Self and Liang (1987)
    
    MEAN MODEL
    
    Parameter Estimates
    -------------------
                            Estimate        SE   95% LCL   95% UCL Z value
    (Intercept)             50.64281   5.90390  39.07138  62.21425  8.5779
    female                   3.40821   5.29036  -6.96070  13.77712  0.6442
    age_echo_std            -5.69541   3.14193 -11.85347   0.46265 -1.8127
    I(n_lmna_vars > 0)TRUE   3.38168   5.73935  -7.86724  14.63059  0.5892
    I(n_oth_vars > 0)TRUE  -13.31161   5.15733 -23.41980  -3.20342 -2.5811
                            Pr(>|Z|)    
    (Intercept)            < 2.2e-16 ***
    female                  0.519426    
    age_echo_std            0.069876 .  
    I(n_lmna_vars > 0)TRUE  0.555721    
    I(n_oth_vars > 0)TRUE   0.009849 ** 
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    
    NOTE: Wald tests and CIs are displayed in the above output

---

        fmid id r_star_hat_2
     1:    L  2  0.730043390
     2:    N  5  0.049937665
     3:    N 15  1.742088016
     4:    N 19  0.003407704
     5:    N 21  0.292806987
     6:    N 22  0.200560829
     7:    N 23  1.002810794
     8:    N 25  2.437753831
     9:    N 26  1.125532038
    10:    N 27  4.406180801
    11:    N 28  0.006121703
    12:    N 29  0.033222414
    13:    N 30  2.589934099
    14:    N 31  1.313347943
    15:    N 32  0.028774127
    16:    N 33  0.125852263
    17:    N 34  0.003365751
    18:    N 35  1.338057278
    19:    N 36  0.091152927
    20:    O  1  0.683717864
    21:    O  9  0.133905801
    22:    O 10  1.529241332
    23:    P 20  1.126301426
    24:    P 22  0.083385519
    25:    P 23  5.625341970
    26:    S  5  0.088331867
    27:    S  6  0.624835540
    28:    S  7  0.139160922
    29:    S  8  0.560964832
    30:    S 10  1.761385831
    31:    S 11  1.299945913
        fmid id r_star_hat_2

---

       fmid c_star_hat c_star_hat_df p_c_star_hat
    1:    L  0.7300434             1   0.39286912
    2:    N 16.7379992            18   0.54118244
    3:    O  2.1896029             3   0.53399918
    4:    P  6.8623872             3   0.07641573
    5:    S  4.4798565             6   0.61202782

---

    
    ===LINEAR MIXED MODEL RESULTS===
    DATA: lmna_data
    MEAN MODEL: lvedd_z ~ female + age_echo_std + I(n_lmna_vars > 0) + I(n_oth_vars > 
        0)
    VARIANCE PARAMETER GROUPS: ~1
    
    FAMILIES USED: 4
    SUBJECTS USED: 31
    PROBANDS: 4
    FAMILY SIZE DISTRIBUTION:
     3 4 6 18
             
     1 1 1  1
    CONVERGENCE ACHIEVED AT -2 LL = 102.1452
    EVALUATIONS:
    function gradient 
          18       18 
    MESSAGE: CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
    MAX ABSOLUTE ELEMENT OF LL GRADIENT (g) AT SOLUTION: 2.386636e-04
    NEGATIVE LL HESSIAN (-H) CHARACTERISTICS AT SOLUTION:
       SMALLEST EIGENVALUE: 9.318933e-01
       RECIPROCAL CONDITION NUMBER: 3.982013e-02
    SCALED LL GRADIENT (-g' * H^-1 * g) CRITERION AT SOLUTION: 2.856608e-08
    
    VARIANCE PARAMETERS
    
    Parameter Estimates
    -------------------
          Estimate      SE
    h2_a   0.13429 0.52085
    sigma  1.60994 0.22302
    
    Likelihood Ratio Tests
    ----------------------
           Ho      Max |g| -g' * H^-1 * g Min lambda(-H) 1 / kappa(-H)   LR X^2 DF
     h2_a = 0 1.346326e-04   1.159326e-09   1.190739e+00  5.623017e-02 0.074573  1
     Pr(> X^2)
        0.3924
    
    NOTE: P-values are calculated from a 50:50 mixture of chi-square(0)
          and chi-square(1) per Self and Liang (1987)
    
    MEAN MODEL
    
    Parameter Estimates
    -------------------
                           Estimate       SE  95% LCL  95% UCL Z value Pr(>|Z|)   
    (Intercept)             2.24162  0.80966  0.65472  3.82851  2.7686  0.00563 **
    female                 -0.29486  0.67298 -1.61389  1.02417 -0.4381  0.66129   
    age_echo_std            0.67347  0.45334 -0.21506  1.56199  1.4856  0.13739   
    I(n_lmna_vars > 0)TRUE -0.82794  0.65748 -2.11658  0.46070 -1.2593  0.20794   
    I(n_oth_vars > 0)TRUE   1.17017  0.68393 -0.17030  2.51064  1.7110  0.08709 . 
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    
    NOTE: Wald tests and CIs are displayed in the above output

---

        fmid id r_star_hat_2
     1:    N  5  0.007399985
     2:    N 15  4.933747014
     3:    N 19  0.579627627
     4:    N 21  1.955515388
     5:    N 22  1.316808203
     6:    N 23  0.855979016
     7:    N 25  0.842785037
     8:    N 26  0.545983577
     9:    N 28  0.348419829
    10:    N 29  0.196420835
    11:    N 30  0.407111335
    12:    N 31  0.172471369
    13:    N 32  0.004259401
    14:    N 33  0.067660699
    15:    N 34  0.027169011
    16:    N 35  0.073200680
    17:    N 36  0.419129517
    18:    O  9  3.958837335
    19:    O 10  0.293656364
    20:    P 20  0.560162025
    21:    P 22  1.035805896
    22:    P 23  6.279545269
    23:    S  5  0.462379576
    24:    S  7  1.252586816
    25:    S  8  0.299198265
    26:    S 10  0.068550003
    27:    S 11  0.181230595
        fmid id r_star_hat_2

---

       fmid c_star_hat c_star_hat_df p_c_star_hat
    1:    N  12.200361            17   0.78785561
    2:    O   4.247641             2   0.11957394
    3:    P   8.078991             3   0.04440709
    4:    S   2.473123             5   0.78053719

