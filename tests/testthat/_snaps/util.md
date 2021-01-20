# Printed output is correct

    Code
      print_ests(theta_hat, V_theta_hat)
    Output
      
      Parameter Estimates
      -------------------
                Estimate       SE  95% LCL  95% UCL Z value Pr(>|Z|)  
      n_1_0.9    1.28155  1.00000 -0.67841  3.24152  1.2816     0.20  
      n_2_0.95   3.28971  2.00000 -0.63022  7.20964  1.6449     0.10  
      n_3_0.975  5.87989  3.00000  0.00000 11.75978  1.9600     0.05 .
      n_4_0.995 10.30332  4.00000  2.46346 18.14317  2.5758     0.01 *
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      print_ests(theta_hat, V_theta_hat, trans = exp)
    Output
      
      Parameter Estimates
      -------------------
      Transformation: exp 
                Estimate (tr)         SE 95% LCL (tr) 95% UCL (tr) Z value Pr(>|Z|)  
      n_1_0.9      3.6022e+00 1.0000e+00   5.0742e-01   2.5572e+01  1.2816     0.20  
      n_2_0.95     2.6835e+01 2.0000e+00   5.3247e-01   1.3524e+03  1.6449     0.10  
      n_3_0.975    3.5777e+02 3.0000e+00   1.0000e+00   1.2800e+05  1.9600     0.05 .
      n_4_0.995    2.9831e+04 4.0000e+00   1.1745e+01   7.5767e+07  2.5758     0.01 *
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

