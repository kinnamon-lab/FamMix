#include <TMB.hpp>
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(phi);
  DATA_IVECTOR(f_pr_idxs);
  DATA_IVECTOR(f_sizes);
  DATA_MATRIX(f_pops);
  // Parameters
  PARAMETER_VECTOR(betas);
  PARAMETER_VECTOR(h2_a);
  PARAMETER_VECTOR(sigma);
  // Calculate population mean
  vector<Type> X_beta = X * betas;
  // Calculate population mean residuals and declare vector to hold Pearson-type
  // population mean residuals
  vector<Type> r_pop = y - X_beta;
  // Calculate population-specific variance parameters for each family
  vector<Type> h2_a_fam = f_pops * h2_a;
  vector<Type> sigma_fam = f_pops * sigma;
  // Loglikelihood for independent families
  Type nll = Type(0);
  int fidx = 0;
  for (int i = 0; i < f_sizes.size(); i++) {
    // Get zero-based index of family's proband in r_pop
    int f_pr_idx = f_pr_idxs(i) - 1;
    // Additive polygenic effect covariance matrix for family
    matrix<Type> D_fam =
      Type(2.0) * h2_a_fam(i) * sigma_fam(i) * sigma_fam(i) *
      phi.block(fidx, fidx, f_sizes(i), f_sizes(i));
    // Environmental effect covariance matrix for family
    matrix<Type> R_fam(f_sizes(i), f_sizes(i));
    R_fam
      .setZero()
      .diagonal()
      .setConstant((Type(1.0) - h2_a_fam(i)) * sigma_fam(i) * sigma_fam(i));
    // Covariance matrix for family
    matrix<Type> Sigma_fam = D_fam + R_fam;
    // Add family contribution to loglikelihood. MVNORM returns -log density for
    // multivariate normal
    nll += MVNORM(Sigma_fam)(r_pop.segment(fidx, f_sizes(i))) +
      dnorm(r_pop(f_pr_idx), Type(0.0), sigma_fam(i), 1L);
    // Update family index counter
    fidx += f_sizes(i);
  }
  // REPORT fitted values and population residuals to R for checks
  REPORT(X_beta);
  REPORT(r_pop);
  // Return negative loglikehood
  return(nll);
}
