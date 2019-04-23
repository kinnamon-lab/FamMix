#include <TMB.hpp>
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR(y_pr);
  DATA_MATRIX(X_pr);
  DATA_SPARSE_MATRIX(phi);
  DATA_IVECTOR(f_sizes);
  // Parameters
  PARAMETER_VECTOR(betas);
  PARAMETER(h2_g);
  PARAMETER(sigma2);
  // Loglikelihood for families
  // Center Y so that it has a mean of zero
  vector<Type> y_cent = y - X * betas;
  Type ll = Type(0);
  int fidx = 0;
  for (int i = 0; i < f_sizes.size(); i++) {
    matrix<Type> Sigma_fam = Type(2.0) * h2_g * sigma2 *
      phi.block(fidx, fidx, f_sizes(i), f_sizes(i));
    Sigma_fam.diagonal().setConstant(sigma2);
    // MVNORM returns -log density for multivariate normal with mean zero and
    // covariance matrix Sigma_fam evaluated at y_cent for this family
    ll -= MVNORM(Sigma_fam)(y_cent.segment(fidx, f_sizes(i)));
    fidx += f_sizes(i);
  }
  // Subtract marginal loglikelihood for (iid) probands
  vector<Type> y_pr_cent = y_pr - X_pr * betas;
  ll -= dnorm(y_pr_cent, Type(0.0), sqrt(sigma2), 1).sum();
  // Return negative loglikehood
  return(-ll);
}
