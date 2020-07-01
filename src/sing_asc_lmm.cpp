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
  PARAMETER_VECTOR(sigma2);
  // Calculate population mean residuals
  vector<Type> r_pop = y - X * betas;
  // Declare vector of residuals from distribution conditional on proband
  vector<Type> r_cond(y.size() - f_pr_idxs.size());
  // Declare vector of Cholesky-scaled residuals from distribution conditional
  // on proband
  vector<Type> r_cholesky(y.size() - f_pr_idxs.size());
  // Calculate population-specific variance parameters for each family
  vector<Type> h2_a_fam = f_pops * h2_a;
  vector<Type> sigma2_fam = f_pops * sigma2;
  // Loglikelihood for independent families
  Type ll = Type(0);
  int fidx = 0;
  for (int i = 0; i < f_sizes.size(); i++) {
    // Additive polygenic effect covariance matrix for family
    matrix<Type> D_fam = Type(2.0) * h2_a_fam(i) * sigma2_fam(i) *
      phi.block(fidx, fidx, f_sizes(i), f_sizes(i));
    // Environmental effect covariance matrix for family
    matrix<Type> R_fam(f_sizes(i), f_sizes(i));
    R_fam.setZero().diagonal()
      .setConstant((Type(1.0) - h2_a_fam(i)) * sigma2_fam(i));
    // Covariance matrix for family
    matrix<Type> Sigma_fam = D_fam + R_fam;
    // Create proband extractor matrix to create conditional distribution
    // matrices. Note that 1 needs to be subtracted from 1-based proband index
    // passed from R to obtain 0-based C++ index
    int f_pr_idx = f_pr_idxs(i) - 1;
    int f_pr_idx_fam = f_pr_idx - fidx;
    int f_size_non_pr = f_sizes(i) - 1;
    matrix<Type> pr_extractor(f_size_non_pr, f_sizes(i));
    pr_extractor.setZero();
    pr_extractor.topLeftCorner(f_pr_idx_fam, f_pr_idx_fam).setIdentity();
    pr_extractor
      .bottomRightCorner(f_size_non_pr - f_pr_idx_fam,
                         f_size_non_pr - f_pr_idx_fam)
      .setIdentity();
    // Intermediate declarations are used to help provide type clarity for
    // automatic differentiation
    vector<Type> r_pop_fam = r_pop.segment(fidx, f_sizes(i));
    vector<Type> Sigma_fam_pr_col = Sigma_fam.col(f_pr_idx_fam);
    vector<Type> r_pop_fam_no_pr = pr_extractor * r_pop_fam;
    vector<Type> Sigma_fam_pr_col_no_pr = pr_extractor * Sigma_fam_pr_col;
    matrix<Type> Sigma_fam_no_pr = pr_extractor * Sigma_fam *
      pr_extractor.transpose();
    matrix<Type> Sigma_fam_no_pr_cor = Sigma_fam_pr_col_no_pr.matrix() *
      Sigma_fam_pr_col_no_pr.matrix().transpose() /
      Sigma_fam(f_pr_idx_fam, f_pr_idx_fam);
    vector<Type> r_cond_fam = r_pop_fam_no_pr -
      Sigma_fam_pr_col_no_pr * r_pop(f_pr_idx) /
      Sigma_fam(f_pr_idx_fam, f_pr_idx_fam);
    matrix<Type> Sigma_cond_fam = Sigma_fam_no_pr - Sigma_fam_no_pr_cor;
    matrix<Type> Sigma_cond_fam_chol = Sigma_cond_fam.llt().matrixL();
    // MVNORM returns -log density for multivariate normal conditional on the
    // proband
    ll -= MVNORM(Sigma_cond_fam)(r_cond_fam);
    // Update residuals
    r_cond.segment(fidx - i, f_size_non_pr) = r_cond_fam;
    r_cholesky.segment(fidx - i, f_size_non_pr) =
      atomic::matinv(Sigma_cond_fam_chol) * r_cond_fam;
    // Update family index counter
    fidx += f_sizes(i);
  }
  // Return residuals to R
  REPORT(r_pop);
  REPORT(r_cond);
  REPORT(r_cholesky);
  // Return negative loglikehood
  return(-ll);
}
