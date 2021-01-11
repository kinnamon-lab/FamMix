#' Nonsegregating *LMNA* family data
#'
#' A dataset containing idiopathic DCM pedigrees with nonsegregating
#' *LMNA* variants as well as other potentially relevant variants along
#' with quantitative endophenotype measurements. Note that there were no twins
#' in this data set.
#'
#' @format A [`data.table`] with 88 rows and 13 columns:
#' * __family_ID__: Family ID
#' * __proband__: 1 = Yes, 0 = No
#' * __individual_ID__: Individual identifier (unique within family)
#' * __paternal_ID__: Individual_ID of person's father
#' * __maternal_ID__: Individual_ID of person's mother
#' * __sex__: 1 = Male, 2 = Female
#' * __n_lmna_vars__: # potentially relevant \emph{LMNA} variants
#' * __n_oth_vars__: # potentially relevant non-\emph{LMNA} variants
#' * __grade__: Phenotypic severity grade, 0 - 4 (see Cowan et al. 2018
#'     reference)
#' * __age_echo_yrs__: Age at diagnostic echo, in years
#' * __lvedd_z__: Left ventricular end-diastolic dimension z-score
#'     (standardized for height and sex)
#' * __lvef__: Left ventricular ejection fraction, \%
#' * __mzpair__: Same value for individuals who are MZ twins; \code{NA}
#'     otherwise
#' * __dzpair__: Same value for individuals who are DZ twins; \code{NA}
#'     otherwise
#'
#'
#' @source <https://bit.ly/lmna_nonseg>
#'
#' @references
#'
#' Cowan JR, Kinnamon DD, Morales A, Salyer L, Nickerson DA, Hershberger RE.
#'   Multigenic Disease and Bilineal Inheritance in Dilated Cardiomyopathy
#'   Is Illustrated in Nonsegregating LMNA Pedigrees.
#'   *Circ Genom Precis Med*. 2018 Jul;11(7):e002038.
#'   <https://doi.org/10.1161/CIRCGEN.117.002038>
"lmna_nonseg"
