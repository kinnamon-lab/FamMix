#' Nonsegregating \emph{LMNA} family data
#'
#' A dataset containing idiopathic DCM pedigrees with nonsegregating
#' \emph{LMNA} variants as well as other potentially relevant variants along
#' with quantitative endophenotype measurements. Note that there were no twins
#' in this data set.
#'
#' @format A \code{\link[data.table]{data.table}} with 88 rows and 13 columns:
#' \describe{
#'   \item{family_ID}{Family ID}
#'   \item{proband}{1 = Yes, 0 = No}
#'   \item{individual_ID}{Individual identifier (unique within family)}
#'   \item{paternal_ID}{Individual_ID of person's father}
#'   \item{maternal_ID}{Individual_ID of person's mother}
#'   \item{sex}{1 = Male, 2 = Female}
#'   \item{n_lmna_vars}{# potentially relevant \emph{LMNA} variants}
#'   \item{n_oth_vars}{# potentially relevant non-\emph{LMNA} variants}
#'   \item{grade}{Phenotypic severity grade, 0 - 4 (see Cowan et al. 2018
#'     reference)}
#'   \item{age_echo_yrs}{Age at diagnostic echo, in years}
#'   \item{lvedd_z}{Left ventricular end-diastolic dimension z-score
#'     (standardized for height and sex)}
#'   \item{lvef}{Left ventricular ejection fraction, \%}
#'   \item{mzpair}{Same value for individuals who are MZ twins; \code{NA}
#'     otherwise}
#'   \item{dzpair}{Same value for individuals who are DZ twins; \code{NA}
#'     otherwise}
#' }
#'
#' @source \url{https://bit.ly/lmna_nonseg}
#'
#' @references
#'
#' Cowan JR, Kinnamon DD, Morales A, Salyer L, Nickerson DA, Hershberger RE.
#'     Multigenic Disease and Bilineal Inheritance in Dilated Cardiomyopathy
#'     Is Illustrated in Nonsegregating LMNA Pedigrees.
#'     \emph{Circ Genom Precis Med}. 2018 Jul;11(7):e002038.
#'     doi: 10.1161/CIRCGEN.117.002038.
"lmna_nonseg"
