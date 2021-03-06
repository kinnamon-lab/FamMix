% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FamLMMFit.R
\name{FamLMMFit}
\alias{FamLMMFit}
\title{A \code{FamModelFit} subclass specialized for linear mixed models}
\description{
A specialized subclass of the \code{\link{FamModelFit}} abstract class for
linear mixed model results. Methods of this subclass can be used to print
model results in a nice format, access model objects, perform inference on
general contrasts, and perform likelihood ratio tests on variance
parameters. Due to reliance on the internal structures of other objects,
the object generator function is not exported, and instances of the
\code{FamLMMFit} class can be initizalized \emph{only from within functions in the
\code{FamModel} namespace}.
}
\section{Super class}{
\code{\link[FamModel:FamModelFit]{FamModel::FamModelFit}} -> \code{FamLMMFit}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{FamLMMFit$new()}}
\item \href{#method-get_formula}{\code{FamLMMFit$get_formula()}}
\item \href{#method-get_model_res}{\code{FamLMMFit$get_model_res()}}
\item \href{#method-get_objfun}{\code{FamLMMFit$get_objfun()}}
\item \href{#method-get_objfun_par_names}{\code{FamLMMFit$get_objfun_par_names()}}
\item \href{#method-print}{\code{FamLMMFit$print()}}
\item \href{#method-get_h2_a_lrts}{\code{FamLMMFit$get_h2_a_lrts()}}
\item \href{#method-clone}{\code{FamLMMFit$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="FamModel" data-topic="FamModelFit" data-id="contrast">}\href{../../FamModel/html/FamModelFit.html#method-contrast}{\code{FamModel::FamModelFit$contrast()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="FamModel" data-topic="FamModelFit" data-id="get_V_theta_hat">}\href{../../FamModel/html/FamModelFit.html#method-get_V_theta_hat}{\code{FamModel::FamModelFit$get_V_theta_hat()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="FamModel" data-topic="FamModelFit" data-id="get_data">}\href{../../FamModel/html/FamModelFit.html#method-get_data}{\code{FamModel::FamModelFit$get_data()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="FamModel" data-topic="FamModelFit" data-id="get_optres">}\href{../../FamModel/html/FamModelFit.html#method-get_optres}{\code{FamModel::FamModelFit$get_optres()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="FamModel" data-topic="FamModelFit" data-id="get_theta_hat">}\href{../../FamModel/html/FamModelFit.html#method-get_theta_hat}{\code{FamModel::FamModelFit$get_theta_hat()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Constructs a new instance of this class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$new(data, formula, objfun, optres)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{The \code{\link{FamData}} object that produced this model
fit.}

\item{\code{formula}}{The \code{\link[Formula:Formula]{Formula::Formula}} object that describes
the model.}

\item{\code{objfun}}{A \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}} object that has already been run
through an optimizer. This object is expected to have a "type"
attribute of "LMM-SA".}

\item{\code{optres}}{A \code{list} from the optimizer function detailing results.
This list is expected to have an attribute "opttype" that contains a
string identifying the function used for optimization.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Due to their reliance on the internal structures of other
objects, instances of the \code{FamLMMFit} class can be initizalized \emph{only
from within functions in the \code{FamModel} namespace}.
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_formula"></a>}}
\if{latex}{\out{\hypertarget{method-get_formula}{}}}
\subsection{Method \code{get_formula()}}{
Returns the \code{\link[Formula:Formula]{Formula::Formula}} object used to construct
the instance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$get_formula()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_model_res"></a>}}
\if{latex}{\out{\hypertarget{method-get_model_res}{}}}
\subsection{Method \code{get_model_res()}}{
Returns the original data with fitted values, residuals,
and diagnostics.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$get_model_res(all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{all}}{If \code{FALSE} (default), include only observations that
contributed to the model fit. If \code{TRUE}, include all observations in
the \code{\link{FamData}} object \code{data} member. This can be useful for
generating a data set to refit the model without certain outliers.
Fitted values, residuals, and diagnostics columns will have \code{NA}
values for observations that did not contribute to the model fit.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
See \code{vignette("linear_mixed_models")} for a detailed
description of each residual and diagnostic, suggestions for its use,
and relevant references.
}

\subsection{Returns}{
A \code{\link{data.table}} containing the \code{data} member of the
\code{\link{FamData}} object used to fit this model with the following
additional columns (notation defined in
\code{vignette("linear_mixed_models")}):
\itemize{
\item \code{mu_hat}: The fitted population mean, \eqn{\hat{\mu}_{ij}}.
\item \code{eta_hat}: The mean conditional upon the data in the proband,
\eqn{\hat{\eta}_{ij}}. These are defined for family members
(\eqn{j \neq j_i}) only and are \code{NA} in the proband.
\item \code{r_c_hat}: The estimated Cholesky residuals,
\eqn{\hat{r}_{\mathrm{c}, ij}}. These are defined for family members
(\eqn{j \neq j_i}) only and are \code{NA} in the proband.
\item \code{r_s_hat}: The estimated Pearson-type conditional residuals,
\eqn{\hat{r}_{\mathrm{s}, ij}}. These are defined for family members
(\eqn{j \neq j_i}) only and are \code{NA} in the proband.
\item \code{c_star_hat}: Family-level chi-square goodness of fit statistic,
\eqn{\hat{c}^{*}_{i}}. The same value is provided for each family
member.
\item \code{c_star_hat_df}: Degrees of freedom for \eqn{\hat{c}^{*}_{i}}, which
is equal to the number of non-probands in the family. The same value
is provided for each family member.
\item \code{p_c_star_hat}: Probability of obtaining a deviate as or more extreme
than  \eqn{\hat{c}^{*}_{i}}, \eqn{\hat{p}_{\hat{c}^{*}_{i}}}. The
same value is provided for each family member.
\item \code{r_star_hat}: Individual-level goodness of fit statistic,
\eqn{\hat{r}^{*}_{ij}}. These are defined for family members
(\eqn{j \neq j_i}) only and are \code{NA} in the proband.
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_objfun"></a>}}
\if{latex}{\out{\hypertarget{method-get_objfun}{}}}
\subsection{Method \code{get_objfun()}}{
Returns \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}} object used to construct the
instance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$get_objfun()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_objfun_par_names"></a>}}
\if{latex}{\out{\hypertarget{method-get_objfun_par_names}{}}}
\subsection{Method \code{get_objfun_par_names()}}{
Returns a vector of names derived from the original formula
that map to each element of \code{objfun$par}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$get_objfun_par_names()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Prints contents of \code{FamLMMFit} objects with nice formatting.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Arguments passed on to \code{\link[=print_ests]{print_ests()}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_h2_a_lrts"></a>}}
\if{latex}{\out{\hypertarget{method-get_h2_a_lrts}{}}}
\subsection{Method \code{get_h2_a_lrts()}}{
Perform likelihood ratio test(s) for no polygenic effect(s)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$get_h2_a_lrts(print = TRUE, use_cached = TRUE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{print}}{If \code{TRUE} (default), prints the likelihood ratio test
results in a nice format.}

\item{\code{use_cached}}{If \code{TRUE} (default), use cached likelihood ratio
test results if any exist. Otherwise, new likelihood ratio tests will
be run, and these results will replace any in the cache.}

\item{\code{...}}{Additional parameters to pass to the \code{control} list for
\code{\link[=optim]{optim()}} with \code{method = "L-BGFS-B"}. Note that \code{parscale} and
\code{fnscale} cannot be modified.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Tests that each narrow-sense heritability parameter is zero
individually using a likelihood ratio test. See
\code{vignette("linear_mixed_models")} for a detailed description of the
methods used. By default, results are cached in the \code{\link{FamLMMFit}}
object after the first call to avoid redundant optimizations on
subsequent calls.
}

\subsection{Returns}{
A \code{\link{data.table}} containing likelihood ratio test results,
invisibly.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FamLMMFit$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
