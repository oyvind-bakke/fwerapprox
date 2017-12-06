#' The fwerapprox package
#'
#' (1) Computes score test statistics for testing whether each of a large number of
#' coefficients (typically corresponding to genetic markers) in a GLM is zero in
#' presence of a smaller number of covariates (typically environmental covariates),
#' and provides estimates of correlations between the test
#' statistics, and (2) computes Glaz–Johnson-type intersection approximations of the
#' asymptotic multivariate normal distribution of the test statistics. See Halle et al.,
#' "Efficient and powerful familywise error control in genome-wide association studies
#' using generalised linear models".
#'
#' Use \code{library(help = "fwerapprox")} to get a list of functions and data sets.
#'
#' @section Functions:
#' \describe{
#' \item{\code{scorestatcorr}}{Computes score test statistics for testing whether each of a large number of
#' coefficients (typically corresponding to genetic markers) in a GLM with canonical
#' link is zero
#' in presence of a smaller number of covariates (typically environmental
#' covariates), and provides estimates of correlations
#' between the test statistics.}
#'
#' \item{\code{gamma2}}{Computes a second-order Glaz--Johnson approximation to a multivariate
#' standard normal probability with a given correlation matrix. Gives a
#' familywise error rate level bound in multiple testing for a given local
#' (per-hypothesis) significance level.}
#'
#' \item{\code{gamma_k}}{Computes an order \eqn{k} Glaz--Johnson approximation to a multivariate
#' standard normal probability with a given correlation matrix. Gives a familywise error rate level
#' bound in multiple testing for a given local (per-hypothesis) significance level.}
#' }
#' Use \code{?scorestatcorr}, \code{?gamma2}, \code{?gamma_k} for more information.
#' @docType package
#' @import stats
#' @name fwerapprox
NULL
