#' Data: Activity
#'
#' Artificial data. Approximately normally distributed covariate measuring activity of
#' 1000 individuals.
"activity"

#' Data: Age category
#'
#' Artificial data. Five-level factor giving age category of 1000 individuals
"agecategory"

#' Data: Sex
#'
#' Artificial data. Two-level factor giving sex of 1000 individuals
"sex"

#' Data: Genetic markers
#'
#' Artificial data. Matrix containing markers (0, 1, 2, corresponding to genotypes aa, aA, AA)
#' at 2000 markers (columns) for 1000 individuals (rows).
#'
#' The frequencies of 0, 1 and 2 are approximately 0.81, 0.18 and 0.01, respectively, at each
#' marker. Sample correlation beween neighbouring markers is approximately 0.6.
"xg"

#' Data: Normal response
#'
#' Artificial data. Normally distributed response for 1000 individuals, with some effect of
#' covariates \code{activity}, \code{agecategory} and \code{sex}, and of the first genetic
#' marker (first column of \code{xg}).
"y_normal"

#' Data: Bernoulli response
#'
#' Artificial data. Bernoulli response for 1000 individual, with effect of the first genetic
#' marker (first column of \code{xg}).
"y_logistic"
