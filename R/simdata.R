#' @name simdata
#'
#' @title Simulated Multi-trait Fine-Mapping Data Used in Tutorial
#'
#' @description Simulated fine-mapping data set used to illustrate
#'   mvSuSiE in the tutorial. The data set includes genotype and
#'   phenotype data for 574 samples, 1,001 genetic markers and 20
#'   traits. The traits were simulated from the mvSuSiE model with
#'   coefficients \code{simdata$B} and residual \code{simdata$par$V}.
#'   This is a simulation with three causal genetic variants at
#'   positions 255, 335 and 493; that is, these are the only genetic
#'   variants witih nonzero coefficients.
#' 
#' @docType data
#' 
#' @format \code{simdata} is a list with the following elements:
#' 
#' \describe{
#' 
#'   \item{raw$X}{The matrix of simulated genotypes.}
#' 
#'   \item{raw$Y}{The matrix of simulated traits.}
#'
#'   \item{Btrue}{The coefficients used to simulate the data.}
#'
#'   \item{par$V}{The residual covariance matrix used to simulated the data.}
#'
#'   \item{par$U}{The collection of covariance matrices specifying the
#'     mvsusie prior.}
#'
#'   \item{par$w}{The weights associated with the covariance matrices.}
#'
#'   \item{sumstats$n}{The sample size.}
#'
#'   \item{sumstats$LD}{The LD computed from \code{raw$X}.}
#'
#'   \item{sumstats$bhat}{The least-squares effect estimates from the
#'     single-marker association tests.}
#'
#'   \item{sumstats$sehat}{The standard errors of the least-squares effect
#'     estimates.}}
#' 
#' @keywords data
#'
NULL
