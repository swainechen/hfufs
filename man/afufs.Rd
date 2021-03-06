% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afufs.R
\name{afufs}
\alias{afufs}
\title{Calculate Fu's Fs for arbitrary parameters - asymptotic approximation}
\usage{
afufs(n, k, theta)
}
\arguments{
\item{n}{The number of total sequences/individuals}

\item{k}{The number of unique alleles}

\item{theta}{The average pairwise nucleotide divergence}
}
\value{
Fu's Fs
}
\description{
Returns Fu's Fs statistic.
}
\details{
Fu's Fs is a population genetics statistic that is useful for detecting
loci that are responsible for population expansion, for example. Fu's Fs
can be formulated as a calculation that involves Stirling numbers of the
first kind. These can get large very quickly and exceed the floating point
range for modern genomic data sets.

`afufs` calculates Fu's Fs using a single term asymptotic approximation
This contrasts with `hfufs` which uses a Stirling number approximator for
each term. `afufs` is therefore simultaneously more accurate and faster.
}
\examples{
n <- 100
k <- 30
theta <- 12.345
afufs(n, k, theta)
# -0.7368616

}
\references{
Chen, S.L. and Temme, N. (2020) Computing sums of Stirling numbers of the
first kind: Application to population genetics statistics. In preparation.
}
