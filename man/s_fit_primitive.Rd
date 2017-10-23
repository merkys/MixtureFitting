\name{s_fit_primitive}
\alias{s_fit_primitive}
\docType{methods}
\title{Estimate Student's t distribution parameters using Batch Approximation Algorithm.}
\description{
    Estimates parameters for univariate Student's t distribution
    parameters using Batch Approximation Algorithm, according to Fig. 2
    of Aeschliman et al. (2010).
}
\usage{
    s_fit_primitive( x )
}
\arguments{
    \item{x}{data vector}
}
\value{
    Vector \code{c( mu, k, ni )}, where
    mu is the center,
    k is the concentration and
    ni is the degrees of freedom of the distribution.
}
\references{
    Aeschliman, C.; Park, J. & Kak, A. C. A
    Novel Parameter Estimation Algorithm for the Multivariate t-Distribution
    and Its Application to Computer Vision
    European Conference on Computer Vision 2010, 2010
    \url{https://engineering.purdue.edu/RVL/Publications/Aeschliman2010ANovel.pdf}
}
\author{Andrius Merkys}
