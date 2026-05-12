\name{snmm_init_vector}
\alias{snmm_init_vector}
\docType{methods}
\title{Estimate Skew-Normal Mixture parameters using Expectation Maximization.}
\description{
    Estimate an initialization vector for skew-normal mixture fitting via Expectation Maximization.
    Estimation method follows Lin et al. (2007) suggestion to use k-means clustering for initial cluster assignment.
    Calculation of moments and dzeta and sigma parameters are done according to Equation 3 of Lin et al. (2007), but lambda is estimated from Equation 18c of Arnold et al. (1993).
}
\usage{
    snmm_init_vector( x, n = 1 )
}
\arguments{
    \item{x}{data vector}
    \item{n}{number of mixture components}
}
\value{
    Parameter vector of 4*\emph{n} parameters, where \emph{n} is
    number of mixture components. Structure of p vector is
    p = c( omega1, omega2, ..., omega\emph{n}, dzeta1, dzeta2, ..., dzeta\emph{n}, sigma1, sigma2, ..., sigma\emph{n}, lambda1, lambda2, ..., lambda\emph{n} ),
    where omega\emph{i} is the proportion of \emph{i}-th component,
    dzeta\emph{i} is the location of \emph{i}-th component,
    sigma\emph{i} is the scale of \emph{i}-th component and
    lambda\emph{i} is the shape of \emph{i}-th component.
}
\references{
    Arnold et al.
    The nontruncated marginal of a truncated bivariate normal distribution,
    Psychometrika 58, pages 471--488 (1993)

    Lin et al.
    Finite mixture modelling using the skew normal distribution,
    Statistica Sinica 17 (2007), 909--927
    \url{https://www3.stat.sinica.edu.tw/statistica/oldpdf/A17n35.pdf}
}
\author{Andrius Merkys}
