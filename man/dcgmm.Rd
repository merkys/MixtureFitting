\name{dcgmm}
\alias{dcgmm}
\docType{methods}
\title{Density of The Cauchy-Gaussian Distribution}
\description{
    Density function for the Cauchy-Gaussian distribution, according to
    Eqn. 2 of Swami (2000).
}
\usage{
    dcgmm( x, p )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 5*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, e1, e2, ..., e\emph{n}, mu1, mu2, ..., mu\emph{n}, gamma1, gamma2, ..., gamma\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        e\emph{i} is the proportion of Cauchy subcomponent of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        gamma\emph{i} is the Cauchy concentration of \emph{i}-th component and
        sigma\emph{i} is the Gaussian standard deviation of \emph{i}-th component.
    }
}
\value{A vector.}
\references{
    Swami, A.
    Non-Gaussian mixture models for detection and estimation in heavy-tailed noise
    2000 IEEE International Conference on Acoustics, Speech, and Signal Processing. Proceedings (Cat. No.00CH37100), 2000, 6, 3802-3805
}
\author{Andrius Merkys}
