\name{gmm_intersections}
\alias{gmm_intersections}
\docType{methods}
\title{Intersections of Two Gaussian Distributions}
\description{
    Finds intersections of two Gaussian distributions by finding roots of a quadratic equation.
}
\usage{
    gmm_intersections( p )
}
\arguments{
    \item{p}{
        parameter vector of 6 parameters. Structure of p vector is
        p = c( A1, A2, mu1, mu2, sigma1, sigma2 ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the location of \emph{i}-th component,
        sigma\emph{i} is the scale of \emph{i}-th component.
    }
}
\value{
    A vector of x values of intersections (zero, one or two).
    Returns NaN if both distributions are identical.
}
\author{Andrius Merkys}
