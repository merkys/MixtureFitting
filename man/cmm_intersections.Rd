\name{cmm_intersections}
\alias{cmm_intersections}
\docType{methods}
\title{Intersections of Two Cauchy Distributions}
\description{
    Finds intersections of two Cauchy distributions by finding roots of a quadratic equation.
}
\usage{
    cmm_intersections( p )
}
\arguments{
    \item{p}{
        parameter vector of 6 parameters. Structure of p vector is
        p = c( A1, A2, mu1, mu2, gamma1, gamma2 ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the location of \emph{i}-th component,
        gamma\emph{i} is the Cauchy scale of \emph{i}-th component.
    }
}
\value{
    A vector of x values of intersections (zero, one or two).
    Returns NaN if both distributions are identical.
}
\author{Andrius Merkys}
