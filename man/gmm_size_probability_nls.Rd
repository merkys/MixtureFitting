\name{gmm_size_probability_nls}
\alias{gmm_size_probability_nls}
\docType{methods}
\title{The Gaussian Mixture Distribution}
\description{
    Calculates the posterior probability of a Gaussian mixture with \emph{n} components.
    Internally, it bins the data vector and calls nls() to optimize the mixture fit.
    Returns the list of the same form as received from optim().
}
\usage{
    gmm_size_probability_nls( x, n, bins = 100, trace = FALSE )
}
\arguments{
    \item{x}{data vector}
    \item{n}{number of mixture components}
    \item{bins}{number of bins}
    \item{trace}{should debug trace be printed?}
}
\value{
    List of the same form as received from optim().
}
\author{Andrius Merkys}
