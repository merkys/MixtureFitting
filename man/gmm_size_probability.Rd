\name{gmm_size_probability}
\alias{gmm_size_probability}
\docType{methods}
\title{The Gaussian Mixture Distribution}
\description{
    Calculates the posterior probability of a Gaussian mixture with \emph{n} components.
    Internally, it attempts to maximize log-likelihood of data by calling optim() and returns the list as received from optim().
}
\usage{
    gmm_size_probability( x, n, method = "SANN" )
}
\arguments{
    \item{x}{data vector}
    \item{n}{number of mixture components}
    \item{method}{optimization method passed to optim()}
}
\value{
    List representing the converged optim() run.
}
\author{Andrius Merkys}
