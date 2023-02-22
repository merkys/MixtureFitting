\name{bic}
\alias{bic}
\docType{methods}
\title{Bayesian Information Criterion (BIC)}
\description{
    Calculates Bayesian Information Criterion (BIC) for any type of mixture model.
    Log-likelihood function has to be provided.
}
\usage{
    bic( x, p, llf )
}
\arguments{
    \item{x}{data vector}
    \item{p}{vector of mixture model parameters}
    \item{llf}{function calculating log-likelihood, called as llf( x, p )}
}
\value{
    Bayesian Information Criterion value.
}
\author{Andrius Merkys}
