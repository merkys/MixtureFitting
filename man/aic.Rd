\name{aic}
\alias{aic}
\docType{methods}
\title{Akaike Information Criterion (AIC)}
\description{
    Calculates Bayesian Information Criterion (AIC) for any type of mixture model.
    Log-likelihood function has to be provided.
}
\usage{
    aic( x, p, llf )
}
\arguments{
    \item{x}{data vector}
    \item{p}{vector of mixture model parameters}
    \item{llf}{function calculating log-likelihood, called as llf( x, p )}
}
\value{
    Akaike Information Criterion value.
}
\author{Andrius Merkys}
