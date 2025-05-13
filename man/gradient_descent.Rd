\name{gradient_descent}
\alias{gradient_descent}
\docType{methods}
\title{Gradient Descent}
\description{
    Simple implementation of gradient descent method.
    Given a derivative function, it follows its decrease until convergence criterion is met.
}
\usage{
    gradient_descent( gradfn, start, gamma = 0.1, ..., epsilon = 0.01 )
}
\arguments{
    \item{gradfn}{derivative function}
    \item{start}{starting value}
    \item{gamma}{learning rate}
    \item{...}{additional arguments passed to derivative function}
    \item{epsilon}{convergence threshold for absolute squared difference}
}
\value{log-likelihood}
\author{Andrius Merkys}
