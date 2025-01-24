% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_optim_logit.R
\name{my_optim_logit}
\alias{my_optim_logit}
\title{One-step optimization for logistic regression model}
\usage{
my_optim_logit(index, x, y, n, d, n_cand = 10)
}
\arguments{
\item{index}{variable index}

\item{x}{data matrix}

\item{y}{binary response variable (zero-one)}

\item{n}{sample size}

\item{d}{number of random selection dimension}

\item{n_cand}{number of candidate beta coefficient for warm start}
}
\description{
One-step optimization for logistic regression model
}
