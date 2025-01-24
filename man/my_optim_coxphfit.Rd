% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_optim_coxphfit.R
\name{my_optim_coxphfit}
\alias{my_optim_coxphfit}
\title{One-step optimization for coxph model}
\usage{
my_optim_coxphfit(
  index,
  x,
  y,
  control_opt = list(strata = NULL, init = NULL, offset = NULL, eps = 1e-09, maxiter =
    20, weights = NULL, method = "efron", rownames = NULL)
)
}
\arguments{
\item{index}{variable index}

\item{x}{data matrix}

\item{y}{'Surv' class variable}

\item{control_opt}{list variable for hyper-parameters in coxph.fit. Similar to coxph.control()}
}
\description{
One-step optimization for coxph model
}
