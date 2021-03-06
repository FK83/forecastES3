% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/murphy_procs_Jan22.R
\name{murphy_VaRES}
\alias{murphy_VaRES}
\alias{murphy_VaRES_diff}
\title{Murphy Diagrams for VaR and Expected Shortfall}
\usage{
murphy_VaRES(TSfc, alpha = 0.025, plot_nr = 2, labels = NULL,
  grid = NULL, legend_loc = NULL, murphy_cols = NULL, cex_gen = 1.6,
  ...)

murphy_VaRES_diff(TSfc, alpha = 0.025, plot_nr = 2, grid = NULL,
  labels = NULL, legend_loc = NULL, HAC_lags = 9, conf_level = 0.95,
  cex_gen = 1.6, ...)
}
\arguments{
\item{TSfc}{Matrix containing forecasts and realizations. First two columns contain VaR/ES forecasts from method A; cols 3-4 contain VaR/ES forecasts from method B; col 5 contains realizations.}

\item{alpha}{scalar, level of VaR/ES}

\item{plot_nr}{Either 1 (first elementary score) or 2 (second elementary score)}

\item{labels}{Vector of labels (set to NULL to omit labels)}

\item{grid}{Either NULL (to use grid at all data points) or numeric vector of customized grid points.}

\item{legend_loc}{Location of legend}

\item{murphy_cols}{Colors to be used for first plot}

\item{cex_gen}{Font size}

\item{...}{other plotting parameters}

\item{HAC_lags}{Truncation lag for Newey-West variance estimator (needed for confidence bands)}

\item{conf_level}{Level of confidence bands}
}
\description{
Murphy Diagrams for VaR and Expected Shortfall
}
