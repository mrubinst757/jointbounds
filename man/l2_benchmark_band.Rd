\name{l2_benchmark_band}
\alias{l2_benchmark_band}
\title{Simultaneous multiplier intervals for calibration benchmarks}
\usage{
l2_benchmark_band(
  benchmarks, estimator = "aipw", B = 1000L, conf_level = 0.95,
  multiplier = c("rademacher", "normal"), seed = 1L
)
}
\arguments{
  \item{benchmarks}{Output from \code{l2_aipw_benchmarks()} with
  \code{variance = TRUE}.}
  \item{estimator}{The AIPW benchmark estimator.}
  \item{B}{Number of multiplier draws.}
  \item{conf_level}{Simultaneous confidence level.}
  \item{multiplier}{Rademacher or standard-normal multipliers.}
  \item{seed}{Random seed.}
}
\value{Benchmark estimates with simultaneous intervals, influence scores, the
critical value, and multiplier draws.}
\description{Forms a standardized supremum multiplier interval jointly over
all paired covariate-omission AIPW benchmarks.}
