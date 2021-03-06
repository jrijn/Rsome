% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effplot.R
\name{eff.plot}
\alias{eff.plot}
\title{Plot concentration curve and linear regression.}
\usage{
eff.plot(tablename, logConcentration = "logConcentration",
  groupColumn = "Target", meanCq = "meanCq", plotWidth = 10,
  plotHeight = 10)
}
\arguments{
\item{tablename:}{a relex object with an added column of log2([cDNA]) called logConcentration}

\item{groupColumn:}{defaults to "Target". Defines the grouping column which will define coloring.}

\item{logConcentration:}{defaults to "logConcentration". Defines the column containing log2([cDNA]) values.}

\item{meanCq:}{defaults to "meanCq". Defines the column containing meanCq values.}

\item{plotWidth:}{width of output pdf file}

\item{plotHeight:}{height of output pfd file}
}
\description{
This function plots the untransformed Cq values for each target against the log2([cDNA]).
Then it performs linear regression and plots the R^2 and y ~ x formula.
}
\examples{

cq <- Rsome::cqimport(cqfile)
mc <- Rsome::mcimport(
  cqimport = cq, 
  meltderivative = meltfile)
  
re <- Rsome::relex(
  cq, 
  household = "Gapdh",
  SDcutoff = 1,
  Cqcutoff = 35)

library(tidyr)
library(dplyr)

eff <- separate(re, "Sample", c("Exp", "Condition", "Concentration"), sep = "_")
eff$Concentration <- as.numeric(eff$Concentration)
eff <- eff \%>\%
  filter(Concentration > 0) \%>\%
  mutate(logConcentration = log2(Concentration))

p1 <- eff.plot(
  eff,
  logConcentration = "logConcentration",
  groupColumn = "Target",
  meanCq = "meanCq")
  
}
\keyword{expression,}
\keyword{plot}
\keyword{qPCR,}
\keyword{relative}
