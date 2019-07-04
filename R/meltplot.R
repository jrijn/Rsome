#' Plot Relative Expression Values
#'
#' This function plots the meltcurves for each target. 
#' 
#' @param tablename: an mcimport object
#' @param groupColumn: defaults to "Sample". Defines the grouping column which will define coloring.
#' @param plotWidth: width of output pdf file
#' @param plotHeight: height of output pfd file
#' @keywords qPCR, relative expression, plot
#' @export mc.plot
#' @examples
#' 
#' setwd("parent_directory")
#' filename <- "example_file_name"
#' meltderivative <- "example_file_name"
#' 
#' df <- Rsome::cqimport(filename)
#' p1 <- Rsome::cq.plot(df)
#' mc <- Rsome::mcimport(cqimport = df, meltderivative = meltderivative)
#' p2 <- Rsome::mc.plot(mc)
#' 

mc.plot <- function(tablename, groupColumn = "Sample", plotWidth = 10, plotHeight = 10){
  
  library(ggplot2)
  library(dplyr)
  
  #Plot meltcurves
  p1 <- ggplot(tablename, aes(x=Temperature, y=Fluorescence))+
    #scale_y_continuous(limits=c(0,NA))+
    geom_line(aes(colour=!!sym(groupColumn)))+
    facet_wrap(~ Target, scales = "free")+
    theme_classic()
  #geom_hline(yintercept=0, size=.5)
  
  return(p1)
}
