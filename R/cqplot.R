#' Plot Relative Expression Values
#'
#' This function plots the untransformed Cq values for each target. 
#' 
#' @param tablename: a cqimport object
#' @param groupColumn: defaults to "Sample". Defines the grouping column which will define coloring.
#' @param plotWidth: width of output pdf file
#' @param plotHeight: height of output pfd file
#' @keywords qPCR, relative expression, plot
#' @export cq.plot
#' @examples
#' 
#' df <- cqimport(filename)
#' p1 <- cq.plot(df)
#' 

cq.plot <- function(tablename, groupColumn = "Sample", plotWidth = 10, plotHeight = 10){
  
  #Load required packages
  
  library(ggplot2)
  library(dplyr)
  
  #Remove -RT values from the dataframe
  
  samples <- df %>%
    filter(!grepl('-RT', Sample))
  
  #Then plot the Cq values, faceting is set to the target gene.
  
  p1 <- ggplot(samples, aes(x=Sample, y=Cq))+
    scale_y_continuous(limits=c(12.5,NA))+
    geom_dotplot(
      shape = 23,
      binaxis = "y",
      stackdir = "center")+
    facet_wrap(
      ~ Target)+
    theme_classic()+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(
      yintercept=30,
      color="red",
      size=.5)
  
  return(p1)
}
