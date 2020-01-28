#' Plot Relative Expression Values
#'
#' This function plots the relative expression values for each target normalized to the household gene. 
#' 
#' @param tablename: a relex object
#' @param groupColumn: defaults to "Sample". Defines the grouping column which will define coloring.
#' @param plotWidth: width of output pdf file
#' @param plotHeight: height of output pfd file
#' @keywords qPCR, relative expression, plot
#' @export re.plot
#' @examples
#' 
#' df <- cqimport(filename)
#' re <- relex(df, household = "Gapdh")
#' p1 <- re.plot(re)
#' 

re.plot <- function(tablename, groupColumn = "Sample", plotWidth = 10, plotHeight = 10){
  
  library(ggplot2)
  import::from(dplyr, sym)
  
  if(dir.exists('./output') == FALSE){
    dir.create('./output')
  }
  
  p1 <- ggplot(tablename, aes(x=Sample, y=rel.expression), group = !!sym(groupColumn))+
    scale_y_continuous(limits=c(0,NA))+
    geom_bar(
      aes(fill=!!sym(groupColumn)), 
      stat="identity", 
      position="dodge")+
    facet_wrap(
      ~ Target, 
      scales = "free")+
    theme_classic()+
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1))+
    geom_hline(
      yintercept=0,
      size=.5)
  
  return(p1)
}
