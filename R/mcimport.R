#' Import BioRad qPCR Data
#'
#' This function imports a meltcurve derivative .txt file from CFX manager and returns a dataframe.
#' 
#' First simply use the file name (without .txt) as input, assuming subfolder /input as source.
#' 
#' @param cqimport: A cqimport object. 
#' @param meltderivative: This is the file containing the meltcurve derivative data.
#' @keywords qPCR, import
#' @export mcimport
#' @examples
#' #The .txt files should be in: "parent_directory/input/example_file_name.txt"
#' 
#' setwd("parent_directory")
#' filename <- "example_file_name"
#' meltderivative <- "example_file_name"
#' df <- cqimport(filename)
#' mc <- mcimport(cqimport = df, meltderivative = meltderivative)

mcimport <- function(cqimport = df, meltderivative = meltderivative){
 
  #First load the required packages
  
  library(tidyr)
  library(dplyr)
  
  # Meltcurve analysis
  ## Import data
  
  filename <- paste("input/", meltderivative, ".txt", sep = "")
  
  mc <- read.table(
    filename, 
    header = TRUE, 
    dec = "."
  )
  
  ## Format data
  
  dat <- mc %>% 
    gather("Well", "Fluorescence", -Temperature) %>%
    separate("Well", c("Y", "X"), remove = FALSE, sep = 1) %>%
    mutate(X = as.numeric(X))
  
  melt <- cqimport %>% 
    separate("Well", c("Y", "X"), remove = FALSE, sep = 1) %>%
    dplyr::select(X, Y, Target, Sample) %>%
    mutate(X = as.numeric(X)) %>%
    left_join(dat, ., by=c("X", "Y"))
  
  melt <- melt %>%
    filter(Sample != "NA")
  
  return(melt)
}
