#' Import BioRad qPCR Data
#'
#' This function imports a .txt file from CFX manager and returns a dataframe.
#' 
#' First simply use the file name (without .txt) as input, assuming subfolder /input as source.
#' The Cq is set to 40 in those wells were the threshold was not reached.
#' 
#' @param tablename: This is the filename of the .txt file containing Cq values as copied from the file manager. 
#' @param dropunnecessary: Defaults to TURE. If TRUE, returns only the columns "Well", "Target", "Sample", "Cq".
#' @keywords qPCR, import
#' @export cqimport
#' @examples
#' #The .txt files should be in: "parent_directory/input/example_file_name.txt"
#' 
#' setwd("parent_directory")
#' filename <- "example_file_name"
#' df <- cqimport(filename)

cqimport <- function(tablename, dropunnecessary = TRUE){
  
  #First compile the .txt file name, assuming subfolder /input as source.
  #Then import the table with headings, assuming decimal separator ",".
 
  if(grepl("exampleData", tablename) == FALSE){
    filename <- paste("input/", tablename, ".txt", sep = "")
  }
  
  df <- read.table(
    filename, 
    header=TRUE, 
    na.strings=c("", "NA"), 
    sep = "\t", 
    check.names = TRUE, 
    dec = ","
  )
  
  #If the column "Cq" does not contain numeric values, import again with decimal separator "."
  
  cqclass <- sapply(df, class)
  
  if(cqclass["Cq"] != "numeric") {
    
    df <- read.table(
      filename, 
      header=TRUE, 
      na.strings=c("", "NA"), 
      sep = "\t", 
      check.names = TRUE, 
      dec = "."
    )
    
  }
  
  #Now clean up the unnecessary columns and drop empty wells.
  
  if(dropunnecessary == TRUE){
    
    import::from(dplyr, select, filter)
    import::from(magrittr, "%>%")
    
    df <- df %>%
      select(Well, Target, Sample, Cq) %>%
      filter(Target != "NaN")
    
  }
  
  #Finally, set the Cq to 40 in those wells were the threshold was not reached.
  
  df$Cq[is.nan(df$Cq)] <- 40
  
  #Return the resulting dataframe!
  
  return(df)
}
