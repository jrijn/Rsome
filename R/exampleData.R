#' Import BioRad qPCR Data
#'
#' This function retrieves the path for the example file.
#' 
#' @keywords example
#' @export exampleData

exampleData <- function(){
  exampleData <- system.file("extdata", "exampleData.txt", package = "Rsome")
  return(exampleData)
}
