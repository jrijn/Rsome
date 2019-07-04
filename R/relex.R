#' Calculate The 2^dCt
#'
#' This function calculates the relative expression from a cqimport object. 
#' Furthermore, it creates separate .csv files of the relative expression summary, 
#' the samples which don't meet the set criteria, and the -RT values.
#' 
#' @param tablename: a cqimport object
#' @param household: defaults to "Gapdh". Accepts a string to define the household gene.
#' @param SDcutoff: Defines the cutoff value for the standard deviation of the mean Cq. Default = 1.
#' @param Cqcutoff: Defines the cutoff value for the mean Cq. Default = 35.
#' @keywords qPCR, relative expression
#' @export relex
#' @examples
#' 
#' df <- cqimport(filename)
#' rel.expression <- relex(df, household = "Gapdh")
#' 

relex <- function(tablename, household = "Gapdh", SDcutoff = 1, Cqcutoff = 35){
  
  library(dplyr)
  import::from(tidyr, separate, gather)
  import::from(magrittr, "%>%")
  
  meanCq <- tablename %>%
    group_by(Target, Sample) %>%	
    summarise(meanCq = mean(Cq), SD = sd(Cq), n = n())
  
  rejectCq <- meanCq %>%
    filter(SD > SDcutoff | meanCq > Cqcutoff)

  warning <- rejectCq %>%
    filter(!grepl('-RT', Sample))

  RT <- meanCq %>%
    filter(grepl('-RT', Sample))

  if (length(warning$Sample) > 0){
    print("Some samples do not meet SD or cycle number criteria:")
    print(warning)

  }
  
  #Check if the folder ./output exists. If not, create it.
  
  if(dir.exists('./output') == FALSE){
    dir.create('./output')
  }
  
  write.csv(rejectCq, file="output/rejectCq.csv")
  write.csv(RT, file="output/RTcontrols.csv")
  
  samples	<- unique(meanCq$Sample)
  dCq <- data.frame()
  
  for (i in samples){
    samp <- filter(meanCq, Sample == i)
    hh <- filter(samp, Target == household)[,which(colnames(samp)=="meanCq")]
    dSamp <- samp %>%
      mutate(dCq = meanCq - unlist(hh)) %>%
      mutate(rel.expression = 2^-dCq)
    
    dCq <- bind_rows(dCq, dSamp)
    dCq <- dCq %>%
      filter(!grepl('-RT', Sample))
  }
  
  write.csv(dCq, file = "output/dCq_results.csv")
  
  return(dCq)
}