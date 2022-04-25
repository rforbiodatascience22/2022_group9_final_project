# Define project functions ------------------------------------------------

change_format <- function(ID){
  string <- str_split(ID, pattern = ".\\d\\dTCGA")[[1]]
  string <- str_c("TCGA-", string)[1]
  return(string)
}
