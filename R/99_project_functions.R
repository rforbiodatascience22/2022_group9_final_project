# Define project functions ------------------------------------------------

change_format <- function(ID){
  string <- str_split(ID, pattern = ".\\d\\dTCGA")[[1]]
  string <- str_c("ATCG-", string)[1]
  return(string)
}


count_na_func <- function(x) sum(is.na(x))