# Define project functions ------------------------------------------------

change_format <- function(ID){
  string <- str_split(ID, pattern = ".\\d\\dTCGA")[[1]]
  string <- str_c("TCGA-", string)[1]
  return(string)
}


count_na_func <- function(x) sum(is.na(x))

subtype_glm <- function(subtype, BC_data){
  result <- BC_data %>% 
    select(subtype, starts_with("NP")) %>% 
    pivot_longer(cols = -subtype,
                 names_to = "proteome",
                 values_to = "expr_level") %>% 
    group_by(proteome) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(mdl = map(data, ~glm(as.formula(paste(subtype, " ~ expr_level", sep = "")),
                                data = .,
                                family = binomial(link = "logit"))),
           mdl_tidy = map(mdl, tidy, conf.int = TRUE)) %>% 
    unnest(mdl_tidy) %>% 
    filter(term != "(Intercept)") %>% 
    mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                     p.value >= 0.05 ~ "insignificant"),
           identified_as = as.factor(identified_as))
  return(result)
}


new_theme <- theme(legend.position = "none",
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text = element_text(colour = "black"),
                   axis.line = element_line(colour = "black"),
                   legend.key=element_blank(),
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 8),
                   text = element_text(family = "Avenir",
                                       size = 12))