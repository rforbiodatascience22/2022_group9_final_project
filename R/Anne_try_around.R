#analysis

BC_data_clean_aug %>% 
  select(c("PAM50 mRNA",starts_with("NP"))) %>% 
  colnames()

#^the data we want to do a linear Regression on
#Next up: respond numbers (PAM50). convert into 1's and 0's
#then follow week 6.



BC_data_clean_aug %>% 
  select("PAM50 mRNA") %>% 
  mutate(respons = case_when("Basal-like", respons=1))


BC_data_clean_aug$`PAM50 mRNA`=="Basal-like"

BC_data_clean_aug$Subtype


# Modification of column names in proteomes so they are comparable with patients_clean_aug
adjusted_names <- proteomes_clean %>% 
  select(-c("RefSeq_accession_number","gene_symbol","gene_name")) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes_clean)[4:86] <- adjusted_names


# Creating new data sets with columns consisting of dublicates and too little data removed:
proteomes_clean <- proteomes_clean %>% 
  select(unique(colnames(.))) %>%                          # Removing duplicates
  select(-c("gene_symbol","gene_name")) %>%                # Removing unnecessary describtions of protein
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10) %>%                               # Removing columns consisting of more than 10% NAs
  pivot_longer(cols= -1,                                   # Transposing
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)


proteomes_clean %>%
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10)

proteomes_clean_aug %>% 
  select(unique(colnames(.))) %>%           # removing duplicates
  select(-c("gene_symbol","gene_name")) %>% # removing unnecessary describtions of protein
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10)                    # removing columns consisting of more than 10% NAs


proteomes_clean_aug$gene_name

colnames(proteomes_clean) %>% 
  filter(colnames(proteomes_clean) == colnames(proteomes_clean))

colnames(proteomes_clean)[duplicated(colnames(proteomes_clean_aug))]


# Only selecting the proteomes that have XX% of data
proteomes_clean %>%
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10)







# Creating new tibble that is a transposed and reduced version of proteomes_clean_aug
Gene_Expresion_clean_aug <- proteomes_clean_aug %>%
  select(-c(2,3,13,71,77)) %>%  #Deletes gene_symbol, gene_name and the 3 duplicates 
  pivot_longer(cols= -1,
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)



# Transposed data
Gene_Expresion <- proteomes_raw[,c(1,4:86)] %>% 
  pivot_longer(cols= -1) %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)
view(Gene_Expresion)


# Part of data we wish to merge by
substr(patients_raw$`Complete TCGA ID`,6,nchar(patients_raw$`Complete TCGA ID`)) #ALL unique (105)
gsub("TCGA-",substr(Gene_Expresion$`TCGA ID`,0,7))                              #80/83 are unique (can't do anything about it)                          

#Which ones are there several of?
table(substr(Gene_Expresion$`TCGA ID`,0,7))

