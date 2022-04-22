
# Analyse data
dim(proteomes_raw)
view(proteomes_raw)
view(PAM50_raw)
view(patients_raw)
length(unique(proteomes_raw$RefSeq_accession_number))


#transposed data
Gene_Expresion <- proteomes_raw[,c(1,4:86)] %>% 
  pivot_longer(cols= -1) %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value")

# TRASH ------------------------------------------------------------------------


#Edit gene names out of Proteomes:
Gene_Expresion <-proteomes_raw[,c(1,4:86)] %>% 
  column_to_rownames(var = "RefSeq_accession_number")
view(Gene_Expresion)

#transposed 
Gene_Expresion <- proteomes_raw[,c(1,4:86)] %>% 
  pivot_longer(cols= -1) %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value")



view(proteomes_raw)

Gene_Expresion %>% 
  pivot_longer(cols = -people, names_to = 'People') %>% 
  pivot_wider(names_from = people, values_from = value)

#Transpose


patients_raw %>% colnames(Gene_Expresion)(year,month)



test <- full_join(patients_raw,
                  Gene_Expresion,
                  by = c("Complete TCGA ID" = "TCGA_ID"),
                  suffix = c(".x", ".y")) %>% head()

#gather and rather to transpose



proteomes_raw$


view(Gene_Expresion)

column_to_rownames(proteomes_raw[,c(1,4:86)], var = "RefSeq_accession_number")  %>% head()

proteomes_raw$RefSeq_accession_number

has_rownames(Gene_Expresion)
t(proteomes_raw[,4:86])
view(Gene_Expresion)

length(Gene_Expresion)
length(proteomes_raw[,1])

add_rownames(Gene_Expresion) <- c(proteomes_raw[,1])

c(proteomes_raw[,1])
  
length(colnames(proteomes_raw[,4:86]))
#12,553 x 86

proteomes_raw
#12,553 x 86

proteomes_raw[]
  
colnames(proteomes_raw)

proteomes_raw$RefSeq_accession_number

patients_raw$`Complete TCGA ID`
proteomes_long$TCGA_ID

test <- full_join(patients_raw,
                  proteomes_long,
                  by = c("Complete TCGA ID" = "TCGA_ID")) %>%
  


