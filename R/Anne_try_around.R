
# We fix the 3 replicates --> takes mean
# Which ones are dublicated?
table(colnames(proteomes)) > 1
proteomes["ATCG-AO-A12D"]

proteomes %>% 
  filter("ATCG-AO-A12D")

substr(colnames(proteomes)[c(11,69,75)])

colnames(proteomes)[c(13,71,77)] <-

  sub_replace(colnames(proteomes)[c(13,71,77)],"2")
  
proteomes %>% 
  select("ATCG-AO-A12D")

proteomes %>% 
  mutate("ATCG-AO-A12D" =  )
proteomes[,c(4,13)]

proteomes %>% 
  select("ATCG-AO-A12D")

proteomes[c(13,71,77)]

proteomes %>% 
  select(-c(13,71,77)) %>% 
  dim()

dim(proteomes)
 
#* "ATCG-AO-A12D" at locations 2 and 11.
#* "ATCG-C8-A131" at locations 3 and 69.
#* "ATCG-AO-A12B" at locations 4 and 75.

table(colnames(proteomes)) > 1

proteomes$`ATCG-C8-A131`
proteomes$`ATCG-AO-A12D`
proteomes$`ATCG-AO-A12B`
proteomes$`ATCG-AO-A12B`

select(table(colnames(proteomes))[,] )

colnames(proteomes)[count(table(colnames(proteomes)))[[ ]] > 1]
table(colnames(proteomes))[[2]]



proteomes$

proteomes %>% 
  select(starts_with())

proteomes %>% 
  C8-A131

proteomes %>% 
  select(4:86)

# Adds "ID_short" to join easier
patients <- patients_raw %>% 
  mutate("ID_short" = substr(patients_raw$`Complete TCGA ID`,6,nchar(patients_raw$`Complete TCGA ID`))) #ALL unique (105)
Gene_Expresion <- Gene_Expresion %>% 
  mutate("ID_short" = substr(Gene_Expresion$`TCGA ID`,0,7))                                             #80/83 are unique (TAKE MEAN)

#Which ones are there several of in proteomes?
#table(substr(Gene_Expresion$`TCGA ID`,0,7))



patients$`Complete TCGA ID`
Gene_Expresion$`TCGA ID`
%>% 
  select(-c(`ID_short`))
#view(BC_Data)

#correlation and how to take non unique into account
Gene_Expresion[,1:200] %>% 
  select(-c(`TCGA ID`)) %>%
  cor(use="complete.obs")

#dim(Gene_Expresion)

#str(Gene_Expresion)

#cor(proteomes_raw$`C8-A131.32TCGA`,proteomes_raw$`C8-A131.01TCGA`,use="complete.obs")
#transpose t() if you want to check correlation between patients --> to check the 3 not unique in proteomes




## old -----------------------------------------------------------------------

# Transposed data
Gene_Expresion <- proteomes_raw[,c(1,4:86)] %>% 
  pivot_longer(cols= -1) %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)

# Part of data we wish to merge by
substr(patients_raw$`Complete TCGA ID`,6,nchar(patients_raw$`Complete TCGA ID`)) #ALL unique (105)
gsub("TCGA-",substr(Gene_Expresion$`TCGA ID`,0,7))                              #80/83 are unique (can't do anything about it)                          

#Which ones are there several of?
table(substr(Gene_Expresion$`TCGA ID`,0,7))


# Adds "ID_short" to make merge easier
patients <- patients_raw %>% 
  mutate("ID_short" = substr(patients_raw$`Complete TCGA ID`,6,nchar(patients_raw$`Complete TCGA ID`))) 

Gene_Expresion <- Gene_Expresion %>% 
  mutate("ID_short" = substr(Gene_Expresion$`TCGA ID`,0,7)) 

#### Wants to short this down and edit Gene_Expresion to be like patients_raw

# Merge data
test <- full_join(patients,
                  Gene_Expresion,
                  by = "ID_short")
view(test)



# TRASH ------------------------------------------------------------------------


# Analyse data
dim(proteomes_raw)
view(proteomes_raw)
view(PAM50_raw)
view(patients_raw)
length(unique(proteomes_raw$RefSeq_accession_number))



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
  


