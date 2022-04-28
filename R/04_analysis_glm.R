# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
<<<<<<< HEAD
proteomes_clean_aug <- read_csv(file = "data/03_proteomes_clean_aug.csv")
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
=======
BC_data_clean_aug <- read_csv(file = "data/03_BC_data_clean_aug.csv")
>>>>>>> 9f24866f1c91c060df444e040a654955095a0794


# Wrangle data ------------------------------------------------------------
# Add subtype to Gene Expression data
<<<<<<< HEAD
Expresion_Subtype <- left_join(proteomes_clean_aug,
                                    BC_data_clean_aug,
                                    by = c("TCGA ID" = "Complete TCGA ID")) %>%
  column_to_rownames(var = "TCGA ID") %>% 
  select(c(1:12553), "Subtype")
=======
Expresion_Subtype <- BC_data_clean_aug %>%
  column_to_rownames(var = "Complete TCGA ID") %>% 
  select(-c(1:29)) %>% 
  relocate(Subtype)
>>>>>>> 9f24866f1c91c060df444e040a654955095a0794


# Create long table
Expresion_Subtype_long <- Expresion_Subtype %>% 
  pivot_longer(cols = -Subtype,
               names_to = "protein",
               values_to = "log2_expr_level")

# Create nested tibble
Expresion_Subtype_long_nested <- Expresion_Subtype_long %>% 
  group_by(protein) %>% 
  nest() %>% 
  ungroup()


<<<<<<< HEAD
# Model data
#my_data_clean_aug %>% ...
=======
# Model data ---------------------------------------------------------------
# Build glm model
Expresion_Subtype_long_nested <- Expresion_Subtype_long_nested %>% 
  mutate(mdl = map(data, ~glm(Subtype ~ log2_expr_level,
                              data = .x)))
>>>>>>> 9f24866f1c91c060df444e040a654955095a0794

# Clean up model objects
Expresion_Subtype_long_nested <- Expresion_Subtype_long_nested %>% 
  mutate(tidy = map(mdl, broom::tidy, conf.int = TRUE)) %>%
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>%  # Remove the (Intercept)-rows
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",    # Add an indicator variable
                                   p.value >= 0.05 ~ "insignificant")) 

<<<<<<< HEAD
# Visualise data ----------------------------------------------------------
#my_data_clean_aug %>% ...
=======
>>>>>>> 9f24866f1c91c060df444e040a654955095a0794

# Visualise data -----------------------------------------------------------
Expresion_Subtype_long_nested %>% 
  filter(identified_as == "significant") %>% 
  ggplot(aes(x = estimate,
             y = fct_reorder(protein, desc(estimate)),
             colour = identified_as)) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     height = 0.2)) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom") +
  labs(y = "") 
### We need to think of another way of visualization

<<<<<<< HEAD
# Write data --------------------------------------------------------------
#write_tsv(...)
#ggsave(...)
=======
# Write data ---------------------------------------------------------------
write_csv(Expresion_Subtype_long_nested, file = "data/04_glm_result.csv")
ggsave(...)
>>>>>>> 9f24866f1c91c060df444e040a654955095a0794
