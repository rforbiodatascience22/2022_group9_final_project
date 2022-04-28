# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
proteomes_clean_aug <- read_csv(file = "data/03_proteomes_clean_aug.csv")
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
BC_data_clean_aug <- read_csv(file = "data/03_BC_data_clean_aug.csv")


# Wrangle data ------------------------------------------------------------

# Wrangling data so that it can be used for glm
BC_data_clean_aug_test <- BC_data_clean_aug %>%
  select(c(Luminal_A,
           Luminal_B,
           HER2_enriched,
           Basal_like,
           starts_with("NP"))) %>% 
  pivot_longer(cols = -c(Luminal_A,
                         Luminal_B,
                         HER2_enriched,
                         Basal_like),
               names_to = "proteom",
               values_to = "expr_level") %>% 
  group_by(proteom) %>% 
  nest() %>% 
  ungroup()

# Making four different logistic models for each of the subtypes
# Luminal A --------------------------------------------------------------------
BC_data_clean_aug_LuminalA = BC_data_clean_aug_test %>%
  mutate(mdl_LuminalA = map(data, ~glm(Luminal_A ~ expr_level,
                              data = .,
                              family = binomial(link = "logit"))),
         mdl_LuminalA_tidy = map(mdl_LuminalA, tidy, conf.int = TRUE)) %>% 
  unnest(mdl_LuminalA_tidy) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(identified_as_LuminalA = case_when(p.value < 0.05 ~ "significant",    # Add an indicator variable
                                   p.value >= 0.05 ~ "insignificant"),
         identified_as_LuminalA = as.factor(identified_as_LuminalA))
  

# Density plot of Luminal A
BC_data_clean_aug_LuminalA %>%
  ggplot(aes(estimate,
             fill = identified_as_LuminalA)) +
  geom_density(alpha = 0.6) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom")
# Missing labels etc.


# Luminal B --------------------------------------------------------------------
BC_data_clean_aug_LuminalB = BC_data_clean_aug_test %>%
  mutate(mdl_LuminalB = map(data, ~glm(Luminal_B ~ expr_level,
                                       data = .x,
                                       family = binomial(link = "logit"))))

# Adding model information
BC_data_clean_aug_LuminalB = BC_data_clean_aug_LuminalB %>%
  mutate(mdl_LuminalB_tidy = map(mdl_LuminalB, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_LuminalB_tidy) %>% 
  filter(term != "(Intercept)") %>%  
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",  
                                   p.value >= 0.05 ~ "insignificant"))

# HER2 --------------------------------------------------------------------
BC_data_clean_aug_HER2 = BC_data_clean_aug_test %>%
  mutate(mdl_HER2 = map(data, ~glm(HER2_enriched ~ expr_level,
                                       data = .x,
                                       family = binomial(link = "logit"))))

# Adding model information
BC_data_clean_aug_HER2 = BC_data_clean_aug_HER2 %>%
  mutate(mdl_HER2_tidy = map(mdl_HER2, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_HER2_tidy) %>% 
  filter(term != "(Intercept)") %>%  
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",  
                                   p.value >= 0.05 ~ "insignificant"))

# Basal-like -------------------------------------------------------------------
BC_data_clean_aug_basal = BC_data_clean_aug_test %>%
  mutate(mdl_basal = map(data, ~glm(Basal_like ~ expr_level,
                                       data = .x,
                                       family = binomial(link = "logit"))))

# Adding model information
BC_data_clean_aug_basal = BC_data_clean_aug_basal %>%
  mutate(mdl_basal_tidy = map(mdl_basal, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_basal_tidy) %>% 
  filter(term != "(Intercept)") %>%  
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",  
                                   p.value >= 0.05 ~ "insignificant"))


# TRIED TO PUT EVERYTHING IN ONE, BUT THINK IT IS TOO LONG
#BC_data_clean_aug_test_model = BC_data_clean_aug_test %>%
#  mutate(mdl_LuminalA = map(data, ~glm(Luminal_A ~ expr_level,
#                                       data = .x,
#                                       family = binomial(link = "logit"))),
#         mdl_LuminalB = map(data, ~glm(Luminal_B ~ expr_level,
#                                       data = .x,
#                                       family = binomial(link = "logit"))),
#         mdl_HER2 = map(data, ~glm(HER2_enriched ~ expr_level,
#                                   data = .x,
#                                   family = binomial(link = "logit"))),
#         mdl_Basal = map(data, ~glm(Basal_like ~ expr_level,
#                                    data = .x,
#                                    family = binomial(link = "logit"))))


# Adding model information
BC_data_clean_aug_test_model = BC_data_clean_aug_test_model %>%
  mutate(mdl_LuminalA_tidy = map(mdl_LuminalA, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_LuminalA_tidy)
         #mdl_LuminalB_tidy = map(mdl_LuminalB, ~tidy(.x, conf.int = TRUE)),
         #mdl_HER2_tidy = map(mdl_HER2, ~tidy(.x, conf.int = TRUE)),
         #mdl_Basal_tidy = map(mdl_Basal, ~tidy(.x, conf.int = TRUE))) %>% 
         #mdl_LuminalB_tidy,
         #mdl_HER2_tidy,
         #mdl_Basal_tidy)



gene_expr_data_long_nested = gene_expr_data_long_nested %>%
  mutate(mdl_tidy = map(mdl, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_tidy)
gene_expr_data_long_nested

#BC_data_clean_aug_test_model = BC_data_clean_aug_test %>%
#  mutate(mdl = map(data, ~glm(response ~ expr_lvl,
#                              data = subset(BC_data_clean_aug_test == 0),
#                              family = binomial(link = "logit"))))


# Add subtype to Gene Expression data
Expresion_Subtype <- BC_data_clean_aug %>%
  column_to_rownames(var = "Complete TCGA ID") %>% 
  select(-c(1:29)) %>% 
  relocate(Subtype)




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


# Model data
#my_data_clean_aug %>% ...

# Model data ---------------------------------------------------------------
# Build glm model
Expresion_Subtype_long_nested <- Expresion_Subtype_long_nested %>% 
  mutate(mdl = map(data, ~glm(Subtype ~ log2_expr_level,
                              data = .x)))

# Clean up model objects
Expresion_Subtype_long_nested <- Expresion_Subtype_long_nested %>% 
  mutate(tidy = map(mdl, broom::tidy, conf.int = TRUE)) %>%
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>%  # Remove the (Intercept)-rows
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",    # Add an indicator variable
                                   p.value >= 0.05 ~ "insignificant")) 

# Visualise data ----------------------------------------------------------
#my_data_clean_aug %>% ...

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


# Write data ---------------------------------------------------------------
write_csv(Expresion_Subtype_long_nested, file = "data/04_glm_result.csv")
ggsave(...)