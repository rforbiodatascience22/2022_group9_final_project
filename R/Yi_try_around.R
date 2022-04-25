# Modeling part: How to find the significant proteins?
# Method 1: Use glm() with the nested gene expression table---------------------

# Load data
Gene_Expresion <- read_csv(file = "data/01_gene_expression.csv")

# Create long table
BC_data_long <- Gene_Expresion %>% 
  pivot_longer(cols = -Subtypes,
               names_to = "gene",
               values_to = "log2_expr_level")