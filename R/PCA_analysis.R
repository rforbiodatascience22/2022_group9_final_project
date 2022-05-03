### From the BC_data I prepare to do some PCA visualization----
### Now we use the BC_clean data to ruduce the dimensions 
### 

library(tidyverse)
library(broom)  
library(cowplot)
library(patchwork) # required to arrange plots side-by-side
library(ggthemes) # for colorblind color scale
library(dplyr)

BC_PCA <- BC_data_clean %>% 
  select(29:33) %>% # as there are too much demensions I decide to have 
  # the first 4 dimensions 
  scale() %>% 
  prcomp() ## run PCA

view(BC_PCA)
view(head(BC_PCA$x)) ## check the standard division contributions

### add the subtypes to the PCA data and visualize it.
BC_PCA_plt <- data.frame(BC_PCA$x,
                         subtype = BC_data_clean$`PAM50 mRNA`)
head(BC_PCA_plt)

### visualization
ggplot(BC_PCA_plt,aes(x = PC1,
                      y = PC2,
                      col = subtype)) + geom_point() +  
  scale_color_colorblind()

### view the rotation matrix
View(BC_PCA$rotation)

### rotation visualize 
rotation_data <- BC_PCA$rotation

### reduce  the arrorws 
rotation_data <- data.frame(BC_PCA$rotation,
                            variable = row.names(BC_PCA$rotation))

# define a pleasing arrow style
arrow_style <- arrow(
  length = unit(0.05, "inches"),
  type = "closed")

# now plot, using geom_segment() for arrows and geom_text() for labels
ggplot(rotation_data) + 
  geom_segment(aes(xend = PC1, 
                   yend = PC2), 
               x = 0,
               y = 0, 
               arrow = arrow_style) + 
  geom_text(aes(x = PC1, 
                y = PC2, 
                label = variable), 
            hjust = 0, 
            size = 3, 
            color = "red") + 
  xlim(-1., 1.25) + 
  ylim(-1., 1.) +
  coord_fixed() # fix aspect ratio to 1:1

percent <- 100*BC_PCA$sdev^2 / sum(BC_PCA$sdev^2)
percent
perc_data <- data.frame(percent = percent, PC = 1:5)
ggplot(perc_data, aes(x = PC, y = percent)) + geom_col()


### as we can see that the demension1 contains over 99% of the information 
### I suggest that we make a only one dimension plot