
#### Variation Partition
# In excel I calculated percentage of SS values for each source of variation. 

# Load necessary library
library(ggplot2)

# Data
ss_values <- c(Age = 8.111, Sex = 16.32, Tissue = 41.93, 
               Age_Sex = 5.63, Age_Tissue = 8.76, 
               Tissue_Sex = 11.71, Age_Tissue_Sex = 7.53)

# Convert the named vector to a data frame for ggplot2
ss_data <- data.frame(
  Source = factor(names(ss_values), levels = names(ss_values)), # Maintain order
  SS = as.numeric(ss_values)
)

# Create the bar plot
ggplot(ss_data, aes(x = Source, y = -SS)) +
  geom_bar(stat = "identity", fill = alpha("#97bebf", 0.7), color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(SS, 2), "%")), vjust = 1.5, color = "black", size = 3.5) +
  labs(title = "Variation Partition",
       subtitle = "Percentage Contribution of Each Source of Variation",
       x = "Source of Variation",
       y = "Percentage Contribution (%)") +
  scale_y_continuous(labels = function(x) paste0(abs(x), "%")) +
  scale_x_discrete(expand = c(0.1, 0.1)) +  # Add space between bars and outer border
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )



