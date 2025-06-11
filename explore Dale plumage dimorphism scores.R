# Plot plumage dimorphism exploration
# 
# 


library(ggplot2)
library(dplyr)


# Better histogram of log dimorphism values
library(ggplot2)

# Create a data frame for ggplot
plot_data <- data.frame(
  log_dimorphism = df$logMaleFemalePlumageDiffAbs
)

# Remove any infinite or NA values
plot_data <- plot_data[is.finite(plot_data$log_dimorphism), , drop = FALSE]

# Create better histogram
ggplot(plot_data, aes(x = log_dimorphism)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Log-Transformed Plumage Dimorphism",
    x = "Log(Male-Female Plumage Difference)",
    y = "Number of Species"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  # Add vertical lines at key points
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 0.2, y = max(table(cut(plot_data$log_dimorphism, 50))) * 0.8, 
           label = "log(1) = 0\n(no dimorphism)", hjust = 0, size = 3)

# Alternative: Base R with better labels
hist(df$logMaleFemalePlumageDiffAbs, 
     breaks = 50,
     main = "Distribution of Log-Transformed Plumage Dimorphism",
     xlab = "Log(Male-Female Plumage Difference)",
     ylab = "Number of Species",
     col = "steelblue",
     border = "black")

# Add reference line at log(1) = 0 (no dimorphism)
abline(v = 0, col = "red", lty = 2, lwd = 2)
text(0.2, max(hist(df$logMaleFemalePlumageDiffAbs, plot=FALSE)$counts) * 0.8, 
     "No dimorphism\n(log = 0)", cex = 0.8)

# For better x-axis labels, you can also create custom breaks
custom_breaks <- c(-2, -1, 0, 1, 2, 3)
custom_labels <- paste0("log(", round(exp(custom_breaks), 2), ")")

hist(df$logMaleFemalePlumageDiffAbs, 
     breaks = 50,
     main = "Distribution of Log-Transformed Plumage Dimorphism",
     xlab = "Log(Male-Female Plumage Difference)",
     ylab = "Number of Species",
     col = "steelblue",
     border = "black",
     xaxt = "n")  # Suppress default x-axis

# Add custom x-axis
axis(1, at = custom_breaks, labels = custom_labels)
abline(v = 0, col = "red", lty = 2, lwd = 2)




#### Scatter plots ----
# Create the plotting dataframe
plot_df <- data.frame(
  male = df$Male_plumage_score_Dale2015,
  female = df$Female_plumage_score_Dale2015,
  dimorphism = df$MaleFemalePlumageDiffAbs
)

# Remove NA values
plot_df <- plot_df[complete.cases(plot_df), ]

# Find the top 5 most frequent dimorphism values
top_5_dimorphism <- names(sort(table(plot_df$dimorphism), decreasing = TRUE))[1:12]
top_5_values <- as.numeric(top_5_dimorphism)

# Create color variable
plot_df$color_group <- "Other"
for(i in 1:12) {
  plot_df$color_group[plot_df$dimorphism == top_5_values[i]] <- 
    paste0("Rank ", i, ": ", round(top_5_values[i], 3))
}

# Count points at each unique coordinate
coord_counts <- plot_df %>%
  group_by(male, female) %>%
  summarise(
    count = n(),
    color_group = first(color_group),
    dimorphism = first(dimorphism),
    .groups = 'drop'
  )

# Create named vector for colors (this fixes the syntax error)
color_names <- c("Other", paste0("Rank ", 1:12, ": ", round(top_5_values, 3)))
color_values <- c("black", "red", "blue", "green", "orange", "purple", rainbow(7))
names(color_values) <- color_names

# Create the plot
ggplot(coord_counts, aes(x = female, y = male)) +
  geom_point(aes(size = count, color = color_group), alpha = 0.7) +
  scale_size_continuous(
    name = "Species\nCount",
    range = c(0.5, 3),
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
    labels = as.character(1:8),
    guide = guide_legend(override.aes = list(alpha = 1, color = "black")),
  ) +
  scale_color_manual(
    name = "Dimorphism Value",
    values = color_values
  ) +
  labs(
    title = "Plumage Dimorphism Patterns in Dale et al. Data",
    subtitle = "Point size = number of species; Colors = top 12 most frequent dimorphism values",
    x = "Female Plumage Score",
    y = "Male Plumage Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right"
  ) +
  # Add diagonal line for reference (equal dimorphism)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  coord_fixed(ratio = 1)


# Create a version showing each rank separately
plot_df_top5 <- plot_df[plot_df$color_group != "Other", ]

ggplot(plot_df_top5, aes(x = female, y = male)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~color_group, ncol = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Top 12 Most Frequent Dimorphism Values - Spatial Distribution")
