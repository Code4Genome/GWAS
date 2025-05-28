library(ggplot2)

# Simulate or load your p-values
# Example: pvals <- read.table("results.txt", header = TRUE)$P
set.seed(123)

# Calculate expected and observed -log10(P)
n <- length(pvals)
expected <- -log10(ppoints(n))
observed <- -log10(sort(pvals))

# Create a data frame
qq_data <- data.frame(Expected = expected, Observed = observed)

# Plot
ggplot(qq_data, aes(x = Expected, y = Observed)) +
  geom_point(size = 1, color = "black") +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Q-Q Plot of P-values",
       x = "Expected -log10(P)",
       y = "Observed -log10(P)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )
