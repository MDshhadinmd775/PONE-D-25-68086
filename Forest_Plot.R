# Load necessary libraries
library(meta)
library(readr)

# Read your CSV file (adjust file path as needed)
data <- read_csv("FPG_forest.csv")  # <-- Replace with your actual file name

# Perform meta-analysis using metacont for change-from-baseline data
meta_result <- metacont(
  n.e = N_Treatment,
  mean.e = Mean_Treatment,
  sd.e = SD_Treatment,
  n.c = N_Control,
  mean.c = Mean_Control,
  sd.c = SD_Control,
  studlab = Study,
  data = data,
  sm = "MD",             # Mean Difference
  method.smd = "Hedges", # Not used here, but in case you switch to SMD
  method.tau = "REML",   # Random-effects model
  hakn = TRUE,
  title = "Meta-analysis of FBG (Change from Baseline)"
)

# Print summary results
summary(meta_result)

# Create forest plot
forest(
  meta_result,
  # sortvar = TE,              # <-- REMOVE this line to preserve CSV order
  xlab = "Mean Difference (FBG change in mmol/L)",
  leftcols = c("studlab", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c"),
  leftlabs = c("Study", "N (T)", "Mean (T)", "SD (T)", "N (C)", "Mean (C)", "SD (C)"),
  rightlabs = c("MD", "95% CI"),
  col.diamond = "blue",
  col.square = "steelblue",       # Updated color as per earlier request
  col.square.lines = "black",
  col.label.right = "black",
  col.label.left = "black",
  print.tau2 = TRUE,
  test.overall = TRUE,
  test.subgroup = FALSE
)