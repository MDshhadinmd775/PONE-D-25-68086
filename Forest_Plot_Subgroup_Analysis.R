# Load necessary libraries
library(meta)
library(readr)

# Read the updated CSV file
data <- read_csv("HbA1c_subgroup_weeks_forest.csv")  # Make sure it now has 'Duration_Group'

# Ensure 'Duration_Group' is treated as a factor
data$Duration_Group <- factor(data$Duration_Group, levels = c(">=12 weeks", "<12 weeks"))

# Run meta-analysis with subgroups
meta_result <- metacont(
  n.e = N_Treatment,
  mean.e = Mean_Treatment,
  sd.e = SD_Treatment,
  n.c = N_Control,
  mean.c = Mean_Control,
  sd.c = SD_Control,
  studlab = Study,
  data = data,
  sm = "MD",
  method.tau = "REML",
  hakn = TRUE,
  byvar = Duration_Group,        # This enables subgroup analysis
  print.byvar = TRUE,            # Show the subgroup names
  title = "Meta-analysis of HbA1c (%) by Treatment Duration"
)

# Print results
summary(meta_result)

# Subgrouped Forest Plot
forest(
  meta_result,
  xlab = "Mean Difference (HbA1c change in %)",
  leftcols = c("studlab", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c"),
  leftlabs = c("Study", "N (T)", "Mean (T)", "SD (T)", "N (C)", "Mean (C)", "SD (C)"),
  rightlabs = c("MD", "95% CI"),
  col.diamond = "blue",
  col.square = "steelblue",
  col.square.lines = "black",
  col.label.right = "black",
  col.label.left = "black",
  print.tau2 = TRUE,
  test.overall = TRUE,
  test.subgroup = TRUE  # <-- this will add p-value for subgroup difference
)
#funnel for subgroup if i want to show
# Subset data by duration group
data_12plus <- subset(data, Duration_Group == ">=12 weeks")
data_less12 <- subset(data, Duration_Group == "<12 weeks")

# Meta-analysis for each subgroup separately
meta_12plus <- metacont(
  n.e = N_Treatment,
  mean.e = Mean_Treatment,
  sd.e = SD_Treatment,
  n.c = N_Control,
  mean.c = Mean_Control,
  sd.c = SD_Control,
  data = data_12plus,
  studlab = Study,
  sm = "MD",
  method.tau = "REML",
  hakn = TRUE
)

meta_less12 <- metacont(
  n.e = N_Treatment,
  mean.e = Mean_Treatment,
  sd.e = SD_Treatment,
  n.c = N_Control,
  mean.c = Mean_Control,
  sd.c = SD_Control,
  data = data_less12,
  studlab = Study,
  sm = "MD",
  method.tau = "REML",
  hakn = TRUE
)