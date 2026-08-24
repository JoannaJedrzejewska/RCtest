# =============================================================================
# Automated check for Table 1 and Table 2
# Output:
#   reproducibility/output/table1_table2_consistency_check.csv
#   Console warning() if the check fails.
# =============================================================================

output_dir <- file.path("reproducibility", "output")

table1_path <- file.path(output_dir, "table1_AR1_raw.csv")
table2_path <- file.path(output_dir, "table2_AR1_raw.csv")

if (!file.exists(table1_path) || !file.exists(table2_path)) {
  stop(
    "Run table1_AR1.R and table2_AR1.R first -- missing: ",
    paste(c(table1_path, table2_path)[!file.exists(c(table1_path, table2_path))],
          collapse = ", ")
  )
}

table1_raw <- read.csv(table1_path, stringsAsFactors = FALSE)
table2_raw <- read.csv(table2_path, stringsAsFactors = FALSE)

alpha <- 0.05

# -----------------------------------------------------------------------------
# Step 1: summarize Table 2's implied direction and strength of evidence.
# -----------------------------------------------------------------------------

p_value_col <- grep("^P_Value", names(table2_raw), value = TRUE)[1]

if (is.na(p_value_col)) stop("Could not locate a P_Value column in table2_AR1_raw.csv")

frac_positive_mean_diff <- mean(table2_raw$Mean_Loss_Diff > 0, na.rm = TRUE)
frac_significant <- mean(table2_raw[[p_value_col]] <= alpha, na.rm = TRUE)
frac_significant_and_positive <- mean(
  table2_raw$Mean_Loss_Diff > 0 & table2_raw[[p_value_col]] <= alpha,
  na.rm = TRUE
)

cat("\n=== TABLE 2 SUMMARY (per-model DM results) ===\n")
cat(sprintf("Models with positive Mean_Loss_Diff (competitor beats AR_1): %.1f%%\n",
            100 * frac_positive_mean_diff))
cat(sprintf("Models significant at alpha=%.2f:                            %.1f%%\n",
            alpha, 100 * frac_significant))
cat(sprintf("Models BOTH positive AND significant:                        %.1f%%\n",
            100 * frac_significant_and_positive))

# -----------------------------------------------------------------------------
# Step 2: extract Table 1's joint WRC/SPA p-values for the MSE loss measure
# (the measure Table 2 is computed on).
# -----------------------------------------------------------------------------

wrc_row <- table1_raw[grepl("^WRC MSE$", table1_raw$Test), ]
spa_row <- table1_raw[grepl("^SPA MSE consistent$", table1_raw$Test), ]

if (nrow(wrc_row) == 0 || nrow(spa_row) == 0) {
  stop("Could not locate 'WRC MSE' / 'SPA MSE consistent' rows in table1_AR1_raw.csv")
}

wrc_p <- wrc_row$P_value[1]
spa_p <- spa_row$P_value[1]

cat("\n=== TABLE 1 SUMMARY (joint tests, MSE) ===\n")
cat(sprintf("WRC p-value: %.4f\n", wrc_p))
cat(sprintf("SPA p-value (consistent): %.4f\n", spa_p))

# -----------------------------------------------------------------------------
# Step 3: the actual consistency check.
#
# Logical rule encoded from the review: if a large majority of models are
# BOTH positive (competitor beats benchmark) AND individually significant,
# then WRC/SPA -- which test the null that the BEST competitor is no better
# than the benchmark -- cannot legitimately return p approx 1. A joint test
# failing to reject while most/all individual comparisons reject is the
# specific contradiction Comment 1 describes. We do not expect WRC/SPA and
# DM p-values to be numerically equal (multiple-testing correction makes
# WRC/SPA more conservative), but the DIRECTION of the conclusion must agree
# once the majority of per-model evidence is one-sided and strong.
# -----------------------------------------------------------------------------

MAJORITY_THRESHOLD <- 0.5
WRC_SPA_SUSPICIOUS_THRESHOLD <- 0.99  # p-values this close to 1 are the
# specific pattern flagged in Comment 1

contradiction_detected <- (frac_significant_and_positive > MAJORITY_THRESHOLD) &&
  (wrc_p >= WRC_SPA_SUSPICIOUS_THRESHOLD || spa_p >= WRC_SPA_SUSPICIOUS_THRESHOLD)

consistency_summary <- data.frame(
  Frac_Models_Positive_And_Significant = round(frac_significant_and_positive, 4),
  WRC_MSE_P_Value = round(wrc_p, 4),
  SPA_MSE_P_Value_Consistent = round(spa_p, 4),
  Contradiction_Detected = contradiction_detected,
  Check_Description = paste0(
    "FAIL if >", MAJORITY_THRESHOLD * 100,
    "% of models are individually significant AND superior, ",
    "while WRC/SPA p-value >= ", WRC_SPA_SUSPICIOUS_THRESHOLD
  ),
  stringsAsFactors = FALSE
)

cat("\n=== CONSISTENCY CHECK RESULT ===\n")
print(consistency_summary, row.names = FALSE)

write.csv(
  consistency_summary,
  file.path(output_dir, "table1_table2_consistency_check.csv"),
  row.names = FALSE
)

if (contradiction_detected) {
  warning(
    "CONTRADICTION STILL PRESENT: ", round(100 * frac_significant_and_positive, 1),
    "% of models are individually significant and superior to the benchmark, ",
    "but WRC/SPA report p >= ", WRC_SPA_SUSPICIOUS_THRESHOLD, ". ",
    "This is exactly the pattern described in Reviewer #1's Major Comment 1. ",
    "Do NOT treat Table 1/Table 2 as resolved until this check passes."
  )
} else {
  message(
    "No contradiction detected under the encoded rule: joint and per-model ",
    "results are directionally consistent."
  )
}