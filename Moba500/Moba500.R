##### Moba500 ####
counts_df <- load_counts("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500_FilteredSampleInfo.csv",
                         remove.sgids = "sgDummy",
                         keep.bc.end = "CTGA",
                         keep.distance = 0)
result <- plot_pretran_histogram(counts_df, percentile.threshold = 0.99, filter.data = TRUE, return.list = TRUE)
print(result$plot)
counts_filtered <- result$filtered_df

result <- plot_barcode_seeding(counts_filtered, return.data = TRUE)
print(result$plot)

counts_filtered_2 <- counts_filtered %>%
  # Filter out unwanted entries
  filter(
    sgid != "sgDummy",                      # remove sgDummy
    bc_end == "CTGA",                       # keep only CTGA
    distance == 0,                          # keep only distance == 0
    (tissue == "preTran" & cell_num >= 2) |
      (tissue == "Blood" & cell_num >= 1) |
      (tissue == "Brain" & cell_num >= 10) |
      (!tissue %in% c("preTran", "Blood", "Brain") & cell_num >= 100)
  ) %>%
  # Create gene column
  mutate(
    gene = ifelse(
      grepl("^(sgnt|sgneo|sgsafe)", sgid, ignore.case = TRUE),
      "Ctrl",
      sgid
    )
  )
write.csv(counts_filtered_2, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500_counts_filtered_final.csv", row.names = FALSE)
pretran_data <- generate_pretran_data(counts_filtered_2)

plot_composition_many(counts_filtered_2, stack.by = "gene", color.by = "gene", mouse.genotype = "NSG", tissues = c("Liver", "Lung", "Brain"))

plot_jitter_ordered(counts_filtered_2,
                    group.by = "sgid",
                    mouse.genotypes = "NSG",
                    tissues = "Liver",
                    y.limits = c(100, 1e7),
                    min.colonies = 1,
                    title = "Colonies in NSG Liver",
                    max.groups = 200)

boot_moba500_bs <- calculate_metastatic_statistics(
  counts.df = counts_filtered_2,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver"),
  statistics = c("seeding", "burden"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  min.colonies.seeding = 10,
  min.colonies.burden = 10,
  debug = T,
  R = 1000
)

results <- dormancy_cutoff_by_genotype_tissue_time_v4(
  counts_filtered_2,
  tissues = c("Liver", "Lung", "Brain"),
  time.points = c("3weeks"),
  min.colonies = 1,
  bw.adjust = 2,
  split.by.tissue = FALSE
)
results$plot_list$BL6_AllTissues
results$plot_list$NSG_AllTissues

plot_dormancy_cutoffs_v2(
  counts.df = counts_filtered_2,
  cutoffs_df = results$cutoffs_df,
  min.colonies = 1,
  value.col = "log2_cell_num",
  time.points = c("3weeks"),
  tissues = c("Liver", "Lung", "Brain"),
  genotypes = unique(counts_filtered_2$mouse_genotype),
  split.by.tissue = FALSE
)
geno_cutoffs <- results$cutoffs_df %>%
  select(mouse_genotype, dormant_peak, dormancy_cutoff)
geno_cutoffs
subset_df <- counts_filtered_2 %>%
  filter(mouse_genotype == "NSG", tissue == "Liver",
         time_point == "3weeks", gene == "Ctrl")

res <- calculate_dormancy_gmm_constrained(
  data = subset_df,
  tag = "gene",
  dormant.mean.constraint = 7.78,
  valley.mean.constraint = 11.05,
  tissue = "Liver",
  time.point = "3weeks",
  plot = TRUE
)

boot_moba500_dorm <- calculate_metastatic_statistics(
  counts.df = counts_filtered_2,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver"),
  statistics = c("dormancy_cutoff"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  dormancy.cutoff = geno_cutoffs,
  debug = T,
  min.colonies.dormancy = 100,
  R = 1000
)
write.csv(boot_moba500_dorm, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500dorm.csv")

# Filter and calculate all transformations
df_filtered <- boot_moba500_dorm %>%
  filter(mouse_genotype == "NSG", tissue == "Liver") %>%
  mutate(
    pseudocount = 1,

    # Base relative dormancy (already exists in your data, but recalculating for consistency)
    relative_dormancy_percent = (dormancy_percent + pseudocount) /
      (ctrl_dormancy_percent + pseudocount),

    # Log2 transformations
    relative_dormancy_log2 = log2(relative_dormancy_percent),
    relative_dormancy_neg_log2 = -log2(relative_dormancy_percent),

    # Natural log transformations
    relative_dormancy_ln = log(relative_dormancy_percent),
    relative_dormancy_neg_ln = -log(relative_dormancy_percent),

    # Log10 transformations
    relative_dormancy_log10 = log10(relative_dormancy_percent),
    relative_dormancy_neg_log10 = -log10(relative_dormancy_percent)
  )

# Select and reshape to long format
df_long <- df_filtered %>%
  select(
    gene,  # Keep gene for reference if needed
    relative_dormancy_percent,
    relative_dormancy_log2,
    relative_dormancy_neg_log2,
    relative_dormancy_ln,
    relative_dormancy_neg_ln,
    relative_dormancy_log10,
    relative_dormancy_neg_log10
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    # Create type categories
    Type = case_when(
      Metric == "relative_dormancy_percent" ~ "Ratio",
      str_detect(Metric, "log2") ~ "Log2",
      str_detect(Metric, "ln") ~ "Natural Log",
      str_detect(Metric, "log10") ~ "Log10",
      TRUE ~ "Other"
    ),
    # Create direction categories
    Direction = case_when(
      str_detect(Metric, "neg") ~ "Negative (Escape)",
      Metric == "relative_dormancy_percent" ~ "Ratio",
      TRUE ~ "Positive (Dormancy)"
    ),
    # Combine for detailed labeling
    Type_Direction = paste(Type, "-", Direction)
  )

# Define facet order
df_long$Metric <- factor(df_long$Metric, levels = c(
  "relative_dormancy_percent",
  "relative_dormancy_log2",
  "relative_dormancy_neg_log2",
  "relative_dormancy_ln",
  "relative_dormancy_neg_ln",
  "relative_dormancy_log10",
  "relative_dormancy_neg_log10"
))

# Create cleaner labels for facets
metric_labels <- c(
  "relative_dormancy_percent" = "Ratio (Raw)",
  "relative_dormancy_log2" = "Log2 FC",
  "relative_dormancy_neg_log2" = "Neg Log2 FC (Escape)",
  "relative_dormancy_ln" = "Ln FC",
  "relative_dormancy_neg_ln" = "Neg Ln FC (Escape)",
  "relative_dormancy_log10" = "Log10 FC",
  "relative_dormancy_neg_log10" = "Neg Log10 FC (Escape)"
)

# Plot
ggplot(df_long, aes(x = Value, fill = Direction, color = Direction)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(~ Metric, scales = "free", ncol = 2,
             labeller = labeller(Metric = metric_labels)) +
  geom_vline(
    data = subset(df_long, Metric != "relative_dormancy_percent"),
    aes(xintercept = 0),
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = c(
    "Ratio" = "#94A3B8",
    "Positive (Dormancy)" = "#8B5CF6",
    "Negative (Escape)" = "#F97316"
  )) +
  scale_color_manual(values = c(
    "Ratio" = "#64748B",
    "Positive (Dormancy)" = "#7C3AED",
    "Negative (Escape)" = "#EA580C"
  )) +
  labs(
    title = "Relative Dormancy Transformations — NSG Liver",
    subtitle = "Comparison of Log2, Natural Log, and Log10 scales",
    x = "Value",
    y = "Density",
    fill = "Direction",
    color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
  )
plot_ranked_log2fc(boot_moba500_dorm,
                   metric = "relative_dormancy",
                   transform = "neg_log2",
                   tissues = "Lung",
                   fc.threshold = 0.5,
                   mouse.genotypes = c("NSG", "BL6"),
                   title = "Dormancy escape (-log2FC of relative dormancy)",
                   facet.by = "mouse_genotype")

boot_moba500_sdp <- calculate_metastatic_statistics(
  counts.df = counts_filtered_2,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver"),
  statistics = c("size_percentile", "supermet", "dormancy_cutoff", "peak_mode"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  dormancy.cutoff = 2^geno_cutoffs,
  debug = T,
  R = 1000
)

# Grid faceting (2x2 layout)
plot_ranked_log2fc(boot_moba500_sdp,
                   metric = "dormancy",
                   tissues = c("Liver"),
                   log2fc.threshold = 0.5,
                   mouse.genotypes = c("NSG", "BL6"),
                   facet.by = c("mouse_genotype"))

metrics <- c("supermet", "dormancy", "peak_mode")
tissues <- c("Liver", "Lung", "Brain")

output_dir <- "RankPlots"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (m in metrics) {
  for (t in tissues) {
    message("Processing: ", m, " - ", t)

    # Identify column name for the metric
    metric_col <- switch(m,
                         "supermet" = "supermet_percent_Log2FC",
                         "peak_mode" = "peak_mode_Log2FC",
                         "dormancy" = "dormancy_percent_Log2FC"
    )

    # Filter once before plotting
    subset_data <- boot_moba500_sdp %>%
      filter(
        tissue == t,
        mouse_genotype %in% c("NSG", "BL6"),
        !is.na(.data[[metric_col]])
      )

    # Skip if there's no usable data
    if (nrow(subset_data) == 0) {
      message("⏩ No valid data for ", m, " - ", t)
      next
    }

    # Try plotting
    p <- tryCatch({
      plot_ranked_log2fc(
        subset_data,
        metric = m,
        facet.by = c("mouse_genotype"),
        facet.scales = "fixed",
        title = paste0(t, " - ", m)
      )
    },
    error = function(e) {
      message("⚠️ Skipping ", m, " - ", t, " (error: ", e$message, ")")
      return(NULL)
    })

    # Skip if still NULL
    if (is.null(p)) next

    # Save the PDF safely
    filename <- file.path(output_dir, paste0(m, "_", t, "_rankplot.pdf"))
    suppressWarnings(
      ggsave(
        filename = filename,
        plot = p,
        width = 12,
        height = 6,
        units = "in"
      )
    )
    message("✅ Saved: ", filename)
  }
}

boot_moba500_merged <- merge(
  boot_moba500_bs,
  boot_moba500_sdp,
  by = c("sgid", "tissue", "time_point", "n_colonies", "mouse_genotype", "gene"),
  all = TRUE  # use all.x = TRUE or all.y = TRUE if you want a left or right join instead
)
write.csv(boot_moba500_merged, "moba500_all_stats.csv")

plot_met_metric_category_comparison(
  met.results = boot_moba500_merged,
  metric.col = "peak_mode_Log2FC",
  comparison.var = "mouse_genotype",
  category1 = "NSG",
  category2 = "BL6",
  filter.vars = list(tissue = "Liver", time_point = "3weeks"),
  title = "Liver: NSG vs BL6 Peak Mode"
)

# Define the metrics to compare
metrics <- c(
  "peak_mode_Log2FC",
  "supermet_percent_Log2FC",
  "dormancy_percent_Log2FC"
)

# Define metric display names for titles
metric_names <- c(
  "peak_mode_Log2FC" = "Peak Mode",
  "supermet_percent_Log2FC" = "Supermet Percent",
  "dormancy_percent_Log2FC" = "Dormancy Percent"
)

# Define tissues
tissues <- c("Liver", "Lung", "Brain")

# Create output directory if it doesn't exist
output_dir <- "metric_comparison_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through all combinations
plots_created <- 0
plots_skipped <- 0

for (tissue in tissues) {
  for (metric in metrics) {

    # Create clean filename
    metric_short <- gsub("_Log2FC", "", metric)
    filename <- paste0(output_dir, "/", tissue, "_NSG_vs_BL6_", metric_short, ".pdf")

    # Get display name for title
    metric_display <- metric_names[metric]

    # Generate the plot with error handling
    cat("Generating:", filename, "... ")

    tryCatch({
      p <- plot_met_metric_category_comparison(
        met.results = boot_moba500_merged,
        metric.col = metric,
        comparison.var = "mouse_genotype",
        category1 = "NSG",
        category2 = "BL6",
        filter.vars = list(tissue = tissue, time_point = "3weeks"),
        title = paste0(tissue, ": NSG vs BL6 ", metric_display)
      )

      # Save as 8x8 inch PDF
      ggsave(
        filename = filename,
        plot = p,
        width = 8,
        height = 8,
        units = "in",
        device = "pdf"
      )

      cat("Done\n")
      plots_created <- plots_created + 1

    }, error = function(e) {
      cat("Skipped (no data or error)\n")
      plots_skipped <- plots_skipped + 1
    })

  }
}

plot_size_percentile_rank(boot_moba500_corrected, percentile = 0.90, filter.genotype = "NSG", filter.tissue = "Liver")

# Define the combinations
genotypes <- c("NSG", "BL6")
tissues <- c("Liver", "Lung", "Brain")

# Create output directory if it doesn't exist
output_dir <- "size_percentile_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through all combinations
for (geno in genotypes) {
  for (tis in tissues) {

    # Create filename
    filename <- paste0(output_dir, "/",
                       "size_p90_", geno, "_", tis, ".pdf")

    # Try to create the plot (in case some combinations don't exist)
    tryCatch({

      # Generate plot
      p <- plot_size_percentile_rank(
        boot_moba500_merged,
        percentile = 0.90,
        filter.genotype = geno,
        filter.tissue = tis
      )

      # Save as 8x8 inch PDF
      ggsave(
        filename = filename,
        plot = p,
        width = 8,
        height = 8,
        units = "in",
        device = "pdf"
      )

      cat("Saved:", filename, "\n")

    }, error = function(e) {
      cat("Skipped:", geno, "-", tis, "- Error:", e$message, "\n")
    })
  }
}

plot_met_burden_rank(boot_moba500_bs, fdr.threshold = 0.05, log2fc.threshold = 1, filter.genotype = "NSG",
                     filter.tissue = "Lung", title = "Bootstrapped burden rank for NSG lung", min.colonies = 10)
# Define genotypes and tissues
genotypes <- c("NSG", "BL6")
tissues <- c("Liver", "Lung", "Brain")

# Loop through all combinations
for (geno in genotypes) {
  for (tissue in tissues) {
    # Create title
    title <- paste0("Bootstrapped burden rank for ", geno, " ", tissue)

    # Generate plot
    p <- plot_met_burden_rank(
      boot_moba500_bs,
      fdr.threshold = 0.05,
      log2fc.threshold = 1,
      filter.genotype = geno,
      filter.tissue = tissue,
      title = title,
      min.colonies = 10
    )

    # Create filename
    filename <- paste0("burden_rank_", tissue, "_", geno, ".pdf")

    # Save as 8x8 PDF
    pdf(filename, width = 8, height = 8)
    print(p)
    dev.off()

    cat("Generated:", filename, "\n")
  }
}

plot_met_burden_category_comparison(
  met.burden.results = boot_moba500_bs,
  comparison.var = "tissue",
  category1 = "Liver",
  category2 = "Lung",
  filter.vars = list(mouse_genotype = "BL6", time_point = "3weeks"),
  title = "BL6: Liver vs Lung"
)

plot_met_seeding_category_comparison(
  met.seeding.results = boot_moba500_bs,
  comparison.var = "tissue",
  category1 = "Lung",
  category2 = "Brain",
  filter.vars = list(mouse_genotype = "BL6", time_point = "3weeks"),
  title = "BL6: Lung vs Brain"
)

plot_met_burden_category_comparison(
  met.burden.results = boot_moba500_bs,
  comparison.var = "mouse_genotype",
  category1 = "NSG",
  category2 = "BL6",
  filter.vars = list(tissue = "Brain", time_point = "3weeks"),
  title = "Brain: NSG vs BL6"
)

# Define tissue priority order
tissues <- c("Liver", "Lung", "Brain")
genotypes <- c("NSG", "BL6")
time_point <- "3weeks"

# Generate all comparison combinations
comparisons <- list()

# 1. Same genotype, different tissues (6 comparisons)
for (geno in genotypes) {
  for (i in 1:(length(tissues) - 1)) {
    for (j in (i + 1):length(tissues)) {
      comparisons[[length(comparisons) + 1]] <- list(
        type = "same_genotype",
        genotype1 = geno,
        genotype2 = geno,
        tissue1 = tissues[i],
        tissue2 = tissues[j],
        title_prefix = geno
      )
    }
  }
}

# 2. Same tissue, different genotypes (3 comparisons)
for (tissue in tissues) {
  comparisons[[length(comparisons) + 1]] <- list(
    type = "same_tissue",
    genotype1 = genotypes[1],
    genotype2 = genotypes[2],
    tissue1 = tissue,
    tissue2 = tissue,
    title_prefix = tissue
  )
}

# Generate and save all plots
for (comp in comparisons) {
  # Determine filter variables based on comparison type
  if (comp$type == "same_genotype") {
    filter_vars_burden <- list(mouse_genotype = comp$genotype1, time_point = time_point)
    filter_vars_seeding <- list(mouse_genotype = comp$genotype1, time_point = time_point)

    # Create burden plot
    burden_title <- paste0(comp$title_prefix, ": ", comp$tissue1, " vs ", comp$tissue2)
    p_burden <- plot_met_burden_category_comparison(
      met.burden.results = boot_moba500_bs,
      comparison.var = "tissue",
      category1 = comp$tissue1,
      category2 = comp$tissue2,
      filter.vars = filter_vars_burden,
      title = burden_title
    )

    # Save burden plot
    filename_burden <- paste0(comp$genotype1, "_", comp$tissue1, "_vs_", comp$tissue2, "_burden.pdf")
    pdf(filename_burden, width = 8, height = 8)
    print(p_burden)
    dev.off()

    # Create seeding plot
    seeding_title <- paste0(comp$title_prefix, ": ", comp$tissue1, " vs ", comp$tissue2)
    p_seeding <- plot_met_seeding_category_comparison(
      met.seeding.results = boot_moba500_bs,
      comparison.var = "tissue",
      category1 = comp$tissue1,
      category2 = comp$tissue2,
      filter.vars = filter_vars_seeding,
      title = seeding_title
    )

    # Save seeding plot
    filename_seeding <- paste0(comp$genotype1, "_", comp$tissue1, "_vs_", comp$tissue2, "_seeding.pdf")
    pdf(filename_seeding, width = 8, height = 8)
    print(p_seeding)
    dev.off()

  } else {  # same_tissue
    filter_vars_burden <- list(tissue = comp$tissue1, time_point = time_point)
    filter_vars_seeding <- list(tissue = comp$tissue1, time_point = time_point)

    # Create burden plot
    burden_title <- paste0(comp$title_prefix, ": ", comp$genotype1, " vs ", comp$genotype2)
    p_burden <- plot_met_burden_category_comparison(
      met.burden.results = boot_moba500_bs,
      comparison.var = "mouse_genotype",
      category1 = comp$genotype1,
      category2 = comp$genotype2,
      filter.vars = filter_vars_burden,
      title = burden_title
    )

    # Save burden plot
    filename_burden <- paste0(comp$tissue1, "_", comp$genotype1, "_vs_", comp$genotype2, "_burden.pdf")
    pdf(filename_burden, width = 8, height = 8)
    print(p_burden)
    dev.off()

    # Create seeding plot
    seeding_title <- paste0(comp$title_prefix, ": ", comp$genotype1, " vs ", comp$genotype2)
    p_seeding <- plot_met_seeding_category_comparison(
      met.seeding.results = boot_moba500_bs,
      comparison.var = "mouse_genotype",
      category1 = comp$genotype1,
      category2 = comp$genotype2,
      filter.vars = filter_vars_seeding,
      title = seeding_title
    )

    # Save seeding plot
    filename_seeding <- paste0(comp$tissue1, "_", comp$genotype1, "_vs_", comp$genotype2, "_seeding.pdf")
    pdf(filename_seeding, width = 8, height = 8)
    print(p_seeding)
    dev.off()
  }

  cat("Generated plots for:", comp$title_prefix, "\n")
}

# Purple: Outlier metastatic suppressors
# points significant in one measurement (Log2FC > 1, FDR < 0.05) but differing by at least 1 Log2FC in the opposite direction
plot_burden_seeding_comparison(
  combined.results = boot_moba500_corrected,
  filter.vars = list(tissue = "Brain", mouse_genotype = "NSG", time_point = "3weeks"),
  title = "NSG Brain: Burden vs Seeding"
)

# Define tissues and genotypes to loop through
tissues <- c("Liver", "Lung", "Brain")
genotypes <- c("NSG", "BL6")

# Create output directory if it doesn't exist
output_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/burden_vs_stats"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through all combinations
for (tissue in tissues) {
  for (genotype in genotypes) {

    # Create informative title
    plot_title <- paste0(genotype, " ", tissue, ": Burden vs Seeding")

    # Create filename (replace spaces with underscores, make lowercase)
    filename <- paste0(output_dir, "/",
                       tolower(genotype), "_",
                       tolower(tissue),
                       "_burden_vs_seeding.pdf")

    # Generate the plot
    cat("Generating plot for", genotype, tissue, "...\n")

    tryCatch({
      p <- plot_burden_seeding_comparison(
        combined.results = boot_moba500_bs,
        filter.vars = list(
          tissue = tissue,
          mouse_genotype = genotype,
          time_point = "3weeks"
        ),
        title = plot_title
      )

      # Save as 8x8 inch PDF
      ggsave(
        filename = filename,
        plot = p,
        width = 8,
        height = 8,
        units = "in",
        device = "pdf"
      )

      cat("  Saved:", filename, "\n")

    }, error = function(e) {
      cat("  Error generating plot for", genotype, tissue, ":", e$message, "\n")
    })
  }
}

cat("\nAll plots generated successfully!\n")
cat("Plots saved in:", output_dir, "\n")

# merge new dormancy with moba500 stats
# Columns to remove
cols_to_remove <- c(
  "dormancy_percent", "expanding_percent", "dormant_mode", "expanding_mode",
  "n_dormant", "n_proliferative", "mean_dormant", "mean_proliferative",
  "cutoff_value", "constrained", "tag_type", "dormancy_percent_ratio", "dormancy_percent_Log2FC"
)

# 1. Remove those columns from boot_moba500_merged
boot_moba500_merged_clean <- boot_moba500_merged[, !(names(boot_moba500_merged) %in% cols_to_remove)]

# 2. Merge with boot_moba500_dorm
boot_moba500_merged2 <- merge(
  boot_moba500_merged_clean,
  boot_moba500_dorm,
  by = c("sgid", "tissue", "time_point", "n_colonies", "mouse_genotype", "gene"),
  all.x = TRUE
)

# 3. Optional: check structure
str(boot_moba500_merged2)
write.csv(boot_moba500_merged2, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500newdorm_allstats.csv")

plot_met_metric_category_comparison(
  met.results = boot_moba500_merged2,
  metric.col = "relative_dormancy_escape",
  comparison.var = "mouse_genotype",
  category1 = "NSG",
  category2 = "BL6",
  filter.vars = list(tissue = "Lung", time_point = "3weeks"),
  title = "Lung: NSG vs BL6 Relative Dormancy Escape",
  divergence.threshold = 1.5,
  show.diagonal.lines = FALSE,
  use.log2 = F
)

# Define parameters
x.metrics <- c("burden_value", "seeding_value", "size_rel_cell_num_p90", "peak_mode_ratio")
tissues <- c("Liver", "Lung", "Brain")
genotypes <- c("BL6", "NSG")

# Loop through all combinations and save PDFs
for (x in x.metrics) {
  for (tissue in tissues) {
    for (geno in genotypes) {
      filename <- paste0("scatter_", x, "_vs_relative_dormancy_", tissue, "_", geno, ".pdf")
      tryCatch({
        pdf(filename, width = 10, height = 10)

        p <- plot_met_metric_scatter(
          met.results = boot_moba500_merged2,
          x.metric = x,
          y.metric = "relative_dormancy_percent",
          tissue = tissue,
          mouse_genotype = geno,
          time_point = "3weeks",
          label.genes = "top",
          n.labels = 45,
          show.diagonal = FALSE,
          highlight.significance = TRUE,
          use.log2 = FALSE
        )
        print(p)  # CRITICAL: explicitly print the ggplot object

        dev.off()
        message("Saved: ", filename)
      }, error = function(e) {
        message("Skipping ", filename, " due to error: ", conditionMessage(e))
        if (length(dev.list()) > 0) {
          dev.off()
        }
      })
    }
  }
}

plot_statistic_comparison(
  combined.results = boot_moba500_corrected,
  stat.x = "seeding", stat.y = "burden",
  filter.vars = list(tissue = "Liver", mouse_genotype = "NSG", time_point = "3weeks"),
  title = "NSG Liver: Burden vs Seeding"
)

# Define the parameter combinations
stat_y <- "burden"
stat_x_options <- c("seeding", "size", "dormancy_escape", "peak_mode", "supermet_percent")
tissues <- c("Liver", "Lung", "Brain")
genotypes <- c("NSG", "BL6")
time_point <- "3weeks"

# Create output directory if it doesn't exist
output_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/burden_comparison_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through all combinations
for (genotype in genotypes) {
  for (tissue in tissues) {
    for (stat_x in stat_x_options) {

      # Create a clean filename
      filename <- paste0(
        output_dir, "/",
        stat_y, "_vs_", stat_x, "_",
        genotype, "_", tissue, "_", time_point,
        ".pdf"
      )

      # Create a title for the plot
      plot_title <- paste0(
        stat_y, " vs ", stat_x, " | ",
        genotype, " ", tissue, " (", time_point, ")"
      )

      # Generate the plot
      tryCatch({
        p <- plot_statistic_comparison(
          combined.results = boot_moba500_corrected,
          stat.x = stat_x,
          stat.y = stat_y,
          filter.vars = list(
            tissue = tissue,
            mouse_genotype = genotype,
            time_point = time_point
          ),
          title = plot_title
        )

        # Save as 10x10 inch PDF
        pdf(filename, width = 10, height = 10)
        print(p)
        dev.off()

        cat("Saved:", filename, "\n")

      }, error = function(e) {
        cat("Error with", genotype, tissue, stat_x, ":", e$message, "\n")
      })
    }
  }
}

cat("\nAll plots saved to:", output_dir, "\n")

# Compare seeding vs relative dormancy with significance highlighting
plot_met_metric_scatter(
  met.results = boot_moba500_merged2,
  x.metric = "size_rel_cell_num_p90",
  y.metric = "relative_dormancy_escape",
  tissue = c("Lung", "Liver"),
  mouse_genotype = "NSG",
  time_point = "3weeks",
  label.genes = "divergent",
  n.labels = 5,
  show.diagonal = F,
  highlight.significance = TRUE,
  use.log2 = F,
  show.regression = T,
  regression.fdr.filter = T
)

# Make radar plots: show only genes with all values and only show normalized 0-100%
plot_radar_comparison(boot_moba500_merged2, gene = "Tsc2", tissues = c("Liver"), normalize.scale = T, time.points = "3weeks", dormancy.metric = "relative_dormancy_escape")

boot_moba500_re <- calculate_metastatic_statistics(
  counts.df = counts_filtered_2,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver"),
  statistics = c("supermet", "peak_mode"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  dormancy.cutoff = geno_cutoffs,
  debug = T,
  R = 1000
)

# Step 1: Remove the incorrectly calculated columns from boot_moba500_merged2
boot_moba500_merged2_cleaned <- boot_moba500_merged2 %>%
  dplyr::select(
    -peak_mode,
    -peak_density,
    -n_super,
    -n_regular,
    -n_total,
    -supermet_percent,
    -supermet_peak_mode,
    -supermet_percent_ratio,
    -supermet_percent_Log2FC,
    -peak_mode_ratio,
    -peak_mode_Log2FC
  )

# Step 2: Join the corrected data from boot_moba500_re
boot_moba500_corrected <- boot_moba500_merged2_cleaned %>%
  dplyr::left_join(
    boot_moba500_re,
    by = c("sgid", "tissue", "time_point", "mouse_genotype")
  )
boot_moba500_corrected <- boot_moba500_corrected %>%
  dplyr::mutate(
    gene = dplyr::coalesce(gene.x, gene.y),
    n_colonies = dplyr::coalesce(n_colonies.x, n_colonies.y)
  ) %>%
  dplyr::select(-gene.x, -gene.y, -n_colonies.x, -n_colonies.y)
boot_moba500_corrected$ctrl_dormancy_escape = 1/boot_moba500_corrected$ctrl_dormancy_percent
# Make radar plots: show only genes with all values and only show normalized 0-100%
plot_radar_comparison(boot_moba500_corrected, gene = "Tsc2", tissues = c("Liver"), normalize.scale = T, time.points = "3weeks", dormancy.metric = "relative_dormancy_escape")

# Normalize data for radar plotting
radar_ready_data <- normalize_for_radar(
  boot_moba500_corrected,
  genes = c("Tsc2"),
  tissues = c("Liver"),
  time.points = "3weeks",
  genotypes = c("BL6", "NSG"),
  dormancy.metric = "relative_dormancy_escape",
  check.ctrl.mean = TRUE
)

# Now plot
plot_radar_comparison(
  radar_ready_data,
  gene = "Tsc2",
  tissues = c("Liver"),
  time.points = "3weeks",
  dormancy.metric = "relative_dormancy_escape"
)

# Define your genes
genes_to_plot <- c("Tsc1", "Pten", "Tsc2", "Gpatch8", "Mga", "Gata6", "Crebbp", "Zmiz1", "Zeb1", "Csnk1a1", "Grhpr", "Hnf4a")

# Loop through each gene
for (gene in genes_to_plot) {

  # Normalize data for this gene
  radar_ready_data <- normalize_for_radar(
    boot_moba500_corrected,
    genes = gene,
    tissues = c("Liver"),
    time.points = "3weeks",
    genotypes = c("BL6", "NSG"),
    dormancy.metric = "relative_dormancy_escape",
    check.ctrl.mean = TRUE
  )

  # Create output filename
  output_file <- paste0("radar_nsg_bl6_", gene, "_Liver_3weeks.pdf")

  # Plot and save
  plot_radar_comparison(
    radar_ready_data,
    gene = gene,
    genotype1 = "NSG",
    genotype2 = "BL6",
    comparison.mode = "genotype",
    tissues = c("Liver"),
    time.points = "3weeks",
    dormancy.metric = "relative_dormancy_escape",
    export.pdf = TRUE,
    output.file = output_file,
    width = 8,
    height = 8
  )

  cat("Saved:", output_file, "\n")
}
write.csv(boot_moba500_corrected, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500_updatedstats.csv", row.names = FALSE)

# Compare mutation stats
moba500_stats <- read.csv("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500_updatedstats.csv", header = TRUE)
sclc_alteration_freq <- readxl::read_excel("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/SCLC-gene_alteration_frequency_finalpool.xlsm")
colnames(sclc_alteration_freq) <- c("Gene_Symbol", "OQL_Line", "Num_Samples_Altered", "Proportion_Samples_Altered")
sclc_enrichment <- read.csv("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/dependencies enriched in Small Cell Lung Cancer.csv", header = T)
sclc_enrichment <- subset(sclc_enrichment, Type != "compound")

# Default binning with priority hierarchy
result <- plot_tstat_bins(
  moba_stats.df = moba500_stats,
  mutation.df = sclc_alteration_freq,
  dependency.df = sclc_enrichment,
  bin.breaks = c(0, 3, 4.5, 6, 7.5, Inf),
  label.top.n = 10
)

print(result$plot)
print(result$dataset_summary)

moba500_body_burden <- calculate_whole_body_burden(
  counts.df = counts_filtered_2,
  tag = "sgid",
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  blood.tissue = "Blood",
  control.label = "Ctrl",
  R = 1000,
  min.colonies.burden = 10,
  debug = FALSE
)

result <- compare_mutation_body_burden(
  mutation.df = sclc_alteration_freq,
  body.burden.df = moba500_body_burden,
  label.top.n = 5,
  mouse.genotype = "BL6",
  add.correlation = F,
  add.regression = F
)
result$plot

result <- compare_dependency_body_burden(
  dependency.df = sclc_enrichment,
  body.burden.df = moba500_body_burden,
  mouse.genotype = "NSG",
  burden.column = "whole_body_burden_Log2FC",
  use.facet.break = TRUE,
  label.top.n = 5,
  add.correlation = F,
  add.regression = F,
  label.neg.limit = 10
)

print(result$plot)
View(pretran_data$pretran_cell_num)
boot_moba500_newseed <- calculate_metastatic_statistics(
  counts.df = counts_filtered_2,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver"),
  statistics = c("seeding"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  min.colonies.seeding = 10,
  pretran.cutoff.percentile = 5,
  debug = T,
  R = 1000
)

plot_met_seed_rank(boot_moba500_newseed, fdr.threshold = 0.05, log2fc.threshold = 1, filter.genotype = "BL6",
                   filter.tissue = "Liver", min.colonies = 10)
# Define combinations
genotypes <- c("BL6", "NSG")
tissues <- c("Liver", "Lung", "Brain")

# Loop through each combination
for (g in genotypes) {
  for (t in tissues) {
    # Construct file name
    file_name <- paste0("seed_rank_", g, "_", t, ".pdf")

    # Create PDF device
    pdf(file_name, width = 8, height = 8)

    # Run plotting function with explicit print
    print(plot_met_seed_rank(
      boot_moba500_newseed,
      fdr.threshold = 0.05,
      log2fc.threshold = 1,
      filter.genotype = g,
      filter.tissue = t,
      min.colonies = 10
    ))

    # Close PDF device
    dev.off()

    # Optional: Print progress
    message("Saved: ", file_name)
  }
}

plot_whole_body_burden_rank(moba500_body_burden, fdr.threshold = 0.05, log2fc.threshold = 1, filter.genotype = "NSG")
# size_Log2FC_p90 burden_Log2FC seeding_log2FC seeding_value burden_value size_rel_cell_num_p90 relative_dormancy_escape relative_dormancy_neg_log2FC
plot_met_metric_scatter(
  met.results = boot_moba500_corrected,
  x.metric = "burden_value",
  y.metric = "relative_dormancy_escape",
  tissue = c("Lung"),
  mouse.genotype = "NSG",
  time.point = "3weeks",
  label.genes = "divergent",
  n.labels = 5,
  show.diagonal = F,
  highlight.significance = TRUE,
  use.log2 = T,
  show.regression = T,
  regression.fdr.filter = T
)

# Define parameter sets
x.metrics <- c("burden_value", "seeding_value", "size_rel_cell_num_p90")
tissues <- c("Liver", "Lung")
log2.options <- c(FALSE, TRUE)

# Loop to generate and save plots
for (x in x.metrics) {
  for (t in tissues) {
    for (use_log in log2.options) {
      p <- plot_met_metric_scatter(
        met.results = boot_moba500_corrected,
        x.metric = x,
        y.metric = "relative_dormancy_escape",
        tissue = t,
        mouse.genotype = "NSG",
        time.point = "3weeks",
        label.genes = "divergent",
        n.labels = 5,
        show.diagonal = FALSE,
        highlight.significance = TRUE,
        use.log2 = use_log,
        show.regression = TRUE,
        regression.fdr.filter = TRUE
      )

      # Construct a descriptive filename
      fname <- sprintf("%s_vs_relative_dormancy_escape_%s_log%s_nsg.pdf",
                       x, t, ifelse(use_log, "2T", "2F"))

      # Save 9x9 inch plot
      ggsave(filename = fname, plot = p, width = 9, height = 9, units = "in", dpi = 300)
    }
  }
}

suppressor_genes <- c("Tsc1", "Tsc2", "Pten", "Gpatch8", "Zeb1", "Mga", "Crebbp", "Gata6", "Csnk1a1")
rac1_genes <- c("Rac1", "Nckap1", "Brk1", "Arpc4", "Arpc2", "Actr3", "Actr2")

plot_met_rank(boot_moba500_corrected,
              genes.highlight.increased = suppressor_genes,
              genes.highlight.decreased = rac1_genes,
              metric = "seeding",
              #percentile = 0.90,
              filter.tissue = "Liver",
              filter.timepoint = "3weeks",
              filter.genotype = "BL6")
plot_ranked_log2fc(boot_moba500_corrected,
                   genes.highlight.increased = suppressor_genes,
                   genes.highlight.decreased = rac1_genes,
                   metric = "relative_dormancy",
                   transform = "neg_log2",
                   tissues =  "Liver",
                   time.points = "3weeks",
                   mouse.genotypes = "NSG")
test_results <- calculate_metastatic_statistics(
  counts.df = counts_filtered_2,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver"),
  statistics = c("size_percentile"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  debug = T,
  R = 1000,
  return.list = F
)
plot_comparison_bars(boot_moba500_corrected,
                     genes = rac1_genes,
                     metric = "burden",
                     value.type = c("fc"),
                     compare.by = "tissue",
                     tissue1 = "Liver",
                     tissue2 = "Lung",
                     filter.genotype = "BL6",
                     filter.timepoint = "3weeks",
                     x.limits = c(0,1))
# Supermet requires filtering out control samples

#### Part 1: Scatter plots — seeding vs peak_mode and seeding vs size ####
setwd("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500")
boot_moba500_corrected <- read.csv("Moba500_updatedstats_final.csv", header = TRUE)

## Seeding vs Peak Mode — NSG Liver 3weeks
p1 <- plot_statistic_comparison(
  combined.results = boot_moba500_corrected,
  stat.x = "seeding",
  stat.y = "peak_mode",
  filter.vars = list(tissue = "Liver", mouse_genotype = "NSG", time_point = "3weeks"),
  title = "NSG Liver: Seeding vs Peak Mode"
)
ggsave("seeding_vs_peak_mode_NSG_Liver_3weeks.pdf", p1, width = 10, height = 10)

## Seeding vs Size (P90) — NSG Liver 3weeks
p2 <- plot_statistic_comparison(
  combined.results = boot_moba500_corrected,
  stat.x = "seeding",
  stat.y = "size",
  filter.vars = list(tissue = "Liver", mouse_genotype = "NSG", time_point = "3weeks"),
  title = "NSG Liver: Seeding vs Size (P90)"
)
ggsave("seeding_vs_size_NSG_Liver_3weeks.pdf", p2, width = 10, height = 10)


#### Part 2: Multi-percentile size computation ####
# Compute size at each percentile directly from counts data,
# then correlate with seeding for each percentile.
counts_filtered_2 <- read.csv("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba500/Moba500_counts_filtered_final.csv", header = TRUE)
library(tidyverse)
library(pheatmap)

percentiles <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)

size_by_pctile <- counts_filtered_2 %>%
  filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
  group_by(gene) %>%
  filter(n() >= 10) %>%
  reframe(
    n_colonies = n(),
    percentile = percentiles * 100,
    size_value = quantile(cell_num, probs = percentiles, na.rm = TRUE)
  ) %>%
  pivot_wider(
    names_from = percentile,
    values_from = size_value,
    names_prefix = "size_p"
  )

# Compute the same for Ctrl to get relative (Log2FC) values
ctrl_sizes <- counts_filtered_2 %>%
  filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks", gene == "Ctrl") %>%
  pull(cell_num)

ctrl_pctiles <- sapply(percentiles, function(p) quantile(ctrl_sizes, probs = p, na.rm = TRUE))
names(ctrl_pctiles) <- paste0("size_p", percentiles * 100)

# Log2FC relative to Ctrl at each percentile
size_log2fc <- size_by_pctile %>%
  mutate(across(
    starts_with("size_p"),
    ~ log2(.x / ctrl_pctiles[cur_column()]),
    .names = "{.col}_Log2FC"
  ))

# --- Merge with seeding data ---
seeding_df <- boot_moba500_corrected %>%
  filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
  select(gene, seeding_log2FC)

merged <- size_log2fc %>%
  inner_join(seeding_df, by = "gene") %>%
  filter(gene != "Ctrl")


#### Part 3A: R² of seeding vs size at each percentile (summary plot) ####

r2_by_pctile <- tibble(
  percentile = percentiles * 100,
  r_squared = sapply(paste0("size_p", percentiles * 100, "_Log2FC"), function(col) {
    fit <- lm(merged[[col]] ~ merged$seeding_log2FC)
    summary(fit)$r.squared
  })
)

p3 <- ggplot(r2_by_pctile, aes(x = percentile, y = r_squared)) +
  geom_line(color = "#7C3AED", linewidth = 1) +
  geom_point(color = "#7C3AED", size = 3) +
  geom_text(aes(label = sprintf("%.2f", r_squared)), vjust = -1, size = 3.5) +
  scale_x_continuous(breaks = percentiles * 100) +
  labs(
    title = "R² of Seeding vs Size across Percentiles — NSG Liver",
    x = "Size Percentile",
    y = "R² (Seeding Log2FC ~ Size Log2FC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave("seeding_vs_size_R2_by_percentile_NSG_Liver.pdf", p3, width = 10, height = 6)


r2_both <- map_dfr(c("NSG", "BL6"), function(geno) {
  sub <- counts_filtered_2 %>%
    filter(mouse_genotype == geno, tissue == "Liver", time_point == "3weeks")

  ctrl_vals <- sub %>% filter(gene == "Ctrl") %>% pull(cell_num)
  ctrl_q <- sapply(percentiles, function(p) quantile(ctrl_vals, probs = p, na.rm = TRUE))

  size_wide <- sub %>%
    group_by(gene) %>%
    filter(n() >= 10) %>%
    reframe(
      percentile = percentiles * 100,
      size_log2fc = log2(quantile(cell_num, probs = percentiles, na.rm = TRUE) / ctrl_q)
    ) %>%
    pivot_wider(names_from = percentile, values_from = size_log2fc, names_prefix = "p")

  seed <- boot_moba500_corrected %>%
    filter(mouse_genotype == geno, tissue == "Liver", time_point == "3weeks") %>%
    select(gene, seeding_log2FC)

  m <- inner_join(size_wide, seed, by = "gene") %>% filter(gene != "Ctrl")

  tibble(
    genotype = geno,
    percentile = percentiles * 100,
    r_squared = sapply(paste0("p", percentiles * 100), function(col) {
      summary(lm(m[[col]] ~ m$seeding_log2FC))$r.squared
    })
  )
})

ggplot(r2_both, aes(x = percentile, y = r_squared, color = genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = percentiles * 100) +
  labs(
    title = "R² of Seeding vs Size across Percentiles — Liver",
    x = "Size Percentile",
    y = "R²"
  ) +
  theme_minimal(base_size = 13)
ggsave("seeding_vs_size_R2_by_percentile_Liver.pdf", width = 10, height = 6)

#### preTran correlation ####
# Merge preTran cell numbers with seeding results
pretran_df <- pretran_data$pretran_cell_num
colnames(pretran_df)  # should be sgid, gene, cell_num

seeding_nsg_liver <- boot_moba500_corrected %>%
  filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
  select(sgid, gene, seeding_log2FC, seeding_value, n_colonies)

pretran_seed <- inner_join(seeding_nsg_liver, pretran_df, by = "sgid") %>%
  filter(gene.x != "Ctrl")

r2 <- summary(lm(seeding_log2FC ~ log2(cell_num), data = pretran_seed))$r.squared

ggplot(pretran_seed, aes(x = log2(cell_num), y = seeding_log2FC)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
  ggrepel::geom_text_repel(
    data = pretran_seed %>%
      filter(abs(seeding_log2FC) > quantile(abs(seeding_log2FC), 0.95, na.rm = TRUE)),
    aes(label = sgid),
    size = 3, max.overlaps = 20
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("R² = %.2f", r2), size = 5, fontface = "bold") +
  labs(
    title = "preTran Representation vs Seeding — NSG Liver",
    x = "log2(preTran Cell Number)",
    y = "Seeding Log2FC"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("pretran_vs_seeding_NSG_Liver.pdf", width = 10, height = 10)

#### Seeding stability across colony size cutoffs — NSG Liver ####
# Define cutoffs (log-spaced to cover 1–1000 without 1000 iterations)
cutoffs <- c(1, 2, 5, 10, 25, 50, 75, 100, 150, 200, 300, 500, 750, 1000)

# Get all Liver colonies for NSG (no cell_num filter yet)
# Need counts_filtered (before the >=100 cutoff) or reload from counts_df
liver_all <- counts_filtered %>%
  filter(
    sgid != "sgDummy", bc_end == "CTGA", distance == 0,
    mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks"
  ) %>%
  mutate(gene = ifelse(grepl("^(sgnt|sgneo|sgsafe)", sgid, ignore.case = TRUE), "Ctrl", sgid))

# preTran cell numbers
pretran_df <- pretran_data$pretran_cell_num  # data.frame with sgid, gene, cell_num

# Compute seeding at each cutoff
seeding_by_cutoff <- map_dfr(cutoffs, function(min_cells) {
  # Count colonies above cutoff per sgid
  colonies <- liver_all %>%
    filter(cell_num >= min_cells) %>%
    group_by(sgid, gene) %>%
    summarise(n_colonies = n(), .groups = "drop")

  # Normalize by preTran
  merged <- colonies %>%
    inner_join(pretran_df %>% select(sgid, pretran = cell_num), by = "sgid") %>%
    mutate(seeding = n_colonies / pretran)

  # Ctrl mean seeding
  ctrl_seeding <- merged %>%
    filter(gene == "Ctrl") %>%
    summarise(ctrl_mean = mean(seeding, na.rm = TRUE)) %>%
    pull(ctrl_mean)

  merged %>%
    filter(gene != "Ctrl") %>%
    mutate(
      seeding_log2FC = log2(seeding / ctrl_seeding),
      cutoff = min_cells
    ) %>%
    select(sgid, gene, seeding_log2FC, cutoff)
})

# Pivot wide: genes × cutoffs
seeding_wide <- seeding_by_cutoff %>%
  pivot_wider(names_from = cutoff, values_from = seeding_log2FC, names_prefix = "cut_")

# Compute CV (SD / |mean|) across cutoffs for each gene
seeding_mat <- seeding_wide %>%
  select(starts_with("cut_")) %>%
  as.matrix()
rownames(seeding_mat) <- seeding_wide$sgid

gene_cv <- apply(seeding_mat, 1, function(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)))
seeding_mat <- seeding_mat[order(gene_cv), ]  # stable genes on top

# Clean column names
colnames(seeding_mat) <- gsub("cut_", "", colnames(seeding_mat))

# Cap for visualization
seeding_capped <- pmin(pmax(seeding_mat, -4, na.rm = TRUE), 4, na.rm = TRUE)

# Row annotation: CV
row_annot <- data.frame(CV = gene_cv[rownames(seeding_capped)])
rownames(row_annot) <- rownames(seeding_capped)

pheatmap(
  seeding_capped,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annot,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-4, 4, length.out = 101),
  na_col = "grey90",
  show_rownames = TRUE,
  fontsize_row = 3,
  fontsize_col = 10,
  main = "Seeding Log2FC across Colony Size Cutoffs — NSG Liver\n(sorted by CV, lowest on top)",
  angle_col = 45,
  filename = "seeding_stability_cutoff_heatmap_NSG_Liver.pdf",
  width = 10,
  height = 20
)

#### Mark Crebbp on scatter plots ####
## 1. Seeding vs Peak Mode
p1 <- plot_statistic_comparison(
  combined.results = boot_moba500_corrected,
  stat.x = "seeding", stat.y = "peak_mode",
  filter.vars = list(tissue = "Liver", mouse_genotype = "NSG", time_point = "3weeks"),
  title = "NSG Liver: Seeding vs Peak Mode",
  specific.genes = "Crebbp"
)
ggsave("seeding_vs_peak_mode_NSG_Liver_3weeks.pdf", p1, width = 10, height = 10)

## 2. Seeding vs Size
p2 <- plot_statistic_comparison(
  combined.results = boot_moba500_corrected,
  stat.x = "seeding", stat.y = "size",
  filter.vars = list(tissue = "Liver", mouse_genotype = "NSG", time_point = "3weeks"),
  title = "NSG Liver: Seeding vs Size (P90)",
  specific.genes = "Crebbp"
)
ggsave("seeding_vs_size_NSG_Liver_3weeks.pdf", p2, width = 10, height = 10)


## 3. preTran vs Seeding
crebbp_pretran <- pretran_seed %>% filter(sgid == "Crebbp")

p3 <- ggplot(pretran_seed, aes(x = log2(cell_num), y = seeding_log2FC)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
  geom_point(data = crebbp_pretran, color = "red", size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_label_repel(data = crebbp_pretran, aes(label = "Crebbp"),
                            color = "red", size = 4, fontface = "bold", box.padding = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("R² = %.2f", r2), size = 5, fontface = "bold") +
  labs(
    title = "preTran Representation vs Seeding — NSG Liver",
    x = "log2(preTran Cell Number)",
    y = "Seeding Log2FC"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("pretran_vs_seeding_NSG_Liver.pdf", p3, width = 10, height = 10)


## 4. Heatmap — add Crebbp flag to row annotation
row_annot$Crebbp <- ifelse(rownames(row_annot) == "Crebbp", "Yes", "No")

ann_colors <- list(
  Crebbp = c("Yes" = "red", "No" = "white")
)
pheatmap(
  seeding_capped,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annot,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-4, 4, length.out = 101),
  na_col = "grey90",
  show_rownames = TRUE,
  fontsize_row = 3,
  fontsize_col = 10,
  main = "Seeding Log2FC across Colony Size Cutoffs — NSG Liver\n(sorted by CV, lowest on top)",
  angle_col = 45,
  filename = "seeding_stability_cutoff_heatmap_NSG_Liver.pdf",
  width = 10,
  height = 20
)

#### Cross-mouse colony count consistency — NSG Liver ####

# Replace "mouse_id" below with whatever the column is called
mouse_col <- "mouse_id"  # adjust as needed

# Raw colony counts per sgID per mouse (no preTran normalization)
colony_counts <- counts_filtered_2 %>%
  filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
  group_by(.data[[mouse_col]], sgid, gene) %>%
  summarise(n_colonies = n(), .groups = "drop")

# Pick the 2 mice with the most total colonies
top_mice <- colony_counts %>%
  group_by(.data[[mouse_col]]) %>%
  summarise(total_colonies = sum(n_colonies), .groups = "drop") %>%
  slice_max(total_colonies, n = 2) %>%
  pull(.data[[mouse_col]])

mouse1 <- top_mice[1]
mouse2 <- top_mice[2]

# Pivot to wide: one column per mouse
colony_wide <- colony_counts %>%
  filter(.data[[mouse_col]] %in% c(mouse1, mouse2)) %>%
  pivot_wider(
    id_cols = c(sgid, gene),
    names_from = all_of(mouse_col),
    values_from = n_colonies,
    values_fill = 0
  )

col1 <- as.character(mouse1)
col2 <- as.character(mouse2)

r2 <- summary(lm(colony_wide[[col2]] ~ colony_wide[[col1]]))$r.squared

ggplot(colony_wide, aes(x = .data[[col1]], y = .data[[col2]])) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(
    data = colony_wide %>%
      mutate(diff = abs(.data[[col1]] - .data[[col2]])) %>%
      slice_max(diff, n = 10),
    aes(label = sgid), size = 3, max.overlaps = 20
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("R² = %.2f", r2), size = 5, fontface = "bold") +
  labs(
    title = "Raw Colony Count Consistency — NSG Liver (no preTran normalization)",
    subtitle = paste("Mouse", mouse1, "vs Mouse", mouse2),
    x = paste("Colony count —", mouse1),
    y = paste("Colony count —", mouse2)
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("cross_mouse_colony_consistency_NSG_Liver.pdf", width = 10, height = 10)

# Shared data
colony_wide <- colony_wide %>%
  mutate(
    log2_m1 = log2(.data[[col1]] + 1),
    log2_m2 = log2(.data[[col2]] + 1)
  )

r2_raw <- summary(lm(colony_wide[[col2]] ~ colony_wide[[col1]]))$r.squared
r2_log2 <- summary(lm(log2_m2 ~ log2_m1, data = colony_wide))$r.squared

crebbp_row <- colony_wide %>% filter(sgid == "Crebbp")

# Raw counts
p_raw <- ggplot(colony_wide, aes(x = .data[[col1]], y = .data[[col2]])) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(data = crebbp_row, color = "red", size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_label_repel(data = crebbp_row, aes(label = "Crebbp"),
                            color = "red", fontface = "bold", box.padding = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("R² = %.2f", r2_raw), size = 5, fontface = "bold") +
  labs(
    title = "Raw Colony Counts",
    x = paste("Colonies —", mouse1),
    y = paste("Colonies —", mouse2)
  ) +
  theme_minimal(base_size = 13)

# Log2 transformed
p_log2 <- ggplot(colony_wide, aes(x = log2_m1, y = log2_m2)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(data = crebbp_row, color = "red", size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_label_repel(data = crebbp_row, aes(label = "Crebbp"),
                            color = "red", fontface = "bold", box.padding = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("R² = %.2f", r2_log2), size = 5, fontface = "bold") +
  labs(
    title = "Log2(Colonies + 1)",
    x = paste("log2(colonies + 1) —", mouse1),
    y = paste("log2(colonies + 1) —", mouse2)
  ) +
  theme_minimal(base_size = 13)

# Side by side
combined <- p_raw + p_log2 +
  plot_annotation(
    title = "Raw Colony Count Consistency — NSG Liver (no preTran normalization)",
    subtitle = paste("Mouse", mouse1, "vs Mouse", mouse2),
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

ggsave("cross_mouse_colony_consistency_NSG_Liver_lowest.pdf", combined, width = 16, height = 8)

# Pivot wide: one column per mouse
colony_matrix <- colony_counts %>%
  pivot_wider(
    id_cols = c(sgid, gene),
    names_from = all_of(mouse_col),
    values_from = n_colonies,
    values_fill = 0
  )

# Log2 transform
mouse_cols <- setdiff(colnames(colony_matrix), c("sgid", "gene"))
colony_log2 <- colony_matrix %>%
  mutate(across(all_of(mouse_cols), ~ log2(.x + 1)))

# Pairwise R² excluding double-zeros
pairwise_r2 <- matrix(NA, length(mouse_cols), length(mouse_cols),
                      dimnames = list(mouse_cols, mouse_cols))
for (i in seq_along(mouse_cols)) {
  for (j in seq_along(mouse_cols)) {
    sub <- colony_log2 %>%
      filter(.data[[mouse_cols[i]]] > 0 | .data[[mouse_cols[j]]] > 0)
    pairwise_r2[i, j] <- cor(sub[[mouse_cols[i]]], sub[[mouse_cols[j]]])^2
  }
}

# Plot corrected heatmap
pheatmap(
  pairwise_r2,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("white", "#7C3AED"))(100),
  breaks = seq(0, 1, length.out = 101),
  main = "Pairwise R² of Colony Counts (log2, excluding double-zeros) — NSG Liver",
  fontsize_number = 12,
  filename = "mouse_pairwise_R2_heatmap_NSG_Liver.pdf",
  width = 8, height = 7
)
#specific mouse pair
mouse1 <- "RT2713"
mouse2 <- "RT2711"

# Rebuild wide format for just these two
colony_wide <- colony_counts %>%
  filter(.data[[mouse_col]] %in% c(mouse1, mouse2)) %>%
  pivot_wider(
    id_cols = c(sgid, gene),
    names_from = all_of(mouse_col),
    values_from = n_colonies,
    values_fill = 0
  ) %>%
  mutate(
    log2_m1 = log2(.data[[mouse1]] + 1),
    log2_m2 = log2(.data[[mouse2]] + 1)
  )

#### 1. R² plot — y-axis 0 to 1 ####

p3 <- ggplot(r2_by_pctile, aes(x = percentile, y = r_squared)) +
  geom_line(color = "#7C3AED", linewidth = 1) +
  geom_point(color = "#7C3AED", size = 3) +
  geom_text(aes(label = sprintf("%.2f", r_squared)), vjust = -1, size = 3.5) +
  scale_x_continuous(breaks = percentiles * 100) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "R² of Seeding vs Size across Percentiles — NSG Liver",
    x = "Size Percentile",
    y = "R² (Seeding Log2FC ~ Size Log2FC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("seeding_vs_size_R2_by_percentile_NSG_Liver.pdf", p3, width = 10, height = 6)


#### 2. Scatter plots — smaller dots, tighter theme ####

plot_stat_scatter <- function(data, x.col, y.col,
                              x.label = NULL, y.label = NULL,
                              title = NULL, specific.genes = "Crebbp",
                              point.size = 1.5, point.alpha = 0.4,
                              log2fc.threshold = 1, fdr.x = NULL, fdr.y = NULL,
                              fdr.threshold = 0.05) {

  plot_df <- data %>%
    filter(!is.na(.data[[x.col]]), !is.na(.data[[y.col]])) %>%
    rename(x_val = all_of(x.col), y_val = all_of(y.col))

  if (!is.null(fdr.x) && fdr.x %in% colnames(data)) {
    plot_df$fdr_x <- data[[fdr.x]][match(plot_df$sgid, data$sgid)]
  } else {
    plot_df$fdr_x <- 0
  }

  if (!is.null(fdr.y) && fdr.y %in% colnames(data)) {
    plot_df$fdr_y <- data[[fdr.y]][match(plot_df$sgid, data$sgid)]
  } else {
    plot_df$fdr_y <- 0
  }

  plot_df <- plot_df %>%
    mutate(category = case_when(
      x_val > log2fc.threshold &
        y_val > log2fc.threshold &
        fdr_x < fdr.threshold &
        fdr_y < fdr.threshold ~ "Both Up",

      x_val < -log2fc.threshold &
        y_val < -log2fc.threshold &
        fdr_x < fdr.threshold &
        fdr_y < fdr.threshold ~ "Both Down",

      TRUE ~ "NS"
    ))

  r2 <- summary(lm(y_val ~ x_val, data = plot_df))$r.squared

  p <- ggplot(plot_df, aes(x = x_val, y = y_val)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_smooth(
      method = "lm",
      color = "black",
      linetype = "dashed",
      se = TRUE,
      alpha = 0.15
    ) +
    geom_point(
      aes(color = category),
      size = point.size,
      alpha = point.alpha
    ) +
    scale_color_manual(
      values = c(
        "Both Up" = "#DC143C",
        "Both Down" = "#4169E1",
        "NS" = "grey30"
      )
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      hjust = 1.1,
      vjust = -0.5,
      label = sprintf("italic(R)^2 == %.2f", r2),
      size = 4,
      parse = TRUE
    )

  if (!is.null(specific.genes)) {
    spec <- plot_df %>% filter(sgid %in% specific.genes)

    if (nrow(spec) > 0) {
      p <- p +
        geom_point(
          data = spec,
          aes(x = x_val, y = y_val),
          inherit.aes = FALSE,
          color = "red",
          size = point.size + 2,
          shape = 1,
          stroke = 1.2
        ) +
        ggrepel::geom_label_repel(
          data = spec,
          aes(x = x_val, y = y_val, label = sgid),
          inherit.aes = FALSE,
          color = "red",
          fontface = "bold",
          size = 4,
          fill = "white",
          alpha = 0.85,
          box.padding = 1,
          segment.color = "red"
        )
    }
  }

  p +
    labs(
      title = title,
      x = x.label %||% x.col,
      y = y.label %||% y.col
    ) +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
}

# Filter data
nsg_liver <- boot_moba500_corrected %>%
  filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks", gene != "Ctrl")

# Seeding vs Peak Mode
p_sp <- plot_stat_scatter(nsg_liver,
                          x.col = "seeding_log2FC", y.col = "peak_mode_Log2FC",
                          x.label = "Seeding (Log2FC)", y.label = "Peak Mode (Log2FC)",
                          fdr.x = "seeding_FDR",
                          title = "seeding vs peak_mode | NSG Liver (3weeks)")
ggsave("seeding_vs_peak_mode_NSG_Liver_3weeks.pdf", p_sp, width = 10, height = 10)

# Seeding vs Size
p_ss <- plot_stat_scatter(nsg_liver,
                          x.col = "seeding_log2FC", y.col = "size_Log2FC_p90",
                          x.label = "Seeding (Log2FC)", y.label = "Size (Log2FC)",
                          fdr.x = "seeding_FDR", fdr.y = "size_FDR_pval_p90",
                          title = "seeding vs size | NSG Liver (3weeks)")
ggsave("seeding_vs_size_NSG_Liver_3weeks.pdf", p_ss, width = 10, height = 10)


#### 3. Heatmap — sort by mean seeding instead of CV ####

# Compute mean seeding per gene across cutoffs
gene_mean_seeding <- apply(seeding_mat, 1, mean, na.rm = TRUE)
seeding_mat_sorted <- seeding_mat[order(gene_mean_seeding, decreasing = TRUE), ]
seeding_capped_sorted <- pmin(pmax(seeding_mat_sorted, -4), 4)
# Row annotation with both CV and mean seeding
row_annot_sorted <- data.frame(
  CV = gene_cv[rownames(seeding_capped_sorted)],
  Crebbp = ifelse(rownames(seeding_capped_sorted) == "Crebbp", "Yes", "No")
)
rownames(row_annot_sorted) <- rownames(seeding_capped_sorted)

ann_colors <- list(
  Crebbp = c("Yes" = "red", "No" = "white")
)

pheatmap(
  seeding_capped_sorted,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annot_sorted,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-4, 4, length.out = 101),
  na_col = "grey90",
  show_rownames = TRUE,
  fontsize_row = 3,
  fontsize_col = 10,
  main = "Seeding Log2FC across Colony Size Cutoffs — NSG Liver\n(sorted by mean seeding, highest on top)",
  angle_col = 45,
  filename = "seeding_stability_cutoff_heatmap_NSG_Liver.pdf",
  width = 10,
  height = 20
)
