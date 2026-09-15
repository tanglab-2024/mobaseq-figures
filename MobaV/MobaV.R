#### MobaV QC####
setwd("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal")
counts_df <- load_counts("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/MobaV_FilteredSampleInfo.csv")
# Get one row per sample
sample_data <- counts_df %>%
  distinct(sample_id, spike_rsq, tissue)

# Calculate mean and SD for each tissue
tissue_stats <- sample_data %>%
  group_by(tissue) %>%
  summarise(
    mean_rsq = mean(spike_rsq, na.rm = TRUE),
    sd_rsq = sd(spike_rsq, na.rm = TRUE),
    n = n()
  ) %>%
  mutate(label = sprintf("%.4f ± %.4f", mean_rsq, sd_rsq))

# Create the histogram
p <- ggplot(sample_data, aes(x = spike_rsq, fill = tissue)) +
  geom_histogram(binwidth = 0.01, position = "stack",
                 color = "black", linewidth = 0.3) +
  labs(title = "Histogram of R² of Cell Count Distribution Calculated from Spike-In Regression",
       x = "R²",
       y = "Sample Count",
       fill = "Tissue") +
  theme_bw() +
  theme(legend.position = "right")

# Option 1: Add text annotations on the plot
p + geom_text(data = tissue_stats,
              aes(x = 0.905, y = 40 - seq(0, by = 2.5, length.out = nrow(tissue_stats)),
                  label = paste0(tissue, ": ", label), fill = NULL),
              hjust = 0, size = 3, fontface = "bold")

spike_data <- load_spikein_counts("/Volumes/ruit/Active/AndyXu/2024_MobaV_Fullset/outputs_1/MobaV_SpikeInInfo.csv")
spike_data <- spike_data %>%
  mutate(
    mouse_id = sub("-.*", "", sample_id),
    tissue = sub(".*-", "", sample_id)
  )
result <- plot_spikein_regression(spike_data, color.by = "name", shape.by = "tissue", show.legend = T, y.limits = c(0, .22), show.boxplot = T, show.stats = T)
print(result$plot)

# Mouse seeding histogram
# Create full_bc column
counts_df$full_bc <- paste0(counts_df$sgid, "_", counts_df$barcode)
# Sum cell_num for each full_bc in preTran samples
pretran_bc_cell_counts <- counts_df %>%
  filter(grepl("preTran", sample_id, ignore.case = TRUE)) %>%
  mutate(full_bc = paste0(sgid, "_", barcode)) %>%
  group_by(full_bc) %>%
  summarise(total_cells = sum(cell_num), .groups = "drop")
# Calculate statistics
mean_cells <- mean(pretran_bc_cell_counts$total_cells)
max_cells <- max(pretran_bc_cell_counts$total_cells)
percentile_99 <- quantile(pretran_bc_cell_counts$total_cells, 0.99)
total_full_bc <- nrow(pretran_bc_cell_counts)
total_cells <- sum(pretran_bc_cell_counts$total_cells)
# Calculate mode (most frequent value)
mode_cells <- as.numeric(names(sort(table(pretran_bc_cell_counts$total_cells), decreasing = TRUE)[1]))

# Create histogram
ggplot(pretran_bc_cell_counts, aes(x = total_cells)) +
  geom_histogram(bins = 50, fill = "black", color = "white") +
  geom_vline(xintercept = percentile_99, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = percentile_99, y = Inf,
           label = sprintf("99%% percentile\n%.0f", percentile_99),
           color = "red", vjust = 1.5, hjust = -0.1, size = 3.5) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Total Uni_BC: %s\nTotal Cells: %s\nMean: %.0f\nMode: %.0f\nMax: %.0f",
                           scales::comma(total_full_bc),
                           scales::comma(total_cells),
                           mean_cells,
                           mode_cells,
                           max_cells),
           hjust = 1.1, vjust = 1.5, size = 3.5) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = function(x) paste0(x/1000, "k")) +
  labs(x = "Total Cell Count per Uni_BC (log scale)",
       y = "Number of Uni_BC",
       title = "Distribution of Cell Counts for Uni_BC in preTran Samples") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# Calculate 99th percentile of cell counts in preTran
pretran_bc_cell_counts <- counts_df %>%
  filter(grepl("preTran", sample_id, ignore.case = TRUE)) %>%
  mutate(full_bc = paste0(sgid, "_", barcode)) %>%
  group_by(full_bc) %>%
  summarise(total_cells = sum(cell_num), .groups = "drop")

percentile_99 <- quantile(pretran_bc_cell_counts$total_cells, 0.99)

# Filter to keep only full_bc below 99th percentile
pretran_bcs_filtered <- pretran_bc_cell_counts %>%
  filter(total_cells <= percentile_99) %>%
  pull(full_bc)

# Count number of mice for these filtered full_bc, grouped by genotype and time point
bc_mouse_counts_by_genotype_time <- counts_df %>%
  filter(!is.na(mouse_genotype),
         mouse_genotype != "",
         mouse_genotype %in% c("BL6", "NSG"),
         !is.na(time_point),
         time_point %in% c("1week", "3weeks")) %>%  # Keep only 1week and 3weeks
  mutate(full_bc = paste0(sgid, "_", barcode)) %>%
  filter(full_bc %in% pretran_bcs_filtered) %>%
  group_by(full_bc, mouse_genotype, time_point) %>%
  summarise(n_mice = n_distinct(mouse_id), .groups = "drop") %>%
  group_by(mouse_genotype, time_point, n_mice) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mouse_genotype, time_point) %>%
  mutate(percentage = n / sum(n) * 100)

# Plot with facets
ggplot(bc_mouse_counts_by_genotype_time, aes(x = factor(n_mice), y = percentage)) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(aes(label = sprintf("%.1f%%\n(%d)", percentage, n)),
            vjust = -0.5, size = 2.5) +
  facet_grid(time_point ~ mouse_genotype, scales = "free_y") +
  labs(x = "Seed Mice",
       y = "Percentage of Uni_BC",
       title = "Percentage of Uni_BC per Seed Mice Count by Genotype and Time Point\n(preTran BCs ≤99th percentile in full dataset)") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(0, max(bc_mouse_counts_by_genotype_time$percentage) * 1.15)

#### Moba V ####
counts_V <- load_counts("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/MobaV_FilteredSampleInfo.csv",
                        remove.sgids = "sgDummy",
                        keep.bc.end = c("CTGA", "GGAA"),
                        keep.distance = 0)
counts_V <- load_counts("C:/Users/Andy/Box/AX_data_analysis/Mobaseq/MobaVal/MobaV_FilteredSampleInfo.csv",
                        remove.sgids = "sgDummy",
                        keep.bc.end = c("CTGA", "GGAA"),
                        keep.distance = 0)
counts_V <- counts_V %>%
  mutate(tissue = ifelse(tissue == "Cells", "preTran", tissue))

result <- plot_pretran_histogram(counts_V, percentile.threshold = 0.99, filter.data = TRUE, return.list = TRUE)
print(result$plot)
counts_filtered_V <- result$filtered_df

result <- plot_barcode_seeding(counts_filtered_V, return.data = TRUE)
print(result$plot)
counts_filtered_V <- counts_filtered_V %>%
  mutate(
    tissue = recode(tissue,
                    "WholeBlood" = "Blood",
    )
  )
#### For testing
counts_filtered_V <- counts_filtered_V %>%
  mutate(
    gene = case_when(
      grepl("^(sgnt|sgneo|sgsafe)", sgid, ignore.case = TRUE) ~ "Ctrl",
      grepl("mmu-sg", sgid, ignore.case = TRUE) ~ sub("mmu-sg.*", "", sgid, ignore.case = TRUE),
      TRUE ~ sgid
    )
  )
####
counts_filtered_2V <- counts_filtered_V %>%
  # Filter out unwanted entries
  filter(
    sgid != "sgDummy",                      # remove sgDummy
    bc_end %in% c("CTGA", "GGAA"),          # keep only CTGA or GGAA
    distance == 0,                          # keep only distance == 0
    (tissue == "preTran" & cell_num >= 2) |
      (tissue == "Blood" & cell_num >= 1) |
      (tissue == "Brain" & cell_num >= 10) |
      (tissue == "BoneMarrow" & cell_num >= 5) |        # NEW: bone marrow cutoff
      (!tissue %in% c("preTran", "Blood", "Brain", "BoneMarrow") & cell_num >= 100)
  ) %>%
  # Create gene column
  mutate(
    gene = case_when(
      grepl("^(sgnt|sgneo|sgsafe)", sgid, ignore.case = TRUE) ~ "Ctrl",  # control guides
      grepl("mmu-sg", sgid, ignore.case = TRUE) ~ sub("mmu-sg.*", "", sgid, ignore.case = TRUE),  # e.g. Bag1mmu-sg2 → Bag1
      TRUE ~ sgid  # fallback
    )
  )
counts_filtered_2V %>%
  count(bc_end)

GGAA_df <- counts_filtered_2V %>%
  filter(bc_end == "GGAA")

CTGA_df <- counts_filtered_2V %>%
  filter(bc_end == "CTGA")

plot_composition_many(GGAA_df,
                      stack.by = "sgid",
                      color.by = "sgid",
                      mouse.genotype = "BL6",
                      tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                      time.points = "3weeks",
                      color.palette = "category20c",
                      highlight.top = "Crebbp",
                      ctrl.colors = "black")
plot_composition_many(CTGA_df,
                      stack.by = "sgid",
                      color.by = "sgid",
                      mouse.genotype = "NSG",
                      tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                      time.points = "3weeks",
                      color.palette = "category20c",
                      highlight.top = "Crebbp",
                      ctrl.colors = "black")

# Filter for preTran entries based on sample_id
preTran_df <- CTGA_df %>%
  filter(sample_id == "preTran")

# Compute total cell_num per sgid
sg_summary <- preTran_df %>%
  group_by(gene, sgid) %>%
  summarise(
    total_cellnum = sum(cell_num, na.rm = TRUE),
    .groups = "drop"
  )

# Split into Crebbp and Ctrl subsets
crebbp_df <- sg_summary %>% filter(gene == "Crebbp")
ctrl_df   <- sg_summary %>% filter(gene == "Ctrl")

# Compute all pairwise absolute differences in total cell_num
pairwise_diff <- expand.grid(crebbp_df$sgid, ctrl_df$sgid) %>%
  rename(Crebbp_sgid = Var1, Ctrl_sgid = Var2) %>%
  left_join(crebbp_df %>% select(sgid, total_cellnum), by = c("Crebbp_sgid" = "sgid")) %>%
  rename(Crebbp_total = total_cellnum) %>%
  left_join(ctrl_df %>% select(sgid, total_cellnum), by = c("Ctrl_sgid" = "sgid")) %>%
  rename(Ctrl_total = total_cellnum) %>%
  mutate(diff = abs(Crebbp_total - Ctrl_total))

# Find the most similar pair
best_pair <- pairwise_diff %>%
  arrange(diff) %>%
  slice(1)

# Retrieve total cell_num for those sgids
best_totals <- sg_summary %>%
  filter(sgid %in% c(best_pair$Crebbp_sgid, best_pair$Ctrl_sgid)) %>%
  select(gene, sgid, total_cellnum)

# Print results
print(best_totals)

sgid_totals <- counts_filtered_2V %>%
  filter(sample_id == "preTran", gene %in% c("Crebbp", "Ctrl")) %>%
  group_by(gene, sgid) %>%
  summarise(total_cell_num = sum(cell_num, na.rm = TRUE), .groups = "drop")

print(sgid_totals)

plot_cell_num_jitter_faceted(
  GGAA_df,
  facet.by = "tissue",
  color.by = "sgid",
  mouse.genotypes = "NSG",
  time.points = c("3weeks"),
  filter.facets = c("Liver", "Lung", "Brain", "BoneMarrow"),
  filter.colors = c("sgNT2-MMU", "Crebbpmmu-sg2"),
  title = "Cell Number Seeding Comparison NSG 3weeks: Crebbpmmu-sg2 vs sgNT2-MMU",
  facet.nrow = 2
)

plot_cell_num_jitter_faceted(
  CTGA_df,
  facet.by = "tissue",
  color.by = "sgid",
  mouse.genotypes = "NSG",
  time.points = c("1week", "3weeks"),
  filter.facets = c("Liver"),
  filter.colors = c("sgNT2-MMU", "Crebbpmmu-sg2"),
  title = "Cell Number Seeding Comparison NSG 1-3weeks: Crebbpmmu-sg2 vs sgNT2-MMU",
  facet.nrow = 2
)

# resultsV <- dormancy_cutoff_by_genotype_tissue_time_v4(
#   counts_filtered_2V,
#   tissues = c("Liver", "Lung", "Brain"),
#   time.points = c("3weeks"),
#   min.colonies = 100,
#   bw.adjust = 2,
#   split.by.tissue = FALSE
# )
# resultsV$plot_list$BL6_AllTissues
# resultsV$plot_list$NSG_AllTissues
#
# plot_dormancy_cutoffs_v2(
#   counts.df = counts_filtered_2V,
#   cutoffs_df = resultsV$cutoffs_df,
#   value.col = "log2_cell_num",
#   time.points = c("3weeks"),
#   tissues = c("Liver", "Lung", "Brain", "BM"),
#   genotypes = unique(counts_filtered_2V$mouse_genotype),
#   split.by.tissue = FALSE
# )
# resultsV$cutoffs_df
# geno_cutoffsV <- setNames(
#   resultsV$cutoffs_df$dormancy_cutoff,
#   resultsV$cutoffs_df$mouse_genotype
# )
# 2^geno_cutoffsV

pretran_dataVG <- generate_pretran_data(GGAA_df)

boot_mobaVG <- calculate_metastatic_statistics(
  counts.df = GGAA_df,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile"),
  pretran.cell.num = pretran_dataVG$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVG$pretran_cell_num_gene,
  debug = T,
  R = 1000
)
write.csv(boot_mobaVG, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVG_stats.csv", row.names = FALSE)

sgid_totals <- GGAA_df %>%
  filter(sample_id == "preTran", gene %in% c("Crebbp", "Ctrl")) %>%
  group_by(gene, sgid) %>%
  summarise(total_cell_num = sum(cell_num, na.rm = TRUE), .groups = "drop")

print(sgid_totals)

my_targets <- c("Crebbpmmu-sg1", "Crebbpmmu-sg2", "Crebbpmmu-sg3", "sgSafe28-MMU", "sgSafe29-MMU",
                "sgSafe30-MMU")
gene_map <- c("sgSafe28-MMU" = "Ctrl", "sgSafe29-MMU" = "Ctrl", "sgSafe30-MMU" = "Ctrl",
              "Crebbpmmu-sg1" = "Crebbp", "Crebbpmmu-sg2" = "Crebbp", "Crebbpmmu-sg3" = "Crebbp")

plot_bootstrap_timecourse(boot_mobaVG, my_targets, gene_map, metric = "size_p", plot.type = "point", genotypes = "NSG", tissues = "Liver")

# Define parameters
metrics <- c("seeding", "burden", "size_p")
genotypes <- c("BL6", "NSG")
out_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/GGAA"  # replace with your output directory

# Ensure output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Loop through combinations
for (metric in metrics) {
  for (geno in genotypes) {
    # Generate plot
    p <- plot_bootstrap_timecourse(
      boot_mobaVG,
      my_targets,
      gene_map,
      metric = metric,
      plot.type = "point",
      genotypes = geno,
      tissues = "Liver"
    )

    # Define file name
    file_name <- paste0(out_dir, "/", metric, "_", geno, "_timecourse.pdf")

    # Save plot
    ggsave(
      filename = file_name,
      plot = p,
      width = 11,
      height = 6,
      units = "in"
    )

    message("Saved: ", file_name)
  }
}

plot_met_stat_comparison(
  data = boot_mobaVG,
  comparison.var = "tissue",
  category1 = "Liver",
  category2 = "Lung",
  stat.type = "size",
  plot.by = "sgid",
  use.log2 = F,
  filter.vars = list(mouse_genotype = "NSG", time_point = "3weeks"),
  fdr.threshold = 0.05
)

tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
genotypes <- c("NSG", "BL6")
stat.types <- c("seeding", "burden", "size")
time_point <- "3weeks"

output_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/GGAA/Met_Stat_Plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (stat in stat.types) {

  for (gen in genotypes) {
    # ---- Tissue–Tissue comparisons within each genotype ----
    for (i in 1:(length(tissues) - 1)) {
      for (j in (i + 1):length(tissues)) {
        cat("Generating plot for", stat, "in", gen, ":", tissues[i], "vs", tissues[j], "\n")

        p <- tryCatch({
          plot_met_stat_comparison(
            data = boot_mobaVG,
            comparison.var = "tissue",
            category1 = tissues[i],
            category2 = tissues[j],
            stat.type = stat,
            plot.by = "sgid",
            use.log2 = FALSE,
            size.percentile = 90,
            filter.vars = list(mouse_genotype = gen, time_point = time_point),
            fdr.threshold = 0.05,
            title = paste0(
              toupper(substr(stat, 1, 1)), substr(stat, 2, nchar(stat)),
              " comparison: ", tissues[i], " vs ", tissues[j], " (", gen, ")"
            )
          )
        }, error = function(e) {
          message("Skipping ", stat, " ", gen, " ", tissues[i], " vs ", tissues[j], ": ", e$message)
          NULL
        })

        if (!is.null(p)) {
          pdf_file <- file.path(
            output_dir,
            paste0("GGAA_Met_", stat, "_", gen, "_", tissues[i], "_vs_", tissues[j], "_", time_point, ".pdf")
          )
          ggsave(pdf_file, p, width = 8, height = 8)
        }
      }
    }
  }

  # ---- Genotype–Genotype comparisons within each tissue ----
  for (tissue in tissues) {
    cat("Generating plot for", stat, "in", tissue, "NSG vs BL6\n")

    p <- tryCatch({
      plot_met_stat_comparison(
        data = boot_mobaVG,
        comparison.var = "mouse_genotype",
        category1 = "NSG",
        category2 = "BL6",
        stat.type = stat,
        plot.by = "sgid",
        use.log2 = FALSE,
        size.percentile = 90,
        filter.vars = list(tissue = tissue, time_point = time_point),
        fdr.threshold = 0.05,
        title = paste0(
          toupper(substr(stat, 1, 1)), substr(stat, 2, nchar(stat)),
          " comparison: NSG vs BL6 (", tissue, ")"
        )
      )
    }, error = function(e) {
      message("Skipping ", stat, " NSG vs BL6 ", tissue, ": ", e$message)
      NULL
    })

    if (!is.null(p)) {
      pdf_file <- file.path(
        output_dir,
        paste0("GGAA_Met_", stat, "_", tissue, "_NSG_vs_BL6_", time_point, ".pdf")
      )
      ggsave(pdf_file, p, width = 8, height = 8)
    }
  }
}

pretran_dataVC <- generate_pretran_data(CTGA_df)

boot_mobaVC <- calculate_metastatic_statistics(
  counts.df = CTGA_df,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile"),
  pretran.cell.num = pretran_dataVC$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVC$pretran_cell_num_gene,
  debug = T,
  R = 1000
)
write.csv(boot_mobaVC, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVC_stats.csv", row.names = FALSE)

sgid_totalsC <- CTGA_df %>%
  filter(sample_id == "preTran", gene %in% c("Crebbp", "Ctrl")) %>%
  group_by(gene, sgid) %>%
  summarise(total_cell_num = sum(cell_num, na.rm = TRUE), .groups = "drop")

print(sgid_totalsC)

my_targets <- c("Crebbpmmu-sg1", "Crebbpmmu-sg2", "Crebbpmmu-sg3", "sgSafe4-MMU", "sgNT3-MMU",
                "sgSafe30-MMU")
gene_map <- c("sgSafe4-MMU" = "Ctrl", "sgNT3-MMU" = "Ctrl", "sgSafe30-MMU" = "Ctrl",
              "Crebbpmmu-sg1" = "Crebbp", "Crebbpmmu-sg2" = "Crebbp", "Crebbpmmu-sg3" = "Crebbp")

plot_bootstrap_timecourse(boot_mobaVC, my_targets, gene_map, metric = "size_p", plot.type = "point", genotypes = "BL6", tissues = "Liver")

# Define parameters
metrics <- c("seeding", "burden", "size_p")
genotypes <- c("BL6", "NSG")
out_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/CTGA"  # replace with your output directory

# Ensure output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Loop through combinations
for (metric in metrics) {
  for (geno in genotypes) {
    # Generate plot
    p <- plot_bootstrap_timecourse(
      boot_mobaVC,
      my_targets,
      gene_map,
      metric = metric,
      plot.type = "point",
      genotypes = geno,
      tissues = "Liver"
    )

    # Define file name
    file_name <- paste0(out_dir, "/", metric, "_", geno, "_timecourse.pdf")

    # Save plot
    ggsave(
      filename = file_name,
      plot = p,
      width = 11,
      height = 6,
      units = "in"
    )

    message("Saved: ", file_name)
  }
}
tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
genotypes <- c("NSG", "BL6")
stat.types <- c("seeding", "burden", "size")
time_point <- "3weeks"

output_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/CTGA/Met_Stat_Plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (stat in stat.types) {

  for (gen in genotypes) {
    # ---- Tissue–Tissue comparisons within each genotype ----
    for (i in 1:(length(tissues) - 1)) {
      for (j in (i + 1):length(tissues)) {
        cat("Generating plot for", stat, "in", gen, ":", tissues[i], "vs", tissues[j], "\n")

        p <- tryCatch({
          plot_met_stat_comparison(
            data = boot_mobaVC,
            comparison.var = "tissue",
            category1 = tissues[i],
            category2 = tissues[j],
            stat.type = stat,
            plot.by = "sgid",
            use.log2 = FALSE,
            size.percentile = 90,
            filter.vars = list(mouse_genotype = gen, time_point = time_point),
            fdr.threshold = 0.05,
            title = paste0(
              toupper(substr(stat, 1, 1)), substr(stat, 2, nchar(stat)),
              " comparison: ", tissues[i], " vs ", tissues[j], " (", gen, ")"
            )
          )
        }, error = function(e) {
          message("Skipping ", stat, " ", gen, " ", tissues[i], " vs ", tissues[j], ": ", e$message)
          NULL
        })

        if (!is.null(p)) {
          pdf_file <- file.path(
            output_dir,
            paste0("CTGA_Met_", stat, "_", gen, "_", tissues[i], "_vs_", tissues[j], "_", time_point, ".pdf")
          )
          ggsave(pdf_file, p, width = 8, height = 8)
        }
      }
    }
  }

  # ---- Genotype–Genotype comparisons within each tissue ----
  for (tissue in tissues) {
    cat("Generating plot for", stat, "in", tissue, "NSG vs BL6\n")

    p <- tryCatch({
      plot_met_stat_comparison(
        data = boot_mobaVC,
        comparison.var = "mouse_genotype",
        category1 = "NSG",
        category2 = "BL6",
        stat.type = stat,
        plot.by = "sgid",
        use.log2 = FALSE,
        size.percentile = 90,
        filter.vars = list(tissue = tissue, time_point = time_point),
        fdr.threshold = 0.05,
        title = paste0(
          toupper(substr(stat, 1, 1)), substr(stat, 2, nchar(stat)),
          " comparison: NSG vs BL6 (", tissue, ")"
        )
      )
    }, error = function(e) {
      message("Skipping ", stat, " NSG vs BL6 ", tissue, ": ", e$message)
      NULL
    })

    if (!is.null(p)) {
      pdf_file <- file.path(
        output_dir,
        paste0("CTGA_Met_", stat, "_", tissue, "_NSG_vs_BL6_", time_point, ".pdf")
      )
      ggsave(pdf_file, p, width = 8, height = 8)
    }
  }
}

# Comparison between cell lines
plot_df_comparison(
  df1 = boot_mobaVG,
  df2 = boot_mobaVC,
  use.log2 = T,
  plot.by = "sgid",
  stat.type = "burden",
  filter.vars = list(
    tissue = "Liver",
    mouse_genotype = "NSG",
    time_point = "1week"
  ),
  df1.name = "RP48 GGAA",
  df2.name = "RP116 CTGA",
  fdr.threshold = 0.05,
)

# Define tissues, genotypes, and stat types
tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
genotypes <- c("NSG", "BL6")
stat_types <- c("burden", "seeding", "size")

# Define time points per tissue
time_points <- list(
  Liver = c("1week", "3weeks"),
  Lung = "3weeks",
  Brain = "3weeks",
  BoneMarrow = "3weeks"
)

# Define your output directory
output_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/cell_line_comparison"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Loop through all combinations
for (t in tissues) {
  for (g in genotypes) {
    for (s in stat_types) {
      for (tp in time_points[[t]]) {

        # Create subfolder for tissue and time point
        subfolder <- file.path(output_dir, t, tp)
        if (!dir.exists(subfolder)) dir.create(subfolder, recursive = TRUE)

        message(paste("Generating plot:", s, "-", t, "-", g, "-", tp))

        # Try generating and saving the plot safely
        tryCatch({
          p <- plot_df_comparison(
            df1 = boot_mobaVG,
            df2 = boot_mobaVC,
            use.log2 = F,
            plot.by = "sgid",
            stat.type = s,
            filter.vars = list(
              tissue = t,
              mouse_genotype = g,
              time_point = tp
            ),
            df1.name = "RP48 GGAA",
            df2.name = "RP116 CTGA",
            fdr.threshold = 0.05
          )

          # Only save if a plot is returned
          if (!is.null(p)) {
            # Construct descriptive filename
            filename <- paste0(
              tolower(t), "_",
              tp, "_",
              s, "_",
              tolower(g), "_cell_comparison.pdf"
            )

            filepath <- file.path(subfolder, filename)

            # Save the plot as 8x8 PDF
            ggsave(
              filename = filepath,
              plot = p,
              width = 8,
              height = 8
            )

            message(paste("✅ Saved:", filepath))
          } else {
            warning(paste("⚠️ Skipping (no plot returned):", s, "-", t, "-", g, "-", tp))
          }

        }, error = function(e) {
          warning(paste("❌ Failed:", s, "-", t, "-", g, "-", tp, "\nReason:", e$message))
        })
      }
    }
  }
}

#### Count CTCs per mL blood ####
library(readxl)
counting_table_MobaV <- read_excel("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/CountingTable.xlsx")

plot_ctc_counts(
  counts.df = counts_filtered_V,
  counting.table = counting_table_MobaV,
  output.path = "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/ctc_count.pdf"
)

#### Chord plots ####
# Function to create chord diagram for one genotype
library(circlize)

dissem <- plot_chord_by_genotype(counts_filtered_V, "BL6", mouse.ids = "73328",
                                 genes = c("Crebbp", "Ctrl"), project = "Moba V")
dissem <- plot_chord_by_genotype(counts_filtered_V, "BL6",
                                 genes = c("Brk1", "Ctrl"), project = "Moba V")
plot_dissemination_heatmap(dissem)

plot_chord_by_genotype(counts_filtered_V, "NSG", mouse.ids = "73201",
                       genes = c("Crebbp", "Ctrl"), project = "Moba V")
dissem_nsg <- plot_chord_by_genotype(counts_filtered_V, "NSG",
                                     genes = c("Crebbp", "Ctrl"), project = "Moba V")
plot_dissemination_heatmap(dissem_nsg)
# Computes and plots in one call, returns df invisibly
dissem_bl6_all <- compute_dissemination(counts_filtered_V, "BL6", project = "Moba V")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/dissem_heatmap2_BL6.pdf",
       dissem_bl6_all$plot, width = 15, height = 12)
dissem_nsg_all <- compute_dissemination(counts_filtered_V, "NSG", project = "Moba V")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/dissem_heatmap2_NSG.pdf",
       dissem_nsg_all$plot, width = 15, height = 12)

# chord plot for mouse genotype comparison
plot_chord_by_genotype(counts_filtered_V, genotypes = c("BL6", "NSG"),
                       genes = "Crebbp", project = "Moba V",
                       save.pdf = "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/chord_BL6_vs_NSG_Crebbp.pdf")

plot_chord_by_genotype(counts_filtered_V, genotypes = c("BL6", "NSG"),
                       genes = "Ctrl", project = "Moba V",
                       save.pdf = "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/chord_BL6_vs_NSG_Ctrl.pdf")

# Count liver colonies per mouse for each genotype
liver_counts <- counts_filtered_V %>%
  filter(tissue == "Liver") %>%
  mutate(uid_bc = paste0(mouse_id, "_", full_bc)) %>%
  group_by(mouse_id, mouse_genotype) %>%
  summarise(n_colonies = n_distinct(uid_bc), .groups = "drop") %>%
  arrange(mouse_genotype, desc(n_colonies))

# Total per genotype
liver_counts %>%
  group_by(mouse_genotype) %>%
  summarise(total_colonies = sum(n_colonies), n_mice = n())

nsg_subset <- c("73200", "73207")

counts_downsampled <- counts_filtered_V %>%
  filter(
    (mouse_genotype == "BL6") |
      (mouse_genotype == "NSG" & mouse_id %in% nsg_subset)
  )

# Verify
counts_downsampled %>%
  filter(tissue == "Liver") %>%
  mutate(uid_bc = paste0(mouse_id, "_", full_bc)) %>%
  group_by(mouse_genotype) %>%
  summarise(n_colonies = n_distinct(uid_bc))

# Generate liver-focused chord plots
plot_chord_by_genotype(counts_downsampled, genotypes = c("BL6", "NSG"),
                       genes = "Crebbp", project = "Moba V",
                       focal.tissue = "Liver",
                       save.pdf = "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/chord_liver_BL6_vs_NSG_Crebbp.pdf")

plot_chord_by_genotype(counts_downsampled, genotypes = c("BL6", "NSG"),
                       genes = "Ctrl", project = "Moba V",
                       focal.tissue = "Liver",
                       save.pdf = "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/chord_liver_BL6_vs_NSG_Ctrl.pdf")

# Heatmaps: BL6 vs NSG for Crebbp
result_crebbp <- compute_dissemination(counts_filtered_V, genotypes = c("BL6", "NSG"),
                                       genes = "Crebbp", project = "Moba V")

# Heatmaps: BL6 vs NSG for Ctrl
result_ctrl <- compute_dissemination(counts_filtered_V, genotypes = c("BL6", "NSG"),
                                     genes = "Ctrl", project = "Moba V")

# Save
ggsave("dissem_Crebbp_BL6_vs_NSG.pdf", result_crebbp$plot, width = 12, height = 6)
ggsave("dissem_Ctrl_BL6_vs_NSG.pdf", result_ctrl$plot, width = 12, height = 6)

# updating heatmaps
# Compute for BL6 vs NSG, Crebbp and Ctrl separately
dissem_crebbp <- compute_dissemination(counts_filtered_V, genotypes = c("BL6", "NSG"),
                                       genes = "Crebbp", project = "Moba V", plot = FALSE)$data

dissem_ctrl <- compute_dissemination(counts_filtered_V, genotypes = c("BL6", "NSG"),
                                     genes = "Ctrl", project = "Moba V", plot = FALSE)$data

genes_of_interest <- c("Ctrl", "Tsc1", "Pten", "Tsc2", "Gpatch8", "Mga",
                       "Gata6", "Crebbp", "Zmiz1", "Zeb1", "Csnk1a1", "Grhpr", "Hnf4a")

path <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/"

# Compute separately
dissem_bl6 <- compute_dissemination(counts_filtered_V, "BL6",
                                    genes = genes_of_interest,
                                    project = "Moba V", plot = FALSE)$data

dissem_nsg <- compute_dissemination(counts_filtered_V, "NSG",
                                    genes = genes_of_interest,
                                    project = "Moba V", plot = FALSE)$data

# BL6 new format
p1 <- plot_dissemination_heatmap(dissem_bl6, blood.in.diagonal = TRUE)
ggsave(paste0(path, "dissem_BL6_selected_new.pdf"), p1, width = 20, height = 12)

# BL6 old format
p2 <- plot_dissemination_heatmap(dissem_bl6, blood.in.diagonal = FALSE, normalize.to = NULL)
ggsave(paste0(path, "dissem_BL6_selected_old.pdf"), p2, width = 20, height = 12)

# NSG new format
p3 <- plot_dissemination_heatmap(dissem_nsg, blood.in.diagonal = TRUE)
ggsave(paste0(path, "dissem_NSG_selected_new.pdf"), p3, width = 20, height = 12)

# NSG old format
p4 <- plot_dissemination_heatmap(dissem_nsg, blood.in.diagonal = FALSE, normalize.to = NULL)
ggsave(paste0(path, "dissem_NSG_selected_old.pdf"), p4, width = 20, height = 12)

# dissemination scoring
path <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/"

td_bl6 <- plot_tissue_destinations(counts_filtered_V, "BL6", project = "Moba V")
ggsave(paste0(path, "tissue_dest_BL6.pdf"), td_bl6$plot, width = 18, height = 6)

td_nsg <- plot_tissue_destinations(counts_filtered_V, "NSG", project = "Moba V")
ggsave(paste0(path, "tissue_dest_NSG.pdf"), td_nsg$plot, width = 18, height = 6)

#### Mga Analysis - MobaV (RP48 GGAA & RP116 CTGA) ####
setwd("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal")

# Mirrors the Crebbp analysis workflow for Mga

library(mobaseq)
library(dplyr)
library(ggplot2)
library(circlize)

# 1. Load & filter data (same as existing pipeline)

counts_V <- load_counts(
  "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/MobaV_FilteredSampleInfo.csv",
  remove.sgids = "sgDummy",
  keep.bc.end = c("CTGA", "GGAA"),
  keep.distance = 0
)
counts_V <- counts_V %>%
  mutate(tissue = ifelse(tissue == "Cells", "preTran", tissue))

result <- plot_pretran_histogram(counts_V, percentile.threshold = 0.99, filter.data = TRUE, return.list = TRUE)
counts_filtered_V <- result$filtered_df

result <- plot_barcode_seeding(counts_filtered_V, return.data = TRUE)
print(result$plot)

counts_filtered_V <- counts_filtered_V %>%
  mutate(tissue = recode(tissue, "WholeBlood" = "Blood"))

# Create gene column
counts_filtered_V <- counts_filtered_V %>%
  mutate(
    gene = case_when(
      grepl("^(sgnt|sgneo|sgsafe)", sgid, ignore.case = TRUE) ~ "Ctrl",
      grepl("mmu-sg", sgid, ignore.case = TRUE) ~ sub("mmu-sg.*", "", sgid, ignore.case = TRUE),
      TRUE ~ sgid
    )
  )

# Apply tissue-specific cell_num cutoffs
counts_filtered_2V <- counts_filtered_V %>%
  filter(
    sgid != "sgDummy",
    bc_end %in% c("CTGA", "GGAA"),
    distance == 0,
    (tissue == "preTran" & cell_num >= 2) |
      (tissue == "Blood" & cell_num >= 1) |
      (tissue == "Brain" & cell_num >= 10) |
      (tissue == "BoneMarrow" & cell_num >= 5) |
      (!tissue %in% c("preTran", "Blood", "Brain", "BoneMarrow") & cell_num >= 100)
  ) %>%
  mutate(
    gene = case_when(
      grepl("^(sgnt|sgneo|sgsafe)", sgid, ignore.case = TRUE) ~ "Ctrl",
      grepl("mmu-sg", sgid, ignore.case = TRUE) ~ sub("mmu-sg.*", "", sgid, ignore.case = TRUE),
      TRUE ~ sgid
    )
  )

# Split by cell line
GGAA_df <- counts_filtered_2V %>% filter(bc_end == "GGAA")
CTGA_df <- counts_filtered_2V %>% filter(bc_end == "CTGA")

# 2. Discover Mga sgIDs present in data
mga_sgids_GGAA <- GGAA_df %>%
  filter(gene == "Mga") %>%
  distinct(sgid) %>%
  pull(sgid) %>%
  sort()
cat("Mga sgIDs in GGAA (RP48):\n"); print(mga_sgids_GGAA)

mga_sgids_CTGA <- CTGA_df %>%
  filter(gene == "Mga") %>%
  distinct(sgid) %>%
  pull(sgid) %>%
  sort()
cat("Mga sgIDs in CTGA (RP116):\n"); print(mga_sgids_CTGA)

# Identify ctrl sgIDs in each cell line (pick 3 matched controls)
ctrl_sgids_GGAA <- GGAA_df %>%
  filter(gene == "Ctrl") %>%
  distinct(sgid) %>%
  pull(sgid) %>%
  sort()
cat("Ctrl sgIDs in GGAA:\n"); print(ctrl_sgids_GGAA)

ctrl_sgids_CTGA <- CTGA_df %>%
  filter(gene == "Ctrl") %>%
  distinct(sgid) %>%
  pull(sgid) %>%
  sort()
cat("Ctrl sgIDs in CTGA:\n"); print(ctrl_sgids_CTGA)

# 3. Stacked composition plots - Mga highlighted
# RP48 (GGAA) - BL6
plot_composition_many(GGAA_df,
                      stack.by = "sgid",
                      color.by = "sgid",
                      mouse.genotype = "BL6",
                      tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                      time.points = "3weeks",
                      color.palette = "category20c",
                      highlight.top = "Mga",
                      ctrl.colors = "black")

# RP48 (GGAA) - NSG
plot_composition_many(GGAA_df,
                      stack.by = "sgid",
                      color.by = "sgid",
                      mouse.genotype = "NSG",
                      tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                      time.points = "3weeks",
                      color.palette = "category20c",
                      highlight.top = "Mga",
                      ctrl.colors = "black")

# RP116 (CTGA) - BL6
plot_composition_many(CTGA_df,
                      stack.by = "sgid",
                      color.by = "sgid",
                      mouse.genotype = "BL6",
                      tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                      time.points = "3weeks",
                      color.palette = "category20c",
                      highlight.top = "Mga",
                      ctrl.colors = "black")

# RP116 (CTGA) - NSG
plot_composition_many(CTGA_df,
                      stack.by = "sgid",
                      color.by = "sgid",
                      mouse.genotype = "NSG",
                      tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                      time.points = "3weeks",
                      color.palette = "category20c",
                      highlight.top = "Mga",
                      ctrl.colors = "black")

out_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/Mga"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Merge cell lines
merged_df <- bind_rows(
  GGAA_df %>% mutate(cell_line = "RP48_GGAA"),
  CTGA_df %>% mutate(cell_line = "RP116_CTGA")
)

for (geno in c("BL6", "NSG")) {
  p <- plot_composition_many(merged_df,
                             stack.by = "sgid",
                             color.by = "sgid",
                             mouse.genotype = geno,
                             tissues = c("Liver", "Lung", "Brain", "BoneMarrow"),
                             time.points = "3weeks",
                             color.palette = "category20c",
                             highlight.top = "Crebbp",
                             ctrl.colors = "black")

  ggsave(file.path(out_dir, paste0("composition_Crebbp_merged_", geno, ".pdf")),
         p, width = 12, height = 8)
  message("Saved: merged ", geno)
}
# 4. Recompute bootstrap stats with peak_mode and supermet
#    (existing CSVs only have seeding, burden, size_percentile)
# --- RP48 (GGAA) ---
pretran_dataVG <- generate_pretran_data(GGAA_df)

boot_mobaVG_full <- calculate_metastatic_statistics(
  counts.df = GGAA_df,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile", "peak_mode", "supermet"),
  pretran.cell.num = pretran_dataVG$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVG$pretran_cell_num_gene,
  debug = TRUE,
  R = 1000
)
write.csv(boot_mobaVG_full,
          "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVG_stats_full.csv",
          row.names = FALSE)

# --- RP116 (CTGA) ---
pretran_dataVC <- generate_pretran_data(CTGA_df)

boot_mobaVC_full <- calculate_metastatic_statistics(
  counts.df = CTGA_df,
  tag = "sgid",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile", "peak_mode", "supermet"),
  pretran.cell.num = pretran_dataVC$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVC$pretran_cell_num_gene,
  debug = TRUE,
  R = 1000
)
write.csv(boot_mobaVC_full,
          "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVC_stats_full.csv",
          row.names = FALSE)

# 6. Define Mga targets and gene_map
#    UPDATE these once you confirm the exact sgid names from step 2

# --- GGAA (RP48) targets ---
# Mga guides (update if names differ):
mga_targets_GGAA <- mga_sgids_GGAA
# Pick 3 matched controls (using same ones as Crebbp analysis as placeholder):
ctrl_targets_GGAA <- c("sgSafe28-MMU", "sgSafe29-MMU", "sgSafe30-MMU")
my_targets_GGAA <- c(mga_targets_GGAA, ctrl_targets_GGAA)

gene_map_GGAA <- setNames(
  c(rep("Mga", length(mga_targets_GGAA)), rep("Ctrl", length(ctrl_targets_GGAA))),
  my_targets_GGAA
)

# --- CTGA (RP116) targets ---
mga_targets_CTGA <- mga_sgids_CTGA
ctrl_targets_CTGA <- c("sgSafe4-MMU", "sgNT3-MMU", "sgSafe30-MMU")
my_targets_CTGA <- c(mga_targets_CTGA, ctrl_targets_CTGA)

gene_map_CTGA <- setNames(
  c(rep("Mga", length(mga_targets_CTGA)), rep("Ctrl", length(ctrl_targets_CTGA))),
  my_targets_CTGA
)

# 7. Bootstrap timecourse plots - Seeding, Burden, Size
metrics <- c("seeding", "burden", "size_p")
genotypes <- c("BL6", "NSG")

# --- RP48 (GGAA) ---
out_dir_GGAA <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/GGAA/Mga"
if (!dir.exists(out_dir_GGAA)) dir.create(out_dir_GGAA, recursive = TRUE)

for (metric in metrics) {
  for (geno in genotypes) {
    p <- tryCatch({
      plot_bootstrap_timecourse(
        boot_mobaVG_full,
        my_targets_GGAA,
        gene_map_GGAA,
        metric = metric,
        plot.type = "point",
        genotypes = geno,
        tissues = "Liver"
      )
    }, error = function(e) {
      message("Skipping GGAA ", metric, " ", geno, ": ", e$message)
      NULL
    })

    if (!is.null(p)) {
      file_name <- file.path(out_dir_GGAA, paste0("Mga_", metric, "_", geno, "_timecourse.pdf"))
      ggsave(filename = file_name, plot = p, width = 11, height = 6, units = "in")
      message("Saved: ", file_name)
    }
  }
}

# --- RP116 (CTGA) ---
out_dir_CTGA <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/CTGA/Mga"
if (!dir.exists(out_dir_CTGA)) dir.create(out_dir_CTGA, recursive = TRUE)

for (metric in metrics) {
  for (geno in genotypes) {
    p <- tryCatch({
      plot_bootstrap_timecourse(
        boot_mobaVC_full,
        my_targets_CTGA,
        gene_map_CTGA,
        metric = metric,
        plot.type = "point",
        genotypes = geno,
        tissues = "Liver"
      )
    }, error = function(e) {
      message("Skipping CTGA ", metric, " ", geno, ": ", e$message)
      NULL
    })

    if (!is.null(p)) {
      file_name <- file.path(out_dir_CTGA, paste0("Mga_", metric, "_", geno, "_timecourse.pdf"))
      ggsave(filename = file_name, plot = p, width = 11, height = 6, units = "in")
      message("Saved: ", file_name)
    }
  }
}

# 8. Met stat comparison plots (tissue x tissue, genotype x genotype)
tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
genotypes <- c("NSG", "BL6")
stat.types <- c("seeding", "burden", "size")
time_point <- "3weeks"

# --- GGAA ---
output_dir_GGAA <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/Mga/Met_Stat_Plots"
if (!dir.exists(output_dir_GGAA)) dir.create(output_dir_GGAA, recursive = TRUE)

for (stat in stat.types) {
  for (gen in genotypes) {
    # Tissue-tissue comparisons
    for (i in 1:(length(tissues) - 1)) {
      for (j in (i + 1):length(tissues)) {
        p <- tryCatch({
          plot_met_stat_comparison(
            data = boot_mobaVG_full,
            comparison.var = "tissue",
            category1 = tissues[i],
            category2 = tissues[j],
            stat.type = stat,
            plot.by = "sgid",
            use.log2 = FALSE,
            size.percentile = 90,
            filter.vars = list(mouse_genotype = gen, time_point = time_point),
            fdr.threshold = 0.05,
            title = paste0("Mga - ", stat, ": ", tissues[i], " vs ", tissues[j], " (", gen, ")")
          )
        }, error = function(e) {
          message("Skipping GGAA ", stat, " ", gen, " ", tissues[i], " vs ", tissues[j], ": ", e$message)
          NULL
        })
        if (!is.null(p)) {
          ggsave(file.path(output_dir_GGAA,
                           paste0("GGAA_Mga_", stat, "_", gen, "_", tissues[i], "_vs_", tissues[j], "_", time_point, ".pdf")),
                 p, width = 8, height = 8)
        }
      }
    }
  }

  # Genotype-genotype comparisons
  for (tissue in tissues) {
    p <- tryCatch({
      plot_met_stat_comparison(
        data = boot_mobaVG_full,
        comparison.var = "mouse_genotype",
        category1 = "NSG",
        category2 = "BL6",
        stat.type = stat,
        plot.by = "sgid",
        use.log2 = FALSE,
        size.percentile = 90,
        filter.vars = list(tissue = tissue, time_point = time_point),
        fdr.threshold = 0.05,
        title = paste0("Mga - ", stat, ": NSG vs BL6 (", tissue, ")")
      )
    }, error = function(e) {
      message("Skipping GGAA ", stat, " NSG vs BL6 ", tissue, ": ", e$message)
      NULL
    })
    if (!is.null(p)) {
      ggsave(file.path(output_dir_GGAA,
                       paste0("GGAA_Mga_", stat, "_", tissue, "_NSG_vs_BL6_", time_point, ".pdf")),
             p, width = 8, height = 8)
    }
  }
}

# --- CTGA ---
output_dir_CTGA <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/CTGA/Mga/Met_Stat_Plots"
if (!dir.exists(output_dir_CTGA)) dir.create(output_dir_CTGA, recursive = TRUE)

for (stat in stat.types) {
  for (gen in genotypes) {
    for (i in 1:(length(tissues) - 1)) {
      for (j in (i + 1):length(tissues)) {
        p <- tryCatch({
          plot_met_stat_comparison(
            data = boot_mobaVC_full,
            comparison.var = "tissue",
            category1 = tissues[i],
            category2 = tissues[j],
            stat.type = stat,
            plot.by = "sgid",
            use.log2 = FALSE,
            size.percentile = 90,
            filter.vars = list(mouse_genotype = gen, time_point = time_point),
            fdr.threshold = 0.05,
            title = paste0("Mga - ", stat, ": ", tissues[i], " vs ", tissues[j], " (", gen, ")")
          )
        }, error = function(e) {
          message("Skipping CTGA ", stat, " ", gen, " ", tissues[i], " vs ", tissues[j], ": ", e$message)
          NULL
        })
        if (!is.null(p)) {
          ggsave(file.path(output_dir_CTGA,
                           paste0("CTGA_Mga_", stat, "_", gen, "_", tissues[i], "_vs_", tissues[j], "_", time_point, ".pdf")),
                 p, width = 8, height = 8)
        }
      }
    }
  }

  for (tissue in tissues) {
    p <- tryCatch({
      plot_met_stat_comparison(
        data = boot_mobaVC_full,
        comparison.var = "mouse_genotype",
        category1 = "NSG",
        category2 = "BL6",
        stat.type = stat,
        plot.by = "sgid",
        use.log2 = FALSE,
        size.percentile = 90,
        filter.vars = list(tissue = tissue, time_point = time_point),
        fdr.threshold = 0.05,
        title = paste0("Mga - ", stat, ": NSG vs BL6 (", tissue, ")")
      )
    }, error = function(e) {
      message("Skipping CTGA ", stat, " NSG vs BL6 ", tissue, ": ", e$message)
      NULL
    })
    if (!is.null(p)) {
      ggsave(file.path(output_dir_CTGA,
                       paste0("CTGA_Mga_", stat, "_", tissue, "_NSG_vs_BL6_", time_point, ".pdf")),
             p, width = 8, height = 8)
    }
  }
}

# 9. Cell line comparison (RP48 vs RP116) for Mga
output_dir_comp <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/cell_line_comparison/Mga"
if (!dir.exists(output_dir_comp)) dir.create(output_dir_comp, recursive = TRUE)

tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
genotypes <- c("NSG", "BL6")
stat_types <- c("burden", "seeding", "size")
time_points <- list(Liver = c("1week", "3weeks"), Lung = "3weeks", Brain = "3weeks", BoneMarrow = "3weeks")

for (t in tissues) {
  for (g in genotypes) {
    for (s in stat_types) {
      for (tp in time_points[[t]]) {
        subfolder <- file.path(output_dir_comp, t, tp)
        if (!dir.exists(subfolder)) dir.create(subfolder, recursive = TRUE)

        tryCatch({
          p <- plot_df_comparison(
            df1 = boot_mobaVG_full,
            df2 = boot_mobaVC_full,
            use.log2 = FALSE,
            plot.by = "sgid",
            stat.type = s,
            filter.vars = list(tissue = t, mouse_genotype = g, time_point = tp),
            df1.name = "RP48 GGAA",
            df2.name = "RP116 CTGA",
            fdr.threshold = 0.05
          )
          if (!is.null(p)) {
            filename <- paste0("Mga_", tolower(t), "_", tp, "_", s, "_", tolower(g), "_cell_comparison.pdf")
            ggsave(file.path(subfolder, filename), p, width = 8, height = 8)
            message("Saved: ", file.path(subfolder, filename))
          }
        }, error = function(e) {
          message("Failed: Mga ", s, " ", t, " ", g, " ", tp, ": ", e$message)
        })
      }
    }
  }
}

# 10. Chord diagrams & dissemination heatmaps for Mga
dissem_mga_bl6 <- plot_chord_by_genotype(counts_filtered_V, "BL6",
                                         genes = c("Mga", "Ctrl"), project = "Moba V")
plot_dissemination_heatmap(dissem_mga_bl6)

dissem_mga_nsg <- plot_chord_by_genotype(counts_filtered_V, "NSG",
                                         genes = c("Mga", "Ctrl"), project = "Moba V")
plot_dissemination_heatmap(dissem_mga_nsg)

# Save dissemination heatmaps
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/Mga_dissem_heatmap_BL6.pdf",
       dissem_mga_bl6$plot, width = 15, height = 12)
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/chord_plots/Mga_dissem_heatmap_NSG.pdf",
       dissem_mga_nsg$plot, width = 15, height = 12)

# 11. Jitter plots - Mga vs Ctrl cell number comparison
# --- RP48 (GGAA) ---
plot_cell_num_jitter_faceted(
  GGAA_df,
  facet.by = "tissue",
  color.by = "sgid",
  mouse.genotypes = "NSG",
  time.points = c("3weeks"),
  filter.facets = c("Liver", "Lung", "Brain", "BoneMarrow"),
  filter.colors = c("sgNT2-MMU", "Mgammu-sg1"),
  title = "Cell Number Seeding: Mga vs Ctrl (RP48 GGAA, NSG 3weeks)",
  facet.nrow = 2
)

# --- RP116 (CTGA) ---
plot_cell_num_jitter_faceted(
  CTGA_df,
  facet.by = "tissue",
  color.by = "sgid",
  mouse.genotypes = "NSG",
  time.points = c("1week", "3weeks"),
  filter.facets = c("Liver"),
  filter.colors = c("sgNT2-MMU", "Mgammu-sg1"),
  title = "Cell Number Seeding: Mga vs Ctrl (RP116 CTGA, NSG 1-3weeks)",
  facet.nrow = 2
)

#### Mga - Bar Plots, Density Modes, Supermet Jitter ####
# Requires: boot_mobaVG_full, boot_mobaVC_full, GGAA_df, CTGA_df, counts_filtered_V

library(mobaseq)
library(dplyr)
library(ggplot2)

out_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/Mga"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Bootstrap stat bar plots - Seeding, Burden, Size Percentile

genotypes <- c("BL6", "NSG")
tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
time_point <- "3weeks"

# --- Helper: generate bar plots for one cell line ---
generate_bar_plots <- function(boot_data, cell_line_name, out_subdir) {
  if (!dir.exists(out_subdir)) dir.create(out_subdir, recursive = TRUE)

  bar_configs <- list(
    list(
      stat = "seeding",
      value.col = "seeding_log2FC",
      ci.lower.col = "seeding_log2_CI_lower",
      ci.upper.col = "seeding_log2_CI_upper",
      fdr.col = "seeding_FDR",
      ylab = "Met. Seeding (Log2FC)"
    ),
    list(
      stat = "burden",
      value.col = "burden_Log2FC",
      ci.lower.col = "burden_Log2_CI_lower",
      ci.upper.col = "burden_Log2_CI_upper",
      fdr.col = "burden_FDR",
      ylab = "Met. Burden (Log2FC)"
    ),
    list(
      stat = "size_p90",
      value.col = "size_Log2FC_p90",
      ci.lower.col = "size_Log2_CI_lower_p90",
      ci.upper.col = "size_Log2_CI_upper_p90",
      fdr.col = "size_FDR_pval_p90",
      ylab = "Colony Size 90th Percentile (Log2FC)"
    )
  )

  for (geno in genotypes) {
    for (tis in tissues) {
      subset_data <- boot_data %>%
        filter(mouse_genotype == geno, tissue == tis, time_point == time_point)

      if (nrow(subset_data) == 0) next

      for (cfg in bar_configs) {
        # Check columns exist and have non-NA data
        if (!cfg$value.col %in% colnames(subset_data)) next
        plot_subset <- subset_data %>% filter(!is.na(.data[[cfg$value.col]]))
        if (nrow(plot_subset) == 0) next

        p <- tryCatch({
          plot_stat_bar(
            data = plot_subset,
            value.col = cfg$value.col,
            ci.lower.col = cfg$ci.lower.col,
            ci.upper.col = cfg$ci.upper.col,
            fdr.col = cfg$fdr.col,
            fdr.threshold = 0.05,
            log2fc.threshold = 1,
            half.ci = TRUE,
            title = paste0(cell_line_name, " - ", cfg$stat, " - ", geno, " ", tis, " ", time_point),
            ylab = cfg$ylab
          )
        }, error = function(e) {
          message("Skipping ", cell_line_name, " ", cfg$stat, " ", geno, " ", tis, ": ", e$message)
          NULL
        })

        if (!is.null(p)) {
          fname <- file.path(out_subdir,
                             paste0("Mga_bar_", cfg$stat, "_", geno, "_", tis, "_", time_point, ".pdf"))
          ggsave(fname, p, width = 10, height = 6)
          message("Saved: ", fname)
        }
      }
    }
  }
}

# RP48 (GGAA)
generate_bar_plots(boot_mobaVG_full, "RP48 GGAA",
                   file.path(out_dir, "bar_plots/GGAA"))

# RP116 (CTGA)
generate_bar_plots(boot_mobaVC_full, "RP116 CTGA",
                   file.path(out_dir, "bar_plots/CTGA"))

# 2. Density mode plots (peak mode) for Mga vs Ctrl
density_dir <- file.path(out_dir, "density_modes")
if (!dir.exists(density_dir)) dir.create(density_dir, recursive = TRUE)

cell_lines <- list(
  list(df = GGAA_df, name = "RP48_GGAA"),
  list(df = CTGA_df, name = "RP116_CTGA")
)

for (cl in cell_lines) {
  for (geno in genotypes) {
    for (tis in tissues) {
      p <- tryCatch({
        # Filter to specific genotype first
        geno_df <- cl$df %>% filter(mouse_genotype == geno)

        plot_density_modes(
          geno_df,
          genes = c("Ctrl", "Mga"),
          tissues = tis,
          time.points = "3weeks",
          xlim = c(1e1, 1e6)
        )
      }, error = function(e) {
        message("Skipping density ", cl$name, " ", geno, " ", tis, ": ", e$message)
        NULL
      })

      if (!is.null(p)) {
        fname <- file.path(density_dir,
                           paste0("Mga_density_", cl$name, "_", geno, "_", tis, "_3weeks.pdf"))
        ggsave(fname, p, width = 10, height = 6)
        message("Saved: ", fname)
      }
    }
  }
}

# 3. Supermet jitter plots - split by time point, cell line, mouse genotype
supermet_dir <- file.path(out_dir, "supermet_jitter")
if (!dir.exists(supermet_dir)) dir.create(supermet_dir, recursive = TRUE)

time_points_all <- c("1week", "3weeks")

for (cl in cell_lines) {
  for (geno in genotypes) {
    for (tis in tissues) {
      p <- tryCatch({
        geno_df <- cl$df %>% filter(mouse_genotype == geno)

        plot_supermet_jitter(
          geno_df,
          genes = c("Ctrl", "Mga"),
          time.points = time_points_all,
          target.tissue = tis
        )
      }, error = function(e) {
        message("Skipping supermet ", cl$name, " ", geno, " ", tis, ": ", e$message)
        NULL
      })

      if (!is.null(p)) {
        fname <- file.path(supermet_dir,
                           paste0("Mga_supermet_", cl$name, "_", geno, "_", tis, ".pdf"))
        ggsave(fname, p, width = 10, height = 6)
        message("Saved: ", fname)
      }
    }
  }
}

# 4. Combined merged cell line versions (both GGAA + CTGA)
merged_df <- bind_rows(
  GGAA_df %>% mutate(cell_line = "RP48_GGAA"),
  CTGA_df %>% mutate(cell_line = "RP116_CTGA")
)

merged_density_dir <- file.path(out_dir, "density_modes/merged")
merged_supermet_dir <- file.path(out_dir, "supermet_jitter/merged")
if (!dir.exists(merged_density_dir)) dir.create(merged_density_dir, recursive = TRUE)
if (!dir.exists(merged_supermet_dir)) dir.create(merged_supermet_dir, recursive = TRUE)

for (geno in genotypes) {
  geno_merged <- merged_df %>% filter(mouse_genotype == geno)

  for (tis in tissues) {
    # Density
    p_dens <- tryCatch({
      plot_density_modes(
        geno_merged,
        genes = c("Ctrl", "Mga"),
        tissues = tis,
        time.points = "3weeks",
        xlim = c(1e1, 1e6)
      )
    }, error = function(e) {
      message("Skipping merged density ", geno, " ", tis, ": ", e$message)
      NULL
    })

    if (!is.null(p_dens)) {
      ggsave(file.path(merged_density_dir,
                       paste0("Mga_density_merged_", geno, "_", tis, "_3weeks.pdf")),
             p_dens, width = 10, height = 6)
    }

    # Supermet
    p_super <- tryCatch({
      plot_supermet_jitter(
        geno_merged,
        genes = c("Ctrl", "Mga"),
        time.points = time_points_all,
        target.tissue = tis
      )
    }, error = function(e) {
      message("Skipping merged supermet ", geno, " ", tis, ": ", e$message)
      NULL
    })

    if (!is.null(p_super)) {
      ggsave(file.path(merged_supermet_dir,
                       paste0("Mga_supermet_merged_", geno, "_", tis, ".pdf")),
             p_super, width = 10, height = 6)
    }
  }
}

#### Gpatch8 Brain Analysis - Gene-Level Bootstrap & Peak Mode ####
# Gene-level = sgRNAs combined into one "Gpatch8" gene
# Split by cell line: RP48 (GGAA) and RP116 (CTGA)
out_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/Gpatch8"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 0. Verify Gpatch8 sgIDs exist
gpatch8_sgids_GGAA <- GGAA_df %>%
  filter(gene == "Gpatch8") %>%
  distinct(sgid) %>%
  pull(sgid) %>%
  sort()
cat("Gpatch8 sgIDs in GGAA (RP48):\n"); print(gpatch8_sgids_GGAA)

gpatch8_sgids_CTGA <- CTGA_df %>%
  filter(gene == "Gpatch8") %>%
  distinct(sgid) %>%
  pull(sgid) %>%
  sort()
cat("Gpatch8 sgIDs in CTGA (RP116):\n"); print(gpatch8_sgids_CTGA)

# 1. Recompute bootstrap stats at GENE level (tag = "gene")
#    This aggregates all Gpatch8 sgRNAs into one entry

# --- RP48 (GGAA) ---
pretran_dataVG <- generate_pretran_data(GGAA_df)

boot_mobaVG_gene <- calculate_metastatic_statistics(
  counts.df = GGAA_df,
  tag = "gene",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain"),
  statistics = c("seeding", "burden", "size_percentile"),
  pretran.cell.num = pretran_dataVG$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVG$pretran_cell_num_gene,
  debug = TRUE,
  R = 1000
)

# --- RP116 (CTGA) ---
pretran_dataVC <- generate_pretran_data(CTGA_df)

boot_mobaVC_gene <- calculate_metastatic_statistics(
  counts.df = CTGA_df,
  tag = "gene",
  genotypes = c("NSG"),
  tissues = c("Brain"),
  statistics = c("seeding", "burden", "size_percentile"),
  pretran.cell.num = pretran_dataVC$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVC$pretran_cell_num_gene,
  debug = TRUE,
  R = 1000
)

# 2. Bootstrap timecourse plots - Gene level (Gpatch8 vs Ctrl)
#    Seeding, Burden, Size Percentile for Brain
metrics <- c("seeding", "burden", "size_p")
genotypes <- c("BL6", "NSG")
gene_targets <- c("Ctrl", "Gpatch8")

tc_dir <- file.path(out_dir, "timecourse")
if (!dir.exists(tc_dir)) dir.create(tc_dir, recursive = TRUE)

# --- RP48 (GGAA) ---
for (metric in metrics) {
  for (geno in genotypes) {
    p <- tryCatch({
      plot_bootstrap_timecourse(
        results.df = boot_mobaVG_gene,
        targets = gene_targets,
        metric = metric,
        p = 90,
        plot.type = "point",
        genotypes = geno,
        tissues = "Brain",
        title = paste0("RP48 GGAA - Gpatch8 Brain - ", metric, " (", geno, ")")
      )
    }, error = function(e) {
      message("Skipping GGAA gene ", metric, " ", geno, ": ", e$message)
      NULL
    })

    if (!is.null(p)) {
      fname <- file.path(tc_dir, paste0("Gpatch8_gene_", metric, "_GGAA_", geno, "_Brain.pdf"))
      ggsave(fname, plot = p, width = 8, height = 5)
      message("Saved: ", fname)
    }
  }
}

# 3. Peak mode density plots - Gpatch8 vs Ctrl in Brain
density_dir <- file.path(out_dir, "density_modes")
if (!dir.exists(density_dir)) dir.create(density_dir, recursive = TRUE)

cell_lines <- list(
  list(df = GGAA_df, name = "RP48_GGAA"),
  list(df = CTGA_df, name = "RP116_CTGA")
)

for (cl in cell_lines) {
  for (geno in genotypes) {
    geno_df <- cl$df %>% filter(mouse_genotype == geno)

    p <- tryCatch({
      plot_density_modes(
        geno_df,
        genes = c("Ctrl", "Gpatch8"),
        tissues = "Brain",
        time.points = "3weeks",
        xlim = c(1e1, 1e6)
      )
    }, error = function(e) {
      message("Skipping density ", cl$name, " ", geno, ": ", e$message)
      NULL
    })

    if (!is.null(p)) {
      fname <- file.path(density_dir,
                         paste0("Gpatch8_density_", cl$name, "_", geno, "_Brain_3weeks.pdf"))
      ggsave(fname, p, width = 10, height = 6)
      message("Saved: ", fname)
    }
  }
}

#### Gene-Level Met Stat Comparison Plots ####
# All tissues, both cell lines, NSG vs BL6 and tissue vs tissue
#### 1. Recompute gene-level bootstrap stats (all tissues) ####

# --- RP48 (GGAA) ---
pretran_dataVG <- generate_pretran_data(GGAA_df)

boot_mobaVG_gene_all <- calculate_metastatic_statistics(
  counts.df = GGAA_df,
  tag = "gene",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile"),
  pretran.cell.num = pretran_dataVG$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVG$pretran_cell_num_gene,
  debug = TRUE,
  R = 1000
)

write.csv(boot_mobaVG_gene_all,
          "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVG_gene_stats.csv",
          row.names = FALSE)

# --- RP116 (CTGA) ---
pretran_dataVC <- generate_pretran_data(CTGA_df)

boot_mobaVC_gene_all <- calculate_metastatic_statistics(
  counts.df = CTGA_df,
  tag = "gene",
  genotypes = c("BL6", "NSG"),
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile"),
  pretran.cell.num = pretran_dataVC$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataVC$pretran_cell_num_gene,
  debug = TRUE,
  R = 1000
)

write.csv(boot_mobaVC_gene_all,
          "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVC_gene_stats.csv",
          row.names = FALSE)

#### 2. Generate NSG vs BL6 comparison plots per tissue ####

tissues <- c("Liver", "Lung", "Brain", "BoneMarrow")
stat.types <- c("seeding", "burden", "size")
time_points <- list(
  Liver = c("1week", "3weeks"),
  Lung = "3weeks",
  Brain = "3weeks",
  BoneMarrow = "3weeks"
)

cell_lines <- list(
  list(data = boot_mobaVG_gene_all, name = "GGAA", label = "RP48 GGAA"),
  list(data = boot_mobaVC_gene_all, name = "CTGA", label = "RP116 CTGA")
)

for (cl in cell_lines) {
  out_dir <- paste0("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/",
                    cl$name, "/Gene_Met_Stat_Plots")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  for (stat in stat.types) {
    #### Genotype comparisons (NSG vs BL6) per tissue ####
    for (tis in tissues) {
      for (tp in time_points[[tis]]) {
        cat("Generating", cl$name, "gene-level:", stat, tis, "NSG vs BL6", tp, "\n")

        p <- tryCatch({
          plot_met_stat_comparison(
            data = cl$data,
            comparison.var = "mouse_genotype",
            category1 = "NSG",
            category2 = "BL6",
            stat.type = stat,
            plot.by = "gene",
            use.log2 = FALSE,
            size.percentile = 90,
            filter.vars = list(tissue = tis, time_point = tp),
            fdr.threshold = 0.05,
            title = paste0(cl$label, " Gene-Level ", stat,
                           ": NSG vs BL6 (", tis, " ", tp, ")")
          )
        }, error = function(e) {
          message("  Skipping: ", e$message)
          NULL
        })

        if (!is.null(p)) {
          fname <- file.path(out_dir,
                             paste0(cl$name, "_gene_", stat, "_", tis,
                                    "_NSG_vs_BL6_", tp, ".pdf"))
          ggsave(fname, p, width = 8, height = 8)
        }
      }
    }

    #### Tissue-tissue comparisons per genotype ####
    for (gen in c("NSG", "BL6")) {
      for (i in 1:(length(tissues) - 1)) {
        for (j in (i + 1):length(tissues)) {
          # Use the latest shared time point
          shared_tp <- intersect(time_points[[tissues[i]]], time_points[[tissues[j]]])
          for (tp in shared_tp) {
            cat("Generating", cl$name, "gene-level:", stat, gen,
                tissues[i], "vs", tissues[j], tp, "\n")

            p <- tryCatch({
              plot_met_stat_comparison(
                data = cl$data,
                comparison.var = "tissue",
                category1 = tissues[i],
                category2 = tissues[j],
                stat.type = stat,
                plot.by = "gene",
                use.log2 = FALSE,
                size.percentile = 90,
                filter.vars = list(mouse_genotype = gen, time_point = tp),
                fdr.threshold = 0.05,
                title = paste0(cl$label, " Gene-Level ", stat,
                               ": ", tissues[i], " vs ", tissues[j],
                               " (", gen, " ", tp, ")")
              )
            }, error = function(e) {
              message("  Skipping: ", e$message)
              NULL
            })

            if (!is.null(p)) {
              fname <- file.path(out_dir,
                                 paste0(cl$name, "_gene_", stat, "_", gen, "_",
                                        tissues[i], "_vs_", tissues[j], "_", tp, ".pdf"))
              ggsave(fname, p, width = 8, height = 8)
            }
          }
        }
      }
    }
  }
}

#### 3. Cell line comparison (GGAA vs CTGA) at gene level ####

comp_dir <- "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/cell_line_comparison/gene_level"
if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

for (tis in tissues) {
  for (gen in c("NSG", "BL6")) {
    for (s in stat.types) {
      for (tp in time_points[[tis]]) {
        subfolder <- file.path(comp_dir, tis, tp)
        if (!dir.exists(subfolder)) dir.create(subfolder, recursive = TRUE)

        p <- tryCatch({
          plot_df_comparison(
            df1 = boot_mobaVG_gene_all,
            df2 = boot_mobaVC_gene_all,
            stat.type = s,
            plot.by = "gene",
            use.log2 = FALSE,
            filter.vars = list(tissue = tis, mouse_genotype = gen, time_point = tp),
            df1.name = "RP48 GGAA",
            df2.name = "RP116 CTGA",
            fdr.threshold = 0.05
          )
        }, error = function(e) {
          message("Skipping cell line comp ", s, " ", tis, " ", gen, " ", tp, ": ", e$message)
          NULL
        })

        if (!is.null(p)) {
          fname <- file.path(subfolder,
                             paste0("gene_", tolower(tis), "_", tp, "_", s, "_",
                                    tolower(gen), "_cell_comparison.pdf"))
          ggsave(fname, p, width = 8, height = 8)
          message("Saved: ", fname)
        }
      }
    }
  }
}


#### MobaV: Seeding vs preTran at sgRNA level ####

# Load MobaV data — adjust paths as needed
mobav_counts <- load_counts("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/MobaV_FilteredSampleInfo.csv",
                        remove.sgids = "sgDummy",
                        keep.bc.end = c("CTGA", "GGAA"),
                        keep.distance = 0)
mobav_boot_ctga <- read.csv("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVC_stats_full.csv")
mobav_boot_ggaa <- read.csv("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/MobaVal/mobaVG_stats_full.csv")

# Genes to highlight
highlight_genes <- c("Crebbp", "Tsc1", "Tsc2", "Pten", "Mga")

# Function to build one plot
plot_pretran_vs_seeding_sgrna <- function(counts_df, boot_df, bc_end_filter, cell_line_label) {

  # preTran is "Cells" in MobaV
  pretran_df <- counts_df %>%
    filter(tissue == "Cells", sgid != "sgDummy", distance == 0) %>%
    group_by(sgid) %>%
    summarise(pretran_cell_num = sum(cell_num), .groups = "drop")

  # Seeding from bootstrap (already has gene column)
  seeding_df <- boot_df %>%
    filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
    select(sgid, gene, seeding_log2FC)

  # Merge by sgid — use gene from bootstrap
  merged <- inner_join(seeding_df, pretran_df, by = "sgid") %>%
    filter(gene != "Ctrl") %>%
    mutate(log2_pretran = log2(pretran_cell_num + 1))

  # Genes with 2+ sgRNAs
  multi_sg <- merged %>%
    group_by(gene) %>%
    filter(n() >= 2) %>%
    ungroup()

  # Gene centroids for labeling
  centroids <- multi_sg %>%
    group_by(gene) %>%
    summarise(x = mean(log2_pretran), y = mean(seeding_log2FC), .groups = "drop")

  # Within-gene SD
  within_sd <- multi_sg %>%
    group_by(gene) %>%
    summarise(sd = sd(seeding_log2FC, na.rm = TRUE), .groups = "drop")
  cat(sprintf("%s — Median within-gene SD: %.3f\n", cell_line_label, median(within_sd$sd, na.rm = TRUE)))

  # R² of preTran vs seeding
  r2 <- summary(lm(seeding_log2FC ~ log2_pretran, data = merged))$r.squared

  # Available highlight genes
  hl <- highlight_genes[highlight_genes %in% merged$gene]

  ggplot(merged, aes(x = log2_pretran, y = seeding_log2FC)) +
    geom_line(data = multi_sg, aes(group = gene), color = "grey80", linewidth = 0.3) +
    geom_point(alpha = 0.4, size = 2, color = "grey40") +
    geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
    geom_line(data = merged %>% filter(gene %in% hl), aes(group = gene, color = gene), linewidth = 0.8) +
    geom_point(data = merged %>% filter(gene %in% hl), aes(color = gene), size = 3) +
    geom_label_repel(
      data = centroids %>% filter(gene %in% hl),
      aes(x = x, y = y, label = gene, color = gene),
      fontface = "bold", size = 4, box.padding = 1, show.legend = FALSE
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("R² = %.2f", r2), size = 5, fontface = "bold") +
    labs(
      title = cell_line_label,
      x = "log2(preTran Cell Number + 1)",
      y = "Seeding Log2FC",
      color = "Gene"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# Generate both plots
p_ctga <- plot_pretran_vs_seeding_sgrna(mobav_counts, mobav_boot_ctga, "CTGA", "MobaV CTGA — NSG Liver")
p_ggaa <- plot_pretran_vs_seeding_sgrna(mobav_counts, mobav_boot_ggaa, "GGAA", "MobaV GGAA — NSG Liver")

# Side by side
combined <- p_ctga + p_ggaa +
  plot_annotation(
    title = "MobaV: preTran vs Seeding per sgRNA",
    subtitle = "Lines connect sgRNAs targeting the same gene",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12)
    )
  ) +
  plot_layout(guides = "collect")

ggsave("mobav_pretran_vs_seeding_sgRNA_both_lines.pdf", combined, width = 20, height = 10)

# Simple scatter — same style as Moba500 preTran vs seeding plot

plot_pretran_seeding_simple <- function(counts_df, boot_df, cell_line_label) {

  pretran_df <- counts_df %>%
    filter(tissue == "Cells", sgid != "sgDummy", distance == 0) %>%
    group_by(sgid) %>%
    summarise(pretran_cell_num = sum(cell_num), .groups = "drop")

  seeding_df <- boot_df %>%
    filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
    select(sgid, gene, seeding_log2FC)

  merged <- inner_join(seeding_df, pretran_df, by = "sgid") %>%
    filter(gene != "Ctrl") %>%
    mutate(log2_pretran = log2(pretran_cell_num + 1))

  r2 <- summary(lm(seeding_log2FC ~ log2_pretran, data = merged))$r.squared
  crebbp_row <- merged %>% filter(gene == "Crebbp")

  ggplot(merged, aes(x = log2_pretran, y = seeding_log2FC)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", color = "#7C3AED", se = TRUE) +
    geom_point(data = crebbp_row, color = "red", size = 4, shape = 1, stroke = 1.5) +
    ggrepel::geom_label_repel(data = crebbp_row, aes(label = "Crebbp"),
                              color = "red", fontface = "bold", box.padding = 1) +
    ggrepel::geom_text_repel(
      data = merged %>%
        filter(abs(seeding_log2FC) > quantile(abs(seeding_log2FC), 0.95, na.rm = TRUE)),
      aes(label = sgid), size = 3, max.overlaps = 20
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("R² = %.2f", r2), size = 5, fontface = "bold") +
    labs(
      title = paste("preTran vs Seeding —", cell_line_label),
      x = "log2(preTran Cell Number + 1)",
      y = "Seeding Log2FC"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p_ctga_simple <- plot_pretran_seeding_simple(mobav_counts, mobav_boot_ctga, "MobaV CTGA — NSG Liver")
p_ggaa_simple <- plot_pretran_seeding_simple(mobav_counts, mobav_boot_ggaa, "MobaV GGAA — NSG Liver")

combined_simple <- p_ctga_simple + p_ggaa_simple +
  plot_annotation(
    title = "MobaV: preTran Representation vs Seeding (per sgRNA)",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
  )

ggsave("mobav_pretran_vs_seeding_simple_both_lines.pdf", combined_simple, width = 20, height = 10)

plot_pretran_seeding_rainbow <- function(counts_df, boot_df, cell_line_label) {

  pretran_df <- counts_df %>%
    filter(tissue == "Cells", sgid != "sgDummy", distance == 0) %>%
    group_by(sgid) %>%
    summarise(pretran_cell_num = sum(cell_num), .groups = "drop")

  seeding_df <- boot_df %>%
    filter(mouse_genotype == "NSG", tissue == "Liver", time_point == "3weeks") %>%
    select(sgid, gene, seeding_log2FC)

  merged <- inner_join(seeding_df, pretran_df, by = "sgid") %>%
    filter(gene != "Ctrl") %>%
    mutate(log2_pretran = log2(pretran_cell_num + 1))

  r2 <- summary(lm(seeding_log2FC ~ log2_pretran, data = merged))$r.squared
  n_genes <- n_distinct(merged$gene)

  ggplot(merged, aes(x = log2_pretran, y = seeding_log2FC, color = gene)) +
    geom_line(aes(group = gene), linewidth = 0.5, alpha = 0.5) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(aes(group = 1), method = "lm", color = "black", se = TRUE, linetype = "dashed") +
    scale_color_manual(values = rainbow(n_genes)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("R² = %.2f", r2), size = 5, fontface = "bold") +
    labs(
      title = cell_line_label,
      x = "log2(preTran Cell Number + 1)",
      y = "Seeding Log2FC",
      color = "Gene"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
}

p_ctga_rainbow <- plot_pretran_seeding_rainbow(mobav_counts, mobav_boot_ctga, "MobaV CTGA — NSG Liver")
p_ggaa_rainbow <- plot_pretran_seeding_rainbow(mobav_counts, mobav_boot_ggaa, "MobaV GGAA — NSG Liver")

combined_rainbow <- p_ctga_rainbow + p_ggaa_rainbow +
  plot_annotation(
    title = "MobaV: preTran vs Seeding per sgRNA (colored by gene)",
    subtitle = "Lines connect sgRNAs targeting the same gene",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12)
    )
  ) +
  plot_layout(guides = "collect")

ggsave("mobav_pretran_vs_seeding_rainbow_both_lines.pdf", combined_rainbow, width = 22, height = 10)
