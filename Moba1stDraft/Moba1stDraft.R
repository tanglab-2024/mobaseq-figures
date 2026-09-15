#### Moba1stdraft ####
counts_df_1 <- load_counts("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/Moba1st_FilteredSampleInfo.csv",
                           remove.sgids = "sgDummy",
                           keep.distance = 0)
result <- plot_pretran_histogram(counts_df_1, percentile.threshold = 0.99, filter.data = TRUE, return.list = TRUE)
print(result$plot)
counts_filtered_1 <- result$filtered_df

result <- plot_barcode_seeding(counts_filtered_1, return.data = TRUE)
print(result$plot)

# Filter counts_df to remove these barcodes
counts_filtered_1_2 <- counts_filtered_1 %>%
  filter(sgid %in% c("sgSafe23", "sgSafe28", "sgSafe34", "sgSafe38", "sgSafe6",
                     "sgNf1b-1", "sgNf1b-2", "sgNf1b-3", "sgNf1b-4", "sgNf1b-5")) %>%
  mutate(gene = case_when(
    grepl("Safe", sgid, ignore.case = TRUE) ~ "Ctrl",
    grepl("Nf1b", sgid, ignore.case = TRUE) ~ "Nf1b",
    TRUE ~ NA_character_
  ))

# 2 days
counts_df_f_2days <- counts_filtered_1_2 %>%
  filter((time_point == "2days" | is.na(time_point)) & tissue != "Blood")
p1 <- plot_composition(counts_df_f_2days, stack.by = "gene", color.by = "sgid", title = "2 days")
print(p1)
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/2days_composition_plot.pdf",
       plot = p1, width = 10, height = 7, dpi = 300)

# 1 week
counts_df_f_1wk <- counts_filtered_1_2 %>%
  filter((time_point == "1week" | is.na(time_point)) & tissue != "Blood")
p2 <- plot_composition(counts_df_f_1wk, stack.by = "gene", color.by = "sgid", title = "1 week")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/1wk_composition_plot.pdf",
       plot = p2, width = 10, height = 7, dpi = 300)

# 2 weeks
counts_df_f_2wk <- counts_filtered_1_2 %>%
  filter((time_point == "2weeks" | is.na(time_point)) & tissue != "Blood")
p3 <- plot_composition(counts_df_f_2wk, stack.by = "gene", color.by = "sgid", title = "2 weeks")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/2wk_composition_plot.pdf",
       plot = p3, width = 10, height = 7, dpi = 300)

# 3 weeks
counts_df_f_3wk <- counts_filtered_1_2 %>%
  filter((time_point == "3weeks" | is.na(time_point)) & tissue != "Blood")
p4 <- plot_composition(counts_df_f_3wk, stack.by = "gene", color.by = "sgid", title = "3 weeks")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/3wk_composition_plot.pdf",
       plot = p4, width = 10, height = 7, dpi = 300)

# Plot seeding comparison for 3 weeks for sgid sgSafe-23 vs sgNf1b-4 for liver brain and lung with all pairs on the same plot
# (3 facets, brain liver lung, then each facet contains the jitter plots ranked by time point 2days, 1week, 2weeks, 3weeks)
counts_filtered_12 <- counts_filtered_1_2 %>%
  filter(
    (tissue == "preTran" & cell_num >= 2) |
      (tissue == "Blood" & cell_num >= 1) |
      (tissue == "Brain" & cell_num >= 10) |
      (!tissue %in% c("preTran", "Blood", "Brain") & cell_num >= 100)
  )
plot_tumor_colonies(counts_filtered_12, tissues = c("Liver", "Lung", "Brain"))

plot_cell_num_jitter_faceted(
  counts_filtered_12,
  facet.by = "tissue",
  color.by = "sgid",
  time.points = c("1week", "2weeks", "3weeks"),
  filter.facets = c("Liver", "Lung"),
  filter.colors = c("sgSafe23", "sgNf1b-4"),
  title = "Cell Number Seeding Comparison: sgSafe-23 vs sgNf1b-4",
  facet.nrow = 1
)

# Filter to 3 weeks only
count_df_3weeks <- counts_filtered_12 %>%
  filter(time_point == "3weeks")

# Plot by gene at 3 weeks
plot_boxplot_stats_faceted(
  count_df_3weeks,
  facet.by = "tissue",
  group.by = "sgid",
  metric = "total",
  label.points = F,
  exclude.facets = c("Blood"),
  facet.scales = "free_y",
  facet.nrow = 1,
  title = "Total Metastatic Burden at 3 Weeks by sgid"
)

# Or plot by sgid at 3 weeks
plot_boxplot_stats_faceted(
  count_df_3weeks,
  facet.by = "tissue",
  group.by = "gene",
  metric = "count",
  label.points = T,
  exclude.facets = c("Blood"),
  facet.scales = "free_y",
  facet.nrow = 1,
  title = "Number of Metastatic Seeds at 3 Weeks by gene"
)

# Regenerate composition plots using tissue filtering
# 2 days
counts_df_f_2days_t <- counts_filtered_12 %>%
  filter((time_point == "2days" | is.na(time_point)) & tissue != "Blood")
p1 <- plot_composition(counts_df_f_2days_t, stack.by = "gene", color.by = "sgid", title = "2 days")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/2days_composition_plot_filt.pdf",
       plot = p1, width = 10, height = 7, dpi = 300)

# 1 week
counts_df_f_1wk_t <- counts_filtered_12 %>%
  filter((time_point == "1week" | is.na(time_point)) & tissue != "Blood")
p2 <- plot_composition(counts_df_f_1wk_t, stack.by = "gene", color.by = "sgid", title = "1 week")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/1wk_composition_plot_filt.pdf",
       plot = p2, width = 10, height = 7, dpi = 300)

# 2 weeks
counts_df_f_2wk_t <- counts_filtered_12 %>%
  filter((time_point == "2weeks" | is.na(time_point)) & tissue != "Blood")
p3 <- plot_composition(counts_df_f_2wk_t, stack.by = "gene", color.by = "sgid", title = "2 weeks")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/2wk_composition_plot_filt.pdf",
       plot = p3, width = 10, height = 7, dpi = 300)

# 3 weeks
counts_df_f_3wk_t <- counts_filtered_12 %>%
  filter((time_point == "3weeks" | is.na(time_point)) & tissue != "Blood")
p4 <- plot_composition(counts_df_f_3wk_t, stack.by = "gene", color.by = "sgid", title = "3 weeks")
ggsave("/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2/3wk_composition_plot_filt.pdf",
       plot = p4, width = 10, height = 7, dpi = 300)

pretran_data1 <- generate_pretran_data(counts_filtered_12)

all_results <- boot_all_combinations(counts_filtered_12, pretran_data1$pretran_cell_num,
                                     pretran_data1$pretran_cell_num_gene, R = 1000)
# Create a mapping function
create_gene_map <- function(targets, patterns_map) {
  gene_map <- setNames(rep(NA, length(targets)), targets)

  for (pattern in names(patterns_map)) {
    matches <- grepl(pattern, targets)
    gene_map[matches] <- patterns_map[pattern]
  }

  return(gene_map)
}
# Define targets
sgid_targets <- c("sgSafe23", "sgSafe28", "sgSafe34", "sgSafe38", "sgSafe6",
                  "sgNf1b-1", "sgNf1b-2", "sgNf1b-3", "sgNf1b-4", "sgNf1b-5")

# Use it
sgid_gene_map <- create_gene_map(
  sgid_targets,
  patterns_map = c("sgSafe" = "Ctrl", "sgNf1b" = "Nf1b")
)

sgid_gene_map

gene_targets <- c("Ctrl", "Nf1b")
gene_gene_map <- c("Ctrl" = "Ctrl", "Nf1b" = "Nf1b")

# Step 2: Generate all 4 plots at once
plots <- plot_all_bootstrap_combinations(
  all.results = all_results,
  targets.sgid = sgid_targets,
  targets.gene = gene_targets,
  gene.mapping.sgid = sgid_gene_map,
  gene.mapping.gene = gene_gene_map,
  plot.type = "point",
  exclude.patterns = "2days",
  fdr.threshold = 0.05
)

# Step 3: View individual plots
plots$met_seed_sgid
plots$met_seed_gene
plots$met_burden_sgid
plots$met_burden_gene
save_all_bootstrap_plots(plots, "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2", format = "pdf")
plot_density_modes(counts_filtered_12,
                   genes = c("Ctrl", "Nf1b"),
                   tissues = "Brain",
                   time.points = c("3weeks"),
                   xlim = c(1e1, 1e6))

plot_supermet_jitter(counts_filtered_12,
                     genes = c("Ctrl", "Nf1b"),
                     time.points = c("1week", "2weeks", "3weeks"),
                     target.tissue = "Liver")

size_percentiles <- boot_size_percentile(
  counts.df = counts_filtered_12,
  tag = "both",
  percentiles = c(.1, .2, .3, .4, .5, .6, .7, .8, .9),
  tissues = c("Liver", "Lung", "Brain"),
  R = 1000
)

plot_size_percentiles(
  size_percentiles,
  tissues = c("Liver", "Lung", "Brain"),
  time.points = c("2days", "1week", "2weeks", "3weeks"),
  output.dir = "/Users/ruit/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/V2",
  save.format = "pdf",
  width = 12,
  height = 6
)

# Analyze all targets
plot_dormancy_analysis(
  data = counts_filtered_12,
  tissues = c("Liver"),
  time.points = c("2days", "1week", "2weeks", "3weeks"),
  include.ctrl = T,
  output.dir = "/Users/tanglab/Library/CloudStorage/Box-Box/AX_data_analysis/Mobaseq/Moba1st/test",
  width = 10,
  height = 7
)

results1 <- dormancy_cutoff_by_genotype_tissue_time_v4(
  counts_filtered_12,
  tissues = c("Liver", "Lung", "Brain"),
  time.points = c("3weeks"),
  min.colonies = 1,
  bw.adjust = 2,
  split.by.tissue = FALSE
)
results1$plot_list$NSG_AllTissues

plot_dormancy_cutoffs_v2(
  counts.df = counts_filtered_12,
  cutoffs_df = results1$cutoffs_df,
  value.col = "log2_cell_num",
  time.points = c("3weeks"),
  tissues = c("Liver", "Lung", "Brain"),
  genotypes = unique(counts_filtered_12$mouse_genotype),
  split.by.tissue = FALSE,
  min.colonies = 1
)
geno_cutoffs1 <- results1$cutoffs_df %>%
  select(mouse_genotype, dormant_peak, dormancy_cutoff)
geno_cutoffs1
subset_df <- counts_filtered_12 %>%
  filter(mouse_genotype == "NSG", tissue == "Liver",
         time_point == "3weeks", gene == "Nf1b")

res <- calculate_dormancy_gmm_constrained(
  data = subset_df,
  tag = "gene",
  dormant.mean.constraint = 9.035225,
  valley.mean.constraint = 11.63531,
  tissue = "Liver",
  time.point = "3weeks",
  plot = TRUE
)

# Analyze all targets
results <- plot_dormancy_analysis_v2(
  data = counts_filtered_12,
  tissues = c("Liver", "Lung", "Brain"),
  genotypes = c("NSG"),
  time.points = c("3weeks"),
  include.ctrl = TRUE,
  use.constraints = TRUE,
  tag.type = "gene",
  facet.by = c("tissue", "gene"),
  min.colonies = 100
)

results_df_1C <- calculate_metastatic_statistics(
  counts_filtered_12,
  tissues = c("Liver", "Lung"),
  statistics = c("seeding", "burden", "size_percentile", "supermet", "dormancy_cutoff", "peak_mode"),
  pretran.cell.num = pretran_data1$pretran_cell_num,
  pretran.cell.num.gene = pretran_data1$pretran_cell_num_gene,
  R = 1000,
  dormancy.cutoff = geno_cutoffs1,
  debug = T
)

plot_radar_comparison(results_df_1C, gene = "Nf1b", tissues = c("Liver", "Lung"), time.points = "3weeks", dormancy.metric = "relative_dormancy_escape")

plot_density_modes(count_df_filtered2, tissues = "Liver", time.points = "3weeks", genes = c("Ctrl"), xlim = c(1e2, 1e6))

# Downsample sgSafe data colony number by 10% to see how peak mode changes (within 10% of ground truth) until peak mode is inaccurate
results <- analyze_colony_sampling(
  count_df_filtered2,
  tissue = "Liver",
  time.point = "3weeks",
  gene = "Ctrl",
  n.resamples = 10000,
  percentages = c(seq(90, 10, -10), 5, 1, 0.5),
  seed = 123
)

plot_sampling_results(results, tolerance = 25, crop.outliers = T, percentile.cutoff = 0.3)
