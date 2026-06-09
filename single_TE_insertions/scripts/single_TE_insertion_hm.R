#### Cut&Tag heatmaps on TE insertions. ####
# publication: nr...
# 2026-06-02
# Artem Ilin (modified by Thomas)
# Data from ngs_te_mapper2 ran on all C&T FASTQ files was used
# Raw data at: "GEO ACCESSION NR"
#---------------------------------------------------------------

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(stringr)
library(ggrastr)
library(writexl)  
library(GenomeInfoDb)

#### General info ####
# MTD = KD across germline development
# BamTOsk = KD from right after the GSC stage
# TOsk = KD from ~stage 2B
# BamMut = KD driven by nanos across early stage germline development

#### setting up folder paths ####
base_dir <- paste0(getwd(), "/single_TE_insertions")
input_url <- "http://ftp.genome-ftp.mbg.au.dk/public/PAN/Data_piDev_paper/cutntag_single_insertions_processed_data/"
output_dir <- paste0(base_dir, "/output")
if (!dir.exists(output_dir)) dir.create(output_dir)
scripts_dir <- paste0(base_dir, "/scripts/")

#### setting up logging and outputting a few info####
.log_path <- file.path(output_dir, "heatmap_creation.log")
.log_con  <- file(.log_path, open = "wt")
sink(.log_con, split = TRUE)
sink(.log_con, type = "message")

#### Load inputs ####
load_from_ftp <- function(filename, url = input_url) {
  temp <- tempfile()
  download.file(paste0(url, filename), destfile = temp)
  load(temp, envir = .GlobalEnv, verbose = T)
}
load_from_ftp("bw_avg_split.RData")
load_from_ftp("mtd_bam_common_all_info.RData")

if (exists("bammutbw") && exists("mtdbw") && exists("mtd_bam_common_all_info")) {
  cat("Input files loaded successfully\n")
} else {
  cat("ERROR: one or more input files failed to load\n")
}

# Functions for heatmap generation
source(paste0(scripts_dir, "cat_hm_helper.R"))


cat(sprintf("=== junction_catalog_creation.R  [%s] ===\n\n", Sys.time()))
mtd_bam_common_df <- as.data.frame(mtd_bam_common_all_info)
cat(sprintf(paste("mtd_bam_common_all_info table was loaded into workspace \n This contains a total of: \n", as.character(length(mtd_bam_common_df$seqnames)), "single TE insertions \n\n")))
smpls <- names(bammutbw)
cat(sprintf(paste("In addition the following bw data frames were loaded into workspace \n bammutbw  mtdbw \n Containing cutntag data for bam mutant knockdown and full ovary mtd knockdown, respectively, for the following samples: \n", 
                  paste(smpls, collapse =" \n "))))
cat(sprintf("\n Each samplerepresent the average CPM+1 values from of 2 biological replicates "))


# Functions for heatmap generation
source(paste0(scripts_dir, "/cat_hm_helper.R"))

###############################################
# USER-CONFIGURABLE PARAMETERS
###############################################

# 1. Toogle whether to plot subset TEs or all
SUB_TEs <- TRUE
if (SUB_TEs == TRUE){
  SOMATIC <- TRUE
  EARLY <- TRUE
  LATE <- TRUE #set which to include in heatmap. final figure has somatic, early and late in main and het only in sub
}

# 2. Filtering for euchromatic reads

EUCH_FILT <- TRUE # set false to include all insertions

if (EUCH_FILT == TRUE) {
  MIN_NEXT_TE <- 1500 # set min distance to neighboring TE insertion
  MIN_TE_COV <- 0 # set min TE coverage, only possible for TEs found in nanopore reference! see MM for more info on this
} else { EUCH_FILT <- FALSE}


# 3. Filtering parameters for H3K4me3
APPLY_H3K4_FILTER <- TRUE # Set to FALSE to skip filtering
H3K4_THRESHOLD <- 2       # H3K4 avg threshold in control samples
H3K4_MAX_FC <- 2          # Minimum fold change (piwi/white) (if high H3K4 in control and also no fold change, discard)

# Filtering mode: "union" = filter if fails in ANY dataset, "intersection" = filter only if fails in ALL datasets
FILTER_MODE <- "union"     # Options: "union" or "intersection"

# 5. Include H3K27me3 data in heatmaps
INCLUDE_H3K27 <- FALSE     # Set to TRUE to include H3K27me3 heatmaps

# 6. set true if export table
export_table <- FALSE

###############################################

cat(sprintf("─────── PARAMETERS SET ───────\n\n"))

cat(paste("Plot subset TE insertions:  ", SUB_TEs, "\n"))

cat(paste("including only euchromatic reads:  ", EUCH_FILT, "\n"))
if(EUCH_FILT == TRUE){
  cat(paste("minimal distance to neighboring TE:   ", MIN_NEXT_TE, "\n"))
  cat(paste("minimal insertion TE coverage over TE consensus :   ", MIN_TE_COV, "\n"))
  cat(paste("apply H3K4me3 filter:   ", APPLY_H3K4_FILTER, "\n"))
  cat(paste("filter mode:   ", FILTER_MODE, "\n"))
  cat(paste("include H3K27:   ", INCLUDE_H3K27, "\n\n\n"))
}

#### Heatmaps based on TE insertions subsets ####
#small correction: 

names(mcols(mtd_bam_common_all_info))[names(mcols(mtd_bam_common_all_info)) == "name"] <- "te"
mcols(mtd_bam_common_all_info)$name <- names(mtd_bam_common_all_info)

mtd_bam <- mtd_bam_common_all_info

# removing NA 
mtd_bam_na <- mtd_bam[(is.na(mtd_bam$lenfrac))]
mtd_bam <- mtd_bam[!(is.na(mtd_bam$lenfrac))]


###### Prepare the subset ######

cat(sprintf("─────── plotting TE families ───────\n\n"))

if (SUB_TEs == TRUE){
  te_keep <- character(0)
  if(SOMATIC == TRUE){
    somatic_tes <- c("17.6", "412", "Tirant")
    te_keep <- c(te_keep, somatic_tes)
    cat(paste("somatic TEs: ", paste(somatic_tes, collapse = " "), "\n"))
  }
  if(EARLY == TRUE){
    early_tes <- c("mdg3", "copia", "accord")
    te_keep <- c(te_keep, early_tes)
    cat(paste("Early TEs: ", paste(early_tes, collapse = " "), "\n"))
  }
  if(LATE == TRUE) {
    late_tes <- c("3S18", "Burdock", "blood")
    te_keep <- c(te_keep, late_tes)
    cat(paste("Late TEs: ", paste(late_tes, collapse = " "), "\n\n"))
  }
  tesub <- mtd_bam[mtd_bam$te %in% te_keep]
  df <- as.data.frame(tesub)
} else {
  tesub <- mtd_bam
  df <- as.data.frame(mtd_bam)
}

# Filtering TE insertions for euchromatic and long insertions (IF IN REFERENCE) size in not known for non-reference. 
if(EUCH_FILT  == TRUE){
  tesub <- tesub[tesub$dist2nearTE > MIN_NEXT_TE] # filter away those that are less than cut distance to neighboring TE
  #tesub <- tesub[!(tesub$source == "reference" & tesub$lenfrac < MIN_TE_COV)] # filter away TEs that are reference, and have a TE coverage that is less than MIN_TE_COV
} else if (SUB_TEs == TRUE) {
  tesub <- mtd_bam_common_premod[mtd_bam_common_premod$te %in% te_keep]
} else {
  tesub <- mtd_bam_common_all_info
}



###### Prepare the matrix ######
te_mat_bammtd_sub <- lapply(c(bammutbw, mtdbw), make_norm_mat, gr = tesub)
cat(paste("Total TE insertions plotted in heatmap: ", as.character(length(as.data.frame(te_mat_bammtd_sub$BamMut_piwi_GSC_H3K27)[,1]))))

if(export_table == TRUE){
  te_mat_bammtd_DF <- as.data.frame(te_mat_bammtd_sub)
  te_mat_bammtd_DF$TE_insertion <- rownames(te_mat_bammtd_DF)
  colnames(te_mat_bammtd_DF) <- sub("MTD","fullOv", colnames(te_mat_bammtd_DF))
  write_xlsx(te_mat_bammtd_DF, paste0(output_dir,"/CutnTag_bam_Table_for_heatmap_single_TE_insertions.xlsx"))
  mtd_bam_common_df$TE_insertion <- rownames(mtd_bam_common_df)
  write_xlsx(mtd_bam_common_df, paste0(output_dir, "/CutnTag_single_TE_insertions_table.xlsx"))
} 

###### Graphic parameters for heatmaps ######
group_levels <- c("Somatic", "Germline early", "Germline late")

group_cols <- c(
  "Somatic"        = "#787878",
  "Germline early" = "#E69F00",
  "Germline late"  = "#56B4E9"
)

col_fun_k9  <- colorRamp2(c(0, 6.00), c("gray95", "#D50000"))
col_fun_k4  <- colorRamp2(c(0, 6.00), c("gray95", "#0072B2"))
col_fun_k27 <- colorRamp2(c(0, 6.00), c("gray95", "#cb79a6"))
col_fun_lfc <- colorRamp2(c(-2.75, 0, 2.75), c("#1042B9", "gray95", "#D55E00"))

hm_axis_name    <- c("-500 bp", "insert", "+500 bp")
hm_axis_name_gp <- gpar(fontsize = 7)
hm_title_gp     <- gpar(fontsize = 7)


###### Data extraction (BEFORE log2 transformation) ######
# Select marks based on H3K27 toggle
if (INCLUDE_H3K27) {
  mtd_TEss <- te_mat_bammtd_sub[grepl("MTD.*H3K(27|9|4)", names(te_mat_bammtd_sub))]
  bam_TEss <- te_mat_bammtd_sub[grepl("BamMut.*H3K(27|9|4)", names(te_mat_bammtd_sub))]
} else {
  mtd_TEss <- te_mat_bammtd_sub[grepl("MTD.*H3K[49]", names(te_mat_bammtd_sub))]
  bam_TEss <- te_mat_bammtd_sub[grepl("BamMut.*H3K[49]", names(te_mat_bammtd_sub))]
}

###### Apply unified H3K4 filtering across both datasets ######
if (APPLY_H3K4_FILTER) {
  cat("\n=== Applying H3K4me3 filtering ===\n")
  
  # Filter MTD (on RAW data, before log2 transformation)
  cat("\nMTD dataset:\n")
  mtd_filter_result <- filter_te_by_h3k4(
    mtd_TEss, 
    threshold = H3K4_THRESHOLD,
    max_fc = H3K4_MAX_FC,
    white_pattern = "MTD_white.*H3K4",
    piwi_pattern = "MTD_piwi.*H3K4"
  )
  
  # Filter BamMut (on RAW data, before log2 transformation)
  cat("\nBamMut dataset:\n")
  bam_filter_result <- filter_te_by_h3k4(
    bam_TEss,
    threshold = H3K4_THRESHOLD,
    max_fc = H3K4_MAX_FC,
    white_pattern = "BamMut_white.*H3K4",
    piwi_pattern = "BamMut_piwi.*H3K4"
  )
  
  # Combine filtering results
  combined_filter <- combine_filter_results(
    mtd_filter_result, 
    bam_filter_result, 
    mode = FILTER_MODE
  )
  # keep combined filtered names as vector
  combined_filtered_TEs <- combined_filter$filtered_names
  
  # Save combined filtered names
  write.table(
    data.frame(filtered_te = combined_filter$filtered_names),
    file = paste0(output_dir, "/combined_filtered_elements.txt"),
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("\nCombined filtered elements saved to: combined_filtered_elements.txt\n")
  
  # NOW apply log2 transformation AFTER filtering decisions are made
  mtd_TEss <- lapply(mtd_TEss, \(x) log2(x + 1))
  bam_TEss <- lapply(bam_TEss, \(x) log2(x + 1))
  
  # Prepare data with COMBINED filtering (same TEs filtered in both)
  mtd_prep <- prep_te_order_2(
    mtd_TEss, 
    "MTD_white_512c_H3K9", 
    tesub,
    keep_idx = combined_filter$keep_idx
  )
  
  bam_prep <- prep_te_order_2(
    bam_TEss,
    "BamMut_white_GSC_H3K9",
    tesub,
    keep_idx = combined_filter$keep_idx
  )
  
  # Create detailed summary showing filtering for each dataset
  combined_summary <- data.frame(
    te_name = names(combined_filter$keep_idx),
    mtd_white_h3k4 = mtd_filter_result$white_upstream_avg,
    mtd_fold_change = mtd_filter_result$fold_change,
    mtd_filtered = !mtd_filter_result$keep_idx,
    bam_white_h3k4 = bam_filter_result$white_upstream_avg,
    bam_fold_change = bam_filter_result$fold_change,
    bam_filtered = !bam_filter_result$keep_idx,
    combined_filtered = !combined_filter$keep_idx
  )
  
  write.csv(
    combined_summary,
    file = paste0(output_dir, "/combined_filter_summary.csv"),
    row.names = FALSE
  )
  
  cat("Combined filter summary saved to: combined_filter_summary.csv\n\n")
  
} else {
  cat("\n=== H3K4me3 filtering disabled ===\n\n")
  
  # Apply log2 transformation
  mtd_TEss <- lapply(mtd_TEss, \(x) log2(x + 1))
  bam_TEss <- lapply(bam_TEss, \(x) log2(x + 1))
  
  # Prepare data without filtering
  mtd_prep <- prep_te_order_2(
    mtd_TEss, 
    "MTD_white_512c_H3K9", 
    tesub
  )
  
  bam_prep <- prep_te_order_2(
    bam_TEss, 
    "BamMut_white_GSC_H3K9", 
    tesub
  )
}

a <- combined_summary %>% filter(combined_filtered == FALSE)



l <- a$te_name
mtd_bam_common_df$tename <- rownames(mtd_bam_common_df)
k <- mtd_bam_common_df %>% filter(tename %in% l)


cat("Number of Reference vs non-reference TEs plotted: \n")
summary(k$source == "nonref")

###### MTD heatmap ######
# Build heatmaps with configured parameters
mtd_cpm_hts <- build_cpm_hts(
  mats = mtd_prep$mats,
  te_group = mtd_prep$te_group,
  show_row_names_ref = "MTD_white_512c_H3K9",
  show_legend_regex = "MTD_white",
  ylim_k9 = c(-0.15, 2.75),
  ylim_k4 = c(-0.15, 2.05),
  ylim_k27 = c(-0.15, 3.5),
  include_k27 = INCLUDE_H3K27
)

# LFC2 heatmap
mtd_lfc_mats <- make_logfc_mats(mtd_prep$mats)
mtd_lfc_hts  <- build_lfc_hts(mtd_lfc_mats, mtd_prep$te_group)

# Create heatmap layout based on K27 inclusion
file_path <- paste0(output_dir, "/te_grouped_k4_k9_mtd.pdf")
if (INCLUDE_H3K27) {
  pdf(file_path, width = 12, height = 10)  # Wider for more heatmaps
  ht_mtd <- draw(
    mtd_prep$row_annot +
      mtd_cpm_hts[[1]] + mtd_cpm_hts[[2]] +  # K27
      mtd_cpm_hts[[3]] + mtd_cpm_hts[[4]] +  # K9
      mtd_cpm_hts[[5]] + mtd_cpm_hts[[6]],   # K4
    merge_legends = FALSE
  )
} else {
  pdf(file_path, width = 8, height = 10)
  ht_mtd <- draw(
    mtd_prep$row_annot +
      mtd_cpm_hts[[1]] + mtd_cpm_hts[[2]] +  # K9
      mtd_cpm_hts[[3]] + mtd_cpm_hts[[4]],   # K4
    merge_legends = FALSE
  )
}
for (nm in names(mtd_cpm_hts)) add_insert_line_all_slices(ht_mtd, nm, x = 0.5)
dev.off()

cat("MTD heatmap saved to:", file_path, "\n\n")

###### BamMut heatmap ######
# Build heatmaps with configured parameters
bam_cpm_hts <- build_cpm_hts(
  mats = bam_prep$mats,
  te_group = bam_prep$te_group,
  show_row_names_ref = "BamMut_white_GSC_H3K9",
  show_legend_regex = "BamMut_white",
  ylim_k9 = c(-0.15, 3.15),
  ylim_k4 = c(-0.15, 2.05),
  ylim_k27 = c(-0.15, 3.5),
  include_k27 = INCLUDE_H3K27
)

# LFC2 heatmap
bam_lfc_mats <- make_logfc_mats(bam_prep$mats)
bam_lfc_hts  <- build_lfc_hts(bam_lfc_mats, bam_prep$te_group)

# Create heatmap layout based on K27 inclusion
file_path <- paste0(output_dir, "/te_grouped_k4_k9_bam.pdf")
if (INCLUDE_H3K27) {
  pdf(file_path, width = 12, height = 10)  # Wider for more heatmaps
  ht_bam <- draw(
    bam_prep$row_annot +
      bam_cpm_hts[[1]] + bam_cpm_hts[[2]] +  # K27
      bam_cpm_hts[[3]] + bam_cpm_hts[[4]] +  # K9
      bam_cpm_hts[[5]] + bam_cpm_hts[[6]],   # K4
    merge_legends = FALSE
  )
} else {
  pdf(file_path, width = 8, height = 10)
  ht_bam <- draw(
    bam_prep$row_annot +
      bam_cpm_hts[[1]] + bam_cpm_hts[[2]] +  # K9
      bam_cpm_hts[[3]] + bam_cpm_hts[[4]],   # K4
    merge_legends = FALSE
  )
}
for (nm in names(bam_cpm_hts)) add_insert_line_all_slices(ht_bam, nm, x = 0.5)
dev.off()

cat("BamMut heatmap saved to: \n", file_path, "\n\n")


#-----------------------------------------------------------------------------------------------------------------------


# making scatterplots:
# one is absolute values of K9 ctrl vs. K4 KD
# the other is fold change KD vs. control. 
# this is only for subset data points, so SUB_TEs must be set to TRUE: 

if(SUB_TEs == TRUE) {

# creating long format of te_mat_bammtd_sub 
mat_to_long <- function(mat, mat_name) {
  parts <- strsplit(mat_name, "_")[[1]]
  as.data.frame(mat) |>
    mutate(insert = rownames(mat)) |>
    pivot_longer(
      cols      = -insert,
      names_to  = "bin",
      values_to = "signal"
    ) |>
    mutate(
      region          = if_else(startsWith(bin, "u"), "upstream", "downstream"),
      experiment_type = parts[1],
      knockdown       = parts[2],
      tissue          = parts[3],
      histone_mark    = parts[4]
    )
}

cat("Converting matrices to long format...\n")
df_long <- bind_rows(
  mapply(mat_to_long, te_mat_bammtd_sub, names(te_mat_bammtd_sub), SIMPLIFY = FALSE)
)


df_long <- df_long |>
  mutate(
    family = sub("_[0-9]+$", "", insert),
    family_highlight = case_when(
      family %in% early_tes   ~ "early",
      family %in% late_tes    ~ "late",
      family %in% somatic_tes ~ "somatic",
      TRUE                          ~ "other"
    )
  )

# Filter away TEs that were filtered away in heatmap if APPLY_H3K4_FILTER is set to TRUE: 
if (APPLY_H3K4_FILTER) { 
  df_long <- df_long %>% filter(!(insert %in% unique(c(combined_filtered_TEs))))
}

highlight_all <- c(early_tes, late_tes, somatic_tes)

source(paste0(scripts_dir,"/aesthetics_functions.R"))
sty <- get_plot_style("paper_smaller")
plot_theme <-   theme(plot.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(color = sty$axis_line_c, linewidth = sty$line_w),
                      axis.ticks = element_line(color = sty$line_c, linewidth = sty$line_w),
                      axis.text.x = element_text(size = sty$axis_txt, family= sty$font),
                      axis.title.x = element_text(size = sty$title_txt, family= sty$font),
                      axis.text.y = element_text(size = sty$axis_txt, family= sty$font), 
                      plot.title = element_text(size = sty$title_txt, family = sty$font), 
                      axis.title.y= element_text(size = sty$title_txt, family = sty$font), 
                      legend.text = element_text(size = sty$legend_txt, family = sty$font), 
                      legend.title = element_text(size = sty$legend_txt, family = sty$font))



# making scatterplots:
# one is absolute values of K9 ctrl vs. K4 KD
# the other is fold change KD vs. control. 
df_scatter <- df_long %>% 
  filter(!(histone_mark == "H3K27")) %>% 
  filter(region == "upstream")

# summarizing signal upstream of insertion site
df_scatter <- df_scatter  %>%
  group_by(across(-c(bin, signal))) %>%
  summarise(
    upstr_signal = sum(signal, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(desc(upstr_signal))

# pivot to wide format for making scatterplots
df_scatter <- df_scatter %>% 
  select(
    insert, tissue, family, family_highlight,
    knockdown, histone_mark,
    upstr_signal
  ) %>%
  mutate(condition = paste(knockdown, histone_mark, sep = "_")) %>%
  tidyr::pivot_wider(
    id_cols     = c(insert, tissue, family, family_highlight),
    names_from  = condition,
    values_from = upstr_signal
  )

scat_plots <- list()
scat_plots_FC <- list()

tissues <- c("GSC","512c")

for (i in 1:length(tissues)){
  plot_df <- df_scatter %>% filter(tissue == tissues[[i]])
  plot_df <- plot_df %>%
    mutate(
      log2FC_K4 = log2((piwi_H3K4 + 10) / (white_H3K4 + 10)),
      log2FC_K9 = log2((piwi_H3K9 + 10) / (white_H3K9 + 10)))
  p <- ggplot(plot_df, aes(x = white_H3K9+1, y = piwi_H3K4+1)) +
    geom_point(shape= sty$dot_shape, size = sty$gdot_size, alpha = sty$dot_alpha) +
    geom_point(data=subset(plot_df, family_highlight == "early"), shape=sty$dot_shape, size=sty$dot_size, alpha=sty$dot_alpha, color="#E69F00") +
    geom_point(data=subset(plot_df, family_highlight == "late"), shape=sty$dot_shape, size=sty$dot_size, alpha=sty$dot_alpha, color="#56B4E9") +
    scale_x_continuous(trans = "log10",limits = c(1,1000)) +
    scale_y_continuous(trans = "log10", limits = c(1,1000)) +
    ggtitle(tissues[[i]]) +
    plot_theme
  scat_plots[[i]] <- p
  file <- paste0(output_dir, "/", tissues[[i]], "scatter_all.pdf")
  ggsave(filename = file, plot = p, width = sty$plot_width, height = sty$plot_height, units="cm", useDingbats=FALSE)
  
  # fold change all
  xlim <- c(-5,5)
  ylim <- c(-5,5)
  p <- ggplot(plot_df, aes(x = log2FC_K9, y = log2FC_K4)) +
    geom_point(shape= sty$dot_shape, size = sty$gdot_size, alpha = sty$dot_alpha) +
    geom_point(data=subset(plot_df, family_highlight == "early"), shape=sty$dot_shape, size=sty$dot_size, alpha=sty$dot_alpha, color="#E69F00") +
    geom_point(data=subset(plot_df, family_highlight == "late"), shape=sty$dot_shape, size=sty$dot_size, alpha=sty$dot_alpha, color="#56B4E9") +
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    ggtitle(tissues[[i]])+
    geom_hline(yintercept=0, linewidth=sty$line_w, color=sty$line_c)+
    geom_vline(xintercept=0, linewidth=sty$line_w, color=sty$line_c)+
    plot_theme
  scat_plots_FC[[i]] <- p
  file <- paste0(output_dir, "/",  tissues[[i]], "scatter_FC_all.pdf")
  ggsave(filename = file, plot = p, width = sty$plot_width, height = sty$plot_height, units="cm", useDingbats=FALSE)
  
  # computing linear regression to get R-squared value
  r2_df <- plot_df %>%
    mutate(
      x = log10(white_H3K9 + 1),
      y = log10(piwi_H3K4 + 1)
    ) %>%
    group_by(family_highlight) %>%
    summarise(
      r2 = summary(lm(y ~ x))$r.squared,
      .groups = "drop"
    )

  ylim <- c(1 ,400)
  xlim <- c(1,400)
  p <- ggplot(plot_df, aes(x = white_H3K9+1, y = piwi_H3K4+1)) +
    geom_point(shape= sty$dot_shape, size = sty$dot_size+0.7, alpha = sty$dot_alpha) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linewidth = sty$line_w, linetype = "dashed") + #this adds linear regression  fitting to the log10-transformed axes. geom_smooth() fits one regression per facet
    facet_wrap(~ family_highlight, ncol = 3) +
    scale_x_continuous(trans = "log10",limits = xlim, oob = scales::oob_keep) +
    scale_y_continuous(trans = "log10", limits = ylim, oob = scales::oob_keep) +
    geom_text(data = r2_df,
              aes(x = Inf, y = Inf, label = paste0("R² = ", round(r2, 2))),
              hjust = 1.1,
              vjust = 1.3,
              inherit.aes = FALSE,
              size = 3
    ) +
    ggtitle(tissues[[i]]) +
    plot_theme +
    annotation_logticks(sides="bl", outside = TRUE, linewidth = sty$line_w, colour = sty$line_c, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) + # nice ticks for log scale axis
    coord_cartesian(clip = "off") # the ticks needs clipping to be disabled otherwise they will be invisible. 
  scat_plots[[i+2]] <- p
  file <- paste0(output_dir, "/",  tissues[[i]], "scatter_sub.pdf")
  ggsave(filename = file, plot = p, width = sty$plot_width*2.5, height = sty$plot_height, units="cm", useDingbats=FALSE)
  
  # fold change
  xlim <- c(-5,5)
  ylim <- c(-5,5)
  
  p <- ggplot(plot_df, aes(x = log2FC_K9, y = log2FC_K4)) +
    geom_point(shape= sty$dot_shape, size = sty$dot_size+0.7, alpha = sty$dot_alpha) +
    facet_wrap(~family_highlight) +
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    ggtitle(tissues[[i]])+
    geom_hline(yintercept=0, linewidth=sty$line_w, color=sty$line_c)+
    geom_vline(xintercept=0, linewidth=sty$line_w, color=sty$line_c)+
    plot_theme
  scat_plots_FC[[i+2]] <- p
  file <- paste0(output_dir, "/", tissues[[i]], "scatter_FC_sub.pdf")
  ggsave(filename = file, plot = p, width = sty$plot_width*2, height = sty$plot_height*2, units="cm", useDingbats=FALSE)
}

} else {
  cat("Scatterplots were not created because subset TEs were not selected!")
}

sink(type = "message")
sink()
close(.log_con)
message("log written to: ", .log_path)
