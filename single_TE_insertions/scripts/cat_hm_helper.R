# Updated helper functions with new features

make_logfc_mats <- function(mats) {
  lfc <- lapply(c(2, 1), function(i) log2((mats[[i]] + 1) / (mats[[i + 2]] + 1)))
  names(lfc) <- paste(sub("_piwi", "", names(mats)[c(2, 1)]), "log2FC", sep = "_")
  lfc
}

build_lfc_hts <- function(lfc_mats, te_group, groups_in_metagene = NULL) {
  line_cols <- group_cols[levels(te_group)]
  
  lapply(names(lfc_mats), function(nm) {
    # For metagene: filter groups if specified
    if (!is.null(groups_in_metagene)) {
      metagene_idx <- te_group %in% groups_in_metagene
      metagene_mat <- lfc_mats[[nm]][metagene_idx, , drop = FALSE]
      metagene_group <- te_group[metagene_idx]
      metagene_cols <- group_cols[levels(metagene_group)]
    } else {
      metagene_mat <- lfc_mats[[nm]]
      metagene_group <- te_group
      metagene_cols <- line_cols
    }
    
    EnrichedHeatmap(
      lfc_mats[[nm]],
      name = nm,
      col = col_fun_lfc,
      split = te_group,
      pos_line = FALSE,
      row_order = NULL,
      row_title = NULL,
      row_names_side = "left",
      show_row_names = FALSE,
      column_title = nm,
      column_title_gp = hm_title_gp,
      axis_name = hm_axis_name,
      axis_name_gp = hm_axis_name_gp,
      heatmap_legend_param = list(title = "log2FC(Piwi/White KD)"),
      show_heatmap_legend = grepl("K4", nm),
      top_annotation = HeatmapAnnotation(
        enrich = anno_enriched(
          mat = metagene_mat,
          gp = gpar(col = metagene_cols, lwd = 1.2),
          axis_param = list(side = "left", facing = "inside")
        )
      )
    )
  }) %>% setNames(names(lfc_mats))
}

add_insert_line_all_slices <- function(ht_drawn, hm_name, x = 0.5, lwd = 0.8) {
  rs <- ComplexHeatmap::row_order(ht_drawn)
  for (s in seq_along(rs)) {
    ComplexHeatmap::decorate_heatmap_body(hm_name, slice = s, {
      grid::grid.lines(
        x = grid::unit(c(x, x), "npc"),
        y = grid::unit(c(0, 1), "npc"),
        gp = grid::gpar(lwd = lwd, lty = 2)
      )
    })
  }
}

build_cpm_hts <- function(mats, te_group, show_row_names_ref, show_legend_regex, 
                          ylim_k9, ylim_k4, ylim_k27 = NULL, show_row_names = TRUE, 
                          groups_in_metagene = NULL, include_k27 = FALSE) {
  
  line_cols <- group_cols[levels(te_group)]
  
  # Determine order based on whether K27 is included
  if (include_k27) {
    # Order: K27_white, K27_piwi, K9_white, K9_piwi, K4_white, K4_piwi
    order_names <- names(mats)[grepl("H3K27|H3K9|H3K4", names(mats))]
    # Sort to get white first, then piwi, for each mark
    k27_names <- order_names[grepl("H3K27", order_names)]
    k9_names <- order_names[grepl("H3K9", order_names)]
    k4_names <- order_names[grepl("H3K4", order_names)]
    
    order_names <- c(
      k27_names[grepl("white", k27_names)],
      k27_names[grepl("piwi", k27_names)],
      k9_names[grepl("white", k9_names)],
      k9_names[grepl("piwi", k9_names)],
      k4_names[grepl("white", k4_names)],
      k4_names[grepl("piwi", k4_names)]
    )
  } else {
    # Original order: K9_white, K9_piwi, K4_white, K4_piwi
    order_names <- names(mats)[c(4, 2, 3, 1)]
  }
  
  hts <- lapply(order_names, function(nm) {
    # Determine color function
    mat_col_fun <- if (grepl("H3K4", nm)) {
      col_fun_k4
    } else if (grepl("H3K27", nm)) {
      col_fun_k27
    } else {
      col_fun_k9
    }
    
    # Determine y-axis limits
    ylim_use <- if (grepl("H3K4", nm)) {
      ylim_k4
    } else if (grepl("H3K27", nm)) {
      if (is.null(ylim_k27)) c(-0.15, 2.75) else ylim_k27
    } else {
      ylim_k9
    }
    
    # Determine if row names should be shown for this heatmap
    show_names_this <- show_row_names && (nm == show_row_names_ref)
    
    # For metagene: filter groups if specified
    if (!is.null(groups_in_metagene)) {
      metagene_idx <- te_group %in% groups_in_metagene
      metagene_mat <- mats[[nm]][metagene_idx, , drop = FALSE]
      metagene_group <- te_group[metagene_idx]
      metagene_cols <- group_cols[levels(metagene_group)]
    } else {
      metagene_mat <- mats[[nm]]
      metagene_group <- te_group
      metagene_cols <- line_cols
    }
    
    EnrichedHeatmap(
      mats[[nm]],
      name = nm,
      col = mat_col_fun,
      split = te_group,
      pos_line = FALSE,
      row_order = NULL,
      row_title = NULL,
      row_names_side = "left",
      show_row_names = show_names_this,
      row_names_gp = gpar(fontsize = 6),
      column_title = nm,
      column_title_gp = hm_title_gp,
      axis_name = hm_axis_name,
      axis_name_gp = hm_axis_name_gp,
      heatmap_legend_param = list(title = "log2(CPM+1)"),
      show_heatmap_legend = grepl(show_legend_regex, nm),
      top_annotation = HeatmapAnnotation(
        enrich = anno_enriched(
          mat = metagene_mat,
          gp = gpar(col = metagene_cols, lwd = 1.2),
          ylim = ylim_use,
          axis_param = list(side = "left", facing = "inside")
        )
      )
    )
  })
  
  setNames(hts, order_names)
}

make_norm_mat <- function(bw, gr) {
  m <- normalizeToMatrix(
    bw, gr,
    value_column = "score",
    w = 50,
    extend = c(500, 500),
    target_ratio = 0,
    include_target = FALSE,
    mean_mode = "w0",
    background = NA,
    smooth = FALSE
  )
  m[complete.cases(m), , drop = FALSE]
}

# NEW FUNCTION: Filter TE insertions based on H3K4me3 criteria
# Now supports multiple datasets and returns union/intersection
filter_te_by_h3k4 <- function(mat_list, threshold = 5, max_fc = 1.5, 
                              white_pattern = "white.*H3K4", 
                              piwi_pattern = "piwi.*H3K4") {
  
  # Find the white and piwi H3K4 matrices
  white_names <- names(mat_list)[grepl(white_pattern, names(mat_list))]
  piwi_names <- names(mat_list)[grepl(piwi_pattern, names(mat_list))]
  
  if (length(white_names) == 0 || length(piwi_names) == 0) {
    stop("Could not find white or piwi H3K4 matrices with the given patterns")
  }
  
  # Use the first match if multiple
  white_mat <- mat_list[[white_names[1]]]
  piwi_mat <- mat_list[[piwi_names[1]]]
  
  # Calculate upstream averages (columns 1-10)
  white_upstream_avg <- rowMeans(white_mat[, 1:10, drop = FALSE])
  piwi_upstream_avg <- rowMeans(piwi_mat[, 1:10, drop = FALSE])
  
  # Calculate fold change
  fold_change <- piwi_upstream_avg / (white_upstream_avg + 0.001) # add small value to avoid division by zero
  
  # Identify elements to filter out
  # High H3K4 in white (> threshold) AND small change upon knockdown (FC < max_fc)
  filter_out <- (white_upstream_avg > threshold) & (abs(fold_change - 1) < (max_fc - 1))
  
  filtered_names <- rownames(white_mat)[filter_out]
  
  # Create keep_idx with names
  keep_idx <- !filter_out
  names(keep_idx) <- rownames(white_mat)
  
  cat("  Dataset:", white_names[1], "\n")
  cat("  White H3K4 upstream average > ", threshold, "\n")
  cat("  Fold change (piwi/white) between ", 1/max_fc, " and ", max_fc, "\n")
  cat("  Filtered out ", sum(filter_out), " out of ", length(filter_out), " elements\n")
  
  return(list(
    keep_idx = keep_idx,
    filtered_names = filtered_names,
    white_upstream_avg = white_upstream_avg,
    fold_change = fold_change,
    dataset_name = white_names[1]
  ))
}

# NEW FUNCTION: Combine filtering results from multiple datasets
combine_filter_results <- function(..., mode = "union") {
  filter_results <- list(...)
  
  if (mode == "union") {
    # Filter out elements that fail in ANY dataset
    all_filtered_names <- unique(unlist(lapply(filter_results, function(x) x$filtered_names)))
    cat("\n=== UNION filtering mode ===\n")
    cat("Total unique elements filtered across all datasets:", length(all_filtered_names), "\n")
  } else if (mode == "intersection") {
    # Filter out elements that fail in ALL datasets
    all_filtered_lists <- lapply(filter_results, function(x) x$filtered_names)
    all_filtered_names <- Reduce(intersect, all_filtered_lists)
    cat("\n=== INTERSECTION filtering mode ===\n")
    cat("Elements filtered in ALL datasets:", length(all_filtered_names), "\n")
  } else {
    stop("mode must be 'union' or 'intersection'")
  }
  
  # Get all element names from the first result (they should all have the same elements)
  all_element_names <- names(filter_results[[1]]$white_upstream_avg)
  
  # Create combined keep_idx
  combined_keep_idx <- !all_element_names %in% all_filtered_names
  names(combined_keep_idx) <- all_element_names
  
  return(list(
    keep_idx = combined_keep_idx,
    filtered_names = all_filtered_names,
    individual_results = filter_results,
    mode = mode
  ))
}

prep_te_order_2 <- function(mat_list, ref_name, tesub_gr, keep_idx = NULL) {
  ref_mat <- mat_list[[ref_name]]
  
  # Apply filtering if provided
  if (!is.null(keep_idx)) {
    ref_mat <- ref_mat[keep_idx, , drop = FALSE]
    mat_list <- lapply(mat_list, function(m) m[keep_idx, , drop = FALSE])
  }
  
  te_family <- sub("_[^_]*$", "", rownames(ref_mat))
  
  te_group <- factor(
    case_when(
      te_family %in% c("17.6", "412", "Tirant") ~ "Somatic",
      te_family %in% c("mdg3", "accord", "copia") ~ "Germline early",
      te_family %in% c("3S18", "Burdock", "blood") ~ "Germline late",
      TRUE ~ "Heterogenous"
    ),
    levels = group_levels
  )
  
  te_subgroup <- factor(te_family)
  k9_score <- rowMeans(ref_mat[, 1:10, drop = FALSE])
  
  tesub_meta <- tibble(
    insert = mcols(tesub_gr)$name,
    source = mcols(tesub_gr)$source
  )
  
  ord_df <- tibble(
    idx = seq_len(nrow(ref_mat)),
    Group = te_group,
    Subgroup = te_subgroup,
    k9 = k9_score,
    insert = names(k9_score)
  ) |>
    left_join(tesub_meta, by = "insert") |>
    arrange(Group, Subgroup, desc(k9))
  
  rowo <- ord_df$idx
  
  mats_ord <- lapply(mat_list, \(x) x[rowo, , drop = FALSE])
  te_group_ord <- te_group[rowo]
  te_subgroup_ord <- te_subgroup[rowo]
  
  sub_cols <- setNames(
    colorRampPalette(brewer.pal(8, "Set3"))(nlevels(te_subgroup_ord)),
    levels(te_subgroup_ord)
  )
  
  group_lab <- block_center_labels(te_group_ord)
  subgroup_lab <- block_center_labels(te_subgroup_ord)
  
  row_annot <- rowAnnotation(
    GroupLab = anno_text(group_lab, rot = 90, just = "center", gp = gpar(fontsize = 8)),
    Group = te_group_ord,
    Subgroup = te_subgroup_ord,
    SubgroupLab = anno_text(subgroup_lab, gp = gpar(fontsize = 7)),
    col = list(Group = group_cols, Subgroup = sub_cols),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  list(mats = mats_ord, te_group = te_group_ord, row_annot = row_annot)
}

# Helper function for block labels (assuming this exists in your original code)
# If not, here's a simple implementation:
block_center_labels <- function(groups) {
  rle_groups <- rle(as.character(groups))
  positions <- cumsum(rle_groups$lengths) - rle_groups$lengths/2
  labels <- rep("", length(groups))
  labels[round(positions)] <- rle_groups$values
  labels
}