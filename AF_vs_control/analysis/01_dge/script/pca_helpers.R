# pca_helpers.R

# Function to Extract Specific Loadings with Gene Names
extract_specific_loadings_with_names <- function(pca_obj, gene_list, annot) {
  missing_genes <- setdiff(gene_list, rownames(pca_obj$rotation))
  if (length(missing_genes) > 0) {
    warning("Some genes in the list were not found in PCA loadings: ", paste(missing_genes, collapse = ", "))
  }
  
  specific_loadings <- pca_obj$rotation[gene_list, , drop = FALSE]
  specific_loadings_df <- as.data.frame(specific_loadings)
  specific_loadings_df$ENSEMBL <- rownames(specific_loadings_df)  # Add ENSEMBL IDs
  
  specific_loadings_df <- merge(specific_loadings_df, annot[, c("ENSEMBL", "GENENAME")], by = "ENSEMBL", all.x = TRUE)
  return(specific_loadings_df)
}

# Function to Create PCA Plot with Loadings, Optional Labels, and Probability Ellipse
create_publication_pca_plot_with_loadings <- function(pca_df, loadings_df, x, y, color_var, colors, title_text = "", show_labels = TRUE, scaling_factor = 150, add_ellipse = TRUE) {
  # Create PCA plot with sample points
  pca_plot <- ggplot(pca_df, aes_string(x = x, y = y, color = color_var)) +
    geom_point(size = 5, alpha = 0.9, stroke = 1.2) +  
    scale_color_manual(values = colors) +              
    theme_minimal(base_size = 14) +                    
    theme(
      axis.title = element_text(size = 16, face = "bold"),  
      axis.text = element_text(size = 14),                 
      legend.text = element_text(size = 12),               
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "right",                           
      plot.title = element_text(size = 18, hjust = 0.5)    
    ) +
    labs(title = title_text, x = x, y = y)
  
  # Add probability ellipses if requested
  if (add_ellipse) {
    pca_plot <- pca_plot +
      stat_ellipse(aes_string(fill = color_var), type = "norm", alpha = 0.3, geom = "polygon") +
      scale_fill_manual(values = colors)  # Use same colors for ellipse fill
  }
  
  # Add loadings as arrows (vectors) if loadings are present and scaling factor is greater than 0
  if (!is.null(loadings_df) && nrow(loadings_df) > 0) {
    # Apply scaling factor to loadings
    loadings_df[, c(x, y)] <- loadings_df[, c(x, y)] * scaling_factor
    
    pca_plot <- pca_plot +
      geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = !!sym(x), yend = !!sym(y)),
                   arrow = arrow(length = unit(0.3, "cm")), color = "darkgrey", size = 0.8)
    
    if (show_labels) {
      pca_plot <- pca_plot +
        geom_text(data = loadings_df, aes_string(x = x, y = y, label = "GENENAME"),
                  color = "black", size = 4, vjust = -0.5, nudge_x = 0.2, nudge_y = 0.2, fontface = "bold")
    }
  }
  
  # Restore grid lines in the plot by adding back default settings for `panel.grid`
  pca_plot <- pca_plot + theme(panel.grid.major = element_line(color = "grey90"),
                               panel.grid.minor = element_line(color = "grey95"))
  return(pca_plot)
}

# Function to Generate PCA Plots and Save with/without Labels
generate_and_save_pca_plots <- function(pca_df, loadings_df, dimensions, colors, base_path = "../output/", add_ellipse = TRUE) {
  for (dim in dimensions) {
    x_dim <- dim[1]
    y_dim <- dim[2]
    
    pca_plot_with_labels <- create_publication_pca_plot_with_loadings(
      pca_df = pca_df,
      loadings_df = loadings_df,
      x = x_dim,
      y = y_dim,
      color_var = "Region",
      colors = colors,
      title_text = paste("PCA Plot of Healthy Horse Transcriptome with Loadings:", x_dim, "vs", y_dim),
      show_labels = TRUE,
      add_ellipse = add_ellipse  # Add ellipses based on function argument
    )
    
    pca_plot_without_labels <- create_publication_pca_plot_with_loadings(
      pca_df = pca_df,
      loadings_df = loadings_df,
      x = x_dim,
      y = y_dim,
      color_var = "Region",
      colors = colors,
      title_text = paste("PCA Plot of Healthy Horse Transcriptome with Loadings (No Labels):", x_dim, "vs", y_dim),
      show_labels = FALSE,
      add_ellipse = add_ellipse  # Add ellipses based on function argument
    )
    
    # Save each plot separately
    ggsave(filename = paste0(base_path, "PCA_Healthy_Horses_Loadings_", x_dim, "_", y_dim, "_With_Labels.png"),
           plot = pca_plot_with_labels, width = 8, height = 6, dpi = 600)
    ggsave(filename = paste0(base_path, "PCA_Healthy_Horses_Loadings_", x_dim, "_", y_dim, "_No_Labels.png"),
           plot = pca_plot_without_labels, width = 8, height = 6, dpi = 600)
    
    # Print the plots for visualization
    print(pca_plot_with_labels)
    print(pca_plot_without_labels)
  }
}
