# Function to create a customized volcano plot
create_custom_volcano_plot <- function(df, logFC_col = "logFC", pvalue_col = "P.Value", adj_pvalue_col = "adj.P.Val", 
                                       contrast_name = "", fc_cutoff = 0, pvalue_cutoff = 0.05, 
                                       upregulated_color = "#8B1A1A",  # firebrick3 for upregulated genes
                                       downregulated_color = "#00004D",  # navy for downregulated genes
                                       nonsig_color = "grey70",  # Light grey for non-significant genes
                                       plot_width = 6, plot_height = 6, plot_resolution = 600, 
                                       save_plot = FALSE, output_path = "", show_labels = TRUE) {
  
  # Define significance categories (removing logFC filter)
  df$Significance <- ifelse(df[[adj_pvalue_col]] < pvalue_cutoff & df[[logFC_col]] > 0, "Upregulated",
                            ifelse(df[[adj_pvalue_col]] < pvalue_cutoff & df[[logFC_col]] < 0, "Downregulated", 
                                   "Not significant"))
  
  # Base volcano plot with fully filled circles and NO rim
  base_plot <- ggplot(df, aes(x = !!sym(logFC_col), y = -log10(!!sym(adj_pvalue_col)), color = Significance)) + 
    geom_point(size = 2, alpha = 0.6, shape = 16) +  # Fully filled circles
    scale_color_manual(values = c("Upregulated" = upregulated_color, 
                                  "Downregulated" = downregulated_color, 
                                  "Not significant" = nonsig_color)) + 
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "grey") +  # Only horizontal cutoff
    labs(title = contrast_name, x = "Log2 Fold Change", y = "-Log10 FDR") + 
    theme_minimal(base_size = 14) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold")
    )
  
  # Create version with labels if requested
  volcano_plot_with_labels <- base_plot
  if (show_labels) {
    volcano_plot_with_labels <- base_plot +
      geom_text_repel(data = df %>% filter(Significance != "Not significant"), 
                      aes(x = !!sym(logFC_col), y = -log10(!!sym(adj_pvalue_col)), label = GENENAME),
                      size = 3, max.overlaps = 10, segment.color = "black", segment.size = 0.4, 
                      box.padding = 0.2, point.padding = 0.2, color = "black")  # Ensures black text
    
  }
  
  # Save plots if requested
  if (save_plot) {
    ggsave(filename = paste0(output_path, "volcano_plot_", contrast_name, "_No_Labels.png"), 
           plot = base_plot, width = plot_width, height = plot_height, dpi = plot_resolution)
    ggsave(filename = paste0(output_path, "volcano_plot_", contrast_name, "_With_Labels.png"), 
           plot = volcano_plot_with_labels, width = plot_width, height = plot_height, dpi = plot_resolution)
  }
  
  # Return the plots
  return(list(No_Labels = base_plot, With_Labels = volcano_plot_with_labels))
}
