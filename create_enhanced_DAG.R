#' Create an enhanced DAG plot with better styling
#' 7/23/2025 - added ability to set different variables, use those variable names throughout fxn
#'
#' @param phylopath_result Result from phylopath analysis or path to RDS file
#' @param output_file Path for saving the figure
#' @param title Title for the plot
create_enhanced_dag <- function(phylopath_result = NULL,
                                output_file = "Outputs/PhylopathPlots/enhanced_dag.pdf",
                                title = "Phylogenetic Path Analysis",
                                female_song_var = "FemaleSong_Agg01", coop_breeding_var = "HighConfidence_Coop", territoriality_var = "TerritorialityWeakVsStrong", mass_var = "logMass_AVONET") {
  
  library(ggplot2)
  library(ggrepel)
  
  # Load the result if needed
  if (is.null(phylopath_result) && file.exists("phylopath_full_dataset_result.rds")) {
    phylopath_result <- readRDS("phylopath_full_dataset_result.rds")
  }
  
  # Extract coefficients
  if (!is.null(phylopath_result$result)) {
    result <- phylopath_result$result
    cond_avg <- phylopath::average(result, cut_off = 2, avg_method = "conditional")
    coef_matrix <- cond_avg$coef
  } else {
    stop("Could not extract phylopath result")
  }
  
  # Define node positions and labels - spread out top nodes more
  nodes <- data.frame(
    name = c(mass_var, territoriality_var, 
             coop_breeding_var, female_song_var),
    label = c("Body Mass", "Territoriality", 
              "Cooperative\nBreeding", "Female Song"),
    x = c(5, 5, 1, 9),  # Moved CB left and FS right
    y = c(1, 5, 9, 9),
    stringsAsFactors = FALSE
  )
  
  # Create edge list from coefficient matrix
  edges <- data.frame()
  for (i in 1:nrow(coef_matrix)) {
    for (j in 1:ncol(coef_matrix)) {
      if (coef_matrix[i,j] != 0 && !is.na(coef_matrix[i,j])) {
        from_node <- nodes[nodes$name == rownames(coef_matrix)[i], ]
        to_node <- nodes[nodes$name == colnames(coef_matrix)[j], ]
        
        edges <- rbind(edges, data.frame(
          from = rownames(coef_matrix)[i],
          to = colnames(coef_matrix)[j],
          from_x = from_node$x,
          from_y = from_node$y,
          to_x = to_node$x,
          to_y = to_node$y,
          coefficient = coef_matrix[i,j],
          abs_coef = abs(coef_matrix[i,j]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Adjust end positions to account for node size
  node_radius <- 0.7  # Increased to ensure arrows don't overlap nodes
  edges$dx <- edges$to_x - edges$from_x
  edges$dy <- edges$to_y - edges$from_y
  edges$dist <- sqrt(edges$dx^2 + edges$dy^2)
  
  # For straight edges, adjust normally
  edges$to_x_adj <- edges$to_x - (edges$dx/edges$dist) * node_radius
  edges$to_y_adj <- edges$to_y - (edges$dy/edges$dist) * node_radius
  edges$from_x_adj <- edges$from_x + (edges$dx/edges$dist) * node_radius
  edges$from_y_adj <- edges$from_y + (edges$dy/edges$dist) * node_radius
  
  # For curved edges, we need to account for the curve
  # We'll adjust this after identifying bidirectional edges
  
  # Calculate midpoints
  edges$mid_x <- (edges$from_x + edges$to_x) / 2
  edges$mid_y <- (edges$from_y + edges$to_y) / 2
  
  # Identify bidirectional edges
  edges$edge_id <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "_")
  edge_counts <- table(edges$edge_id)
  edges$is_bidirectional <- edges$edge_id %in% names(edge_counts[edge_counts > 1])
  
  # For bidirectional edges, offset the label positions
  edges$label_offset <- 0
  for (id in unique(edges$edge_id[edges$is_bidirectional])) {
    idx <- which(edges$edge_id == id)
    if (length(idx) == 2) {
      # Offset labels perpendicular to the edge
      # Determine which edge goes which direction to position labels consistently
      if (edges$from_y[idx[1]] == edges$from_y[idx[2]]) {
        # Horizontal bidirectional edge (CB <-> FS)
        # The CB->FS arrow curves upward, so put its label on the upper curve
        # The FS->CB arrow curves downward, so put its label on the lower curve
        if (edges$from[idx[1]] == coop_breeding_var) {
          edges$label_offset[idx[1]] <- 0.6   # On upper curve for CB->FS (0.56)
          edges$label_offset[idx[2]] <- -0.6  # On lower curve for FS->CB (0.39)
        } else {
          edges$label_offset[idx[1]] <- -0.6  # On lower curve for FS->CB (0.39)
          edges$label_offset[idx[2]] <- 0.6   # On upper curve for CB->FS (0.56)
        }
      } else {
        # Other bidirectional edges
        edges$label_offset[idx[1]] <- 0.3
        edges$label_offset[idx[2]] <- -0.3
      }
    }
  }
  
  # Calculate label positions with offset
  edges$perp_x <- -edges$dy / edges$dist
  edges$perp_y <- edges$dx / edges$dist
  
  # For bidirectional edges, ensure perpendicular vectors point the same direction
  for (id in unique(edges$edge_id[edges$is_bidirectional])) {
    idx <- which(edges$edge_id == id)
    if (length(idx) == 2) {
      # If perpendicular vectors point in opposite directions, flip one
      if (edges$perp_y[idx[1]] * edges$perp_y[idx[2]] < 0) {
        edges$perp_y[idx[2]] <- -edges$perp_y[idx[2]]
        edges$perp_x[idx[2]] <- -edges$perp_x[idx[2]]
      }
    }
  }
  
  edges$label_x <- edges$mid_x + edges$perp_x * edges$label_offset
  edges$label_y <- edges$mid_y + edges$perp_y * edges$label_offset
  
  # No need for special handling - the perpendicular offset now works correctly
  
  # Create the plot
  p <- ggplot()
  
  # Separate bidirectional and unidirectional edges
  bidirectional_edges <- edges[edges$is_bidirectional, ]
  unidirectional_edges <- edges[!edges$is_bidirectional, ]
  
  # Draw straight arrows for unidirectional paths
  if (nrow(unidirectional_edges) > 0) {
    p <- p + geom_segment(data = unidirectional_edges,
                          aes(x = from_x_adj, y = from_y_adj, 
                              xend = to_x_adj, yend = to_y_adj,
                              linewidth = abs_coef,
                              #alpha = abs_coef,
                              color = coefficient > 0),
                          arrow = arrow(length = unit(0.4, "cm"),  # Larger arrow heads
                                        type = "closed",
                                        ends = "last"))
  }
  
  # Draw curved arrows for bidirectional paths
  if (nrow(bidirectional_edges) > 0) {
    # Process each unique bidirectional pair
    for (edge_id in unique(bidirectional_edges$edge_id)) {
      pair_edges <- bidirectional_edges[bidirectional_edges$edge_id == edge_id, ]
      if (nrow(pair_edges) == 2) {
        # Find which edge is CB->FS and which is FS->CB
        cb_to_fs_idx <- which(pair_edges$from == coop_breeding_var)
        fs_to_cb_idx <- which(pair_edges$from == female_song_var)
        
        # Draw CB->FS arrow curving upward
        if (length(cb_to_fs_idx) > 0) {
          p <- p + geom_curve(data = pair_edges[cb_to_fs_idx, ],
                              aes(x = from_x_adj, y = from_y_adj, 
                                  xend = to_x_adj, yend = to_y_adj,
                                  linewidth = abs_coef,
                                  #alpha = abs_coef,
                                  color = coefficient > 0),
                              curvature = -0.2, 
                              arrow = arrow(length = unit(0.4, "cm"),
                                            type = "closed",
                                            ends = "last"),
                              ncp = 30)
        }
        
        # Draw FS->CB arrow curving downward
        if (length(fs_to_cb_idx) > 0) {
          p <- p + geom_curve(data = pair_edges[fs_to_cb_idx, ],
                              aes(x = from_x_adj, y = from_y_adj, 
                                  xend = to_x_adj, yend = to_y_adj,
                                  linewidth = abs_coef,
                                  #alpha = abs_coef,
                                  color = coefficient > 0),
                              curvature = -0.2,  
                              arrow = arrow(length = unit(0.4, "cm"), 
                                            type = "closed",
                                            ends = "last"),
                              ncp = 30)
        }
      }
    }
  }
  
  # Continue with the rest of the plot
  p <- p +
    # Add edge labels separately for better control
    # First add labels for unidirectional edges
    geom_label(data = edges[!edges$is_bidirectional, ],
               aes(x = mid_x, y = mid_y, 
                   label = sprintf("%.2f", coefficient)),
               size = 5,
               fill = "white",
               label.padding = unit(0.12, "lines"),
               label.size = 0.4) +
    # Then add labels for bidirectional edges with proper offset
    geom_label(data = edges[edges$is_bidirectional, ],
               aes(x = label_x, y = label_y, 
                   label = sprintf("%.2f", coefficient)),
               size = 5,
               fill = "white",
               label.padding = unit(0.12, "lines"),
               label.size = 0.4) +
    # Draw nodes - larger circles
    geom_point(data = nodes,
               aes(x = x, y = y),
               size = 40,  # Increased from 25
               color = "#2E86AB") +
    # Add node labels - smaller text
    geom_text(data = nodes,
              aes(x = x, y = y, label = label),
              size = 4.9,  # Reduced from 4
              color = "white"#, fontface = "bold"
              ) +
    # Styling
    scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2E86AB"),
                       guide = "none") +
    labs(title = title,
         subtitle = bquote(paste("     n = ", .(phylopath_result$nSpecies), " species | Conditional averaging (", Delta, "CICc < 2)"))) +
    theme_void() +
    theme(plot.title = element_text(size = 19, face = "bold", hjust = 0), #hjust was 0.5
          plot.subtitle = element_text(size = 18, hjust = 0), #hjust was 0.5
          plot.margin = margin(20, 10, 20, 10)) +
    coord_fixed(ratio = 1, xlim = c(0, 10), ylim = c(0, 10))
  
  # Save the plot
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, p, width = 8, height = 8, dpi = 600)
    #ggsave(gsub(".pdf", ".png", output_file), p, 
    #       width = 8, height = 8, dpi = 300)
  }
  
  return(p)
}
