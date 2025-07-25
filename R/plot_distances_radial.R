#' Radial Distance Visualization with Collision Avoidance
#'
#' Creates a radial plot with automatic label placement to prevent overlaps between
#' nodes and text labels.
#'
#' @param distance_result Data.frame from calculate_nearest_distances()
#' @param reference_type Name of the reference cell type (center node)
#' @param id_col Name of ID column (default: "barcode")
#' @param scale_radius Scaling factor for layout (default: 1)
#' @param show_labels Whether to show distance labels (default: TRUE)
#' @param palette Color palette name (default: "Set2")
#' @param label_padding Radial padding for labels (default: 0.15)
#' @param center_label_expansion Center expansion for labels (default: 1.5)
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom scales rescale
#' @examples
#' plot_radial_distance(distance_results,id_col = "Newbarcode",
#'                      reference_type = "Macrophages",label_padding = 0.3,
#'                      show_labels = TRUE,palette = "Dark2")

plot_radial_distance <- function(distance_result,
                                 reference_type,
                                 id_col = "barcode",
                                 scale_radius = 1,
                                 show_labels = TRUE,
                                 palette = "Set2",
                                 label_padding = 0.15,
                                 center_label_expansion = 1.5) {

  # --- Input Validation ---
  if (!is.data.frame(distance_result)) {
    stop("distance_result must be a data.frame")
  }

  if (!id_col %in% colnames(distance_result)) {
    stop(sprintf("ID column '%s' not found", id_col))
  }

  target_types <- setdiff(colnames(distance_result), id_col)
  if (length(target_types) == 0) {
    stop("No target types found in distance_result")
  }

  # --- Prepare Node Data ---
  mean_dist <- distance_result %>%
    select(-all_of(id_col)) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(everything(),
                 names_to = "target",
                 values_to = "distance") %>%
    mutate(
      radius = rescale(distance, to = c(0.2, 0.8)) * scale_radius
    )

  # --- Calculate Coordinates ---
  n_targets <- nrow(mean_dist)
  angle <- seq(0, 2*pi*(1-1/n_targets), length.out = n_targets)

  nodes <- data.frame(
    type = c(reference_type, mean_dist$target),
    x = c(0, mean_dist$radius * cos(angle)),
    y = c(0, mean_dist$radius * sin(angle)),
    size = c(10, rep(8, n_targets)),
    is_ref = c(TRUE, rep(FALSE, n_targets))
  ) %>%
    mutate(
      label_angle = atan2(y, x),
      hjust = 0.5,  # 中心标签居中对齐
      vjust = 0.5,
      # 为中心标签添加额外偏移
      label_offset = ifelse(is_ref, center_label_expansion, 1),
      label_x = x * (1 + label_padding * label_offset),
      label_y = y * (1 + label_padding * label_offset),
      # 对非中心节点保持原有智能定位
      label_x = ifelse(!is_ref, x * (1 + label_padding + 0.1 * abs(cos(label_angle))), label_x),
      label_y = ifelse(!is_ref, y * (1 + label_padding + 0.1 * abs(sin(label_angle))), label_y),
      # 单独设置中心标签的对齐方式
      hjust = ifelse(is_ref, 0.5, ifelse(x > 0, 0, 1)),
      vjust = ifelse(is_ref, 0.5, ifelse(y > 0, 0, 1))
    )

  edges <- data.frame(
    x_from = 0,
    y_from = 0,
    x_to = mean_dist$radius * cos(angle),
    y_to = mean_dist$radius * sin(angle),
    distance = mean_dist$distance
  ) %>%
    mutate(
      # Position distance labels at 60% of edge length
      label_x = x_to * 0.6,
      label_y = y_to * 0.6,
      # Adjust for quadrant
      label_x = ifelse(x_to > 0, label_x * 1.05, label_x * 0.95),
      label_y = ifelse(y_to > 0, label_y * 1.05, label_y * 0.95)
    )

  # --- Create Plot ---
  p <- ggplot() +
    # Edges
    geom_segment(
      data = edges,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      color = "grey80",
      linewidth = 0.8
    ) +
    # Nodes
    geom_point(
      data = nodes,
      aes(x, y, size = size, fill = type),
      shape = 21,
      color = "white",
      stroke = 1.2
    ) +
    # Node labels with collision avoidance
    geom_text(
      data = nodes,
      aes(x = label_x, y = label_y,
          label = type,
          hjust = hjust,
          vjust = vjust,
          fontface = ifelse(is_ref, "bold", "plain")),
      size = 4.5,
      color = "black"
    ) +
    scale_size_identity() +
    scale_fill_brewer(palette = palette) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1,1,1,1), "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 20))
    ) +
    labs(title = paste("Radial Layout:", reference_type, "as Center"))

  # Add distance labels if requested
  if (show_labels) {
    p <- p +
      geom_label(
        data = edges,
        aes(x = label_x, y = label_y, label = sprintf("%.1f", distance)),
        color = "firebrick",
        fill = alpha("white", 0.8),
        size = 3.2,
        label.padding = unit(0.15, "lines"),
        label.size = 0
      )
  }

  return(p)
}
