% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize_network.R
\name{visualize_spatial_network}
\alias{visualize_spatial_network}
\title{Visualize spatial relationships between cell types}
\usage{
visualize_spatial_network(
  spatial_data,
  sample,
  reference_type,
  target_type,
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  type_col = "Epi_strom",
  color_palette = c(Macrophage = "#90ee90", Epithelial_cells_A = "#377EB8"),
  alpha = 0.7
)
}
\arguments{
\item{spatial_data}{Spatial coordinates data frame}

\item{sample}{Sample name in the spatial transcriptome data}

\item{reference_type}{Reference cell type}

\item{target_type}{Target cell type}

\item{x_col}{Column name for x-coordinates}

\item{y_col}{Column name for y-coordinates}

\item{type_col}{Column name for cell type information}

\item{color_palette}{Named vector of colors for cell types}

\item{alpha}{Transparency level for points and lines}
}
\value{
A ggplot object showing the spatial relationships
}
\description{
Visualize spatial relationships between cell types
}
\examples{
visualize_spatial_network(posi, sample="SP8", reference_type="Macrophage",
                          target_type="Epithelial_cells_A",
                          x_col = "pxl_row_in_fullres",
                          y_col = "pxl_col_in_fullres",
                          type_col = "celltype_ABCDepi",
                          color_palette = c("Macrophage" = "#90ee90",
                         "Epithelial_cells_A" = "#377EB8"))
}
