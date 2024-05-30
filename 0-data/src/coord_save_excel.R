# Load necessary libraries
library(ggplot2)
library(openxlsx)

# Define a function to save ggplot layer data to an Excel file
coord_save_excel <- function(plot, coordPath) {
  # Get the name of the ggplot object
  plot_name <- deparse(substitute(plot))
  
  # Reorder the characters to form the filename (e.g., "e4" -> "4e")
  filename_pattern <- paste0(substr(plot_name, 2, nchar(plot_name)), substr(plot_name, 1, 1))
  
  # Build the ggplot object to extract the data
  plot_data <- ggplot_build(plot)$data
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Write each layer's data to a separate sheet
  for (i in seq_along(plot_data)) {
    layer_name <- paste("Layer", i)
    addWorksheet(wb, layer_name)
    writeData(wb, layer_name, plot_data[[i]])
  }
  
  # Construct the filename based on the pattern
  filename <- paste0(filename_pattern, ".xlsx")
  
  # Save the workbook to a file
  saveWorkbook(wb, file.path(coordPath, filename), overwrite = TRUE)
  
  # Notify the user
  cat("All layers' plot coordinates have been exported to", filename, "\n")
}


