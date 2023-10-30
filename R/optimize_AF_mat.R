#' Title
#'
#' @param sce sce file
#' @param sensitivity_threshold vector of sensitivity thresholds for rows (variants) and columns (cells)
#'
#' @return
#' @export
#'
#' @examples
optimize_matrix <- function(sce, sensitivity_threshold = c(0.01,0.001)) {
  iteration <- 0  # Initialize the iteration counter
  previous_matrix <- NULL  # Initialize a variable to store the previous matrix
  
  input_matrix<- data.frame(sce@assays@data$AF_mask)
  colnames(input_matrix)<-rownames(sce@colData)
  rownames(input_matrix)<-names(sce@rowRanges)
  
  while (TRUE) {
    iteration <- iteration + 1  # Increment the iteration counter
    
    # Print the current iteration number and dimensions of the input matrix
    cat("Iteration:", iteration, "\n")
    cat("Input Matrix Dimensions:", dim(input_matrix), "\n")
    cat("Cell sensitivity threshold", sensitivity_threshold[2], "\n")
    cat("Variant sensitivity threshold", sensitivity_threshold[1], "\n")
    
    # Check if the entire matrix consists of TRUE values
    if (all(input_matrix)) {
      cat("Optimization completed. Matrix is all TRUE.\n")
      break
    }
    
    # Calculate the fraction of FALSE entries in each row and column
    row_frac_false <- rowMeans(!input_matrix)
    col_frac_false <- colMeans(!input_matrix)
    
    # Find the row and column with the maximum fraction of FALSE entries
    max_row <- which.max(row_frac_false)
    max_col <- which.max(col_frac_false)
    
    # If the maximum fraction of FALSE entries is zero, we're done
    if (row_frac_false[max_row] == 0 && col_frac_false[max_col] == 0) {
      cat("Optimization completed. No more FALSE entries can be removed.\n")
      break
    }
    
    # Check if other rows and columns have similar means within the current threshold
    similar_rows <- which(abs(row_frac_false - row_frac_false[max_row]) <= sensitivity_threshold[1])
    similar_cols <- which(abs(col_frac_false - col_frac_false[max_col]) <= sensitivity_threshold[2])
    
    # Create a copy of the current input matrix
    previous_matrix <- input_matrix
    
    # Remove the similar rows and columns together
    input_matrix <- input_matrix[-similar_rows, -similar_cols]
    
    # Check if the resulting matrix still has at least one row and one column
    if (nrow(input_matrix) == 0 || ncol(input_matrix) == 0) {
      cat("Optimization stopped. Matrix has 0 rows or 0 columns.\n")
      # Restore the previous matrix
      input_matrix <- previous_matrix
      # Set the sensitivity threshold to 0 and continue the loop
      sensitivity_threshold <- c(0,0)
      cat("Sensitivity threshold has been set to 0. Re-entering the loop.\n")
      next
    }
  }
  AF_complete<-t(sce[rownames(input_matrix),colnames(input_matrix)]@assays@data$AF)
  colnames(AF_complete) <-rownames(input_matrix)
  return(AF_complete)
}
