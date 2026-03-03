# functions/helper_functions.R
# Utility functions used across scripts

#' Check if required packages are installed
#' @param packages Character vector of package names
check_packages <- function(packages) {
  missing <- packages[!packages %in% installed.packages()[,"Package"]]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "))
  }
}

#' Safe file reading with error handling
#' @param file_path Path to file
safe_read_csv <- function(file_path) {
  tryCatch({
    read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    warning("Could not read file: ", file_path)
    return(NULL)
  })
}

#' Calculate mutation enrichment between two groups
#' @param group1 Vector of mutation frequencies in group 1
#' @param group2 Vector of mutation frequencies in group 2
calculate_enrichment <- function(group1, group2) {
  fold_change <- mean(group1, na.rm = TRUE) / mean(group2, na.rm = TRUE)
  
  # Fisher's exact test
  contingency <- matrix(
    c(sum(group1 > 0), length(group1) - sum(group1 > 0),
      sum(group2 > 0), length(group2) - sum(group2 > 0)),
    nrow = 2
  )
  
  test_result <- fisher.test(contingency)
  
  list(
    fold_change = fold_change,
    p_value = test_result$p.value,
    significant = test_result$p.value < 0.05
  )
}

#' Format p-values for publication
#' @param p_values Numeric vector of p-values
format_pvalues <- function(p_values) {
  case_when(
    p_values < 0.001 ~ "< 0.001",
    p_values < 0.01 ~ sprintf("%.3f", p_values),
    TRUE ~ sprintf("%.2f", p_values)
  )
}

#' Create publication-ready table
#' @param data Data frame
#' @param caption Table caption
pub_table <- function(data, caption = "") {
  knitr::kable(
    data,
    caption = caption,
    format = "html",
    digits = 2,
    align = "c"
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    )
}

#' Print separator line
#' @param char Character to repeat
#' @param width Width of line
#' @param newline Add newline after?
print_separator <- function(char = "=", width = 80, newline = TRUE) {
  cat(strrep(char, width))
  if (newline) cat("\n")
}

#' Print section header
#' @param title Section title
#' @param width Width of separator
print_section <- function(title, width = 80) {
  print_separator("=", width)
  cat(title, "\n")
  print_separator("=", width)
  cat("\n")
}

