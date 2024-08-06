# Load required libraries
library(genekitr)
library(readxl)
library(dplyr)
library(tibble)

# Define a function to perform analysis and write results
process_sheet <- function(sheet_name, file_path, output_dir) {
  # Load the data
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Prepare gene lists
  up <- data$gene_id[data$log_fold_change > 1.0]
  up <- up[!is.na(up)]
  down <- data$gene_id[data$log_fold_change < -1.0]
  down <- down[!is.na(down)]
  
  # Get GO terms for mouse
  mi_gs <- geneset::getGO(org = 'mouse', ont = 'bp')
  
  # Convert geneset to tibble for easier manipulation
  mi_gs$geneset <- as_tibble(mi_gs$geneset)
  
  # Extract the GO term descriptions
  go_descriptions <- mi_gs$geneset_name
  
  # Filter for signaling pathways
  signaling_terms <- go_descriptions %>%
    filter(grepl("signaling pathway", name, ignore.case = TRUE))
  
  # Extract the relevant GO term IDs for signaling pathways
  signaling_go_ids <- signaling_terms$id
  
  # Filter the geneset to include only those related to signaling pathways
  signaling_gs <- mi_gs$geneset %>%
    filter(bp %in% signaling_go_ids)
  
  # Create updated mi_gs with signaling pathways only
  mi_gs_signaling <- list(
    geneset = signaling_gs,
    geneset_name = signaling_terms,
    organism = mi_gs$organism,
    type = mi_gs$type
  )
  
  # Perform enrichment analysis with all GO terms
  up_go_all <- tryCatch({
    genORA(up, geneset = mi_gs)
  }, error = function(e) {
    tibble()  # Return an empty tibble if there's an error
  })
  
  down_go_all <- tryCatch({
    genORA(down, geneset = mi_gs)
  }, error = function(e) {
    tibble()  # Return an empty tibble if there's an error
  })
  
  # Perform enrichment analysis with signaling pathways only
  up_go_sig <- tryCatch({
    genORA(up, geneset = mi_gs_signaling)
  }, error = function(e) {
    tibble()  # Return an empty tibble if there's an error
  })
  
  down_go_sig <- tryCatch({
    genORA(down, geneset = mi_gs_signaling)
  }, error = function(e) {
    tibble()  # Return an empty tibble if there's an error
  })
  
  # Take top 20 results for visualization
  up_go_all <- head(up_go_all, 20)
  down_go_all <- head(down_go_all, 20)
  up_go_sig <- head(up_go_sig, 20)
  down_go_sig <- head(down_go_sig, 20)
  
  # Conditionally remove the geneID_symbol column if it exists
  if ("geneID_symbol" %in% colnames(up_go_all)) {
    up_go_all <- up_go_all %>% select(-geneID_symbol)
  }
  if ("geneID_symbol" %in% colnames(down_go_all)) {
    down_go_all <- down_go_all %>% select(-geneID_symbol)
  }
  if ("geneID_symbol" %in% colnames(up_go_sig)) {
    up_go_sig <- up_go_sig %>% select(-geneID_symbol)
  }
  if ("geneID_symbol" %in% colnames(down_go_sig)) {
    down_go_sig <- down_go_sig %>% select(-geneID_symbol)
  }
  
  # Write results to CSV files only if they are not empty
  if (nrow(up_go_all) > 0) {
    write.csv(up_go_all, file.path(output_dir_all, paste0(sheet_name, '_up_all.csv')), row.names = FALSE)
  }
  if (nrow(down_go_all) > 0) {
    write.csv(down_go_all, file.path(output_dir_all, paste0(sheet_name, '_down_all.csv')), row.names = FALSE)
  }
  if (nrow(up_go_sig) > 0) {
    write.csv(up_go_sig, file.path(output_dir_sig, paste0(sheet_name, '_up_signaling.csv')), row.names = FALSE)
  }
  if (nrow(down_go_sig) > 0) {
    write.csv(down_go_sig, file.path(output_dir_sig, paste0(sheet_name, '_down_signaling.csv')), row.names = FALSE)
  }
}

# Set file paths
file_path <- './DE-pydeseq2/Copy of DMD-DE_results.xlsx'
output_dir_sig <- '/Users/srivalli/Documents/GitHub/hofmann_dmd/DE-pydeseq2/signalling-pathways'
output_dir_all <- '/Users/srivalli/Documents/GitHub/hofmann_dmd/DE-pydeseq2/all-pathways'

# Get all sheet names
sheet_names <- excel_sheets(file_path)

# Process each sheet
for (sheet in sheet_names) {
  process_sheet(sheet, file_path, output_dir)
}



# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get all sheet names
sheet_names <- excel_sheets(file_path)

# Process each sheet
for (sheet in sheet_names) {
  process_sheet(sheet, file_path, output_dir)
}
