# J. Taroni for ALSF 2019

# options: path to PCL files, path to RNA-seq files, path to QN target, 
# aggregated data directory (output)

#### Libraries -----------------------------------------------------------------

if (!("optparse" %in% installed.packages())) {
  message("Installing optparse library...")
  install.packages("optparse")
}

if (!("preprocessCore" %in% installed.packages())) {
  BiocManager::install("preprocessCore")
}

library(optparse)

#### Functions -----------------------------------------------------------------

`%>%` <- dplyr::`%>%`

# microarray data is "missing" a column name, we need to fill it in
read_in_microarray_data <- function(filename) {
  # read in a given file and set the first column name to "Gene"
  # we're suppressing warnings here because we know the missing column is
  # a problem -- that's why we have this function!
  df <- suppressWarnings(data.table::fread(filename, data.table = FALSE))
  colnames(df)[1] <- "Gene"
  return(df)
}

#### Options -------------------------------------------------------------------

option_list <- list( 
  make_option(opt_str = c("--microarray_path"), 
              type = "character", default = NULL, 
              help = "Path to directory that holds PCL files from refine.bio"),
  make_option(c("--rnaseq_path"),  type = "character", 
              default = NULL, 
              help = "Path to directory that holds lengthScaledTPM.tsv files from refine.bio"),
  make_option(opt_str = c("--qn_target_path"), 
              type = "character", default = NULL, 
              help = "Path to quantile normalization target TSV file from refine.bio"),
  make_option(opt_str = c("--aggregated_path"), 
              type = "character", default = NULL, 
              help = "Path to directory that will hold the output aggregated data")
)

opt <- parse_args(OptionParser(option_list = option_list))
microarray_path <- opt$microarray_path
rnaseq_path <- opt$rnaseq_path
qn_target_path <- opt$qn_target_path

# create the directory to hold the aggregated data if it doesn't exist
aggregated_path <- opt$aggregated_path
dir.create(aggregated_path, recursive = TRUE, showWarnings = FALSE)

#### Microarray data -----------------------------------------------------------

# full path to all the microarray data we obtain
microarray_files <- list.files(microarray_path, full.names = TRUE)

# read in a list of data.frames to be passed to be joined together
df_list <- lapply(microarray_files, read_in_microarray_data)
names(df_list) <- lapply(df_list, function(x) colnames(x)[2])

# any duplicate names need to get dropped
duplicated_names <- names(df_list)[which(duplicated(df_list))]
df_list[duplicated_names] <- NULL

# remove any miRNA arrays in here that will tank our inner join
df_list[grep("miRNA", names(df_list))] <- NULL

# write to file and remove the list of data.frame for memory considerations
microarray_df <- plyr::join_all(df_list, by = "Gene", type = "inner") %>%
  dplyr::mutate(Gene = sub("_at", "", Gene))
readr::write_tsv(microarray_df, file.path(aggregated_path,
                                          "microarray_zebrafish.tsv"))
rm(df_list)

#### RNA-seq data --------------------------------------------------------------

# read in all the TPM tsv files as individual elements (data.frame) of a list
seq_files <- list.files(rnaseq_path, full.names = TRUE)
df_list <- lapply(seq_files, 
                  function(x) data.table::fread(x, data.table = FALSE))

# if all the genes are the same, we can use do.call + cbind which will be much 
# faster than plyr::join_all
seq_genes <- df_list[[1]]$Gene
genes_equal_list <- lapply(df_list, function(x) all.equal(seq_genes, x$Gene))
if (all(unlist(genes_equal_list))) {
  seq_df <- do.call(cbind, df_list)
  # there will many columns named "Gene" -- dplyr hates this
  # we'll snag first gene column and all actual sample data
  seq_df <- seq_df[, c(1, grep("DRR|ERR|SRR", colnames(seq_df)))] 
} else {
  stop("RNA-seq data.frame has unexpected genes")
}

# write to file, remove list
readr::write_tsv(seq_df, file.path(aggregated_path, "seq_zebrafish.tsv"))
rm(df_list)

#### Filter to common genes ----------------------------------------------------

# filter to genes in both technology matrices
common_genes <- intersect(microarray_df$Gene, seq_df$Gene)
common_seq_df <- seq_df %>% 
  dplyr::filter(Gene %in% common_genes)
common_microarray_df <- microarray_df %>% 
  dplyr::filter(Gene %in% common_genes)

# write to file
readr::write_tsv(common_seq_df, file.path(aggregated_path, 
                                          "common_seq_zebrafish.tsv"))
readr::write_tsv(common_microarray_df, 
                 file.path(aggregated_path, "common_microarray_zebrafish.tsv"))

#### Log transform the data ----------------------------------------------------

# early evaluations showed that log2(x + 1) transformation of seq data improved
# imputation

common_seq_matrix <- common_seq_df %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()
log_seq_matrix <- log2(common_seq_matrix + 1)

as.data.frame(log_seq_matrix) %>%
  tibble::rownames_to_column("Gene") %>%
  readr::write_tsv(file.path(aggregated_path, "log_common_seq_zebrafish.tsv"))

#### Combine technologies ------------------------------------------------------

combined_tech_df <- as.data.frame(log_seq_matrix) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::inner_join(common_microarray_df, by = "Gene")
readr::write_tsv(combined_tech_df, path = file.path(aggregated_path, 
                                                    "both_tech_aggregated.tsv"))

#### Quantile normalize --------------------------------------------------------

# read in target distribution file from refine.bio and make a vector
target_df <- 
  readr::read_tsv(qn_target_path, col_names = FALSE)
target <- as.vector(as.matrix(target_df))

# combine the two technologies and quantile normalize
both_mat <- combined_tech_df %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()
qn_mat <- 
  preprocessCore::normalize.quantiles.use.target(both_mat, target = target)
qn_df <- data.frame(Gene = rownames(both_mat), qn_mat)
colnames(qn_df)[2:ncol(qn_df)] <- colnames(both_mat)
readr::write_tsv(qn_df, file.path(aggregated_path, 
                                  "qn_all_zebrafish.tsv"))
