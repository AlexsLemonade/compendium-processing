# J. Taroni for ALSF CCDL 2018

#### Functions -----------------------------------------------------------------

`%>%` <- dplyr::`%>%`

get_smashable_s3_links <- function(files_df, pattern = "lengthScaledTPM.tsv") {
  # given a data.frame of files, get the s3 links for files matching the
  # pattern specified with the pattern argument
  # error handling
  required_columns <- c("s3_bucket", "s3_key", "is_smashable", "filename")
  all_req_col <- all(required_columns %in% colnames(files_df))
  if (!is.data.frame(files_df) | !all_req_col) {
    # return an empty vector
    warning("\nThe files data.frame is in an unexpected format\n")
    return(NULL)
  } else {
    relevant_files <- files_df %>%
      dplyr::filter(is_smashable == TRUE,
                    grepl(pattern, filename),
                    !is.na(s3_bucket),
                    !is.na(s3_key)) %>%
      dplyr::mutate(s3_url =  paste("https://s3.amazonaws.com",
                                    s3_bucket, s3_key, sep = "/"))
    return(relevant_files$s3_url)
  }
}

get_organism_s3_links <- function(sample_list, organism,
                                  file_pattern) {
  # get the relevant s3 links for specified organism's smashable files that 
  # match the file_pattern argument
  # for zebrafish organism = "DANIO_RERIO"
  # index of the organism's samples
  organism_index <- which(sample_list$results$organism$name == organism)
  # extract results for organisms samples, only
  organism_samples <- sample_list$results$results[organism_index]
  # extract the s3 links for smashable lengthScaledTPM
  get_pattern_links <- function(y) { 
    lapply(y, function(z) get_smashable_s3_links(z, pattern = file_pattern))
  }
  s3_urls <- unique(unlist(
    lapply(lapply(organism_samples, function(x) x$files), 
           get_pattern_links)
  ))
  return(s3_urls)
}

#### RNA-seq samples available from refine.bio ---------------------------------

# initialize list of sample file urls
seq_url <- "https://api.refine.bio/samples/?source_database=SRA&limit=100"
seq_samples <- jsonlite::fromJSON(seq_url)
all_s3_url <- get_organism_s3_links(seq_samples, "DANIO_RERIO",
                                    "lengthScaledTPM.tsv")

# while there's still another page, get more links!
seq_counter <- 1
while (!is.null(seq_samples$`next`)) {
  message(paste0("RNA-seq page ", seq_counter, "..."))
  seq_samples <- jsonlite::fromJSON(seq_samples$`next`)
  s3_urls <- get_organism_s3_links(seq_samples, "DANIO_RERIO", 
                                   "lengthScaledTPM.tsv")
  # only add the links that we don't already have
  new_urls <- setdiff(s3_urls, all_s3_url)
  all_s3_url <- append(all_s3_url, new_urls)
  seq_counter <- seq_counter + 1 
}

# write the urls to file for use with wget -i
write(all_s3_url, "selected_RNAseq_samples.txt")
number_of_seq_samples <- length(all_s3_url)
rm(all_s3_url, seq_samples, seq_url, new_urls, s3_urls, seq_counter)

#### GEO samples from refine.bio -----------------------------------------------

# initialize list of sample file urls
geo_url <- "https://api.refine.bio/samples/?source_database=GEO&limit=100"
geo_samples <- jsonlite::fromJSON(geo_url)
all_s3_url <- get_organism_s3_links(geo_samples, "DANIO_RERIO",
".PCL")

# while there's still another page, get more links!
geo_counter <- 1
while (!is.null(geo_samples$`next`)) {
  message(paste0("GEO page ", geo_counter, "..."))
  geo_samples <- jsonlite::fromJSON(geo_samples$`next`)
  s3_urls <- get_organism_s3_links(geo_samples, "DANIO_RERIO", 
                                   ".PCL")
  # only add the links that we don't already have
  new_urls <- setdiff(s3_urls, all_s3_url)
  all_s3_url <- append(all_s3_url, new_urls)
  geo_counter <- geo_counter + 1
}

# if there are more GEO samples available than RNA-seq samples (which we expect)
# then randomly select same number as RNA-seq samples
if (length(all_s3_url) > number_of_seq_samples) {
  set.seed(123)
  sampled_s3_url <- sample(all_s3_url, number_of_seq_samples)
} else {
  sampled_s3_url <- all_s3_url
}

# write the random selection of s3 links to file
write(sampled_s3_url, "selected_GEO_samples.txt")
