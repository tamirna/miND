library("tidyr")
library("dplyr")
library("tibble")
library("ggplot2")
library("magrittr")
library("readr")
library("stringr")
library("tools")
library("R.utils")
library("readxl")

source("scripts/common.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE, defaults = list(
  includeSequence = 0
))

outFile <- args$outFile
combinedData <- tibble()
combinedDataRPM <- tibble()

# read sample sheet for sampleID
sampleSheetSamples <- read_excel(args$sampleSheetFile, sheet = "Sample Group Matrix") %>%
  drop_na(c(2, 3)) %>% # remove if column 2 or 3 is empty
  select(1:7) # only select the first 7 columns, not the ones with descriptions

sampleSheetSamples$Filename %<>% gsub("'", "", .) %>%
  gsub(".fastq.gz", "", .) %>%
  gsub(".fastq", "", .) %>%
  gsub(".fq.gz", "", .) %>%
  gsub(".fq", "", .) %>%
  gsub(".bam", "", .)

sampleSheetSamples$sampleID[is.na(sampleSheetSamples$sampleID)] <-
  sampleSheetSamples$Filename[is.na(sampleSheetSamples$sampleID)]

sampleSheetSamples$sampleID <- as.character(sampleSheetSamples$sampleID)

sampleSheetSamples %<>%
  select_if(~!all(is.na(.))) # remove columns with NA only

for (i in 4:length(args)) {
  message(paste0("Loading data from ", args[[i]]))

  pathSplit <- str_split(args[[i]], "analysis_sampleMappingStats/", simplify = TRUE)
  curDataSet <- file_path_sans_ext(pathSplit[1,2])
  pathName <- pathSplit[1,1]

  fileContent <- read_delim(args[[i]], delim = "\t") %>%
    separate(col = 1, into = c("ID", "seq"), sep = "#seq:", remove = TRUE)

  sampleID <- sampleSheetSamples %>% filter(Filename == curDataSet)

  if(nrow(sampleID) != 1) message("ERROR: Got more than one samples with the given filename!")

  sampleID <- sampleID$sampleID
  if (length(sampleID) == 0 || is.na(sampleID)) sampleID <- curDataSet

  message(paste0("sampleID: ", sampleID))

  fileContentReads <- fileContent %>%
    select(1, 2, 3)
  colnames(fileContentReads) <- c("ID", "seq", sampleID)
  fileContentReadsRPM.lib <- fileContent %>%
    select(1, 2, 4)
  fileContentReadsRPM.miRNA <- fileContent %>%
    select(1, 2, 5)
  colnames(fileContentReadsRPM.lib) <- c("ID", "seq", sampleID)
  colnames(fileContentReadsRPM.miRNA) <- c("ID", "seq", sampleID)

  fileContentReadsRPM.lib[, 3] <- round(fileContentReadsRPM.lib[, 3], 2)
  fileContentReadsRPM.miRNA[, 3] <- round(fileContentReadsRPM.miRNA[, 3], 2)
  if (length(combinedData) == 0) {
    # first file to analyze will initialize the combinedData tibble
    combinedData <- fileContentReads
    combinedDataRPM.lib <- fileContentReadsRPM.lib
    combinedDataRPM.miRNA <- fileContentReadsRPM.miRNA
  } else {
    # else we add the fileContent to the existing tibble
    combinedData <- full_join(combinedData, fileContentReads, by = c("ID", "seq"))
    combinedDataRPM.lib <- full_join(combinedDataRPM.lib, fileContentReadsRPM.lib, by = c("ID", "seq"))
    combinedDataRPM.miRNA <- full_join(combinedDataRPM.miRNA, fileContentReadsRPM.miRNA, by = c("ID", "seq"))
  }
}

combinedData <- combinedData[, order(colnames(combinedData))] %>%
  select("ID", "seq", everything())
combinedDataRPM.lib <- combinedDataRPM.lib[, order(colnames(combinedDataRPM.lib))] %>%
  select("ID", "seq", everything())
combinedDataRPM.miRNA <- combinedDataRPM.miRNA[, order(colnames(combinedDataRPM.miRNA))] %>%
  select("ID", "seq", everything())

# remove the seq column if it is not requested
if (args$includeSequence == 0) {
  combinedData %<>% select(-c("seq"))
  combinedDataRPM.lib %<>% select(-c("seq"))
  combinedDataRPM.miRNA %<>% select(-c("seq"))
}

combinedData %>%
  arrange(ID) %>%
  column_to_rownames("ID") %>%
  rownames_to_column("ID") %>%
  replace(., is.na(.), 0) %>%
  write_delim(args$outFile, delim = "\t")
combinedDataRPM.lib %>%
  arrange(ID) %>%
  column_to_rownames("ID") %>%
  rownames_to_column("ID") %>%
  replace(., is.na(.), 0) %>%
  write_delim(paste0(file_path_sans_ext(args$outFile), "_rpm_lib.", file_ext(args$outFile)), delim = "\t")
combinedDataRPM.miRNA %>%
  arrange(ID) %>%
  column_to_rownames("ID") %>%
  rownames_to_column("ID") %>%
  replace(., is.na(.), 0) %>%
  write_delim(paste0(file_path_sans_ext(args$outFile), "_rpm_miRNA.", file_ext(args$outFile)), delim = "\t")

combinedDataLong <- combinedDataRPM.miRNA %>% pivot_longer(-ID, names_to = "Sample", values_to = "count")
combinedDataLong %>%
  filter(count > 1) %>%
  write_delim(paste0(file_path_sans_ext(args$outFile), "_rpm_dist.dat"), delim = "\t")
