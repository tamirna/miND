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
library(doParallel)
library(foreach)

source("scripts/common.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

fastqLibSize <- function(filePath) {
  return(as.integer(system2("wc", args = c("-l", filePath, " | awk '{print $1}'"), stdout = TRUE)) / 4)
}

if (!file.exists(args$outFile)) {
  sampleSet <- tibble(
    "sample" = character(),
    "group" = character(),
    "miRNA" = numeric(),
    "tRNA" = numeric(),
    "piRNA" = numeric(),
    "rRNA" = numeric(),
    "lncRNA" = numeric(),
    "mRNA" = numeric(),
    "snRNA" = numeric(),
    "snoRNA" = numeric(),
    "yRNA" = numeric(),
    "scRNA" = numeric(),
    "RNA QC spike-in" = numeric(),
    "spike-in calibrator" = numeric(),
    "other RNA species" = numeric(),
    "unclassified genomic" = numeric(),
    "unmapped" = numeric()
  )

  genomeMappingStats <- tibble(
    "sample" = character(),
    "length" = numeric(),
    "stat" = character(),
    "count" = numeric()
  )

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

  inputFiles <- c()

  doParallel::registerDoParallel(cores=args$cores)

  sampleResults <- foreach::foreach(
    i = 5:length(args),
    .combine = append,
    .inorder = FALSE,
    .packages = c()) %dopar% {
      pathSplit <- str_split(args[[i]], "mapping_bowtie_genome/", simplify = TRUE)
      sampleName <- str_replace(pathSplit[1,2], ".map", "")
      pathName <- pathSplit[1,1]

      sampleSheet <- sampleSheetSamples %>% filter(Filename == sampleName)
      sampleID <- sampleSheet$sampleID
      if (length(sampleID) == 0 || is.na(sampleID)) sampleID <- sampleName

      sampleGroup <- sampleSheet[[3]]
      if (length(sampleGroup) == 0) sampleGroup <- "ungrouped"

      fastqLibSizePathName <- paste0(pathName, "fastq_trimmed/", sampleName, ".fastq.libsize")
      sampleSize <- readLines(fastqLibSizePathName)  %>% as.numeric()

      libs <- list(
        "genome" = paste0(pathName, "mapping_bowtie_genome/", sampleName, ".map"),
        "miRNA" = paste0(pathName, "mapping_bowtie_mirna/", sampleName, ".map"),
        "spikeins" = paste0(pathName, "mapping_bowtie_spikeins/", sampleName, ".map"),
        "RNAcentral" = paste0(pathName, "mapping_bowtie_rnacentral/", sampleName, ".map"),
        "cDNA" = paste0(pathName, "mapping_bowtie_cdna/", sampleName, ".map")
      )
      sampleLine <- tibble(sample = as.character(sampleID), group = as.character(sampleGroup))
      sampleMappings <- NA

      for (j in 1:length(libs)) {
        if(file.exists(libs[[j]]) && file.size(libs[[j]]) > 0) {
          libFileContent <- read_delim(libs[[j]],
                                       delim = "\t",
                                       col_names = c(
                                         "ID",
                                         "strand",
                                         "reference",
                                         "offset",
                                         "sequence",
                                         "quality",
                                         "mappings",
                                         "mismatches"
                                       )
          ) %>%
            separate(col = 1, into = c("ID", "reads"), sep = "_x", remove = TRUE, convert = TRUE) %>%
            distinct(ID, .keep_all = TRUE) # remove duplicate entries for the same ID. these occure if the mapping file does not only report one match per sequence, but multiple which might be the case if we analyze a mapping file that was also needed for other mapping statistics. but here we only care about if the read was mapped at all or not, so only one line per read is of interest
          if(names(libs[j]) == "spikeins"){
            libFileContent %<>%
              add_column(lib = "ngsSpikeins")
            if(nrow(libFileContent[grepl("^#uni", libFileContent$reference),]) > 0) {
              libFileContent[grepl("^#uni", libFileContent$reference),]$lib <- "rnaSpikeins"
            }
          } else {
            libFileContent %<>%
              add_column(lib = names(libs[j]))
          }
          if (length(sampleMappings) == 1) {
            sampleMappings <- libFileContent
          } else {
            sampleMappings %<>% bind_rows(libFileContent)
          }
        }
      }

      mappingSummary <- sampleMappings %>%
        group_by(lib) %>%
        summarise(sum = sum(reads)) %>%
        remove_rownames %>%
        column_to_rownames("lib")

      genome <- if (is.na(mappingSummary["genome", "sum"])) 0 else mappingSummary["genome", "sum"]
      rnaSpikeins <- if (is.na(mappingSummary["rnaSpikeins", "sum"])) 0 else mappingSummary["rnaSpikeins", "sum"]
      ngsSpikeins <- if (is.na(mappingSummary["ngsSpikeins", "sum"])) 0 else mappingSummary["ngsSpikeins", "sum"]
      cDNA <- if (is.na(mappingSummary["cDNA", "sum"])) 0 else mappingSummary["cDNA", "sum"]
      miRNA <- if (is.na(mappingSummary["miRNA", "sum"])) 0 else mappingSummary["miRNA", "sum"]
      RNAcentral <- if (is.na(mappingSummary["RNAcentral", "sum"])) 0 else mappingSummary["RNAcentral", "sum"]

      unmapped <- sampleSize - genome - rnaSpikeins - ngsSpikeins
      unknownRNA <- genome - miRNA - cDNA - RNAcentral

      # split up RNAcentral into more detailed statistics
      rnacentralReads <- sampleMappings %>%
        filter(lib == "RNAcentral")

      rnacentral.tRNA <- rnacentralReads %>%
        filter(grepl("tRNA|(transfer RNA)", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.piRNA <- rnacentralReads %>%
        filter(grepl("piRNA|piR-", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.rRNA <- rnacentralReads %>%
        filter(grepl("rRNA|(ribosomal)", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.lncRNA <- rnacentralReads %>%
        filter(grepl("lncRNA|(long non-coding RNA)|( lnc)", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.snRNA <- rnacentralReads %>%
        filter(grepl("snRNA|(spliceosomal RNA)", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.snoRNA <- rnacentralReads %>%
        filter(grepl("snoRNA|(small nucleolar RNA)", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.yRNA <- rnacentralReads %>%
        filter(grepl("y RNA", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()
      rnacentral.scRNA <- rnacentralReads %>%
        filter(grepl("scRNA", reference, ignore.case = TRUE)) %>%
        summarise(sum = sum(reads)) %>%
        as.numeric()

      # calculate "other rna species"
      rnacentral.other <- RNAcentral -
        rnacentral.tRNA -
        rnacentral.piRNA -
        rnacentral.rRNA -
        rnacentral.lncRNA -
        rnacentral.snRNA -
        rnacentral.snoRNA -
        rnacentral.yRNA -
        rnacentral.scRNA

      sampleLine %<>% add_column(
        "miRNA" = miRNA,
        "mRNA" = cDNA,
        "tRNA" = rnacentral.tRNA,
        "piRNA" = rnacentral.piRNA,
        "rRNA" = rnacentral.rRNA,
        "lncRNA" = rnacentral.lncRNA,
        "snRNA" = rnacentral.snRNA,
        "snoRNA" = rnacentral.snoRNA,
        "yRNA" = rnacentral.yRNA,
        "scRNA" = rnacentral.scRNA,
        "RNA QC spike-in" = rnaSpikeins,
        "spike-in calibrator" = ngsSpikeins,
        "other RNA species" = rnacentral.other,
        "unclassified genomic" = unknownRNA,
        "unmapped" = unmapped
      )

      # create histograms of seq length distributions
      genomeMappedSSD <- read_delim(paste0(pathName, "mapping_bowtie_genome/", sampleName, ".mapped.ssd"),
                                    delim = " ",
                                    col_names = c(
                                      "count",
                                      "length"
                                    ),
                                    col_types = "ii"
      ) %>% add_column(state = "mapped")
      genomeUnMappedSSD <- read_delim(paste0(pathName, "mapping_bowtie_genome/", sampleName, ".unmapped.ssd"),
                                      delim = " ",
                                      col_names = c(
                                        "count",
                                        "length"
                                      ),
                                      col_types = "ii"
      ) %>% add_column(state = "unmapped")

      # one stacked bar chart of mapped and unmapped reads  genomeSSD <- full_join(genomeMappedSSD, genomeUnMappedSSD, by = "length") %>%
      genomeSSD <- full_join(genomeMappedSSD, genomeUnMappedSSD, by = "length") %>%
        rename(
          mapped = count.x,
          unmapped = count.y
        ) %>%
        select(length, mapped, unmapped) %>%
        replace_na(list(0, 0, 0)) %>%
        gather(key = "stat", value = "count", -length) %>%
        add_column("sample" = sampleID) %>%
        select(sample, length, stat, count)

      list(list(sampleName = sampleName, sampleLine = sampleLine, genomeMappingStats = genomeSSD))
    }

  # collect parallelized data as with the unparallelized version
  for(i in 1:length(sampleResults)) {
    sampleSet <- bind_rows(sampleSet, sampleResults[[i]]$sampleLine)
    genomeMappingStats %<>% bind_rows(sampleResults[[i]]$genomeMappingStats)
  }

  sampleSetPerc <- sampleSet
  sampleSetMeta <- sampleSet[, 1:2]
  sampleSetPerc <- round(sampleSetPerc[, 3:length(sampleSetPerc)] / rowSums(sampleSetPerc[, 3:length(sampleSetPerc)]) * 100, 2)
  sampleSetPerc <- bind_cols(sampleSetMeta, sampleSetPerc)

  sampleSet %>%
    arrange(sample) %>%
    column_to_rownames("sample") %>%
    rownames_to_column("sample") %>%
    replace(., is.na(.), 0) %>%
    write_delim(args$outFile, delim = "\t")
  sampleSetPerc %>%
    arrange(sample) %>%
    column_to_rownames("sample") %>%
    rownames_to_column("sample") %>%
    replace(., is.na(.), 0) %>%
    write_delim(paste0(file_path_sans_ext(args$outFile), "_perc.", file_ext(args$outFile)), delim = "\t")

  genomeMappingStats %>%
    arrange(sample) %>%
    replace(., is.na(.), 0) %>%
    write_delim(args$outFileGenomeMapping, delim = "\t")
}

sampleSet <- read_delim(args$outFile, delim = "\t")
sampleSetPerc <- read_delim(paste0(file_path_sans_ext(args$outFile), "_perc.", file_ext(args$outFile)), delim = "\t")

sampleSet$sample <- factor(sampleSet$sample)
sampleSet$group <- factor(sampleSet$group)
sampleSetPerc$sample <- factor(sampleSet$sample)
sampleSetPerc$group <- factor(sampleSet$group)

drawStackedBarCharts <- function(data, groupName, perc) {
  if (groupName != FALSE) {
    groupText <- paste0(" - Group ", groupName)
    filename <- paste0(dirname(args$outFile), "/", groupName)
  } else {
    groupText <- ""
    filename <- paste0(file_path_sans_ext(args$outFile))
  }

  if (perc == TRUE) {
    title <- paste0("Reads composition (percentage)", groupText)
    ylabel <- "Percent"
    filename <- paste0(filename, "_perc")
  } else {
    title <- paste0("Reads composition (absolute)", groupText)
    ylabel <- "Read count"
    filename <- paste0(filename)
  }

  data %>% write_delim(paste0(filename, ".dat"), delim = "\t")
}

groupLevels <- lapply(levels(sampleSet$group), function(groupName) {
  data <- gather(sampleSet, "lib", "reads", 3:length(sampleSet), factor_key = TRUE) %>% filter(group == groupName)
  drawStackedBarCharts(data, groupName, FALSE)

  data <- gather(sampleSetPerc, "lib", "reads", 3:length(sampleSet), factor_key = TRUE) %>% filter(group == groupName)
  drawStackedBarCharts(data, groupName, TRUE)
})
data <- gather(sampleSet, "lib", "reads", 3:length(sampleSet), factor_key = TRUE)
drawStackedBarCharts(data, FALSE, FALSE)
data <- gather(sampleSetPerc, "lib", "reads", 3:length(sampleSet), factor_key = TRUE)
drawStackedBarCharts(data, FALSE, TRUE)
