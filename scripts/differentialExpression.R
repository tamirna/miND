library("tidyr")
library("dplyr")
library("tibble")
library("ggplot2")
library("magrittr")
library("stringr")
library("readr")
library("edgeR")
library("genefilter")
library("ggrepel")
library("readxl")

numTopGenes <- 12
NAspecies <- "miRNA"

source("scripts/common.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE, defaults = list(
  includeSequence = 0
))

readsFile <- args$readsFile
sampleSheetFile <- args$sampleSheetFile
outDir <- args$outDir
alpha <- args$alpha %>% as.numeric

# read in the sample contrast sheet
sampleSheetContrast <- read_excel(sampleSheetFile, sheet = "Contrast Selection") %>% drop_na(c(2, 3)) # remove if column 2 or 3 is empty
sampleSheetSamples <- read_excel(sampleSheetFile, sheet = "Sample Group Matrix") %>%
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

# read rawReads from csv file and add sampleID from the xlsx to use later in the reports
rawReads <- read_delim(readsFile, delim = "\t") %>%
  filter(!grepl("^#", ID)) %>% # remove spike ins that all start with a # so that we don't have them in the PCA and other plots
  column_to_rownames("ID")
rawReadsCID <- rawReads

# each line is one contrast, so we need to do DE for each line:
if (nrow(sampleSheetContrast) > 0) {
  for (contrastLineNum in 1:nrow(sampleSheetContrast)) {
    contrastLine <- sampleSheetContrast[contrastLineNum, ]

    print(paste0(contrastLine$GroupA, "_vs_", contrastLine$GroupB))

    contrastOutDir <- paste0(outDir, '/', contrastLine$GroupA, "_vs_", contrastLine$GroupB, ifelse(!is.na(contrastLine$`Grouping Factor`), paste0("_on_", contrastLine$`Grouping Factor`), ''), "/")
    dir.create(contrastOutDir, showWarnings = TRUE, recursive = TRUE)

    # calculate number of groups in the first group contrast string
    nGroupA <- contrastLine$GroupA %>%
      str_count(pattern = "&") %>%
      max() + 1
    nGroupB <- contrastLine$GroupB %>%
      str_count(pattern = "&") %>%
      max() + 1

    contrastLine %<>%
      separate(col = "GroupA", sep = "&", into = paste0("groupA_", 1:nGroupA), remove = FALSE, convert = TRUE) %>%
      separate(col = "GroupB", sep = "&", into = paste0("groupB_", 1:nGroupB), remove = FALSE, convert = TRUE)

    # load all available samples into a tibble and then start applying the filters for the groups in the for loop
    contrastGroupASamples <- contrastGroupBSamples <- sampleSheetSamples

    # get all samples for groupA
    for (groupAGroupNumber in 1:nGroupA) {
      groupSampleDefinition <- contrastLine[[1, paste0("groupA_", groupAGroupNumber)]] %>% strsplit("#")
      column <- groupSampleDefinition[[1]][1]
      value <- groupSampleDefinition[[1]][2]

      contrastGroupASamples %<>% filter(UQ(as.symbol(column)) == value)
    }
    # get all samples for groupB
    for (groupBGroupNumber in 1:nGroupB) {
      groupSampleDefinition <- contrastLine[[1, paste0("groupB_", groupBGroupNumber)]] %>% strsplit("#")
      column <- groupSampleDefinition[[1]][1]
      value <- groupSampleDefinition[[1]][2]

      if(!is.na(column) && column == "everything" && groupBGroupNumber == 1) {
        # we want everyhing in group B that is not yet in group A
        contrastGroupBSamples %<>% filter(!(sampleID %in% contrastGroupASamples$sampleID))
      } else {
        contrastGroupBSamples %<>% filter(UQ(as.symbol(column)) == value)
      }
    }

    write_delim(contrastGroupASamples, paste0(contrastOutDir, "groupASamples.csv"), delim = "\t")
    write_delim(contrastGroupBSamples, paste0(contrastOutDir, "groupBSamples.csv"), delim = "\t")

    # now that we have the samples of groupA and groupB, we can do the DE
    group1 <- as.character(contrastGroupASamples$sampleID)
    group2 <- as.character(contrastGroupBSamples$sampleID)
    reads <- rawReadsCID %>% select(one_of(c(group1, group2)))

    # remove genes with 0 reads in all samples
    reads <- reads[rowSums(reads) > 0, ]

    reads %>%
      rownames_to_column("genes") %>%
      write_delim(paste0(contrastOutDir, "deRawReads.csv"), delim = "\t")

    groupVector <- rep(0, length(c(group1, group2)))
    groupVector[1:length(group1)] <- contrastLine$GroupA
    groupVector[length(group1) + 1:length(group2)] <- contrastLine$GroupB

    # check if we have a grouping factor set for this contrast
    if(!is.na(contrastLine$`Grouping Factor`) && contrastLine$`Grouping Factor` != "") {
      pairingVector <- c(contrastGroupASamples %>% pull(contrastLine$`Grouping Factor`), contrastGroupBSamples %>% pull(contrastLine$`Grouping Factor`))
    } else {
      pairingVector <- c()
    }

    # write details about this contrast into a file
    contrastLine %>% write_delim(paste0(contrastOutDir, "deContrastDetails.csv"), delim = "\t")

    runDe <- function(reads, groupVector, pairingVector) {
      group <<- factor(groupVector, levels = c(contrastLine$GroupA, contrastLine$GroupB))
      dge <<- DGEList(counts = reads, group = group, genes = row.names(reads))

      if(length(pairingVector) > 0) {
        pairing <<- factor(pairingVector)
        design <- model.matrix(~0 + pairing + group, data = dge$samples)
        contrastVector <- c(rep(0, ncol(design)-1), -1) # make a contrastVector that selects the last column of the design matrix as contrast of interest
      } else {
        # no pairing
        contrastVector <- c(1, -1)
        design <- model.matrix(~0 + group, data = dge$samples)
      }

      # do cpm normalization
      dge <<- calcNormFactors(dge)
      dge <<- estimateCommonDisp(dge)
      dge <<- estimateTagwiseDisp(dge)

      # DE testing:
      dge.fit <<- glmQLFit(dge, design, robust = TRUE)

      # differnetial expression test statistics
      dge.test <<- glmQLFTest(dge.fit, contrast = contrastVector) # use this with glmQLFit

      # get CPM of reads
      reads.cpm <<- as_tibble(cpm(dge), rownames = NA)

      # combine cpms and results into one dataframe
      unfiltered.results <<- bind_cols(dge$genes, reads.cpm, dge.test$table)
    }

    indepFiltering <- function(reads.cpm) {
      # based on code from https://statquest.org/statquest-filtering-genes-with-low-read-counts/
      # implementing methods by:
      ## Robinson MD, McCarthy DJ and Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26, pp. -1.
      ## Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550. doi: 10.1186/s13059-014-0550-8.

      filter <- apply(
        X = reads.cpm, MARGIN = 1,
        FUN = function(data) {
          data[order(rank(data), decreasing = TRUE)[2]]
        }
      )

      # percentage of unfiltered genes (mean of TRUE/FALSE)
      lowerQuantile <<- mean(filter == 0)
      if (lowerQuantile < .95) upperQuantile <<- .95 else upperQuantile <<- 1

      # generate a range of cutoff quantiles that we will try out for p value adjustment
      theta <<- seq(lowerQuantile, upperQuantile, length = 100)

      filtPadj <<- filtered_p(
        filter = filter, test = unfiltered.results$PValue,
        theta = theta, method = "BH"
      )

      min.fdr <<- alpha
      numSigGenes <<- colSums(filtPadj < min.fdr, na.rm = TRUE)
      filter.quantiles <<- quantile(filter, probs = theta)
      lo.fit.theta <<- lowess(numSigGenes ~ theta, f = 1 / 10)

      if (max(numSigGenes) <= 10 && 1 == 2) {
        j <<- 1
        thresh <<- lowerQuantile
      } else {
        residual <<- if (all(numSigGenes == 0)) {
          0
        } else {
          numSigGenes[numSigGenes > 0] - lo.fit.theta$y[numSigGenes > 0]
        }
        thresh <<- max(lo.fit.theta$y) - sqrt(mean(residual^2))
        j <<- if (any(numSigGenes > thresh)) {
          which(numSigGenes > thresh)[1]
        } else {
          1
        }
      }
    }

    runDe(reads, groupVector, pairingVector)
    indepFiltering(reads.cpm)

    filtered.results <- unfiltered.results
    filtered.results$FDR <- filtPadj[, j, drop = TRUE]

    # calculate the percentage of removed read counts
    removed.genes <- filtered.results[is.na(filtered.results$FDR), "genes"]
    removed.rc <- reads[removed.genes, ]
    removed.rcPerc <- sum(colSums(removed.rc)) / sum(colSums(reads))

    if (j == 1 || removed.rcPerc < 0.01) {
      minLibSize <- min(colSums(reads)) / 1000000
      reads.cpm <- cpm(reads)

      # https://f1000research.com/articles/5-1438
      keep <- rowSums(reads.cpm > 10 / minLibSize) >= round(min(c(nrow(contrastGroupASamples), nrow(contrastGroupBSamples))) / 2, 0)
      removed <- reads[!keep, ]
      removed.cpm <- reads.cpm[!keep, ]
      removedStats <- list()
      removedStats$rpmCutoff <- 10 / minLibSize
      removedStats$minLibs <- round(min(c(nrow(contrastGroupASamples), nrow(contrastGroupBSamples))) / 2, 0)
      removedStats$removedRC <- sum(removed)
      removedStats$removeGenes <- sum(keep == FALSE)
      removedStats$removedRCPerc <- removedStats$removedRC / sum(reads)

      removedStats %>%
        as_tibble() %>%
        write_delim(paste0(contrastOutDir, "dePrefilteringStats.csv"), delim = "\t")
      removed.cpm %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        write_delim(paste0(contrastOutDir, "dePrefilteringRemovedGenes.csv"), delim = "\t")

      reads <- reads[keep, ]

      # run DE again with basic filtered data
      if(nrow(reads) > 1) {
        runDe(reads, groupVector, pairingVector)
        indepFiltering(reads.cpm)
      } else {
        unlink(paste0(contrastOutDir, "dePrefilteringStats.csv"))
        unlink(paste0(contrastOutDir, "dePrefilteringRemovedGenes.csv"))
      }
    } else {
      unlink(paste0(contrastOutDir, "dePrefilteringStats.csv"))
      unlink(paste0(contrastOutDir, "dePrefilteringRemovedGenes.csv"))
    }

    filtered.results <- unfiltered.results
    filtered.results$FDR <- filtPadj[, j, drop = TRUE]

    # calculate the percentage of removed read counts
    removed.genes <- filtered.results[is.na(filtered.results$FDR), "genes"]
    removed.rc <- reads[removed.genes, ]
    removed.rcPerc <- sum(colSums(removed.rc)) / sum(colSums(reads))

    tibble(theta, numSigGenes, fitX = lo.fit.theta$x, fitY = lo.fit.theta$y) %>% write_delim(paste0(contrastOutDir, "fdrAdjustmentThreshold.dat"), delim = "\t")
    tibble(removeGenes = length(removed.genes), removedRC = sum(colSums(removed.rc)), removedRCPerc = removed.rcPerc, yIntercept = thresh, xIntercept = lo.fit.theta$x[j]) %>% write_delim(paste0(contrastOutDir, "fdrAdjustmentThreshold.stat"), delim = "\t")

    filtered.results <- filtered.results[!is.na(filtered.results$FDR), ]

    filtered.results$de <- sign(filtered.results$logFC) * (filtered.results$FDR < alpha)
    filtered.results$sig <- as.factor(abs(filtered.results$de))

    # volcano plot
    filtered.results %>% write_delim(paste0(contrastOutDir, "volcanoPlot.dat"), delim = "\t")

    # MA plot
    filtered.results %>% write_delim(paste0(contrastOutDir, "MAPlot.dat"), delim = "\t")

    ## top up and down regulated genes
    topUpGenes <- filtered.results %>%
      filter(logFC > 0 & sig == 1) %>%
      arrange(desc(abs(logFC)))
    upSig <- TRUE
    if (nrow(topUpGenes) == 0) {
      topUpGenes <- filtered.results %>%
        filter(logFC > 0 & sig == 0) %>% # filter for unsignificant results
        arrange(desc(abs(logFC)))
      upSig <- FALSE
    }

    if (nrow(topUpGenes) < numTopGenes) {
      numTopUpGenes <- nrow(topUpGenes)
    } else {
      numTopUpGenes <- numTopGenes
    }
    topUpGenesFC <- topUpGenes[1:numTopUpGenes, c("logFC")]
    topUpGenes <- as.character(topUpGenes[1:numTopUpGenes, c("genes")])

    topDownGenes <- filtered.results %>%
      filter(logFC < 0 & sig == 1) %>%
      arrange(desc(abs(logFC)))
    downSig <- TRUE
    if (nrow(topDownGenes) == 0) {
      topDownGenes <- filtered.results %>%
        filter(logFC < 0 & sig == 0) %>%
        arrange(desc(abs(logFC)))
      downSig <- FALSE
    }

    if (nrow(topDownGenes) < numTopGenes) {
      numTopDownGenes <- nrow(topDownGenes)
    } else {
      numTopDownGenes <- numTopGenes
    }
    topDownGenesFC <- topDownGenes[1:numTopDownGenes, c("logFC")]
    topDownGenes <- as.character(topDownGenes[1:numTopDownGenes, c("genes")])

    upregulated <- tibble(
      gene = factor(),
      sample = double(),
      cpm = double(),
      logCpm = double(),
      logFC = double(),
      group = factor()
    )
    downregulated <- upregulated

    # get the cpm values for each group:
    for (i in levels(group)) {
      groupSamples <- as_tibble(reads.cpm[, group == i])
      # set all genes with a 0 CPM to 1, so that they show up in the boxplot
      groupSamples[groupSamples == 0] <- 1
      groupSamples$genes <- rownames(reads.cpm)

      groupSamplesUp <- groupSamples %>%
        filter(genes %in% topUpGenes) %>%
        add_column(logFC = topUpGenesFC) %>%
        gather(sample, cpm, -genes, -logFC)
      groupSamplesUp %<>% add_column(logCpm = log10(groupSamplesUp$cpm)) %>% add_column(group = i)
      upregulated %<>% rbind(groupSamplesUp)

      groupSamplesDown <- groupSamples %>%
        filter(genes %in% topDownGenes) %>%
        add_column(logFC = topDownGenesFC) %>%
        gather(sample, cpm, -genes, -logFC)
      groupSamplesDown %<>% add_column(logCpm = log10(groupSamplesDown$cpm)) %>% add_column(group = i)
      downregulated %<>% rbind(groupSamplesDown)
    }

    if (nrow(upregulated) > 1) {
      if (upSig == FALSE) {
        nonSigTxt <- ".ns"
      } else {
        nonSigTxt <- ""
      }
      topUpGenes %>%
        as_tibble() %>%
        write_delim(paste0(contrastOutDir, "topUpRegulated", NAspecies, "Names", nonSigTxt, ".dat"), delim = "\t")
      upregulated %>% write_delim(paste0(contrastOutDir, "topUpRegulated", NAspecies, nonSigTxt, ".dat"), delim = "\t")
    }

    if (nrow(downregulated) > 1) {
      if (downSig == FALSE) {
        nonSigTxt <- ".ns"
      } else {
        nonSigTxt <- ""
      }
      topDownGenes %>%
        as_tibble() %>%
        write_delim(paste0(contrastOutDir, "topDownRegulated", NAspecies, "Names", nonSigTxt, ".dat"), delim = "\t")
      downregulated %>% write_delim(paste0(contrastOutDir, "topDownRegulated", NAspecies, nonSigTxt, ".dat"), delim = "\t")
    }
    tt <- topTags(dge.test, n = Inf) # save all results
    write_delim(tt$table, paste0(contrastOutDir, "topTags.csv"), delim = "\t")
    write_delim(filtered.results, paste0(contrastOutDir, "de_ttest_results.csv"), delim = "\t")
  }
} else {
  print("No contrasts defined. Exiting.")
}
