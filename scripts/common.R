# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/ 99% accessibility
cbp1 <- c(
  "#e6194B", # miRNA
  "#3cb44b", # tNRA
  "#ffe119", # piRNA
  "#4363d8", # rRNA
  "#f58231", # lncRNA
  "#42d4f4", # mRNA
  "#f032e6", # snRNA
  "#fabebe", # snoRNA
  "#aaffc3", # yRNA
  "#e6beff", # scRNA
  "#fffac8", # rnaSpikeins
  "#000075", # ngsSpikeins
  "#469990", # other RNA species
  "#9a6324", # unclassified genomic
  "#a9a9a9", # unmapped
  "#a9a9a9", # unused
  "#222222"  # unused
)
cbp2 <- rev(cbp1)
scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = cbp1)
}
scale_colour_discrete <- cbp1

newline_vs_lables <- function(x,...){
  gsub('\\sblocking\\s', '\nblocking ', gsub('\\sversus\\s','\n',x))
}

getCondaVersions <- function(path) {
  library(yaml)
  packages <- tibble(package = character(), version = character())
  envs <- list.files(path, pattern = ".*yml$")
  for (env in envs) {
    filePath <- paste(path, env, sep = "/")
    yamlFile <- read_yaml(filePath)
    for (package in yamlFile$dependencies) {
      packageString <- strsplit(package, "=")
      packageName <- packageString[[1]][1]
      packageVersion <- packageString[[1]][2] %>% str_replace_all('\\.\\*', '')
      packages %<>% add_row(package = packageName, version = packageVersion)
    }
  }
  packages %<>% distinct(package, .keep_all = TRUE) %>% column_to_rownames("package")
  return(packages)
}

# http://michaeljw.com/blog/post/subchunkify/
subchunkify <- function(g, fig_height=7, fig_width=5) {
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')

  sub_chunk <- paste0("
  `","``{r sub_chunk_", floor(runif(1) * 10000000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                      "\n(",
                      g_deparsed
                      , ")()",
                      "\n`","``
  ")

  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}


# https://github.com/crj32/M3C/blob/master/R/tsne.R
# https://www.bioconductor.org/packages/3.7/bioc/html/M3C.html
# John C (2018). M3C: Monte Carlo Consensus Clustering. R package version 1.2.0.
# License 	AGPL-3
tsne <- function(mydata, labels=FALSE, perplex=15, printres=FALSE, seed=FALSE, axistextsize = 18,
                 legendtextsize = 18, dotsize = 5, textlabelsize = 4, legendtitle = 'Group',
                 controlscale = FALSE, scale = 1, low = 'grey', high = 'red',
                 colvec = c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen",
                            "violetred", "orange", "midnightblue", "grey31", "black"),
                 printheight = 20, printwidth = 22, text = FALSE){

  ## basic error handling
  
  if ( controlscale == TRUE && class(labels) %in% c( "character", "factor") && scale %in% c(1,2) ) {
    stop("when categorical labels, use scale=3")
  }
  if ( controlscale == TRUE && class(labels) %in% c( "numeric") && scale %in% c(3) ) {
    stop("when continuous labels, use scale=1 or scale=2")
  }
  if ( controlscale == FALSE && scale %in% c(2,3) ) {
    warning("if your trying to control the scale, please set controlscale=TRUE")
  }
  if (sum(is.na(labels)) > 0 && class(labels) %in% c('character','factor')){
    warning("there is NA values in the labels vector, setting to unknown")
    labels <- as.character(labels)
    labels[is.na(labels)] <- 'Unknown'
  }
  if (sum(is.na(text)) > 0 && class(text) %in% c('character','factor')){
    warning("there is NA values in the text vector, setting to unknown")
    text <- as.character(text)
    text[is.na(text)] <- 'Unknown'
  }

  ##

  message('***t-SNE wrapper function***')
  message('running...')

  if (seed != FALSE){
    set.seed(seed)
  }

  ## combinations

  # K FALSE, labels FALSE, text FALSE
  # K TRUE, labels FALSE (when K TRUE just over ride everything)
  # K FALSE, labels TRUE, text ?
  # require the above, but with text or without so new entry
  # K FALSE, labels FALSE, text TRUE

  ##

  if (labels[1] == FALSE && text[1] == FALSE){

    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix

    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(colour='skyblue', size = dotsize) +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize))+
      scale_colour_manual(values = colvec)

    if (printres == TRUE){
      message('printing t-SNE to current directory...')
      png('TSNE.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }

  }else if (labels[1] != FALSE && text[1] == FALSE){ ##### KEY

    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix

    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")
        #scale_colour_gradient(low="red", high="white")
      }else if (scale == 2){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + #scale_colour_distiller(palette = "Spectral")
          scale_colour_gradient(low=low, high=high)
      }else if (scale == 3){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) +
          scale_colour_manual(values = colvec)
      }
    }else{
      p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize, colour = 'black'),
              axis.text.x = element_text(size = axistextsize, colour = 'black'),
              axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) +
        #guides(colour=guide_legend(title=legendtitle)) +
        labs(colour = legendtitle)
    }

    if (printres == TRUE){
      message('printing tSNE to current directory...')
      png('TSNElabeled.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }

  }else if (labels[1] != FALSE && text[1] != FALSE){ ##### KEY

    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    scores$label <- text

    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
        #scale_colour_gradient(low="red", high="white")
      }else if (scale == 2){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + #scale_colour_distiller(palette = "Spectral")
          scale_colour_gradient(low=low, high=high)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
      }else if (scale == 3){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) +
          scale_colour_manual(values = colvec)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
      }
    }else{
      p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize, colour = 'black'),
              axis.text.x = element_text(size = axistextsize, colour = 'black'),
              axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) +
        #guides(colour=guide_legend(title=legendtitle)) +
        labs(colour = legendtitle)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
    }

    if (printres == TRUE){
      message('printing t-SNE to current directory...')
      png('TSNElabeled.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }

  }else if (labels[1] == FALSE && text[1] != FALSE){

    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y)
    scores$label <- text

    p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) +
      geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) +
      theme_bw() +
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = colvec) + geom_text(vjust="inward",hjust="inward",size=textlabelsize)

    if (printres == TRUE){
      message('printing t-SNE to current directory...')
      png('TSNE.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }

  message('done.')

  return(p)
}

#SpikeIn Lots

spikeInLots <- list(
  Vers2 = c( # concentrations in amol/uL equals the amol added to each finalVolume
    "C" = 0.01,
    "E" = 0.005,
    "H" = 0.1,
    "I" = 50,
    "K" = 10,
    "M" = 2.5,
    "N" = 1.5),
  Vers3.April2021 = c( # concentrations in amol/uL equals the amol added to each finalVolume
    "C" = 0.01,
    "E" = 0.005,
    "H" = 0.075,
    "I" = 20,
    "K" = 5,
    "M" = 1.25,
    "N" = 0.3125),
  Vers3.December2020 = c( # concentrations in amol/uL equals the amol added to each finalVolume
    "C" = 0.01,
    "E" = 0.005,
    "H" = 0.075,
    "I" = 20,
    "K" = 5,
    "M" = 1.25,
    "N" = 0.3125),
  Vers3 = c( # concentrations in amol/uL equals the amol added to each finalVolume
    "miND-01" = 20,
    "miND-02" = 5,
    "miND-03" = 1.25,
    "miND-04" = 0.3125,
    "miND-05" = 0.075,
    "miND-06" = 0.01,
    "miND-07" = 0.005)
)
