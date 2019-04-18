#!/usr/bin/Rscript

sh <- suppressPackageStartupMessages
sh(library(spatstat))
sh(library(optparse))
sh(library(org.Hs.eg.db))
sh(library(AnnotationDbi))
sh(library(ggplot2))
sh(library(reshape2))
sh((library(futile.logger)))
# Functions -----------------------------

ensg_to_rgb <- function(ensgid) {
  #' Assigns a color for each
  #' ENSEMBL-id to make 
  #' plots comparable. Uniqueness is not
  #' guaranteed but large variation between
  #' genes are usually observed

  num <- gsub('ENSG','',ensgid)
  num <- as.numeric(gsub('^0*','',num))
  rr <- (num %% 255 ) / 255
  gg <- ((num / 6) %% 255 ) / 255
  bb <- ((num / 7) %% 255 ) / 255
  rgbval <- rgb(r = rr, g = gg, b = bb)
  return(rgbval)
}

glue <- function(x) {
  return(paste(x,collapse = ' '))
}

# Wrapper for multiple argument parsing -----------------

allowMultipleArgs <- function(){
  
  #' Modify trailing arguments passed such that space
  #' separated arguments to same flag becomes joined by
  #' commas; a format supported by optparse, and which later 
  #' easily can be split into separate parts again
  
  
  oriArgs <- commandArgs(trailingOnly = TRUE)
  flags.pos <- which(sapply(oriArgs, function(x) '-' == substr(x,1,1)))
  newArgs <- c()
  
  if (length(flags.pos) > 1) {
    for (i in 1:(length(flags.pos)-1))
    {
      if ((flags.pos[i] + 1) != flags.pos[i+1]) {
        pos <- c((flags.pos[i]+1):(flags.pos[i+1]-1))
        newArgs <- c(newArgs,oriArgs[flags.pos[i]], paste(oriArgs[pos],collapse=','))
      } else {
        newArgs <- c(newArgs,oriArgs[flags.pos[i]])
      }
    }
  }
  
  if (length(oriArgs) > tail(flags.pos,n=1)) {
    pos <- c((flags.pos[length(flags.pos)]+1):length(oriArgs))
    newArgs <- c(newArgs, oriArgs[tail(flags.pos,n=1)],paste(oriArgs[pos],collapse=','))
  } else {
    newArgs <- c(newArgs, oriArgs[tail(flags.pos,n=1)])
  }
  return(newArgs)
}


splitMultipleArgs <- function(optArgs) {
  
  #' Use in combination with allowMultipleArgs
  #' will split all commaseparated arguments
  #' into individual elements in list
  
  for (i in 1:length(optArgs)) {
    if (grepl(",",optArgs[[i]])) {
      optArgs[[i]] <- unlist(strsplit(optArgs[[i]],','))
    }
  }
  
  return(optArgs)
}

# Parser ------------------------------

parser <- OptionParser()

parser <- add_option(parser,
                     c("-d","--dge_res"),
                     default = NULL,
                     help =glue(c("Name of DGE-analysis results",
                                  "file. Can be used as an alternative",
                                  "to specifying specific genes"))
)

parser <- add_option(parser, c("-g","--genes"),
                               default = NULL,
                               help = glue(c("name of specfifc genes",
                                             "to be visualized.",
                                             "Is case insensitive"))
)


parser <- add_option(parser,
                     c("-c","--count_files"),
                     default = NULL,
                     help = glue(c("path to count files.",
                                   "order should match that",
                                   "of feature file(s).",
                                   "supports globbing."
                                   ))
)

parser <- add_option(parser,
                     c("-f","--feature_files"),
                     default = NULL,
                     help = glue(c("path to feature files.",
                                   "order should match that",
                                   "of count file(s)",
                                   "supports globbing."
                     ))
)

parser <- add_option(parser,
                     c("-p","--polydeg"),
                     default = 5,
                     help = glue(c("degree of polynomial",
                                   "to be used for interpolation",
                                   "if method 'polynomial' is used"
                                   ))
)

parser <- add_option(parser,
                     c("-t","--title"),
                     default = NULL,
                     help = glue(c('title of plot',
                                   'if non specified an automated',
                                   'title will be generated'
                                   ))
)

parser <- add_option(parser,
                     c("-r","--inner_tag"),
                     default = 'tumor',
                     help = glue(c("name of variable to be",
                                   "used as inner region.",
                                   "Default is tumor. Will have",
                                   "negative distance values"))
)

parser <- add_option(parser,
                     c("-x","--outer_tag"),
                     default = 'non',
                     help = glue(c("name of variable to be",
                                   "used as inner region.",
                                   "Default is non as in non-tumor",
                                   "Will have positive distance",
                                   "values."
                                   ))
)

parser <- add_option(parser,
                     c("-m","--method"),
                     default = 'polynomial',
                     help = glue(c("method to use for the",
                                   "curve fitting. Available",
                                   "methods are 'polynomial'",
                                   "(Default) and loess."))
)

parser <- add_option(parser,
                     c("-u","--positive_lfc"),
                     default = F,
                     action = 'store_true',
                     help = glue(c("use flag to depict genes",
                                   "with a positive logFoldChange",
                                   "These will have red values in plot."))
)

parser <- add_option(parser,
                     c("-n","--negative_lfc"),
                     default = F,
                     action = 'store_true',
                     help = glue(c("use flag to depict genes",
                                   "with a negative logFoldChange",
                                   "These will have blue values in plot."))
)

parser <- add_option(parser,
                     c("-s","--loess_span"),
                     default = 0.5,
                     help = glue(c("span to use",
                                   "in loess smoothing."))
)

parser <- add_option(parser,
                     c("-z","--z_transform"),
                     default = F,
                     action = 'store_true',
                     help = glue(c("specify if z-transform of",
                                   "values should be performed.",
                                   "REMAIN TO BE IMPLEMENTED."))
)

parser <- add_option(parser,
                     c("-e","--standard_error"),
                     default = F,
                     action = 'store_true',
                     help = glue(c("Visualize standard error",
                                    "in the resulting plots."))
)

parser <- add_option(parser,
                     c("-q","--no_control"),
                     default = T,
                     action = 'store_false',
                     help = glue(c("Do not control for matching",
                                   "feature files and count files order",
                                   "not recommended to use."))
)

# Variables -----------------

if (interactive()) {
  # example set when running in interactive mode to test performance
  
  st_cnt_files <- c('exdata/st-00001_X1.tsv',
                    'exdata/st-00001_X2.tsv',
                    'exdata/st-00001_X3.tsv')
  
  feat_files <- c('exdata/ft-00001_X1.tsv',
                  'exdata/ft-00001_X2.tsv',
                  'exdata/ft-00001_X3.tsv')

  dge_pth <- NULL
  polydeg <- 5
  loess_span <- 0.5
  outer_tag <- 'non'
  inner_tag <- 'tumor'
  positive_lfc <- T
  negative_lfc <- T
  method = 'loess'
  main_title <- "test"
  z_transform <- F
  use_genes <- c("ERBB2","MZB1")
  use_standard_error <- T
  use_se <- F
  use_control <- TRUE

} else {
  # For CLI usage
  args <- splitMultipleArgs(parse_args(parser, args = allowMultipleArgs())) 
  st_cnt_files <- args$count_files
  feat_files <- args$feature_files
  dge_pth <- args$dge_res
  genes_pth <- args$dge_res
  polydeg <- args$polydeg
  main_title <- args$title
  outer_tag <- args$outer_tag
  inner_tag <- args$inner_tag
  method <- args$method
  loess_span <- args$loess_span
  positive_lfc <- args$positive_lfc
  negative_lfc <- args$negative_lfc
  z_transform <- args$z_transform
  use_genes <- args$genes
  use_se <- args$standard_error
  use_control <- args$no_control
}

# Validate input -----------------------

# check for consitency in number of files
if(!(length(st_cnt_files) == length(feat_files))) {
  flog.error('The number of count matrices does not match
             the number of feature files. Exiting.')
  quit(status = 1)
}

# makes sure count and feature fiels are matched
if (use_control) {
flog.info('Control for matching feature files and count files')
  for (filenum in c(1:length(feat_files))) {
      pat <- "[0-9]{5}_[A-Z][0-9]" # section id pattern
      # grep for section specific pattern in filenames
      frgx <- regexpr(pat,feat_files[filenum],perl =T)
      crgx <- regexpr(pat,st_cnt_files[filenum],perl =T)
      # extract match
      len <- attr(frgx, 'match.length')
      fsec <- substr(feat_files[filenum],frgx[1],frgx[1] + len)
      csec <- substr(st_cnt_files[filenum],crgx[1],crgx[1] + len)
      
      # exit if count and feature file section id does not match  
      if (!(fsec == csec)) {
        log.error('Feature Files and Count Files are not properly matched')
        quit(status =1)
      }
    }
}

flog.info('Matching feature files and count files')
flog.info('Will Initiate analysis')

# Main -----------------------
analysis_type <- ifelse(z_transform,'z_transform','relative.freqs')
keepgenes <- c()
cmap <- c()

# if a dge-file is provided
if (!(is.null(dge_pth))) {
  flog.info('Use Gene-Set and foldChange based on DE-results')
  # load genes to be used in analysis 
  genes <- read.table(dge_pth, sep = ',',
                      header =1,
                      row.names = 1,
                      stringsAsFactors = F)
  
  # if genes with positive lFC should be analyzed
  if (positive_lfc) {
    # get index of up-regulated genes
    pidx <- genes$logFC > 0
    # add up-regulated to genes to set to be analyzed
    keepgenes <- c(keepgenes,rownames(genes)[pidx])
    # create colormap for up-regulated genes. Red spectrum
    print(sprintf('>> Using %d Genes With Positive lFC',sum(pidx)))
    rch <- 255 - ((c(1:sum(pidx))*5) %% 255)
    prgb <- cbind(rch, rep(0,length(rch)),rep(0,length(rch))) / 255
    cmap <- c(cmap,rgb(prgb))
  }
  
  # if genes with negative lFC should be analyzed
  if (negative_lfc) {
    # get index of down-regulated genes
    nidx <- genes$logFC < 0
    print(sprintf('>> Using %d genes with negative lFC',sum(nidx)))
    # add down-regulated genes to set to be analyzed
    keepgenes <- c(keepgenes,rownames(genes)[nidx])
    # create colormap for down-regulated genes. Blue spectrum
    bch <- 255 - ((c(1:sum(nidx))*5) %% 255)
    nrgb <- cbind(rep(0,length(bch)),rep(0,length(bch)),bch) / 255
    cmap <- c(cmap,rgb(nrgb))
  }
  
  if (!(negative_lfc) & !(positive_lfc)) {
    flog.error('Specify at least one set of genes to plot (lFC >0 or lFC <0')
    flog.error('Exiting')
    quit(status = 1)
  }
  
  # name of analysis
  bname <- gsub("\\.tsv","",basename(genes_pth))
  # save result to png-file
  imgpth <- paste(c(dirname(genes_pth),paste(c(bname,analysis_type,"count_by_distance.png"),collapse = '.')),collapse = '/')
  flog.info(sprintf('Will save results to >> %s',imgpth))
    
  genes <- genes[keepgenes,]
  genes <- genes$genes # get ENSEMBL ids
  genes <- genes[!(is.na(genes))] # remove any potential NA

  if (any(grepl('ENSG',genes))) {
    symbol <- mapIds(org.Hs.eg.db,
                     keys=as.character(genes),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
  } else {
    symbol <- genes
  }
  
  dropgenes <- which(is.na(genes) | is.na(symbol))
  symbol <- symbol[-dropgenes]
  genes <- genes[-dropgenes]
      
# if gene symbols are provided
} else if (!(is.null(use_genes))) {
  genes <- use_genes
  
  n_pre <- length(use_genes)
  flog.info('Convert Gene Symbols to ENSEMBL Id')
  genes <-mapIds(org.Hs.eg.db,
          keys=as.character(genes),
          column="ENSEMBL",
          keytype="SYMBOL",
          multiVals="first")
  
  symbol <- use_genes
  genes <- genes[!(is.na(genes))]
  symbol <- symbol[!(is.na(genes))]

  n_post <- length(genes)
  
  flog.info(sprintf("Lost %d/%d genes when converting between gene id's ", n_pre - n_post,n_pre))
    
  grgb <- unlist(lapply(genes,ensg_to_rgb))
  cmap <- c(cmap,grgb)
  
  odir <- paste(c(getwd(),'exprbydist'), collapse = "/")
  if (!(dir.exists(odir))) {
    dir.create(odir)
  }
  bname <- paste(c("custom_geneset",gsub(' |:|-','',as.character(Sys.time()))), collapse = "_")
  imgpth <- paste(c(odir,paste(c(bname,analysis_type),collapse = '.')),collapse = '/')
}

mx <- 60 #  maximum distance from distance from tumor spot
dr <- 0.1 #  distance increment
lims <- seq(1-dr,mx+dr,dr) #  lower limits
meanpnt <- lims + dr/2.0 #  mean point within each interval

# dataframe to hold mean expression values for each distance value
datanames <- c(rev(as.character(-meanpnt[-1])),as.character(meanpnt))
spotcount <- data.frame(matrix(0,
                               nrow = length(datanames)
                               ),
                        row.names = datanames
                        )

data <- data.frame(matrix(0,
                          ncol = length(genes),
                          nrow = length(datanames)
                          ),
                   row.names = datanames)
colnames(spotcount) <- 'n_spots'
colnames(data) <- genes

# iterate over all filenames
flog.info('Read data from files')
for (filenum in c(1:length(st_cnt_files))){
  
  
  st_cnt_pth <- st_cnt_files[filenum]
  feat_pth <- feat_files[filenum]
  
  flog.info(sprintf('Using st-count file >> %s',st_cnt_pth))
  flog.info(sprintf('Usign feature file >> %s',feat_pth))
  
  feat <- read.table(feat_pth, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)
  
  # raw count matrix; n_spots x n_genes
  cnt_raw <- read.table(st_cnt_pth,
                        sep = '\t',
                        header = 1,
                        row.names = 1,
                        stringsAsFactors = F)
  
  # Remove bias due to spot library size; divide row by rowsum
  cnt_raw <- sweep(cnt_raw, MARGIN = 1, FUN = "/", rowSums(cnt_raw))

  # extract spots present in feature and count file
  interspt <- as.character(intersect(rownames(feat),rownames(cnt_raw)))
  cnt_raw <- cnt_raw[interspt,]
  feat <- feat[interspt,]
  
  # create matrix with all as columns genes
  cnt <- data.frame(matrix(0, nrow = nrow(cnt_raw),ncol = length(genes)))
  colnames(cnt) <- genes
  rownames(cnt) <- interspt
  
  # interection of genelist and genes in st-count data
  if (sum(grepl('ENSG',colnames(cnt_raw))) > 0.5 * ncol(cnt_raw)) {
    flog.info("Assuming Count Matrix uses ENSEMBL id's")
    inter <- as.character(intersect(genes,colnames(cnt_raw)))
  } else {
    flog.info("Assuming Count Matrix uses Gene Symbols")
    inter <- as.character(intersect(symbols,colnames(cnt_raw)))
  }
  
  # fill count matrix with scaled st-counts
  cnt[,inter] <- cnt_raw[,inter]
  remove(cnt_raw) # free up memory
  
  # get tumor and non-tumor spot coordinates
  innerspts <- feat[feat$tumor == inner_tag,c('xcoord','ycoord')] # within reference tissue
  outerspts <- feat[feat$tumor == outer_tag,c('xcoord','ycoord')] # outside referene tissue
  
  # create distance matrix (euclidian distance)
  dmat = crossdist.default(X = outerspts[,1], # target spots along dim 1
                           Y = outerspts[,2],
                           x2 = innerspts[,1], # reference spots along dim 2
                           y2 = innerspts[,2])
  
  rownames(dmat) = rownames(outerspts)
  colnames(dmat) = rownames(innerspts)
  
  # get minimum distance to a tumor spot for each non-tumor spot
  mindist_ref <- apply(dmat, 2, min) # minimum distance to nearest target spot for all reference spots
  mindist_trgt <- apply(dmat, 1, min) # minimum distance to nearest reference spot for all target spots
  
  # get expression levels for spots within reference tissue
  for (jj in c(1:length(lims))) {
    # find spots with mindist within I = [lim,lim+dr)
    idx_within_ref <- (mindist_ref >= lims[jj]) & (mindist_ref < lims[jj] + dr) 
    # get count matrix index for spots within I = [lim,lim+dr)
    spt_within_ref <- which(rownames(cnt) %in% colnames(dmat)[idx_within_ref])
    spotcount[length(lims)-jj,] <- spotcount[length(lims)-jj,] + sum(idx_within_ref)
    # if I is non-empty
    if (length(spt_within_ref) > 0) {
      # add mean of observed counts taken over points within I
      data[(length(lims) - jj),] <- data[length(lims) - jj,] + colMeans(cnt[spt_within_ref,]) 
    }
  }
  
  # get expression levels for spots within target tissue
  for (ii in c(1:length(lims))) {
    # find spots with mindist within I = [lim,lim+dr)
    idx_within_trgt <- (mindist_trgt >= lims[ii]) & (mindist_trgt < lims[ii] + dr) 
    # get count matrix index for spots within I = [lim,lim+dr)
    spt_within_trgt <- which(rownames(cnt) %in% rownames(dmat)[idx_within_trgt]) #  
    spotcount[(length(lims) + ii),] <-  spotcount[(length(lims) + ii),] + sum(idx_within_trgt)
    # if I is non-empty
    if (length(spt_within_trgt) > 0) {
      # add mean of observed counts taken over points within I
      data[(length(lims) + ii),] <- data[(length(lims) + ii),] +  colMeans(cnt[spt_within_trgt,]) 
    }
  }

}

# adjust for multiple sections having been used
dataf <- data / length(st_cnt_files)

# get distances where 0 counts are observed
keepdist <- as.vector(which(rowSums(data) > 0, useNames = F))
# if bad gene set
if (length(keepdist) == 0) {
  flog.error('None of the selected genes were observed in any of the sections')
  flog.error('Exiting')
  
}
flog.info(sprintf('%d distances have at least one observed count of selected genes',length(keepdist)))
# remove genes and distances where no observetions are made
dataf <- dataf[keepdist,colSums(data) > 0]
symbol <- symbol[colSums(data) >0]
genes <- genes[colSums(data) >0]
# remove "empty" distances
meanpntf <- meanpnt[keepdist] 
datanamesf <- datanames[keepdist]
spotcount$distance <- as.numeric(rownames(spotcount))
spotcountf <- spotcount[keepdist,]

# perform z-transform
# TODO: Needs to be looked over
if (z_transform) {
  flog.info('Z-transform data')
  dataf <- apply(dataf,2,scale)
}

# Plot Results -------------------------
png(paste(c(imgpth,"count_by_distance.png"),collapse = '.'),
    width = 1000,
    height = 500)

# get plot title
if (is.null(main_title)) {
  plt_title <- bname
} else {
  plt_title <- main_title
}

# if loess is specified
if (method == 'loess') {
  flog.info(sprintf('Using Loess Smoothing with span %f', loess_span))
  smooth_method <- loess
  formula <- y ~ x
  
# default to polynomial if loess is not specified
} else {
  flog.info(sprintf("Using Polynomial of degree %d for Smoothing",polydeg))
  smooth_method <- lm
  formula <-  y ~ poly(x, polydeg, raw=TRUE)  
}

cmx <- apply(dataf,2,max) # max value for each gene 
cmn <- apply(dataf,2,min) # min value for each gene

# transform expression levels to unit interval within each gene
# to allow for easier comparision

dataf <- sweep(dataf, FUN = '-', MARGIN = 2, cmn)
dataf <- sweep(dataf, FUN = '/', MARGIN = 2, (cmx-cmn))

colnames(dataf) <- symbol

# add average expression  for reference to dataframe 
dataf$average <- rowMeans(dataf)
# set color to be black for reference 
cmap <- c(as.character(cmap),'#000000')
# add distance as a column; for melting
dataf$distance <- as.numeric(rownames(dataf))

# reformat data to be compatible with ggplot
data_long <- melt(dataf, id="distance")
data_long$value <- as.numeric(data_long$value)
colnames(data_long)[grepl('variable',colnames(data_long))] <- 'gene'
colnames(data_long)[grepl('value',colnames(data_long))] <- 'expression'

viz <- ggplot(data = data_long, aes(x = distance, y = expression, color = gene)) +
       ggtitle(plt_title) +
       annotate('rect', xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "red", color = NA, alpha = 0.1) + 
       annotate('rect', xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "green", color = NA, alpha = 0.1) +
       geom_smooth(method=smooth_method, se = use_se, alpha = 0.2,formula = formula, span = loess_span) +
       geom_point()  +
       scale_color_manual(values=cmap) +
       xlab('Distance') +
       ylab('Relative Expression') +
       theme(plot.title = element_text(hjust = 0.5))

print(viz)

dev.off()
png(paste(c(imgpth,"spot_distribution.png"),collapse = '.'),
    width = 1000,
    height = 500)

bp <- ggplot(data = spotcountf, aes(x = as.factor(distance))) +
    ggtitle(plt_title) +
    geom_bar(aes(y = n_spots, fill = n_spots),stat = 'identity', color = 'black' ) +
    scale_fill_gradient(low = "#966d2a",
                        high = "#ffc463",
                        space = "Lab",
                        na.value = "grey50",
                        aesthetics = "fill") +
     xlab('Distance') +
     ylab('Number of spots') +
     theme(plot.title = element_text(hjust = 0.5))

print(bp)
dev.off()

