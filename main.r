#!/usr/bin/Rscript

library(spatstat)
library(optparse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Functions -----------------------------

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
                     help =""
)

parser <- add_option(parser, c("-g","--genes",
                               default = NULL,
                               help = ""
                               )
)


parser <- add_option(parser,
                     c("-c","--count_files"),
                     default = NULL,
                     help = ""
)

parser <- add_option(parser,
                     c("-f","--feature_files"),
                     default = NULL,
                     help = ""
)

parser <- add_option(parser,
                     c("-p","--polydeg"),
                     default = 5,
                     help = ""
)

parser <- add_option(parser,
                     c("-t","--title"),
                     default = NULL,
                     help = ""
)

parser <- add_option(parser,
                     c("-r","--inner_tag"),
                     default = 'tumor',
                     help = ""
)

parser <- add_option(parser,
                     c("-x","--outer_tag"),
                     default = 'non',
                     help = ""
)

parser <- add_option(parser,
                     c("-m","--method"),
                     default = 'polynomial',
                     help = ""
)

parser <- add_option(parser,
                     c("-u","--positive_lfc"),
                     default = F,
                     action = 'store_true',
                     help = ""
)

parser <- add_option(parser,
                     c("-n","--negative_lfc"),
                     default = F,
                     action = 'store_true',
                     help = ""
)

parser <- add_option(parser,
                     c("-s","--loess_span"),
                     default = 0.5,
                     help = ""
)

parser <- add_option(parser,
                     c("-z","--z_transform"),
                     default = F,
                     action = 'store_true',
                     help = ""
)

# Variables -----------------

if (interactive()) {
  # example set when running in interactive mode to test performance
  # ajdust paths to make compatible with local version
  
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
}

# Validate input -----------------------

# check for consitency in number of files
if(!(length(st_cnt_files) == length(feat_files))) {
  print('Unequal number of count and feature files')
  quit(status = 1)
}

# makes sure count and feature fiels are matched
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
    print('Not properly mathced files')
    quit(status =1)
  }
}

print('>> Matching feature files and count files')
print('>> Will Initiate analysis')

# Main -----------------------

keepgenes <- c()
cmap <- c()

if (!(is.null(dge_pth))) {
  print('using DGE-genes')
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
    print('>> Specify at least one set of genes to plot (lFC >0 or lFC <0')
    print('>> Exiting')
    quit(status = 1)
  }
    
  genes <- genes[keepgenes,]
  genes <- genes$genes # get ENSEMBL ids
  genes <- genes[!(is.na(genes))] # remove any potential NA

  symbol <- mapIds(org.Hs.eg.db,
                   keys=as.character(genes),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")
    
} else if (!(is.null(use_genes))) {
  genes <- use_genes
  grgb <- replicate(3,runif(length(genes)))
  cmap <- c(cmap,rgb(grgb))
  
  genes <- mapIds(org.Hs.eg.db,
          keys=as.character(genes),
          column="ENSEMBL",
          keytype="SYMBOL",
          multiVals="first")
  
  symbol <- use_genes
}

print('>> using genes ')
print(genes)
mx <- 50 #  maximum distance from distance from tumor spot
dr <- 0.1 #  distance increment
lims <- seq(1-dr,mx+dr,dr) #  lower limits
meanpnt <- lims + dr/2.0 #  mean point within each interval

# dataframe to hold mean expression values for each distance value
datanames <- c(rev(as.character(-meanpnt[-1])),as.character(meanpnt))
data <- data.frame(matrix(0,ncol = length(genes), nrow = length(datanames)),
                   row.names = datanames)
colnames(data) <- genes


# if symbols are not found use ENSEMBL id
symbol[is.na(symbol)] <- genes[is.na(symbol)]

# iterate over all files
for (filenum in c(1:length(st_cnt_files))){
  
  st_cnt_pth <- st_cnt_files[filenum]
  feat_pth <- feat_files[filenum]
  
  feat <- read.table(feat_pth, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)
  
  cnt_raw <- read.table(st_cnt_pth, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)
  # Remove bias due to spot library size
  cnt_raw <- sweep(cnt_raw, MARGIN = 1, FUN = "/", rowSums(cnt_raw))
  # Remove bias due to gene expression levels, using relative frequencies
  #cnt_raw <- sweep(cnt_raw, MARGIN = 2, FUN = '/', colSums(cnt_raw))
  
  # extract spots present in feature and count file
  interspt <- as.character(intersect(rownames(feat),rownames(cnt_raw)))
  cnt_raw <- cnt_raw[interspt,]
  feat <- feat[interspt,]
  
  # create matrix with all as columns genes
  cnt <- data.frame(matrix(0, nrow = nrow(cnt_raw),ncol = length(genes)))
  colnames(cnt) <- genes
  rownames(cnt) <- interspt
  
  # interection of genelist and genes in st-count data
  inter <- as.character(intersect(genes,colnames(cnt_raw)))
  
  # fill count matrix with scaled st-counts
  cnt[,inter] <- cnt_raw[,inter]
  remove(cnt_raw) # to free up space
  
  # get tumor and non-tumor spot coordinates
  innerspts <- feat[feat$tumor == inner_tag,c('xcoord','ycoord')]
  outerspts <- feat[feat$tumor == outer_tag,c('xcoord','ycoord')]
  
  # create distance matrix (euclidian distance)
  dmat = crossdist.default(X = outerspts[,1], Y = outerspts[,2], x2 = innerspts[,1], y2 = innerspts[,2])
  rownames(dmat) = rownames(outerspts)
  colnames(dmat) = rownames(innerspts)
  
  # get minimum distance to a tumor spot for each non-tumor spot
  mindist_ref <- apply(dmat, 2, min)
  mindist_trgt <- apply(dmat, 1, min)
  
  # get expression levels for spots within reference tissue
  for (jj in c(1:length(lims))) {
    # find spots with mindist within I = [lim,lim+dr)
    idx_within_ref <- (mindist_ref >= lims[jj]) & (mindist_ref < lims[jj] + dr) 
    # get count matrix index for spots within I = [lim,lim+dr)
    spt_within_ref <- which(rownames(cnt) %in% colnames(dmat)[idx_within_ref]) #  
    # if I is non-empty
    if (length(spt_within_ref) > 0) {
      # add mean of observed counts taken over points within I
      data[(length(lims) - jj),] <- sweep(data[length(lims) - jj,], FUN = '+', MARGIN = 2, colMeans(cnt[spt_within_ref,])) 
    }
  }
  
  # get expression levels for spots within target tissue
  for (ii in c(1:length(lims))) {
    # find spots with mindist within I = [lim,lim+dr)
    idx_within_trgt <- (mindist_trgt >= lims[ii]) & (mindist_trgt < lims[ii] + dr) 
    # get count matrix index for spots within I = [lim,lim+dr)
    spt_within_trgt <- which(rownames(cnt) %in% rownames(dmat)[idx_within_trgt]) #  
    # if I is non-empty
    if (length(spt_within_trgt) > 0) {
      # add mean of observed counts taken over points within I
      data[(length(lims) + ii),] <- sweep(data[(length(lims) + ii),], FUN = '+', MARGIN = 2, colMeans(cnt[spt_within_trgt,])) 
    }
  }

}
print(data)
# adjust for multiple sections having been used
dataf <- data / length(st_cnt_files)
keepdist <- as.vector(which(rowSums(data) > 0, useNames = F))
keepdist <- c(ifelse(lims[keepdist[1]] > lims[1],keepdist[1]-1,c()),
              keepdist,
              ifelse(as.numeric(datanames[keepdist[length(keepdist)]] )< as.numeric(datanames[length(datanames)]),keepdist[length(keepdist)] + 1, c()))

# remove genes and distances where no observetions are made
dataf <- dataf[keepdist,colSums(data) > 0]
# remove distances where no spots are present
meanpntf <- meanpnt[keepdist]
datanamesf <- datanames[keepdist]
# perform z-transform
if (z_transform) {
  dataf <- apply(dataf,2,scale)
}

# name of analysis
bname <- gsub("\\.fancy\\.tsv","",basename(genes_pth))
print(paste(c(dirname(genes_pth),paste(c(bname,"count_by_distance.png"),collapse = '.')),collapse = '/'))
analysis_type <- ifelse(z_transform,'z_transform','relative.freqs')
# save result to png-file
png(paste(c(dirname(genes_pth),paste(c(bname,analysis_type,"count_by_distance.png"),collapse = '.')),collapse = '/'),
    width = 1000,
    height = 500)

# panel for plot and legend
par(mfrow=c(1,2))
# get plot title
if (is.null(main_title)) {
  plt_title <- gsub('DGE_analysis\\.','',bname)
} else {
  plt_title <- main_title
}

# Plot Results -------------------------

cmx <- apply(dataf,2,max)
cmn <- apply(dataf,2,min)
dataf <- sweep(dataf, FUN = '-', MARGIN = 2, cmn)
dataf <- sweep(dataf, FUN = '/', MARGIN = 2, (cmx-cmn))


xmin <- min(as.numeric(datanamesf))
xmax <- max(as.numeric(datanamesf))
# initiate plot
plot(as.numeric(datanamesf),dataf[,1],
   #  ylim = c(ifelse(z_transform,(min(dataf)-10),0), max(dataf + 1)),
     #ylim = c(0, (max(dataf) + max(dataf)*0.2)),
     ylim = c(-0.15,1.15),
     xlim = c(xmin - 0.5,xmax+0.5),
     xlab = paste(c(outer_tag, 'distance to nearest', inner_tag, 'spot'), collapse = ' '),
     ylab = 'mean expression',
     col = cmap[1],
     main = plt_title,
     cex.main = 0.9,
     pch = 19
     )

# add points for all genes
for (ii in 2:(ncol(dataf))) {
  points(x= as.numeric(datanamesf), dataf[,ii], col = cmap[ii], pch = 19)
}

# add fitted curves for each gene using same color as points
fitted <- list()
if (method == 'loess'){

  for (ii in c(1:ncol(dataf))) {

    fitted[[genes[ii]]] <- loess(y ~ x,
                                 data = data.frame(x = as.numeric(datanamesf), # distance is independent variable
                                                   y = dataf[,ii]), # mean expression is dependent variable
                                 span = loess_span)
    
    lines(datanamesf, predict(fitted[[genes[ii]]], data.frame(x=as.numeric(datanamesf))), col=cmap[ii], lwd = 5)
    
  }
  
} else if (method == 'polynomial') {
  for (ii in c(1:ncol(dataf))) {
    # using polynomial of degree "polydeg" to fit data
    fitted[[genes[ii]]] <- lm(formula = y ~ poly(x, polydeg, raw=TRUE),
                              data = data.frame(x = as.numeric(datanamesf), # distance is independent variable
                                                y = dataf[,ii]) # mean expression is dependent variable
                              )
    
    lines(as.numeric(datanamesf), predict(fitted[[genes[ii]]], data.frame(x=datanamesf)), col=cmap[ii], lwd = 5)
  }
}

# plot legend in different subplot for readability
plot(x=NULL, xlim = c(0,2), ylim = c(0,1), xlab = '', ylab = '')
legend(0,1,legend= c(symbol), fill = cmap, cex = 0.5, bty = 'n', ncol = 5  )
dev.off()

