#!/usr/bin/Rscript

list.of.packages <- c("ggplot2",
                      "reshape2",
                      "BiocGenerics",
                      "spatstat",
                      "datasets",
                      "rpart",
                      "IRanges",
                      "graphics",
                      "methods",
                      "S4Vectors",
                      "stats4",
                      "nlme",
                      "grDevices",
                      "base",
                      "optparse",
                      "spatstat.data",
                      "utils")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
    install.packages(new.packages)
}

bio.packages <- c('AnnotationDbi','org.Hs.eg.db')
new.bio.packages <- bio.packages[!(bio.packages %in% installed.packages()[,"Package"])]

if ( 'AnnotationDbi' %in% new.bio.packages ) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
    BiocManager::install("AnnotationDbi", version = "3.8")
}

if ('org.Hs.eg.db' %in% new.bio.packages) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
}
