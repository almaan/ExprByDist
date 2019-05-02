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

#viz <- ggplot(data = data_long, aes(x = distance, y = expression, color = gene)) +
#       ggtitle(plt_title) +
#       annotate('rect', xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "red", color = NA, alpha = 0.1) + 
#       annotate('rect', xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "green", color = NA, alpha = 0.1) +
#       geom_smooth(method=smooth_method, se = use_se, alpha = 0.2,formula = formula, span = loess_span) +
#       geom_point()  +
#       scale_color_manual(values=cmap) +
#       xlab('Distance') +
#       ylab('Relative Expression') +
#       theme(plot.title = element_text(hjust = 0.5))

#print(viz)

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