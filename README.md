# Gene Expression as a function of distance

A simple CLI based implementation which allows one to plot the gene expression as a function of the distance to a given feature. The current application is implemented with the intentio to be used for
cancer samples. The main plot will on the right hand side (positive) have the expression of genes given as a function of the distance to the _nearest_ tumor spot whilst the left hand side (negative) will
plot the expression of the same gene as a function of the distance to the neares non-tumor spot. To clearly illustrate which values are calculated within the tumor these are plotted over a red
background whilst those outside are plotted over a green background. The smoothed curves are prodcued by fitting a gaussian process (GP) to the observed count data (normalized within each spot as to
be represented in terms of relative frequencies). The predicted mean value is visualized with an error bars indicating +/- 2sd. The fitting is performed using kernlab's *gausspr* with default
settings.

You can specify explicity genes which you want to study using the command "--genes GENE1 GENE2 GENE3", where each of the genes will be plotted in the same graph. To make comparisions between different
runs easier, the genes will have the same color in any run independent of the gene set, the color is given as a function of the ENSEMBL id. As an alternative to specifying genes, one can use the genes
identified as DE from differential expression analysis, then using "--dge\_res path\_to\_dge\_analysis\_results.tsv" as input.

Additional to a list of genes or DE results specifying the genes to be studied, two set of files are **necessary** for the analysis to work; 

1. The feature files containing the annotation of each spot as "tumor" or "non". Look at the files exdata/ft\_00001\_X1.tsv for an example of what data should be included and the structure of it.
2. The expression count matrices, with spots as rows and genes as columns. Look at the files exdata/st\_00001\_X1.tsv for an example of what data should be included and the structure of it.

Multiple sections can be included in the analysis (this is recommended), meaning that several count matrices as feature files should be included (one pair for each section). When such is the case then the feature matrices and count matrices should
be matched w.r.t. to order. Meaning that if the feature file of section X is passed as the n:th feature file then the count
matrix of section X should be passed as the n:th count matrix.

If the naming convetion of five digits followed by an underscore and a section identifier (letter + digit) is used, proper order
will be controlled for unless you disable it with "--no\_control"; use this flag if another naming convention of the sections is used.

## Examples

Three "anonymized" sections with count data and feature files are included in the folder "exdata" which can be used to test the tool. Running the following command

```bash
./main.r --genes "ERBB2" "MZB1" "C3" "CD46"  -f exdata/ft* -c exdata/st*
```

Will  generate the following image

![Image of run 1](https://github.com/almaan/ExprByDist/blob/master/img/custom.png?raw=true)

```bash
./main.r --dge_res exdata/dge_res.tsv  -f exdata/ft* -c exdata/st* --positive_lfc --negative_lfc
```

Will generate the following image

![Image of run 2](https://github.com/almaan/ExprByDist/blob/master/img/dge_res.png?raw=true)

When --dge\_res is used rather than when specifying genes explicitly, thoes genes with a positive logFoldChange will be colored in red whilst those with a negative logFoldChange are colored in blue.

Additional features do exists, use ./main.r --help to get more information regarding these
