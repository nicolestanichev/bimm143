# Class 14: RNASeq mini project
Nicole (PID: A18116280)

- [Background](#background)
  - [Remove zero count genes](#remove-zero-count-genes)
- [DESeq analysis](#deseq-analysis)
- [Data Visualization](#data-visualization)
- [Add Annotation](#add-annotation)
- [Pathway Analysis](#pathway-analysis)
  - [GO terms](#go-terms)
  - [Reactome](#reactome)
- [Save Our Results](#save-our-results)

## Background

Here we work through a complete RNASeq analysis project. The input data
comes from a knock-down experiment of a HDX gene. \## Data Import

Reading the `count()` and `metadata` CSV files

``` r
counts <- read.csv("GSE37704_featurecounts (1).csv", row.names = 1)
metadata <- read.csv("GSE37704_metadata.csv")
```

Check on data structure

``` r
head(counts)
```

                    length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ENSG00000186092    918         0         0         0         0         0
    ENSG00000279928    718         0         0         0         0         0
    ENSG00000279457   1982        23        28        29        29        28
    ENSG00000278566    939         0         0         0         0         0
    ENSG00000273547    939         0         0         0         0         0
    ENSG00000187634   3214       124       123       205       207       212
                    SRR493371
    ENSG00000186092         0
    ENSG00000279928         0
    ENSG00000279457        46
    ENSG00000278566         0
    ENSG00000273547         0
    ENSG00000187634       258

``` r
metadata
```

             id     condition
    1 SRR493366 control_sirna
    2 SRR493367 control_sirna
    3 SRR493368 control_sirna
    4 SRR493369      hoxa1_kd
    5 SRR493370      hoxa1_kd
    6 SRR493371      hoxa1_kd

``` r
head(metadata)
```

             id     condition
    1 SRR493366 control_sirna
    2 SRR493367 control_sirna
    3 SRR493368 control_sirna
    4 SRR493369      hoxa1_kd
    5 SRR493370      hoxa1_kd
    6 SRR493371      hoxa1_kd

Some book-keeping is required as there looks to be a mis-match between
metadata and counts columns

``` r
ncol(counts)
```

    [1] 7

``` r
nrow(metadata)
```

    [1] 6

Looks like we need to get rid of the first “length” column of our
`counts` object.

``` r
cleancounts <- counts[ , -1]
```

``` r
colnames(cleancounts)
```

    [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"

``` r
metadata$id
```

    [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"

``` r
all( colnames(cleancounts) == metadata$id)
```

    [1] TRUE

### Remove zero count genes

There are lots of genes with zero counts. We can remove these from
further analysis

``` r
head(cleancounts)
```

                    SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
    ENSG00000186092         0         0         0         0         0         0
    ENSG00000279928         0         0         0         0         0         0
    ENSG00000279457        23        28        29        29        28        46
    ENSG00000278566         0         0         0         0         0         0
    ENSG00000273547         0         0         0         0         0         0
    ENSG00000187634       124       123       205       207       212       258

``` r
to.keep.inds <- rowSums(cleancounts) > 0
nonzero_counts <- cleancounts[to.keep.inds,]
```

## DESeq analysis

Load the package

``` r
library(DESeq2)
```

Setup DESeq object

``` r
dds <- DESeqDataSetFromMatrix(countData = nonzero_counts,
                              colData = metadata,
                              design = ~condition)
```

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

Run DESeq

``` r
dds <-DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

get results

``` r
res <- results(dds)
head(res)
```

    log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    Wald test p-value: condition hoxa1 kd vs control sirna 
    DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE       stat      pvalue
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ENSG00000279457   29.9136      0.1792571 0.3248215   0.551863 5.81042e-01
    ENSG00000187634  183.2296      0.4264571 0.1402658   3.040350 2.36304e-03
    ENSG00000188976 1651.1881     -0.6927205 0.0548465 -12.630156 1.43993e-36
    ENSG00000187961  209.6379      0.7297556 0.1318599   5.534326 3.12428e-08
    ENSG00000187583   47.2551      0.0405765 0.2718928   0.149237 8.81366e-01
    ENSG00000187642   11.9798      0.5428105 0.5215598   1.040744 2.97994e-01
                           padj
                      <numeric>
    ENSG00000279457 6.86555e-01
    ENSG00000187634 5.15718e-03
    ENSG00000188976 1.76553e-35
    ENSG00000187961 1.13413e-07
    ENSG00000187583 9.19031e-01
    ENSG00000187642 4.03379e-01

## Data Visualization

``` r
library(ggplot2)

ggplot(res) +
  aes(log2FoldChange, -log(padj) ) +
  geom_point()
```

    Warning: Removed 1237 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class14_files/figure-commonmark/unnamed-chunk-17-1.png)

Add threshold lines for fold-change and P-value and color our subset of
genes that make these threshold cut-offs in the plot

``` r
mycols <- rep("gray", nrow(res))
mycols[ abs(res$log2FoldChange) > 2] <- "blue"
mycols[ res$padj > 0.05 ] <- "gray"

ggplot(res) +
  aes(log2FoldChange, -log(padj), color = mycols) +
  geom_point() +
  geom_vline(xintercept = c(-2, 2), color = "red") +
  geom_hline(yintercept = -log(0.05), color = "red") +
  scale_color_identity()
```

    Warning: Removed 1237 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class14_files/figure-commonmark/unnamed-chunk-18-1.png)

## Add Annotation

Add gene symbols and entrez ids

``` r
library(AnnotationDbi)
library(org.Hs.eg.db)
```

``` r
res$symbol <- mapIds(x=org.Hs.eg.db, 
                     keys=row.names(res),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$symbol <- mapIds(x=org.Hs.eg.db, 
                     keys=row.names(res),
                     keytype = "ENSEMBL",
                     column = "ENTREZID" )
```

    'select()' returned 1:many mapping between keys and columns

## Pathway Analysis

Run gage analysis

``` r
library(gage)
library(gageData)
library(pathview)
```

We need a named vector of fold-change values as input for gage

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

    [1]  0.17925708  0.42645712 -0.69272046  0.72975561  0.04057653  0.54281049

``` r
data(kegg.sets.hs)

keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

``` r
head(keggres$less, 2)
```

                                             p.geomean stat.mean p.val q.val
    hsa00232 Caffeine metabolism                    NA       NaN    NA    NA
    hsa00983 Drug metabolism - other enzymes        NA       NaN    NA    NA
                                             set.size exp1
    hsa00232 Caffeine metabolism                    0   NA
    hsa00983 Drug metabolism - other enzymes        0   NA

``` r
pathview(pathway.id = "hsa04110", gene.data=foldchanges)
```

    Warning: None of the genes or compounds mapped to the pathway!
    Argument gene.idtype or cpd.idtype may be wrong.

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/nicolestanichev/Desktop/BIMM143/class15/class14

    Info: Writing image file hsa04110.pathview.png

![](hsa04110.xml)

### GO terms

Same analysis but using GO genesets rather than KEGG

``` r
data(go.sets.hs)
data(go.subs.hs)

# Focus on Biological Process subset of GO
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets)
```

``` r
head(gobpres$less, 4)
```

                                                p.geomean stat.mean p.val q.val
    GO:0000002 mitochondrial genome maintenance        NA       NaN    NA    NA
    GO:0000003 reproduction                            NA       NaN    NA    NA
    GO:0000012 single strand break repair              NA       NaN    NA    NA
    GO:0000018 regulation of DNA recombination         NA       NaN    NA    NA
                                                set.size exp1
    GO:0000002 mitochondrial genome maintenance        0   NA
    GO:0000003 reproduction                            0   NA
    GO:0000012 single strand break repair              0   NA
    GO:0000018 regulation of DNA recombination         0   NA

### Reactome

Lots of folks like the reactome web interface. You can also run this as
an R function but lets look at the website first \<
https://reactome.org/ \>

The website wants a text file with one gene symbol per line of the genes
you want to map to pathways.

``` r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), ]$symbol
head(sig_genes) #res$symbol
```

    ENSG00000187634 ENSG00000188976 ENSG00000187961 ENSG00000188290 ENSG00000187608 
           "148398"         "26155"        "339451"         "57801"          "9636" 
    ENSG00000188157 
           "375790" 

and write out to a file:

``` r
write.table(sig_genes, file="significant_genes.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
```

## Save Our Results

``` r
write.csv(res, file="myresults.csv")
```
