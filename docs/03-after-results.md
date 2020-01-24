# Next Steps

## About this chapter

1. Questions
  - How can I filter out the 'significant' genes?
  - How can I find the functions of these genes?
2. Objectives
  - Filter the genes by _p_ value
  - Find the gene annotations on an external service
3. Keypoints
  - The _p_ value we need must be corrected for the large number of genes

In this chapter we'll look at how to take our results table and get the significantly differentially expressed genes out of it.

## The results data frame

The object we created in the previous chapter `results_data` is already in the general format we need. Im going to load a version with different gene names (from Magnaporthe, not Arabidopsis) so we can work on a more familiar genome.


```r
results_data <- read.csv("sample_data/results_mo.csv")
head(results_data)
```

```
##        gene  baseMean log2FoldChange     lfcSE       stat    pvalue
## 1 MGG_00865 1022.3995      0.6385565 0.7527580  0.8482893 0.3962769
## 2 MGG_08134 1522.9579      0.3549260 0.7857086  0.4517273 0.6514655
## 3 MGG_01588  946.5248      0.6718950 0.7428144  0.9045260 0.3657165
## 4 MGG_13806 2230.1049     -0.5918810 0.6036652 -0.9804789 0.3268498
## 5 MGG_06121 1539.4032     -0.4202588 0.5376360 -0.7816790 0.4344032
## 6 MGG_06504 1384.9674     -0.4856173 0.6907157 -0.7030640 0.4820159
##        padj
## 1 0.7472062
## 2 0.7965762
## 3 0.7472062
## 4 0.7472062
## 5 0.7551582
## 6 0.7551582
```

### Which _p_ value?

Note that two columns in this data.frame have _p_ value information in them - `pvalue` and `padj`. Which is the correct one? We need the adjusted _p_ value in `padj`.

The reason for this is that a separate _p_ value was calculated in a separate test for each gene. As _p_ values have a built in error expectation (IE _p_ = 0.05 means the test will be wrong 5 percent of the time on average) then repeating the test means that using a gene level _p_ value means we get lots of errors. So we need a whole gene set _p_ value rather than a per gene value. Hence `DESeq` adjusts the _p_ value for the whole set of genes.   

## Filtering rows with significant _p_ values

To filter the rows in the data frame we can use `tidyverse` tools like `dplyr` (which we study in a separate training course). Let's keep rows from the `results_data` data frame with a `padj` lower than `0.05`


```r
library(dplyr)
significant_genes <- filter(results_data, padj < 0.05)
significant_genes
```

```
##        gene  baseMean log2FoldChange     lfcSE     stat pvalue  padj
## 1 MGG_12738  552.4170     -1.2700540 1.1960414 1.061881  0.005 0.030
## 2 MGG_01482 1912.4502      1.1010845 0.8063151 1.365576  0.005 0.004
## 3 MGG_17878  988.8704      0.9095213 0.9142498 0.994828  0.002 0.010
```

The `filter()` function works by taking the dataframe and the condition and column to filter with. 

### Filtering UP and DOWN genes

An elaboration of this is to find the up or down genes. To do this we need to build a filter on the `log2FoldChange` column. As the fold changes are encoded in a log scale, up regulated genes will have a positive value, down regulated genes will have a negative value


```r
up_genes <- filter(results_data, padj < 0.05, log2FoldChange > 0)
up_genes
```

```
##        gene  baseMean log2FoldChange     lfcSE     stat pvalue  padj
## 1 MGG_01482 1912.4502      1.1010845 0.8063151 1.365576  0.005 0.004
## 2 MGG_17878  988.8704      0.9095213 0.9142498 0.994828  0.002 0.010
```

```r
down_genes <- filter(results_data, padj < 0.05, log2FoldChange < 0)
down_genes
```

```
##        gene baseMean log2FoldChange    lfcSE     stat pvalue padj
## 1 MGG_12738  552.417      -1.270054 1.196041 1.061881  0.005 0.03
```

Note that if you want to find values that are two fold up or down regulated, then you'll need to change the `log2FoldChange` values to `1` and `-1` (as `log2(2) = 1` and `log2(0.5) = -1`).

You can export each of these tables to files with `write.csv()` as previously.


## Finding gene annotations

A common question is 'which pathways and functional categories do my genes belong to?'. Answering this requires quite an involved process, and doing it entirely in R is out of scope for this 'minimal' RNAseq tutorial. Instead of avoiding the question completely, we'll look at how to achieve a basic annotation using webtools. Specifically, the Ensembl BioMart service.

### BioMart

BioMart is a data warehouse for genomic information that can be queried through a web interface. Not all genome projects provide such a service, but the ones on Ensembl generally do. We'll work with the Magnaporthe one here, available at [https://fungi.ensembl.org/Magnaporthe_oryzae/Info/Index](https://fungi.ensembl.org/Magnaporthe_oryzae/Info/Index)

To access it we need to click `BioMart` from the top menu, and be patient, it can take a little while to load.

Then we need to follow this procedure to get a gene list annotated

  1. From the `Choose Database` drop-down select `Ensembl Fungi Genes`
  2. From the `Choose Dataset` drop down select `Magnaporthe oryzae genes (MG8)`
  3. Select the `Attributes` page and on the `External` tab tick `GO Term Name`, `GO Term Definition` and `KEGG Pathway and Enzyme ID`. These attributes are the things you will retrieve from the BioMart.
  4. Select the `Filters` page and click on the `External` tab, paste in the gene IDs of interest into the box or upload a file.
  5. Click `Results` button at the top
  
After a wait, the screen should fill with the annotations that you asked for. You can save this to a file using the `Export` options at the top.

This is all you need to make an annotated gene list.

## Further Questions

Of course, this isn't all you might want to do with your RNAseq data and gene lists. We've achieved our overall goal of getting a minimal RNAseq analysis done. What happens next will be quite different for every experiment. For example, you might want to look at seeing whether a GO Term or enzymatic pathway is enriched. Pretty much everything will be a separate analysis in itself and will require some design and planning. Please feel free to talk to the bioinformatics team when you find yourself at this stage, we'll be extremely happy to work with you!


