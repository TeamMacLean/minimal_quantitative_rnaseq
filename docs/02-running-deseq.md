# Running `DESeq2`

## About this chapter

1. Questions
  - How do I work out which genes are differentially regulated?
2. Objectives
  - Build a `DESeqDataSet` and `group` factor
  - Run `DESeq` 
3. Keypoints
  - `DESeq2` is a package for estimating differential expression
  - `DESeq2` needs you to describe the experiment in order to work

In this chapter we'll look at how to take our count matrix through `DESeq2` to estimate differential expression of genes. 

## Getting the count matrix and describing the experiment for DESeq2

### The count matrix

The object we created in the previous chapter `raw_counts` is already in the format we need. If you carried straight into this chapter from the last one, then you already have what you need. If not, you can load in the saved version (there's a copy in the sample data) as follows


```r
raw_counts <- readRDS("sample_data/raw_counts.RDS")
head(raw_counts)
```

```
##           control_rep1 control_rep2 control_rep3 treatment_rep1
## AT1G01680          670          784          548           1784
## AT1G07160         1104         1266          976            358
## AT1G07920          703          922          198           1373
## AT1G19250         1865         1654         3207           3533
## AT1G32640         1482         1266         1646           1258
## AT1G35210         1186         1416         1458            834
##           treatment_rep2 treatment_rep3
## AT1G01680           2558            368
## AT1G07160           1186           4436
## AT1G07920           1167           1726
## AT1G19250            703           2427
## AT1G32640           1690           1864
## AT1G35210            594           2684
```

## The 'grouping' object

As R is a very powerful statistical programming language, it can support analysis of some very complicated experimental designs. `DESeq` supports this behaviour and as a result we have to describe our experiment in the appropriate manner.

We need to create a `data.frame` object that states which group each column is in. A `data.frame` is basically an R analogue of an Excel sheet. We just need to work out the right order of sample types in the matrix column.

Our experiment names are in the column names of the count matrix, we can see that with the `colnames()` function.


```r
colnames(raw_counts)
```

```
## [1] "control_rep1"   "control_rep2"   "control_rep3"   "treatment_rep1"
## [5] "treatment_rep2" "treatment_rep3"
```

The controls are all in columns 1 to 3 and the treatments are in columns 4 to 6. To make the groupings we can just type in the sample types in the appropriate order and put them in a column of a `data.frame`. That looks like this


```r
grouping <- data.frame(sample_type =  c("control", "control", "control", "treatment", "treatment", "treatment"))
grouping
```

```
##   sample_type
## 1     control
## 2     control
## 3     control
## 4   treatment
## 5   treatment
## 6   treatment
```


## Running DESeq2

Now we have everything we need to run `DESeq2`. First, we must load in the library.


```r
library(DESeq2)
```

Next, we can prepare the `DESeqDataSet` object that combines all the information `DESeq2` needs to work. We run `DESeqDataSetFromMatrix()` to do this.


```r
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts, 
  colData = grouping, 
  design = ~ sample_type)
```

Here we set the arguments 

  1. `countData` which is the actual data, so gets our `raw_counts`
  2. `colData` which tells the group each data column is in so gets `grouping`
  3. `design` is an R-ish way of describing the experiment design, for a standard exoeriment like this you use the `~` and the name of the `grouping` column
  
Don't worry too much about whether the `design` argument makes sense at this stage, its a bit out of scope to discuss the way R expects experimental designs for now. Follow the pattern you see here until you have a really complex design and have motivation to come back to it.

Finally, we can do the `DESeq` analysis. We have a single function for this and all it needs is our prepared data.


```r
de_seq_analysed <- DESeq(dds)
```

And now we can extract the results with the helpful `results()` function. This needs the `contrast` to be described, basically the column name and the types.

The types are ordered so that the first mentioned is the measurement of interest (ie the `treatment`) and the second is the baseline to which it is compared (here `control`). If you get the two the wrong way round, your up-regulated genes will look down-regulated and vice-versa, so take time to check.
 



```r
results_data <- results(de_seq_analysed, contrast = c('sample_type', 'treatment', 'control'))
head(results_data)
```

```
## log2 fold change (MLE): sample_type treatment vs control 
## Wald test p-value: sample type treatment vs control 
## DataFrame with 6 rows and 6 columns
##                   baseMean     log2FoldChange             lfcSE
##                  <numeric>          <numeric>         <numeric>
## AT1G01680 1022.39953137642  0.638556501008842 0.752757958605879
## AT1G07160 1522.95786602724  0.354926012092298 0.785708629982376
## AT1G07920 946.524794250111  0.671894996505363 0.742814441576441
## AT1G19250 2230.10488838027 -0.591880977166901 0.603665159459885
## AT1G32640 1539.40321179978 -0.420258753438091 0.537635950705001
## AT1G35210 1384.96737134404 -0.485617346461465 0.690715741846472
##                         stat            pvalue              padj
##                    <numeric>         <numeric>         <numeric>
## AT1G01680  0.848289272412955 0.396276890458589 0.747206163566895
## AT1G07160   0.45172726701533 0.651465472435186 0.796576241413894
## AT1G07920  0.904526028168531 0.365716539229272 0.747206163566895
## AT1G19250 -0.980478942492676 0.326849759032812 0.747206163566895
## AT1G32640 -0.781679039296026 0.434403223096372 0.755158226458792
## AT1G35210 -0.703063962554663 0.482015889229017 0.755158226458792
```

We get a lot of information back in this column. We can see in amongst all that the important log fold change estimates and the adjusted p-value. Effectively, our analysis is done, we have our differential expression estimates, though we do need to do more to answer questions of interest. That's what we'll do in the next chapter

## Saving the results

As a final step, we can save the results to a CSV file. As in the earlier chapter we can do this with `write.csv()`


```r
write.csv(results_data, "sample_data/results.csv")
```

