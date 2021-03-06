--- 
title: "A minimal quantitative RNAseq pipeline"
author: "Dan MacLean"
date: "2020-01-30"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: TeamMacLean/minimal_quant_rna_seq
description: "A minimal guide to performing quantitative RNAseq"
---

# About this course

In this short course we'll look at a method for getting quantitative estimates of gene expression from RNAseq data. The course assumes that you will already have performed a read alignment so is _not_ a 'read to results' course. The course is very brief and will show you how to use a perform a common pipeline centered around `DESeq` in `R` and `RStudio`

I acknowledge that there are lots of other programs and methods - this course is _not_ meant to be comprehensive, it is meant to get you being productive. Seek out further advice if you need to run other programs or systems. Do be encouraged though, lots of what you learn here will be applicable to other pipelines for the same job (they all run in a similar manner with similar objects) so this is a good place to start. 

The course is intended to run on your 'local' machine, that is to say, your laptop or desktop computer. In general these machines will be powerful enough for most datasets though the pipeline we will learn can be easily adapted for a high performance computing environment if you need greater computational power.

## Prerequisites

This course assumes that you are a little familiar with the basics of running R and R commands from the R console. You'll need to know the basics of typing in commands and getting output, not much more. 


### R and RStudio

#### Installing R

Follow this link and install the right version for your operating system [https://www.stats.bris.ac.uk/R/](https://www.stats.bris.ac.uk/R/)

#### Installing RStudio

Follow this link and install the right version for your operating system [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)

#### Installing R packages in RStudio.

You'll need the following R packages:

  1. devtools
  2. atacR
  3. DESeq
  
For simplicity, install them in that order.

To install `devtools`:

Start RStudio and use the `Packages` tab in lower right panel. Click the install button (top left of the panel) and enter the package name `devtools`, then click install as in this picture

![Installing Packages](fig/package_install.png)

To install `atacR`:

Type the following into the RStudio console, `devtools::install_github("TeamMacLean/atacr")`

To install `DESeq`:

Type the following into the RStudio console, `BiocManager::install("DESeq")`


Now you are done! Everything is installed ready for you to work with. Next we need to get the sample data

### Sample BAM file and counts files

You'll need this zip file of data: [sample_data.zip](https://github.com/TeamMacLean/minimal_quantitative_rnaseq/blob/master/sample_data/sample_data.zip) which contains a lot of ifles including BAM files and counts. Download it, extract the files and put them into a folder on your machine. I suggest something like `Desktop/align_tut`. This will be the directory we'll work from in the rest of the course.

That's all you need to do the lesson. If you have any problems getting this going, then ask someone in the Bioinformatics Team and we'll help.



