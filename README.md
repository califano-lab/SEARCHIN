# SEARCHIN

##  Systematic Elucidation and Assessment of Regulatory Cell-to-cell Interaction Networks (SEARCHIN)

## Description

This is the tutorial that guides you thorugh the application of the SEARCHIN methdology for the identification of cell-to-cell ligand-receptor interactions responsible of non-cell-autonomous communication mechanisms. 



![Searching Schematics](/images/searchin-schematics.jpg)
Format: ![Alt Text](url)

## Introduction

Cell-to-cell etc., XXX

#### Protein-Protein Interaction (PPI) empirical null model using PrePPI

We used Neo4J Graph Database Community Edition version 3.2.5 to build a network of human protein interaction consisting of 20,259 proteins and 203,952,081 interactions with a total store size of 14.28 GB.

### Using SEARCHIN - Example

#### Pipeline Commands

 The **first step** consists in setting a seed for the random number generator that allows us to reproduce the code in the same way 

```R
rm(list = ls())
set.seed(666)
source("libs/searchin-class.R")	# Load SEARCHING libraries and utilities
```

The **second step** is loading the input files. Here, is the code to load the data from our manuscript that can reproduce the findings. 

```R
.isNeo4jServerOffline <- TRUE	# Flag this as FALSE if you have installed the Neo4j server or if you want to run SEARCHIN with your own data: you need the server with the PPI data running

## Reading TXT file with the list of the ligands ----
tmp <- read_csv("input/putative_ligands_mouse_new_data_jul2019.txt",col_types = c("cccccccc"))
ligands_list <- tmp$GeneSymbol
## Reading CSV file with a table with the modualtors analysis ----
modulators_table <- read_csv("input/cindy-analysis.csv", col_types = "ccd")
## Reading CSV file with a table with the VIPER analysis ----
receptors_table <- read_delim("input/viper-analysis.csv", delim = "\t" , col_types = "ccdddd")
```

The **third and last step** is to run the pipeline as one command-only line of code that performs all the steps in cascade. However, each step can be performed as separate and individual line of code.

```R
my_sin <- SEARCHIN() %>% 
  addLigandslist(ligands_list) %>% 
  addModulatorTable(modulators_table) %>% 
  addReceptorTable(receptors_table) %>% 
  runModel() %>%
  generateRank()
```

By printing the SEARCHIN object, you can see the table content with the ranked ligand-receptor interaction inferred by the algorithm.

##### Output Example

```shell
> my_sin
An object of class SEARCHIN 
------------------------------
>>> Using 69 ligands
>>> Using 214 modulators
>>> Using 873 receptors
>>> Using 1503 PPIs interactions
# A tibble: 543 x 3
   interaction   pvalue  rank
   <chr>          <dbl> <int>
 1 Nme1-Ptprn   0.00508     1
 2 Cdh2-Ryk     0.00588     2
 3 App-Tnfrsf21 0.00674     3
 4 Lgals1-Spn   0.0165      4
 5 Mif-Tgfbr2   0.0165      5
 6 App-Adgrl1   0.0165      6
 7 App-Glrb     0.0182      7
 8 Cxcl5-Cx3cl1 0.0185      8
 9 Ccn2-Itgb5   0.0187      9
10 Ccl3-Cx3cl1  0.0200     10
# … with 533 more rows
------------------------------
```

### Troubleshooting

#### Required Libraries

```R
require(tidyverse) 
require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(RNeo4j)
```

