# SEARCHIN

##  Systematic Elucidation and Assessment of Regulatory Cell-to-cell Interaction Networks (SEARCHIN)

## Description

This tutorial guides you through the application of the SEARCHIN methdology for the identification of cell-to-cell ligand-receptor interactions responsible for mediating non-cell-autonomous communication mechanisms. 

## Schematics of the Methodology

![Searching Schematics](/data/images/searchin-schematics.jpg)

## Introduction

[TODO]

### Using SEARCHIN - Example
To reproduce the manuscript findings you can run: `source(sources/searchin.R)` obtaining the published findings in one step.
Alternatively, you can interactively reproduce the analysis in very few steps. 
Basically, you only have to load the libraries and run one command after having loaded the input data.

#### Pipeline Commands

The *first step* consists in loading three datasets:
1. *ligands_list*, is a list with the IDs of the ligands from the LC-MS experiment
2. *modulators_table*, is a data.frame/tibble having as mandatory columns the modulator IDs and the p-value, from CINDy/MINDy analysis
3. *receptors_table*, is a data.frame/tibble having as mandatory columns the receptor IDs and the p-value, from VIPER analysis

We also need to set a *seed* for the random number generator that allows you to reproduce the results. Specifically, to reproduce the manuscript's findings you don't need a running Neo4j server, because we use a caching system that recover from file PPIs previously queried from a running server. You may need a Neo4j running server to query for PPIs pvalues that are not present in the cache. If you have a Neo4J running server with the PrePPI interaction network loaded, please set the flag `.isNeo4jServerOffline <- TRUE`. For further details, please see [Neo4J Server Setup](#neo4j-server-setup)

Here is the first step:
```R
rm(list = ls())
set.seed(666)
.isNeo4jServerOffline <- TRUE	# Flag this as FALSE if you have installed the Neo4j server or if you want to run SEARCHIN with your own data: you need the server with the PPI data running
source("libs/searchin-class.R")	# Load SEARCHING libraries and utilities

# Loading input data
tmp <- read_csv("input/putative_ligands_mouse_new_data_jul2019.txt",col_types = c("cccccccc"))
ligands_list <- tmp$GeneSymbol
modulators_table <- read_csv("input/cindy-analysis.csv", col_types = "ccd")
receptors_table <- read_delim("input/viper-analysis.csv", delim = "\t" , col_types = "ccdddd")
```
The **second and last** step consists in running the algorithm with the loaded data that presents a table with the ranked interactions as output.
Here is how to run the pipeline as one command-only line of code that performs all the steps in cascade. However, each step can be performed as separate individual line of code.

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


#### Protein-Protein Interaction (PPI) empirical null model using PrePPI

We used Neo4J Graph Database Community Edition version 3.2.5 to build a network of human protein interaction consisting of 20,259 proteins and 203,952,081 interactions with a total store size of 14.28 GB.

##### Neo4j Server Setup

* To run Neo4j Server you need to install the software from here: [Neo4j Server App](https://neo4j.com/download/)
* To install the R package to interface with Neo4j Server, you can follow this guide: [RNeo4j Guide](https://github.com/nicolewhite/RNeo4j)

#### Generation of the VIPER table used as input data

Here is the code to generate the VIPER table: [Code to generate the manuscript's VIPER table](sources/examples/als-tables-generation/viper-table-generator.R)

### Troubleshooting

[TODO]

#### Required Libraries

```R
require(tidyverse) 
require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(RNeo4j)
require(sloop)
require(Qvalue)
require(RobustRankAggreg)

install.packages('BiocManager')
BiocManager::install('tidyverse')
BiocManager::install('org.Hs.eg.db')
BiocManager::install('org.Mm.eg.db')
BiocManager::install('devtools')
devtools::install_github("nicolewhite/RNeo4j")
```

