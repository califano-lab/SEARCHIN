# SEARCHIN

##  Systematic Elucidation and Assessment of Regulatory Cell-to-cell Interaction Networks (SEARCHIN)







## Introduction



#### Protein-Protein Interaction (PPI) empirical null model using PrePPI

We used Neo4J Graph Database Community Edition version 3.2.5 to build a network of human protein interaction consisting of 20259 proteins and 203,952, 081 interactions with a total store size of 14.28 GB.



### Using SEARCHIN - Example

#### Pipeline Commands

```R
  my_sin <- SEARCHIN()
  
  x <- read_csv(filename)
  my_sin <- addLigandslist(my_sin,x$GeneSymbol)
  x <- read_csv("/Users/av2729/Workspace/SEARCHIN/data/input/cindy_table.csv", col_types = "ccd")
  my_sin <- addModulatorTable(my_sin,x)
  x <- read_delim("/Users/av2729/Workspace/SEARCHIN/data/input/marina-analysis-on-membrane-receptors-for-paper.txt", delim = "\t" , col_types = "ccdddd")
  my_sin <- addReceptorTable(my_sin,x)
  
  my_sin <- runModel(my_sin)
  my_sin <- generateRank(my_sin)
```



#### Output Example

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

