# NOJAH
NOt Just Another Heatmap

This interactive web application: NOt Just Another Heatmap (NOJAH) is developed in R with Shiny to provide comprehensive workflows to


1) Perform Genome-Wide Heatmap (GWH) analysis of cancer genomic data sets 
2) Perform Cluster results Cluster (CrC) analysis for upto three genomic platforms  
3) Perform Significance of Cluster (SoC) analysis using a robust bootstrap approach

The goal of this tool is to provide a one stop shop to perform genomic analyses. It helps users extract the top most variable gene sets from Genome-Wide data, estimated the optimal number of clusters, identify a core subset of samples on their own data. This tool is not restricted to expression data but can use genomic data of any type: e.g. Variant Allelic frequencies, methylation, copy number segmentation mean values, etc. 

Individual platform data from expression, methylation, copy number for the same samples (in the same order) in each platform can be integrated in the combined results clustering workflow. Typically, the most variable gene/variant/cpg sites subset from each platfrom is used as input. The clustering results from each platform are combined for a second level clustering analysis using the ConsensusClusterPlus bioconductor package. 

For those interested in heatmap cluster analysis, using the Significane of cluster analysis tab, NOJAH lets you generate a HeatMap with only the minimum knowledge of heatmaps and minimal coding experience. The user may upload their own data or use the example TCGA Breast Cancer RNASeq Expression to generate the HeatMap using any distance or clustering method of their choice. Along with the HeatMap, the column and row dendrograms are displayed individually in separate tabs. If you wish to know the samples within each cluster, this information can be downloaded using the cutree download button inside the row/column dendrogram tab. The novelty about this app is that it is highly reproducible and provides minimum information required to reproduce the heatmaps.

This tool can be accessed using http://bbisr.shinyapps.winship.emory.edu/NOJAH/.


For large datasets, we advise using the command line version of NOJAH. For running NOJAH using RStudio,

#### INSTALLATION
Firstly, you should have the most recent version of R or RStudio.
Next install required packages. Cut and paste what's below in an R session.

``` R
install.packages(c("shiny", "shinythemes", "plotly", "RColorBrewer", "gplots", "gdata",  "plyr",  "DT", "ggplot2",  "dendextend", "colourpicker", "cluster", "reshape", "gridExtra", "DiagrammeR", "changepoint", "rhandsontable", "matrixStats"))
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("ConsensusClusterPlus")
```
You only need to do this once.

Then, you may run NOJAH any time in an R session as follows.
``` R
library(shiny)
library(shinythemes)
library(plotly)
library(RColorBrewer)
library(gplots)
library(gdata)
library(plyr)
library(DT)
require(dendextend)
library(colourpicker)
library(ConsensusClusterPlus)
library(cluster) # for silhouette
library(reshape) # for cast
library(gridExtra)
library(DiagrammeR)
library(changepoint)
library(rhandsontable)
library(ggplot2)

runGitHub("NOJAH", "bbisr-shinyapps")
```

It may take upto several minutes for the webpage to show in your browser depending on your computer.

#### INPUT DATA REQUIREMENTS

Data should be input as a TXT or a CSV file. The first two rows of the data file have information about the patients/specimens and their response/subtype; all remaining rows have gene expression data, one row per gene. In the case of Microarray gene expression data in which there are several probes corresponding to a single gene, a unique identifier would need to be created to separately identify each probe such as, 'Gene 1_p1', 'Gene1_p2' indicating Gene 1 has two probes. The columns represent the different experimental samples. A maximum of up to 10 different sample groups and 6 different gene groups may be used with this tool.

##### DATA FORMAT

1.	The first line of the file contains the gene identifier 'gene_id' (column 1), gene group classification 'Groups' (column 2) followed by the patient IDs e.g. TCGA.01.1A2B, one per column, starting at column 3. Column 1 gene identifier has to be labelled 'gene_id' and column 2 header should be labelled 'Groups' for using this tool. Other titles may cause the program to display errors. 
2.	The second line of the file contains the patient response classification e.g. Fav/Unf for favorable outcome group vs the unfavorable outcome group or Normal/Tumor, etc., in alphabetical order, starting at column 3. The first two columns for this row should be blank. Additional sample information such as subtype should be included under each sample. Whereever data is missing, it should be coded to none. Leaving field blank may cause the program to cause errors. The first column for the additional clinical field includes the label (for example 'Subtype'), but the second column should be blank.
3.	The remaining lines contain gene expression measurements one line per gene, described in the format below.
a) Column_1. This should contain the gene name, for the user's reference. Each gene  name should be unique. When using microarray data, the gene and the probe id can be merged using a delimitor (except '|') to make a unique name.  Delimitors such as >,;:#%&(!)_+ are acceptable.
b) Column_ 2. This should contain the gene group classification e.g. O/U for Over-expressed/Under-expressed or Hyper/Hypo for hypermethylated/hypomethylated in alphabetical order. If only one gene group, use any alphabet e.g. A or na for each row instead. 
c) Remaining Columns. These should contain the expression measurements as numbers. Data inputted should be non-negative. Columns and rows with zero variance should be removed from the data. Rows containing missing expression measurements, should be also be removed from the input data or it will cause the tool to run into errors.

| gene_id           | Groups | GSM9981   | GSM1870  | GSM4618 | GSM7689  | GSM8772 | GSM1121  | GSM1250 | GSM3112  | GSM4987 | GSM1277 |
| -------------     |-------:| ---------:| --------:|--------:|---------:|--------:|---------:|--------:|---------:|--------:|--------:|
|                   |        | MM        | MM       | MM      | MM       | MM      | MUGS     | MUGS    | NPC      | SM      | SM      |
| Subtype           |        | Classical | Classical| Neural  | Pro-neural | Classical | Mesenchymal | Classical | Neural | Neural | Neural|
| YWHAE>210996_s_at | na      | 1.47      |  2.18    | 5.87    |	9.12     |	7.34   | 1.56     |	3       |	7.77     |	3.4    |	1.56   |
| YWHAE>201020_at   | na     | 1.98      |  7.93    | 2.76	  | 9.11     |	8.46   | 0.98     |	5.98    |	8.19     |	8.91   |	5.98   |
| YWHAH>33323_r_at  | na      | 8.02      |  8       | 3.19	  | 11.86    |	6.54   | 8.17     |	2       |	0.99     |	2      |	1.17   |


Example dataset for one gene group (marked na) and four patient groups (MM, MUGS, NPC and, SM). Additonal subtype information is appended above the numeric data. Numeric data starts at row four and column three.

##### Terms of Use

This tool was prepared by members of the Winship Biostatistics and Bioinformatics Shared Resource (BBISR) of Emory University. 
Use of either should properly acknowledge the Winship BBISR in publications, abstracts, presentations, posters, grant proposals, etc. by using the following text	

'Research reported in this publication was supported in part by the Biostatistics and Bioinformatics Shared resource of Winship Cancer Institute of Emory University and NIH/NCI under award number P30CA138292. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. 

Authors- Manali Rupji, dual M.S., Bhakti Dwivedi Ph.D. & Jeanne Kowalski Ph.D.

Maintainer- Manali Rupji 'manali(dot)rupji(at)emory(dot)edu'
