# miRCancer  
A database for microRNA-gene/protein correlation in cancer.  
A web interactive interface for this database is available [here](https://mahshaaban.shinyapps.io/miRCancer/)

# Overview  
miRCancer is database of expression correlation of expression of microRNAs and genes/proteins.
The web interface provides an easy access to the database without having to write any code.
Among other features, users can search, choose, visualize and downlaod the particular subsets of the data.  

# Datasets  
The expression data was obtained from The Cancer Genome Atlas [(TCGA)](https://tcga-data.nci.nih.gov) using [RTCGA](http://bioconductor.org/packages/release/bioc/html/RTCGA.html) R package.
And the targets data were obtained form [TargetScan](http://www.targetscan.org) using [targetscan.Hs.eg.db](http://www.bioconductor.org/packages/devel/data/annotation/html/targetscan.Hs.eg.db.html) [Bioconductor](http://bioconductor.org) package.  

# Workflow
[(TCGA)](https://tcga-data.nci.nih.gov) expression data for microRNAs, mRNA and proteins were obtained using [RTCGA](http://bioconductor.org/packages/release/bioc/html/RTCGA.html) R package. 
The expression correlation were calculated for each microRNA to each feature; gene/protein in each [(TCGA)](https://tcga-data.nci.nih.gov) study.
An SQLite database were build using these corelations. In addition, expression profiles of each miroRNA and
microRNA targets data were obtained from [targetscan.Hs.eg.db](http://www.bioconductor.org/packages/devel/data/annotation/html/targetscan.Hs.eg.db.html) [Bioconductor](http://bioconductor.org) package and added to the database.

# Installation  
Clone or download this repository, navigate to the miRCancer directory and run make to build the database from scratch and
launch a shiny browser application.

```
git clone https://github.com/MahShaaban/miRCancer
cd miRCancer
make
```

Indeed, the same could be achieved through Rstudio and the sqlite db file can be accessed using DBI library in R.

## More  
Citation: <manuscript>  
Usign cRegulome: https://github.com/MahShaaban/cRegulome  