################################################################################################################
# File Name: Build_script.R
# Description: This script calls custom functions to extract, tidy and process TCGA data in the RTCGA package
# Created by: Mahmoud Ahmed <mahmoud.s.fahmy@students.kasralainy.edu.eg>
# Usage: Used with build_script.R to build the miRCancer.db
################################################################################################################

# load/require libraries
# load/require RTCGA and data pacakges
require(RTCGA)
require(RTCGA.miRNASeq)
require(RTCGA.rnaseq)
require(RTCGA.RPPA)

# load/require annotation packages
require(targetscan.Hs.eg.db)
require(org.Hs.eg.db)

# load/require packages for data wrangling
require(dplyr)
require(purrr)
require(reshape2)
require(stringr)
require(data.table)

# load/require packages for creating the db
require(RSQLite)

# sourece the custom functions
source('R/build_functions.R')

# It's sellf-explainatory if you ask me. But just in case I need to use it a week latter :)
# First, a call to mirna_info to get only TCGA studies with all three assays performed.
# Second, a call to cor_write. This is where the magic happen, but it takes a while so hold on.
#         Meanwhile you can read what the function does in build_functions.R. cor_write cheats
#         and write the microRNA expression profiles to files on the way.
# Third, read the cor_write output and write it to a db, a sqlite db.
# Finlly, get the targets data and write them to the db, the same sqlite db.

# get cohorts
cohorts <- mirna_info()$Cohort[-20] %>%
  as.character

# make correlations
dir.create('tmp')
map(cohorts, cor_write)

# tidy and write tables to the db
db <- dbConnect(SQLite(),
                'miRCancer.db')

## cor_rnaseq
cor_rnaseq <- cor_tidy('rnaseq') %>%
  mutate_at(vars(3:36), 
            function(x) x * 100)
mirna <- unique(cor_rnaseq$mirna_base)
dbWriteTable(db, 
             name = 'cor_rnaseq',
             value =  cor_rnaseq, row.names = FALSE)
dbSendQuery(db,
            statement = 'create index idx1 on cor_rnaseq (mirna_base);')
rm(cor_rnaseq)

## cor_rppa
cor_rppa <- cor_tidy('rppa') %>%
  mutate_at(vars(3:35), 
            function(x) x * 100)
mirna <- unique(c(cor_rppa$mirna_base, mirna))

dbWriteTable(db, 
             name = 'cor_rppa',
             value =  cor_rppa,
             row.names = FALSE)
dbSendQuery(db,
            statement = 'create index idx2 on cor_rppa (mirna_base);')
rm(cor_rppa)

## profiles
fls <- list.files('tmp', 
                  pattern = 'profile', 
                  full.names = TRUE)
profiles <- lapply(fls, 
                   read.csv, 
                   stringsAsFactors = FALSE)
names(profiles) <- str_split(fls, 
                             pattern = '/|_',
                             simplify = TRUE)[, 2]
profiles <- bind_rows(profiles,
                      .id = 'study')
dbWriteTable(db, 
             name = 'profiles',
             value =  profiles,
             row.names = FALSE)
dbSendQuery(db, 
            statement = 'create index idx3 on profiles (mirna_base);')
rm(profiles)

# targets
targets <- list()
targets$gene <- get_targets('gene')
targets$ptn <- get_targets('protein')

targets <- bind_rows(targets, .id = 'feature_type') %>%
  filter(mirna_base %in% mirna)

dbWriteTable(db,
             name = 'targets',
             value = targets,
             row.names = FALSE)
dbSendQuery(db, 
            statement = 'create index idx4 on targets (mirna_base);')
rm(targets)

dbDisconnect(db)

# application helpers
## save cohorts to www studies.txt
dir.create('www')
writeLines(cohorts, 'www/studies.txt')
writeLines(mirna, 'www/miRBase.txt')
download.file('https://github.com/rstudio/shiny-examples/blob/master/036-custom-input-control/www/chooser-binding.js?raw=1',
              destfile = 'www/chooser-binding.js')
download.file('https://github.com/rstudio/shiny-examples/blob/master/036-custom-input-control/chooser.R?raw=1',
              destfile = 'www/chooser.R')
download.file('http://www.studyinkorea.go.kr/file/imgpreview.do?filename=510040_logo.jpg&fileStorePath=fileStorePath',
              destfile = 'www/logo.jpg')
