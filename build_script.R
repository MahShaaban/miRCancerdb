tictoc::tic()
library(sqlome)

library(RSQLite)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(purrr)
library(readr)

library(RTCGA)
library(RTCGA.rnaseq)
library(RTCGA.miRNASeq)
library(RTCGA.RPPA)

# get the the TCGA cohorts with all three assays
cohorts <- as.character(sqlome_info()$Cohort)

# open a connection to a db
db <- dbConnect(SQLite(), 'miRCancer.db')

# microRNA gene correlaitons
## make cor_mir table
### calculate correlations
cat('1: Calculating microRNA-gene correlations.\n')

df <- map(cohorts[-20], function(x) {
    # make names of RTCGA data.frames
    mi <- paste(x, 'miRNASeq', sep = '.')
    m <- paste(x, 'rnaseq', sep = '.')

    # read data.frames
    mi <- get(mi)
    m <- get(m)

    # tidy data.frames
    mi <- mirna_tidy(mi)
    m <- mrna_tidy(m)

    # calcualte correlations in a tidy data.table
    corr <- cor_make(mi, m, x, tidy = TRUE)
    corr <- as.data.table(corr)

    # print progress
    cat(paste('microRNA-gene correlation for', x, 'is done.\n'))

    # return tidy data.table
    return(corr)
})

### use reduce to merge data.tables
cat('2: Merging microRNA-gene correlations.\n')
df <- Reduce(function(x, y) merge(x, y, all=TRUE), df)

### write cor_mir table to connection db
cat('3: Writing microRNA-gene correlations.\n')
dbWriteTable(db,
             name = 'cor_mir',
             df,
             overwrite = TRUE)
### making index on cor_mir
dbSendQuery(db,
            statement = 'create index idx1 on cor_mir (mirna_base);',
            overwrite = TRUE)

## make targets for/microRNA-gene mapping
cat('4: Extracting microRNA-gene targets.\n')

targets <- list()
targets$genes <- get_targets(unique(df$mirna_base), 'gene')

# microRNA gene correlaitons
## make cor_rppa table
### calculate correlations
cat('5: Calculating microRNA-protein correlations.\n')

df <- map(cohorts[-20], function(x) {
    # make names of RTCGA data.frames
    mi <- paste(x, 'miRNASeq', sep = '.')
    m <- paste(x, 'RPPA', sep = '.')

    # read data.frames
    mi <- get(mi)
    m <- get(m)

    # tidy data.frames
    mi <- mirna_tidy(mi)
    m <- rppa_tidy(m)

    # calcualte correlations in a tidy data.table
    corr <- cor_make(mi, m, x, tidy = TRUE)
    corr <- as.data.table(corr)

    # print progress
    cat(paste('microRNA-protein correlation for', x, 'is done.\n'))

    # return tidy data.table
    return(corr)
})

### use reduce to merge data.tables
cat('6: Merging microRNA-protein correlations.\n')

df <- Reduce(function(x, y) merge(x, y, all=TRUE), df)

### write cor_rppa table to connection db
cat('7: Writing microRNA-protein correlations.\n')

dbWriteTable(db,
             name = 'cor_rppa',
             df,
             overwrite = TRUE)

### making index on cor_rppa
dbSendQuery(db,
            statement = 'create index idx2 on cor_rppa (mirna_base);',
            overwrite = TRUE)

## make targets for/microRNA-gene mapping
cat('8: Extracting microRNA-gene targets.\n')
targets$protein <- get_targets(unique(df$mirna_base), 'protein')

# merge targets tables
cat('9: Merging microRNA gene and protein targets.\n')

targets <- bind_rows(targets, .id = 'feature_type')

# write targets table
cat('10: Wrtigin microRNA gene and protein targets.\n')

dbWriteTable(db,
             name = 'targets',
             targets,
             overwrite = TRUE)

### making index on cor_mir
dbSendQuery(db,
            statement = 'create index idx3 on targets (mirna_base);',
            overwrite = TRUE)

# write miRNASeq profiles
cat('11: Extracting microRNA profiles.\n')

df <- map(cohorts[-20], function(x) {
    # make names of RTCGA data.frames
    mi <- paste(x, 'miRNASeq', sep = '.')

    # read data.frames
    mi <- get(mi)

    # tidy data.frames
    mi <- mirna_tidy(mi)
    mi <- cbind(mirna_base = rownames(mi), as.data.frame(mi))
    mi <- gather(mi, bcr, count, -mirna_base)[, -2]
    mi <- as.data.table(mi)

    # print progress
    cat(paste('microRNAs profiles for', x, 'is done.\n'))

    # return tidy data.table
    return(mi)
})

### use reduce to merge data.tables
cat('12: Merging microRNA profiles.\n')

names(df) <- cohorts[-20]
df <- bind_rows(df, .id = 'study')

# write profiles table
cat('13: Writing microRNA profiles.\n')

dbWriteTable(db,
             name = 'profiles',
             df,
             overwrite = TRUE)

### making index on profiles
dbSendQuery(db,
            statement = 'create index idx4 on profiles (mirna_base);',
            overwrite = TRUE)

# disconnect from the db file
dbDisconnect(db)
tictoc::toc()
