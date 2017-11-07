library(shiny)
library(DBI)
library(RSQLite)
library(dplyr)
library(dbplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(DT)
library(DiagrammeR)
source('www/chooser.R')

# https://github.com/rstudio/shiny-examples/blob/master/036-custom-input-control/www/chooser-binding.js
# https://github.com/rstudio/shiny-examples/blob/master/036-custom-input-control/chooser.R

studies <- readLines('www/studies.txt')
miRBase <- readLines('www/miRBase.txt')

ui <- navbarPage(title = '',
                 tabPanel('Main',
                          sidebarPanel(tags$h1('miRCancerdb'), 
                                       tags$p('miRCancerdb is a free easy-to-use database of microRNA-gene/protein expression correlation in cancer.
                                              It was built mainly based on', tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'),
                                              'and', tags$a(href='http://www.targetscan.org', 'TargetScan'),'data.'),
                                       tags$hr(),
                                       radioButtons('expression',
                                                    label = 'Choose expression type:',
                                                    choices = c('mRNA', 'Protein'),
                                                    selected = 'mRNA'),
                                       textInput("study", 
                                                 label = "Enter TCGA study code:",
                                                 value = "ACC"),
                                       textInput("mirna",
                                                 label = "Enter the miRBase ID:",
                                                 value = "hsa-let-7b"),
                                       tags$hr(),
                                       radioButtons("type",
                                                    label = "Choose feature type:",
                                                    choices = c('Targets Only', 'ALL'),
                                                    selected = 'ALL'),
                                       numericInput("n_row",
                                                    label = "Choose number of features:",
                                                    value = 5),
                                       radioButtons("dir",
                                                    label = "Choose direction of correlation:", 
                                                    choices = c('Positive', 'Negative', 'Both'), 
                                                    selected = 'Both'),
                                       sliderInput("cor", 
                                                   label = "Choose absolute minimum correlation:",
                                                   min = -1,
                                                   max = 1,
                                                   value = 0,
                                                   step = .1),
                                       submitButton('Submit')),
                          mainPanel(tags$img(src = "logo.jpg",
                                             width = "80px",
                                             height = "80px",
                                             align = "right"),
                                    tags$br(),
                                    tags$br(),
                                    tags$br(),
                                    tags$h2("Getting started"),
                                    tags$p("To get started, first, you need to choose the type of features mRNA or protein, enter the", tags$a(href = 'https://cancergenome.nih.gov'), "study and the microRNA", tags$a(href= 'http://www.mirbase.org', 'miRBase'),"IDs."),
                                    tags$li(tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'), "identifiers like ACC, BLCA or BRCA."),
                                    tags$li(tags$a(href= 'http://www.mirbase.org', 'miRBase'), "microRNA identifiers according miRBase databse like hsa-let-7b."),
                                    tags$br(),
                                    tags$p("miRCancerdb web application provides different types of views:"),
                                    tags$li("Table view: a searchable table of three or more columns."),
                                    tags$li("Dot View: the same information in the tables in a dot graph"),
                                    tags$li("Bar View: works pretty much the same like the dot graph."),
                                    tags$li("Profile View: expression profiles of microRNAs."),
                                    tags$li("Summary View: summaries of distribution of non filtered data."),
                                    tags$p("For more information and other ways to access the database check How-to and about tabs."),
                                    tags$hr(),
                                    tabsetPanel(
                                      tabPanel('Table View',
                                               tags$br(),
                                               dataTableOutput("tbl"),
                                               downloadButton('dl_tbl',
                                                              label = 'Download')),
                                      tabPanel('Dot View',
                                               tags$br(),
                                               downloadButton('dl_dot',
                                                              label = 'Download'),
                                               plotOutput("dot")),
                                      tabPanel('Bar View',
                                               tags$hr(),
                                               plotOutput('bar'),
                                               downloadButton('dl_bar',
                                                              label = 'Download')),
                                      tabPanel('Profile View',
                                               tags$br(),
                                               plotOutput('profile'),
                                               downloadButton('dl_profile',
                                                              label = 'Download')),
                                      tabPanel('Summary View',
                                               tags$br(),
                                               tags$p("This summary is based on the unfiltered data of the selected", tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'),"studies and microRNAs.
                                                      This could help you to filter out the data to the shape and amount that suits you question."),
                                               tags$br(),
                                               tags$ol("Percentages of positive and negative correlations among all genes/proteins for each microRNA in each study."),
                                               plotOutput('summary_pie'),
                                               tags$ol("Numbers of feature targets for each microRNA."),
                                               tags$br(),
                                               tableOutput('summary_targets'),
                                               tags$br(),
                                               tags$ol("Distribution of correlations among all features for each microRNA in each study."),
                                               plotOutput('summary_dist'),
                                               tags$br(),
                                               tags$ol("Distribution of correlations among target features for each microRNA in each study."),
                                               plotOutput('summary_trgt_dist')),
                                      tabPanel('Selection View',
                                               tags$br(),
                                               tags$p("This view helps you select from the available", tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'),"studies and micorRNAs",
                                                      tags$a(href= 'http://www.mirbase.org', 'miRBase')," IDS.
                                                      Select IDs of interest, click submit and copy and past the formated IDs in the corresponding celles on the left panel."),
                                               tags$br(),
                                               tags$b("Select from the available TCGA studeis:"),
                                               tags$br(),
                                               tags$br(),
                                               chooserInput('select_study',
                                                            leftLabel = 'Available Studies',
                                                            rightLabel = 'Selected studies',
                                                            leftChoices = studies,
                                                            rightChoices = c(),
                                                            size = 10,
                                                            multiple = TRUE),
                                               tags$br(),
                                               verbatimTextOutput("selected_studies"),
                                               tags$br(),
                                               tags$b("Select from the available miRBase microRNA IDs:"),
                                               tags$br(),
                                               chooserInput('select_miRBase',
                                                            leftLabel = 'Available microRNAs',
                                                            rightLabel = 'Selected microRNAs',
                                                            leftChoices = miRBase,
                                                            rightChoices = c(),
                                                            size = 10,
                                                            multiple = TRUE),
                                               tags$br(),
                                               verbatimTextOutput("selected_miRBase")
                                      ),
                                      tabPanel('How-to',
                                               tags$br(),
                                               tags$b("How to use the miRCaner effeciently?"),
                                               tags$p("We recomend using miRCaner in following way. First, enter the", tags$a(href= 'http://www.mirbase.org', 'miRBase')," microRNA ID/s and the",
                                               tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'),"study ID/s of interest.
                                                      Then, check the Table and the Summary Views to make sure the amouont and the shape of the data suits answering your question.
                                                      Finally, adjust the filters to customize the data and the figures."),
                                               tags$br(),
                                               tags$b("How to enter multiple", tags$a(href= 'http://www.mirbase.org', 'miRBase'),"or ", tags$a(href = 'https://cancergenome.nih.gov', 'TCGA')," IDs?"),
                                               tags$p("You can enter multible IDs in a comma separated list. For example, three", tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'),
                                                      " studies could be used as following, ACC, BLCA, BRCA.
                                                      Three ",tags$a(href= 'http://www.mirbase.org', 'miRBase')," IDs could be used as following, hsa-let-7b, hsa-let-7c, hsa-mir-134."),
                                               tags$br(),
                                               tags$b("How to get data from miRCancerdb to use locally?"),
                                               tags$p("After choosing and filtering the data on the web app, click download button in the Table View.
                                                      The table will be downloaded in csv, comma separated, format that could be used in R or excel."),
                                               tags$br(),
                                               tags$b("How to obtain the full database of miRCancerdb?"),
                                               tags$p("There are a few ways to obtain the database or build it from scratch and use locally in your favorit data management package.
                                                      First, the scripts to build the databas are available on this github repository including the instruction.
                                                      Second, you can obtain the already built database by downloading this file.
                                                      Finally, the R package, Regulome, provides an interface to build or obtain the data and contains helper functions to work with the database locally.
                                                      In all cases, the database is a sqlite file that could be accessed by RSQLITE and dbplyr packages."))
                                    ))),
                 tabPanel('GitHub',
                          "This repository provides the scripts needed to build the database and the application locally. 
                          Comments, issues and contributions are welcomed.",
                          tags$a(href='https://github.com/MahShaaban/miRCancerdb', 'https://github.com/MahShaaban/miRCancerdb')),
                 tabPanel('About',
                          tags$h2('Overview'),
                          tags$p('miRCancerdb is a web interface for the database miRCancer.db which is an expression correlation of expression of microRNAs and genes/proteins.
                                 The web interface provids an easy access to the database without having to write any code.
                                 Among other features, users can search, choose, visualize and downlaod the particular subset of the data that they are interested in.'),
                          tags$br(),
                          tags$h2('Datasets'),
                          tags$p("The expression data was obtained from The Cancer Genome Atlas (", tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'),
                                 ") using", tags$a(href='https://bioconductor.org/packages/release/bioc/html/RTCGA.html', 'RTCGA'),
                                 "R package. And the targets data were obtained form", tags$a(href='http://www.targetscan.org', 'TargetScan')," using",
                                 tags$a(href = 'https://bioconductor.org/packages/release/data/annotation/html/targetscan.Hs.eg.db.html', 'targetscan.Hs.eg.db'),
                                 "Bioconductor package."),
                          tags$br(),
                          tags$h2('Workflow'),
                          tags$p(tags$a(href = 'https://cancergenome.nih.gov', 'TCGA'), "expression data for microRNAs, mRNA and proteins were
                                 obtained using ", tags$a(href='https://bioconductor.org/packages/release/bioc/html/RTCGA.html', 'RTCGA'),
                                 "R package. The expression correlation were calculated for each microRNA to each feature;
                                 gene/protein in each TCGA study.
                                 An SQLite database were build using these corelations. In addition, expression profiles of each miroRNA and
                                 microRNA targets data were obtained from", tags$a(href = 'https://bioconductor.org/packages/release/data/annotation/html/targetscan.Hs.eg.db.html', 'targetscan.Hs.eg.db'),
                                 tags$a(href='https://bioconductor.org', 'Bioconductor'),"package and added to the database."),
                          grVizOutput('workflow'),
                          tags$br(),
                          tags$h2('Other ways to use miRCancerdb'),
                          tags$ol("To use this app locally, download the", tags$a(href='https://github.com/MahShaaban/miRCancerdb', 'GitHub'),
                                  "repository to build the app on your machine.
                                  Alternatively, the R package", tags$a(href = 'https://github.com/MahShaaban/cRegulome', 'cRegulome'),
                                  "provides an interface to the same database and other functionality.")),
                 tabPanel('Contact us',
                          tags$p('Department of Biochemistry and Convergence Medical Sciences Institute of Health Sciences,'),
                          tags$p('Gyeonsange National University School of Medicine'),
                          tags$p('861 Beongil 15 jinju-daero, jinju, Gyeongnam 660-751,'),
                          tags$p('Republic of Korea'),
                          tags$p('Mob:+82-10-4045-1767'))
)

server <- function(input, output, session) {
  # table view
  ## get data for table view
  output_tbl <- reactive({
    # connect to db
    conn <- dbConnect(
      drv = SQLite(),
      'miRCancer.db'
    )
    
    # unpack input of query options
    tcga_study <- str_split(input$study, ', ', simplify = TRUE) %>% as.character
    mirna <- str_split(input$mirna, ', ', simplify = TRUE) %>% as.character
  
    # make data.frame with correlation data
    ## default table and query options
    tab <- 'cor_mir'
    if(input$expression == 'Protein') 
      tab <- 'cor_rppa'
    
    ## subset table
    dat <- conn %>%
      tbl(tab) %>%
      dplyr::select(mirna_base, feature, tcga_study) %>%
      filter(mirna_base %in% mirna) %>%
      collect %>%
      gather(study, cor, -mirna_base, -feature) %>%
      mutate(cor = as.numeric(cor)/100) %>%
      na.omit %>%
      filter(abs(cor) > input$cor) %>%
      group_by(study, mirna_base) %>%
      arrange(desc(abs(cor))) %>%
      slice(1:input$n_row)
    
    ## make data.frame with targets
    feature_type <- 'gene'
    if(input$expression == 'protein') feature_type <- 'protein'
    trgt <- conn %>%
      tbl('targets') %>%
      filter(mirna_base %in% mirna & feature_type == feature_type) %>%
      collect
    
    # return correlations with custom options
    if(input$type == 'Targets Only')
      dat <-inner_join(trgt, dat) %>% unique
    
    if(input$dir == 'Positive')
      dat <- dat %>% filter(cor > 0)
    
    if(input$dir == 'Negative')
      dat <- dat %>% filter(cor < 0)
    
    # disconnect db
    dbDisconnect(conn)
    
    # return data
    return(dat)
  })
  
  ## render table view
  output$tbl <- renderDataTable({
    output_tbl()
    })
  
  ## download table view
  output$dl_tbl <- downloadHandler(
    filename = 'mirna_cancer.csv',
    content = function(file) {
      write.csv(output_tbl(),
                file,
                row.names = FALSE)
    }
  )
  
  # dot view
  ## use output_tbl and make plot
  output_dot <- reactive({
    output_tbl() %>%
      mutate(color = ifelse(cor < 0, 'Negative', 'Positive'),
             size = abs(cor)) %>%
      ggplot(aes(x = '', y = feature, color = color, size = size)) +
      geom_point(alpha = .8) +
      facet_wrap(study~mirna_base, nrow = 1) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            legend.position = 'top',
            legend.direction = 'horizontal') +
      labs(x = '', y = '', color = 'Direction', size = 'Value')
  })
  
  ## render dot plot
  output$dot <- renderPlot({
    output_dot()
  }, height = function() {(length(unique(output_tbl()$feature)) * 25) + 100})
  
  ## download dot plot 
  output$dl_dot <- downloadHandler(
    filename = paste('dot_view', Sys.Date(), 'png', sep = '.'),
    content = function(file) {
      ggsave(file,
             plot = output_dot(),
             device = "png")
    }
  )
  
  # bar view
  ## use output_tbl and make plot
  output_bar <- reactive({
    output_tbl() %>%
      ggplot(aes(x = feature, y = cor)) +
      geom_col() +
      facet_wrap(mirna_base ~ study, scales = 'free_x') + 
      labs(x = '', y = "Pearson's correlation") +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  

  ## render bar plot
  output$bar <- renderPlot({
    output_bar()
  })
  
  ## download bar plot
  output$dl_bar <- downloadHandler(
    filename = paste('bar_view', Sys.Date(), 'png', sep = '.'),
    content = function(file) {
      ggsave(file,
             plot = output_bar(),
             device = "png")
    }
  )
  
  # profile view
  ## get data for profile view
  output_profile_dat <- reactive({
    # connect to db
    conn <- dbConnect(
      drv = SQLite(),
      './miRCancer.db'
    )
    
    # unpack text input
    tcga_study <- str_split(input$study, ', ', simplify = TRUE) %>% as.character
    mirna <- str_split(input$mirna, ', ', simplify = TRUE) %>% as.character
    
    # subset table data
    profiles <- conn %>%
      tbl('profiles') %>%
      filter(mirna_base %in% mirna & study %in% tcga_study) %>%
      collect
    
    # disconnect db
    dbDisconnect(conn)
    
    # return profiles
    return(profiles)
    })
  
  ## make profile plot
  output_profile <- reactive({
    output_profile_dat() %>%
      ggplot(aes(x = mirna_base, y =log(count))) +
      geom_jitter(width = .1, color = 'gray', alpha = .7) +
      facet_wrap(~study, ncol = 1) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = '', y = 'Expression (log count)')
  })
  
  ## render profile plot
  output$profile <- renderPlot({
    output_profile()
  })
  
  ## download profile plot
  output$dl_profile <- downloadHandler(
    filename = paste('profile_view', Sys.Date(), 'png', sep = '.'),
    content = function(file) {
      ggsave(file,
             plot = output_profile(),
             device = "png")
    }
  )
  
  # summary view
  ## get data for summary view
  output_summary <- reactive({
    # connect to databas
    conn <- dbConnect(
      drv = SQLite(),
      './miRCancer.db'
    )
    
    # unpack text input
    tcga_study <- str_split(input$study, ', ', simplify = TRUE) %>% as.character
    mirna <- str_split(input$mirna, ', ', simplify = TRUE) %>% as.character

    # make data.frame with correlation data
    ## default table and optional query options
    tab <- 'cor_mir'
    
    if(input$expression == 'Protein')
      tab <- 'cor_rppa'
    
    # subset tbl
    dat <- conn %>%
      tbl(tab) %>%
      dplyr::select(mirna_base, feature, tcga_study) %>%
      filter(mirna_base %in% mirna) %>%
      collect %>%
      gather(study, cor, -mirna_base, -feature) %>%
      mutate(cor = as.numeric(cor)/100) %>%
      na.omit %>%
      collect
    
    # disconnect db
    dbDisconnect(conn)
    
    # return data
    return(dat)
  })
  
  ## render summary pie
  output$summary_pie <- renderPlot({
    output_summary() %>%
      group_by(study, mirna_base) %>%
      summarise(Positive = sum(cor > 0) / n(),
                Negative = 1 - Positive) %>%
      gather(Direction, percent, Positive, Negative) %>%
      ggplot(aes(x = '', y = percent, fill = Direction)) +
      geom_bar(width = 1, stat = 'identity') +
      coord_polar('y', start = 0) +
      facet_wrap(mirna_base ~ study) +
      theme_light() +
      theme(legend.position = 'top', legend.direction = 'horizontal') +
      labs(x = '', y = '')
  })
  
  ## get data for summary targets
  output_targets <- reactive({
    # connect to databas
    conn <- dbConnect(
      drv = SQLite(),
      './miRCancer.db'
    )
    
    # unpack text input
    mirna <- str_split(input$mirna, ', ', simplify = TRUE) %>% as.character
    
    # default feature type and query options
    feature_type <- 'gene'
    
    if(input$expression == 'protein')
      feature_type <- 'protein'
    
    # subset targets table
    trgt <- conn %>%
      tbl('targets') %>%
      filter(mirna_base %in% mirna,
             feature_type == feature_type) %>%
      collect
    
    # disconnect db
    dbDisconnect(conn)
    
    # return data
    return(trgt)
  })
  
  ## render target summary table
  output$summary_targets <- renderTable({
    output_targets() %>%
      group_by(mirna_base) %>%
      summarise(Total = n())
  })
  
  ## render distribuition of all genes
  output$summary_dist <- renderPlot({
    output_summary() %>%
      ggplot(aes(x = cor, fill = mirna_base)) +
      geom_density(alpha = .5) +
      facet_grid( ~ study) +
      theme_light() +
      theme(legend.position = 'top', legend.direction = 'horizontal') +
      labs(x = "Pearson's Correlation", y = 'Density', fill = 'microRNA')
  })
  
  ## render distribution of target genes
  output$summary_trgt_dist <- renderPlot({
    inner_join(output_summary(),
               output_targets()) %>%
      ggplot(aes(x = cor, fill = mirna_base)) +
      geom_density(alpha = .5) +
      facet_grid(~ study) +
      theme_light() +
      theme(legend.position = 'top', legend.direction = 'horizontal') +
      labs(x = "Pearson's Correlation", y = 'Density', fill = 'microRNA')
  })
  output_selected_studies <- reactive({
    paste(input$select_study$right, collapse = ', ')
  })
  
  # selection view
  ## studies
  ### formate study text input
  output_selected_studies <- reactive({
    paste(input$select_study$right, collapse = ', ')
  })
  
  ### render text formated selected studies
  output$selected_studies <- renderText(
    paste('Copy and past on the right:',
          output_selected_studies())
  )
  
  ## miRBase
  ### formate miRBase text input
  output_selected_miRBase <- reactive({
    paste(input$select_miRBase$right, collapse = ', ')
  })

  ### render text formated selected miRBase
  output$selected_miRBase <- renderText(
    paste('Copy and past in the side panel:',
          output_selected_miRBase())
  )
  
  # about
  # render workflow
  output$workflow <- renderGrViz(
    grViz("
      digraph main {
          graph[splines=ortho]
          node [shape = box]
          n1 [label='TCGA']
          n2 [label='rnaseq \n & RPPA']
          n3 [label='miRNASeq']
          n4 [label='targets']
          n5 [label='Correlations']
          node [shape = none width = 0 height = 0 label = '']
          {rank=same; p1; p2; p3}
          {rank=same; n4; p4}
          {rank=same; n2; n3}
          
          node [shape=none fontcolor=red]
          v1 [label = 'Table View \n Dot View \n Bar View']
          v2 [label = 'Profile View' color='red']
          v3 [label = 'Summary View']
          v4 [label = 'Selection View']
          {rank=same; n5; v1}
          {rank=same; n2; v2}
          {rank=same; n4; v3}
          {rank=same; n1; v4}
          
          n1 -> p1 [arrowhead = none label='RTCGA']
          p2 -> n2
          p3 -> n3
          p2 -> p1 -> p3 [arrowhead = none]
          
          n2 -> p4 [arrowhead = none label = 'cor']
          n3 -> n4 [arrowhead = none label = 'targetscan.Hs.eg.db'] 
          p4 -> n4 [arrowhead = none]
          p4 -> n5
}")
  
    )
}

shinyApp(ui, server)