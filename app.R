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


ui <- fluidPage(theme = shinytheme("spacelab"),
   #tags$head(tags$style(".shiny-output-error{color: grey;}")), # change color of all error messages to grey
   navbarPage(h4(strong("NOJAH : NOt Just Another Heatmap")), 
              
  
    tabPanel("HomePage", 
             fluidRow(
               column(8, offset = 2, 
                      p("This interactive web application: NOt Just Another Heatmap (NOJAH) is developed in R with Shiny to"),
                      br(),
                      em("1) Perform genome wide analysis of cancer genomic data sets"),
                      br(),
                      em("2) Provide visualization, analysis and download of MMRF CoMMpass Expression, Variant and CNV data with Cluster of Cluster Analysis"),
                      br(),
                      em("3) Perform significance of cluster analysis using a robust bootstrap approach."), 
                      br(),
                      br(),
                      "The goal of this tool is to provide a one stop shop to analyze Genome Wide data or CoMMpass data or perform genomic analysis on their own data.",
                      htmlOutput("ReadMe"), tableOutput("Eg"), htmlOutput("Caption1"), tableOutput("Eg2"), htmlOutput("Caption2")
                      ))),    
    tabPanel("Genome Wide Analysis",
             fluidRow(
               column(2,
                wellPanel(h3("Input GW Data"),
                 selectInput("gw_file1",label= "Select an example dataset or upload your own with 'Load my own GW data.'", 
                             choices = c("Example coMMpass IA9 Expression data"= "GW_Example1", #"Example BRCA Methylation 27 K Data"="GW_Example2", 
                                         "Load my own data" = "load_my_own_gw")),
                 conditionalPanel("input.gw_file1 == 'load_my_own_gw'",
                                  fileInput('gw_file2', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                 conditionalPanel("input.gw_file1 == 'GW_Example1'", downloadButton('download_GW_Ex1', 'Download GW CoMMPass Expression DS'))#,
                # conditionalPanel("input.gw_file1 == 'GW_Example2'", downloadButton('download_GW_Ex2', 'Download GW BRCA Example DataSet'))
                 
                 ),
                br(), br(), 
                wellPanel(h3("Data Subsetting"),
                 checkboxGroupInput("gw_subset","Subset GW data by:",choices = list("Variance"= "VAR","Median Absoute Deviation" = "MAD", "Inter Quartile Range" = "IQR"), selected = c("VAR", "MAD", "IQR")),
                 conditionalPanel(condition = "$.inArray('VAR', input.gw_subset) > -1 & $.inArray('MAD', input.gw_subset) > -1 & $.inArray('IQR', input.gw_subset) > -1",
                                  radioButtons("IMVA_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")), ###
                                  conditionalPanel(condition = "input.IMVA_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("IMVA_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=45)),
                                  conditionalPanel(condition = "input.IMVA_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("IMVA_pInput", label = "Percentile value", min = 0, max = 100, value = 45, step = 5))),
                 conditionalPanel(condition = "$.inArray('VAR', input.gw_subset) > -1 & $.inArray('MAD', input.gw_subset) > -1 & $.inArray('IQR', input.gw_subset) <= -1",
                                  radioButtons("var_mad_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")), ###
                                  conditionalPanel(condition = "input.var_mad_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("var_mad_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=75)),
                                  conditionalPanel(condition = "input.var_mad_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("var_mad_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                 ),
                 conditionalPanel(condition = "$.inArray('VAR', input.gw_subset) > -1 & $.inArray('IQR', input.gw_subset) > -1 & $.inArray('MAD', input.gw_subset) <= -1 ",
                                  radioButtons("var_iqr_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")), ###
                                  conditionalPanel(condition = "input.var_iqr_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("var_iqr_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=75)),
                                  conditionalPanel(condition = "input.var_iqr_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("var_iqr_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                 ),
                 
                 conditionalPanel(condition = "$.inArray('MAD', input.gw_subset) > -1 & $.inArray('IQR', input.gw_subset) > -1 & $.inArray('VAR', input.gw_subset) <= -1",
                                  radioButtons("mad_iqr_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")), ###
                                  conditionalPanel(condition = "input.mad_iqr_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("mad_iqr_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=75)),
                                  conditionalPanel(condition = "input.mad_iqr_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("mad_iqr_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                 ),
                 conditionalPanel(condition = "input.gw_subset=='VAR'",
                        radioButtons("var_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")),
                        conditionalPanel(condition = "input.var_PercenChoice == 'Percentile Slider'",
                                  sliderInput("var_pslider", "Percentile Value:", 
                                              min=0, max=100, value=75)),
                        conditionalPanel(condition = "input.var_PercenChoice == 'Manually Enter Percentile'",
                                  numericInput("var_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                   ),
                  conditionalPanel(condition = "input.gw_subset=='MAD'",
                                  radioButtons("mad_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")),
                                  conditionalPanel(condition = "input.mad_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("mad_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=75)),
                                  conditionalPanel(condition = "input.mad_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("mad_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                   ),
                  conditionalPanel(condition = "input.gw_subset=='IQR'",
                                  radioButtons("iqr_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")),
                                  conditionalPanel(condition = "input.iqr_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("iqr_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=75)),
                                  conditionalPanel(condition = "input.iqr_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("iqr_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                   )
                  
                 
                 ),
                br(), br(), br(), 
                wellPanel(h3("Downloads"),
                          br(), 
                  textInput("fname_subset", "Type the file name you would like to save subset data as :", value = "Subset_data"),
                  downloadButton('downloadSubset', 'Download Subset data'),
                   br(), br(),
                #wellPanel(
                  textInput("fname_HM", "Save HeatMap as :", value = "GW_subset_HM"),
                  downloadButton('downloadHM', 'Download HM'),
                  conditionalPanel("input.conditionedPanels == 2",   
                                 h5(downloadLink('downloadCuttree', 'Download Column clusters after cut-tree'))),
                  conditionalPanel("input.conditionedPanels == 3", 
                                 h5(downloadLink('downloadCuttree2', 'Download Row clusters after cut-tree'))))
                ),               
                column(8, 
                       tags$head(tags$style(type="text/css", "
             #loadmessage {
                                           position: fixed;
                                           bottom: 0px;
                                           right: 0px;
                                           width: 100%;
                                           padding: 5px 0px 5px 0px;
                                           text-align: center;
                                           font-weight: bold;
                                           font-size: 100%;
                                           color: #000000;
                                           background-color: #b8b8b8;
                                           z-index: 105;
                                           }
                                           ")),
                       
                     h3(strong("Measures of Spread")),
                     plotlyOutput("Boxplot"), 
                     h3(strong("Most variable Gene Selection")),
                     h5(em("To see the position of your 'gene of interest', use the 'Choose Option' drop down to the right")),
                     plotlyOutput("GW_Scatter_LH"),
                     htmlOutput("n_selected"),
                     br(),
                     br(),
                     h3(strong("Visualization of Selected Subset")),
                     tabsetPanel(type = "tabs", 
                       tabPanel("HeatMap", plotOutput("GW_subset_heatmap", height = 1200), 
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=1), 
                       tabPanel("Column Dendrogram",plotOutput("plot1", height = 800), htmlOutput("display"), br(), DT::dataTableOutput("df"), value=2), 
                       tabPanel("Row Dendrogram",plotOutput("plot2", height = 800), htmlOutput("display2"), br(), DT::dataTableOutput("df2"),  value =3), 
                       id = "conditionedPanels")
                   
                     ),
               column(2, 
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(),
                    uiOutput(outputId="geneSelector", width=NULL),
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                     
                    conditionalPanel("input.conditionedPanels == 1", 
                       wellPanel(  
                        ########## HeatMap Clustering options ##########
                        h4("Heat Map Options"),
                        selectInput("norm", "Normalization Type",
                             c("Z-Score", "Modified Z-Score", "none")),
                        selectInput("norm2", "Normalize by:",
                             c("row", "col", "both")),
                        sliderInput("inSlider", "Scale Range",
                             min = -10, max = 20, value = c(-2, 2)),
                        sliderInput("inSlider2", "Plot Margin dimensions",
                             min = 0, max = 20, value = c(14, 16)) ),
                        wellPanel(  
                          h4("Clustering Measures"),
                          selectInput("dist", "Distance Method",
                             c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
                          selectInput("hclust", "Agglomerative Linkage Method", 
                             c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid")),
                          conditionalPanel("input.conditionedPanels==1 ", 
                                  radioButtons("clust_byrow", "Row dendrogram", inline = TRUE, c("TRUE", "FALSE")) ,
                          #conditionalPanel("input.conditionedPanels==1", 
                                  radioButtons("clust_bycol", "Col dendrogram", inline = TRUE, c("TRUE", "FALSE")) ,
                          #conditionalPanel("input.conditionedPanels== 1", 
                                  radioButtons("dispRow", "Display Row labels?:", inline = TRUE,c("No", "Yes")),
                                  conditionalPanel("input.dispRow == 'Yes'", sliderInput("size1", "If yes, Row Label font size", min = 0.01, max = 3, value = 0.5)),
                                  radioButtons("dispCol", "Display Col labels?:", inline = TRUE, c("No", "Yes")),
                                  conditionalPanel("input.dispCol == 'Yes'", sliderInput("size2", "If yes, Col Label font size", min = 0.01, max = 3, value = 0.5) ))
                          ), 
                          wellPanel(
                           h4("Heat Map colors"),
                            colourInput("low", "low", "green", returnName = TRUE, palette = "limited", showColour = "background"),
                            colourInput("mid", "mid", "black", returnName = TRUE, palette = "limited", showColour = "background"),
                            colourInput("high", "high", "red", returnName = TRUE, palette = "limited", showColour = "background") )
                         ),
                         conditionalPanel("input.conditionedPanels ==2 | input.conditionedPanels == 3",
                              wellPanel(
                                conditionalPanel("input.conditionedPanels==2", 
                                                 sliderInput("sizeClable", "Adjust Column Label font size", min = 0.01, max = 3, value = 0.8) ),
                                conditionalPanel("input.conditionedPanels==3", 
                                                 sliderInput("sizeRlable", "Adjust Row Label font size", min = 0.01, max = 3, value = 0.42) ),
                                
                                # Column Dendrogram tab
                                conditionalPanel("input.conditionedPanels==2", 
                                                 radioButtons("cutcolden", "Cut Col dendrogram?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE)) ),
                                conditionalPanel("input.conditionedPanels==2 & input.cutcolden == 'TRUE'", 
                                                 numericInput("cuttree", "Cut Col Dendrogram at:", 2)),
                                
                                # Row Dendrogram tab
                                conditionalPanel("input.conditionedPanels==3", 
                                                 radioButtons("cutrowden", "Cut Row dendrogram?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE)) ),
                                conditionalPanel("input.conditionedPanels==3 & input.cutrowden == 'TRUE'", 
                                                 numericInput("cuttree2", "Cut Row Dendrogram at:", 2))
                                
                                )
                            )
                          )
              
             
    )),
    tabPanel("Cluster of Cluster analysis", 
             fluidRow(
               column(2,
                      wellPanel(
                        h3("Input file"),
                        conditionalPanel("input.cPanels1 == 1",
                                         selectInput("Exp_file",label= "Select the example CoMMpass Expression ds or upload your own with 'Load my own'", 
                                                     choices = c("Example ds File"="Exp_example", "Load my own data" = "Exp_load_my_own")),
                                         conditionalPanel("input.Exp_file == 'Exp_load_my_own'",
                                                          fileInput('Exp_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.Exp_file == 'Exp_example'",
                                                          downloadButton('Exp_download', 'Download Expression data'))),
                        conditionalPanel("input.cPanels1 == 2",
                                         selectInput("Variant_file",label= "Select the example CoMMpass Variant ds or upload your own with 'Load my own'", 
                                                     choices = c("Example ds File"="Variant_example", "Load my own data" = "Variant_load_my_own")),
                                         conditionalPanel("input.Variant_file == 'Variant_load_my_own'",
                                                          fileInput('Variant_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.Variant_file == 'Variant_example'",
                                                          downloadButton('Variant_download', 'Download Variant data'))),
                        conditionalPanel("input.cPanels1 == 3",
                                         selectInput("CNV_file",label= "Select the example CoMMpass CNV ds or upload your own with 'Load my own'", 
                                                     choices = c("Example ds File"="CNV_example", "Load my own data" = "CNV_load_my_own")),
                                         conditionalPanel("input.CNV_file == 'CNV_load_my_own'",
                                                          fileInput('CNV_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.CNV_file == 'CNV_example'",
                                                          downloadButton('CNV_download', 'Download CNV data'))),
                        conditionalPanel("input.cPanels1 == 4",
                                         selectInput("coca_file",label= "Select the example CoC analysis ds or upload your own with 'Load my own'", 
                                                     choices = c("Example ds File"="coca_example", "Load my own data" = "coca_load_my_own")),
                                         conditionalPanel("input.coca_file == 'coca_load_my_own'",
                                                          fileInput('coca_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.coca_file == 'coca_example'",
                                                          downloadButton('coca_download', 'Download CoC analysis data')))
                        ),
                       wellPanel(h3("Clustering Options"),
                          conditionalPanel("input.cPanels1 == 1",      
                            selectInput("Exp_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson"),
                            selectInput("Exp_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                            radioButtons("Exp_pItems", "No. of iterations", choices = c(100, 500, 1000), selected = 100)),
                          conditionalPanel("input.cPanels1 == 2",      
                            selectInput("Variant_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "maximum"),
                            selectInput("Variant_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                            radioButtons("Variant_pItems", "No. of iterations", choices = c(100, 500, 1000), selected = 100)),
                          conditionalPanel("input.cPanels1 == 3",      
                            selectInput("CNV_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "manhattan"),
                            selectInput("CNV_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                            radioButtons("CNV_pItems", "No. of iterations", choices = c(100, 500, 1000), selected = 100)),
                          conditionalPanel("input.cPanels1 == 4",
                            checkboxGroupInput("coca_platform","Perform CoC analysis using", choices = list("Expression"= "EXP","Variant" = "PROP", "Copy Number" = "CNV"), selected = c("EXP", "PROP", "CNV")),
                            h6("Select atleast two platforms to run CoC analysis"),
                            selectInput("coca_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean"),
                            selectInput("coca_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                            radioButtons("coca_pItems", "No. of iterations", choices = c(100, 500, 1000), selected = 100))
                          
                      ),
                      conditionalPanel("input.cPanels1 == 4",
                            wellPanel(h4("Clinical markers"),
                               selectInput("clinical",label= "Select a clinical feature from below or add your own using 'Load your own'.", 
                                 choices = c("Features available"="available", "Load my own clinical data" = "load_my_own_c")),
                               conditionalPanel("input.clinical == 'available'",
                                 downloadButton('download_clinical', 'Download clinical ds Example'),
                                     br(), br(),
                               checkboxGroupInput("cli_feature","Choose Clinical feature to be added:",choices = list("High vs Not High Risk"= "HR"), selected = "HR")),
                               conditionalPanel("input.clinical == 'load_my_own_c'",
                                    fileInput('clinical_file', 'Choose file to upload (maximum size 500 MB).', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')))
                      )),
                     wellPanel(h3("Download Results"),
                          h4("Consensus Clustering"),
                          conditionalPanel("input.cPanels1 == 1", 
                                  downloadButton("Exp_cc_dl", label = "Download Expression clusters"),
                                  h6("Download may take a while. Once complete, the result pdf file will automatically open.")),
                          conditionalPanel("input.cPanels1 == 2", 
                                  downloadButton("Variant_cc_dl", label = "Download Variant clusters"), 
                                  h6("Download may take a while. Once complete, the result pdf file will automatically open.")),
                          conditionalPanel("input.cPanels1 == 3", 
                                  downloadButton("CNV_cc_dl", label = "Download CNV clusters"),
                                  h6("Download may take a while. Once complete, the result pdf file will automatically open.")),
                          conditionalPanel("input.cPanels1 == 4", 
                                 downloadButton("coca_cc_dl", label = "Download CoC analysis"),
                                 h6("Download may take a while. Once complete, the result pdf file will automatically open."),
                                 h4("Sample Clusters"),
                                 downloadButton("coca_sil", label = "Download CoC Sample clusters"),
                                 br(),br(),
                                 h4("CoC HeatMap"),
                                 downloadButton('download_coca_HM', 'Download CoC Analysis HM'),
                                 br(),br(),
                                 h4("Cluster Interpretation"),
                                 downloadButton('dl_coca_inter', 'Download Cluster Interpretation')
                         )
                      
                          #br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                          #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      )
                     #)
                   ),
               column(8,
                      tags$head(tags$style(type="text/css", "
             #loadmessage {
                                           position: fixed;
                                           bottom: 0px;
                                           right: 0px;
                                           width: 100%;
                                           padding: 5px 0px 5px 0px;
                                           text-align: center;
                                           font-weight: bold;
                                           font-size: 100%;
                                           color: #000000;
                                           background-color: #b8b8b8;
                                           z-index: 105;
                                           }
                                           ")),
                    
                      tabsetPanel(type = "tabs", 
                                  tabPanel("Expresssion", htmlOutput("com_text1"), plotOutput("Exp_cc", height = 500), br(), htmlOutput("com_text2"), plotOutput("Exp_sil", height = 800), 
                                      conditionalPanel(condition="$('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=1), 
                                  tabPanel("Variant", htmlOutput("com_text12"), plotOutput("Variant_cc", height = 500), br(), htmlOutput("com_text22"), plotOutput("Variant_sil", height = 800), 
                                           conditionalPanel(condition="$('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=2), 
                                  tabPanel("Copy Number", htmlOutput("com_text13"),plotOutput("CNV_cc", height = 500), br(), htmlOutput("com_text23"), plotOutput("CNV_sil", height = 800),  
                                           conditionalPanel(condition="$('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value =3),
                                  tabPanel("Cluster of Clusters (CoC) Analysis", htmlOutput("com_text14"), plotOutput("coca_cc", height = 500), br(),htmlOutput("com_text3"), plotOutput("coca_heatmap", height = 1200), htmlOutput("com_text31"), plotOutput("varmean_bxplot"),
                                           conditionalPanel(condition="$('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value = 4),
                                  id = "cPanels1")
                      
                ),
               column(2,
                      wellPanel(h3("Choose Optimal Number of clusters"),
                                conditionalPanel("input.cPanels1 == 1", numericInput("Exp_opt_k", "Optimal k is ",value = 5)),
                                conditionalPanel("input.cPanels1 == 2", numericInput("Variant_opt_k", "Optimal k is ",value = 2)),
                                conditionalPanel("input.cPanels1 == 3", numericInput("CNV_opt_k", "Optimal k is ",value = 2)),
                                conditionalPanel("input.cPanels1 == 4", numericInput("coca_opt_k", "Optimal k is ",value = 2))
                      ),
                      conditionalPanel("input.cPanels1 == 4", 
                                       br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                      
                                       wellPanel(  
                                         h4("Clustering Measures"),
                                         selectInput("dist_4", "Distance Method",
                                                     c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean"),
                                         selectInput("hclust_4", "Agglomerative Linkage Method", 
                                                     c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                                         radioButtons("clust_byrow_4", "Row dendrogram", inline = TRUE, c("TRUE", "FALSE")) ,
                                         radioButtons("clust_bycol_4", "Col dendrogram", inline = TRUE, c("TRUE", "FALSE")) ,
                                         radioButtons("dispRow_4", "Display Row labels?:", inline = TRUE,c("No", "Yes"), selected = "Yes"),
                                         conditionalPanel("input.dispRow_4 == 'Yes'", sliderInput("size1_4", "If yes, Row Label font size", min = 0.01, max = 3, value = 1.5)),
                                         radioButtons("dispCol_4", "Display Col labels?:", inline = TRUE, c("No", "Yes")),
                                         conditionalPanel("input.dispCol_4 == 'Yes'", sliderInput("size2_4", "If yes, Col Label font size", min = 0.01, max = 3, value = 0.5) )
                                       ),
                                       wellPanel(  
                                         ########## HeatMap Clustering options ##########
                                         h4("Heat Map Options"),
                                         # selectInput("norm_4", "Normalization Type",
                                         #             c("Z-Score", "Modified Z-Score", "none")),
                                         # selectInput("norm2_4", "Normalize by:",
                                         #             c("row", "col", "both")),
                                         # sliderInput("inSlider_4", "Scale Range",
                                         #             min = -10, max = 20, value = c(-2, 2)),
                                         sliderInput("inSlider2_4", "Plot Margin dimensions",
                                                     min = 0, max = 20, value = c(13, 15)) )#,
                                       
                                       #wellPanel(
                                      #   h4("Heat Map colors"), 
                                      #   colourInput("low_4", "low", "cyan", returnName = TRUE, palette = "limited", showColour = "background"),
                                      #   colourInput("mid_4", "mid", "khaki1", returnName = TRUE, palette = "limited", showColour = "background"),
                                      #   colourInput("high_4", "high", "chartreuse", returnName = TRUE, palette = "limited", showColour = "background") )
                      )))                    
                      
             ),
    tabPanel("Significance of Clusters Analysis",
             #tags$iframe(src= "http://shinygispa.winship.emory.edu/CASH", width = 1900, height = 1000)
             fluidRow(
               column(2,
                      wellPanel(h3("Input Data to test significance of clusters "),
                                selectInput("Sig_file1",label= "Select an example dataset or upload your own with 'Load my own GW data.'", 
                                            choices = c("Example Exp Data"="Sig_Example1","Example Meth Data"="Sig_Example2", "Load my own data" = "load_my_own_Sig")),
                                conditionalPanel("input.Sig_file1 == 'load_my_own_Sig'",
                                                 fileInput('Sig_file2', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                                conditionalPanel("input.Sig_file1 == 'Sig_Example1'", downloadButton('download_Sig_Ex1', 'Download Example DataSet')),
				 conditionalPanel("input.Sig_file1 == 'Sig_Example2'", downloadButton('download_Sig_Ex2', 'Download Example DataSet'))

                      ),
                      wellPanel(h4("Select Data Rows and column start"),
                                numericInput("DataR", label = "Data row", min = 1, max = 10, value = 4),
                                numericInput("DataC", label = "Data column", min = 1, max = 10, value = 3)
                                ),
                      conditionalPanel("input.cPanels2 == 2 | input.cPanels2 == 3 |input.cPanels2 == 4", 
                                       wellPanel(
                                         textInput("Sig_fname", "Type the file name you would like to save as", value = "HeatMap"),
                                         downloadButton('Sig_downloadPlots', 'Download HeatMap and dendrograms'),
                                         br(),
                                         conditionalPanel("input.cPanels2 == 3",   
                                                          h5(downloadLink('Sig_downloadCuttree', 'Download Column clusters after cut-tree'))),
                                         conditionalPanel("input.cPanels2 == 4", 
                                                          h5(downloadLink('Sig_downloadCuttree2', 'Download Row clusters after cut-tree')))
                                       )
                      )),
               column(8,
                      tabsetPanel(type = "tabs", 
                                  #tabPanel("ReadMe", htmlOutput("ReadMe"), tableOutput("Eg"), htmlOutput("Caption1"), tableOutput("Eg2"), htmlOutput("Caption2"), htmlOutput("blurp"), value = 1),
                                  tabPanel("HeatMap", plotOutput("Sig_plot", width = 1300, height = 1300 ), value=2), 
                                  tabPanel("Column Dendrogram", plotOutput("Sig_plot1", height= 600, width = 1500), htmlOutput("Sig_display"), br(), DT::dataTableOutput("Sig_df"), htmlOutput("Sig_pv"), htmlOutput("Sig_pvalue"),  value=3), 
                                  tabPanel("Row Dendrogram", plotOutput("Sig_plot2", height = 800, width = 1500), htmlOutput("Sig_display2"), br(), DT::dataTableOutput("Sig_df2"), htmlOutput("Sig_pv2"), htmlOutput("Sig_pvalue2"),  value =4),
                                  id = "cPanels2")
                      ),
               column(2,
                      conditionalPanel("input.cPanels2==2",
                                       wellPanel(  
                                         ########## HeatMap Clustering options ##########
                                         h4("Heat Map Options"),
                                         selectInput("Sig_norm", "Normalization Type",
                                                     c("Z-Score", "Modified Z-Score", "none")),
                                         selectInput("Sig_norm2", "Normalize by:",
                                                     c("row", "col", "both"), selected = "both"),
                                         sliderInput("Sig_inSlider", "Scale Range",
                                                     min = -10, max = 20, value = c(-2, 2)),
                                         conditionalPanel("input.cPanels2==1 | input.cPanels2==2", 
                                                          sliderInput("Sig_inSlider2", "Plot Margin dimensions",
                                                                      min = 0, max = 20, value = c(14, 16)) )
                                       ),
                                       wellPanel(  
                                         h4("Clustering Measures"),
                                         selectInput("Sig_dist", "Distance Method",
                                                     c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean"),
                                         selectInput("Sig_hclust", "Agglomerative Linkage Method", 
                                                     c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "complete"),
                                         conditionalPanel("input.cPanels2==1 | input.cPanels2==2", 
                                                          radioButtons("Sig_clust_byrow", "Row dendrogram", inline = TRUE, c("TRUE", "FALSE")) ),
                                         conditionalPanel("input.cPanels2==1 | input.cPanels2==2", 
                                                          radioButtons("Sig_clust_bycol", "Col dendrogram", inline = TRUE, c("TRUE", "FALSE")) ),
                                         conditionalPanel("input.cPanels2==2", 
                                                          radioButtons("Sig_dispRow", "Display Row labels?:", inline = TRUE,c("No", "Yes")),
                                                          conditionalPanel("input.Sig_dispRow == 'Yes'", sliderInput("Sig_size1", "If yes, Row Label font size", min = 0.01, max = 3, value = 0.5)),
                                                          radioButtons("Sig_dispCol", "Display Col labels?:", inline = TRUE, c("No", "Yes")),
                                                          conditionalPanel("input.Sig_dispCol == 'Yes'", sliderInput("Sig_size2", "If yes, Col Label font size", min = 0.01, max = 3, value = 0.5) ))
                                       ), 
                                       
                                       wellPanel(
                                         h4("Heat Map colors"),
                                         colourInput("Sig_low", "low", "green", returnName = TRUE, palette = "limited", showColour = "background"),
                                         colourInput("Sig_mid", "mid", "black", returnName = TRUE, palette = "limited", showColour = "background"),
                                         colourInput("Sig_high", "high", "red", returnName = TRUE, palette = "limited", showColour = "background") )
                      ) ,
                      
                      conditionalPanel("input.cPanels2==3 | input.cPanels2==4",
                                       wellPanel(
                                         conditionalPanel("input.cPanels2==3", 
                                                          sliderInput("Sig_sizeClable", "Adjust Column Label font size", min = 0.01, max = 3, value = 0.8) ),
                                         conditionalPanel("input.cPanels2==4", 
                                                          sliderInput("Sig_sizeRlable", "Adjust Row Label font size", min = 0.01, max = 3, value = 0.42) ),
                                         
                                         # Column Dendrogram tab
                                         conditionalPanel("input.cPanels2==3", 
                                                          radioButtons("Sig_cutcolden", "Cut Col dendrogram?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE)) ),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE'", 
                                                          numericInput("Sig_cuttree", "Cut Col Dendrogram at:", 2)),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE'",
                                                          radioButtons("Sig_pvalue_cal", "Assess Gene set significance in separation of specimens into 2 clusters?:", c("No" = FALSE, "Yes" = TRUE), inline = TRUE)),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'" , 
                                                          selectInput("Sig_file3", label= "Select a dataset or upload your own with 'Load my own data.'", choices = c("Exp Sampling Data" ="Sig_Exp.Example","Meth Sampling Data" ="Sig_Meth.Example", "Load my own sampling data" = "Sig_load_my_own_s_data"))),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_file3 == 'Sig_load_my_own_s_data'",
                                                          fileInput('Sig_file4', 'Choose file to upload to sample from to estimate significance of separation (Maximum size 100 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))) ,
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                                          numericInput("Sig_n", "Sample size for bootstrap:", 1000)),
                                         #conditionalPanel("input.cPanel2s==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                        #                  numericInput("Sig_n_iter", "No. of iterations for bootstrap:", 1000)),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                                          numericInput("Sig_n_iter", "No. of iterations for bootstrap:", 1000)),
                                         
                                         #conditionalPanel("input.conditionedPanels==3 & input.cutcolden == 'TRUE' & input.pvalue_cal == 'TRUE'", 
                                         #                 radioButtons("histplot1", "Display distribution of permuted p-values in comparison to obs p-value ?:", c( "Yes" = TRUE, "No" = FALSE), inline = TRUE)),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                                          actionButton("Sig_goButton", "Go!")),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                                          p("Click the button to start sampling using bootstrap method for estimating the p-value. A progress indicator will appear shortly (~approx 10 seconds), on top of page indicating the status. Once complete, the p-value will be displayed in the main panel.")),
                                         
                                         # Row Dendrogram tab
                                         conditionalPanel("input.cPanels2==4", 
                                                          radioButtons("Sig_cutrowden", "Cut Row dendrogram?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE)) ),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE'", 
                                                          numericInput("Sig_cuttree2", "Cut Row Dendrogram at:", 2)),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE'",
                                                          radioButtons("Sig_pvalue_cal2", "Assess significance of samples in separation of gene set into 2 clusters (only for more than 2 cluster groups)?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE))),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE' " , 
                                                          selectInput("Sig_file5", label= "Select a dataset or upload your own with 'Load my own data.'", choices = c("Example Exp Sampling Data" ="Sig_Exp.Example2", "Example Meth Sampling Data" ="Sig_Meth.Example2", "Load my own sampling data" = "load_my_own_s_data2"))),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_file5 == 'load_my_own_s_data2'",
                                                          fileInput('Sig_file6', 'Choose file to upload to sample from to estimate significance of separation (Maximum size 75 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))) ,
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE'", 
                                                          numericInput("Sig_n2", "Sample size for bootstrap:", 1000)),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE'", 
                                                          numericInput("Sig_n_iter2", "No. of iterations for bootstrap:", 1000)),
                                         #conditionalPanel("input.conditionedPanels==4 & input.cutrowden == 'TRUE' & input.pvalue_cal2 == 'TRUE'", 
                                         #                 radioButtons("histplot2", "Display distribution of permuted p-values in comparison to obs p-value ?:", c( "Yes" = TRUE, "No" = FALSE), inline = TRUE)),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE'", 
                                                          actionButton("Sig_goButton2", "Go!")),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE'", 
                                                          p("Click the button to start sampling using bootstrap method for estimating the p-value. A progress indicator will appear shortly (~approx 10 seconds), on top of page indicating the status. Once complete, the p-value will be displayed in the main panel."))
                                         
                                       )
                        )
                      ) 
                      )
             ),
    tabPanel("Tutorial",
             tags$iframe(src= "NOJAH_tutorial.pdf", width = 1800, height = 1000)),
    navbarMenu("About Us",
               tabPanel("How to Cite",
               fluidRow(
                 column(8, offset = 2, 
                        "The National Cancer Institute (NCI) requires that publications acknowledge the Winship Cancer Institute CCSG support, and they are tracking compliance. When using this tool to report results in your publication, please include the following statement in the acknowledgment section of your publication(s):",
                        br(),
                        br(),
                        em("Research reported in this publication was supported in part by the Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University and NIH/NCI under award number P30CA138292. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health."),
                        br(),
                        br(),
                        "Authors- Manali Rupji, dual M.S., Bhakti Dwivedi Ph.D. & Jeanne Kowalski Ph.D.",
                        br(),
                        "Maintainer- Manali Rupji 'manali(dot)rupji(at)emory(dot)edu'"
                 ))),
               tabPanel("Contact Us",
               fluidRow(
                 column(8, offset = 2, 
                        p("The Biostatistics and Bioinformatics Shared Resource at Winship Cancer Institute of Emory University"),
                          a(href="https://bbisr.winship.emory.edu/", "https://bbisr.winship.emory.edu/"))
                 )),
               tabPanel("Feedback",
               fluidRow(
                 column(8, offset = 2, 
                       "This App is developed and maintained by Manali Rupji at the Biostatistics and Bioinformatics core, Winship Cancer Institute, Emory University.",
                        br(),
                        br(),
                        "As a Biostatistics and Bioinformatics core, we are actively improving and expanding our NGS analysis services and analysis products. For any questions, comments, or suggestions, please email the developer at mrupji@emory.edu."
                 )))
              )
    ))

options(shiny.maxRequestSize=1000*1024^2)
options(shiny.sanitize.errors = TRUE) 
#options(bitmapType='cairo')

source("zClust.R")
source("modzClust.R")
source("consensus_clustering.R")
source("silhouette_plot.R")
source("coca.R")
source("helper.R")
source("colbars.R")
source("plotMeans.R")
source("bxplot.R")
source("bootstrapfun.R")
source("contingencyfun.R")
source("pvaluefunc.R")

server <- function(input, output) {
  
  output$ReadMe <- renderUI({
    str0 <- paste("&emsp;")
    str00 <- paste("DATA INPUT")
    str1 <- paste("Data should be input as a .txt or .csv file. The first two rows of the data file have information about the patients/specimens and their response/subtype; all remaining rows have gene expression data, one row per gene.  In the case of Microarray gene expression data in which there are several probes corresponding to a single gene, a unique identifier would need to be created to separately identify each probe such as, 'Gene 1_p1', 'Gene1_p2' indicating Gene 1 has two probes.  The columns represent the different experimental samples. A maximum of up to 10 different sample groups and 6 different gene groups may be used with this tool.")
    str2 <- paste("Data Format")
    str3 <- paste("&emsp; 	1.	The first line of the file contains the gene identifier 'gene_id' (column 1), gene group classification 'Groups' (column 2) followed by the patient IDs e.g. TCGA.01.1A2B, one per column, starting at column 3. Column 1 gene identifier has to be labelled 'gene_id' and column 2 header should be labelled 'Groups' for using this tool. Other titles may cause the program to display errors. ")
    str4 <- paste("&emsp; 	2.	The second line of the file contains the patient response classification e.g. Fav/Unf for favorable outcome group vs the unfavorable outcome group or Normal/Tumor, etc., in alphabetical order, starting at column 3. The first two columns for this row should be blank.")
    str5 <- paste("&emsp; 	3.	The remaining lines contain gene expression measurements one line per gene, described in the format below.")
    str6 <- paste("&emsp;&emsp;  a) Column_1. This should contain the gene name, for the user's reference.")
    str7 <- paste("&emsp;&emsp;  b)   Column_ 2.  This should contain the gene group classification e.g. O/U for Over-expressed/Under-expressed or Hyper/Hypo for hypermethylated/hypomethylated in alphabetical order. If only one gene group, use any alphabet e.g. A for each row instead.")
    str8 <- paste("&emsp;&emsp;  c)   Remaining Columns. These should contain the expression measurements as numbers. Data inputted should be non-negative. Columns and rows with zero variance should be removed from the data. Rows containing missing expression measurements, should be also be removed from the input data or it will cause the tool to run into errors." )
    str9 <- paste("NOTE: Clustering is based on scaled data, if the user chooses this option, prior to input into heatmap R function.")
    str10 <- paste("Example format for Data")
    HTML(paste(str0,h4(strong(str00)), str1,str0,h5(strong(str2)),str3,str4,str5, str6, str7,str8, str0, str0, strong(str9), str0,str0,str0,strong(em(str10)), str0,sep = '<br/>'))
  })
  
  output$Eg <- renderTable({
    colna1 <- c("gene_id", "Groups","TCGA.01.98GF", "TCGA.08.U5TD", "TCGA.02.D23F", "TCGA.01.12TD","TCGA.02.AOKO", "TCGA.12.T37D", "TCGA.16.Y2S5", "TCGA.01.KITD")
    colna2 <- c(" ", " ", rep("Normal", 3), rep("Tumor", 5))
    s1 <- c("BRCA1", "over", 1.47, 2.18, 5.87, 7.64, 3.40, 7.77, 5.15, 1.56 )
    s2 <- c("YWHAE", "over", 7.93, 2.76, 9.11, 6.96, 5.98, 8.19, 8.91, 0.98)
    s3 <- c("SFN1", "under",8.02, 8.00, 2.17, 1.12, 3.76, 0.02, 3.67, 9.76 )
    s4 <- c("BRAF", "under", 2.75, 5.99, 3.19, 3.09, 2.00, 0.99, 1.28, 8.17)
    
    d <- rbind.data.frame(colna2, s1, s2, s3, s4)
    colnames(d) <- colna1
    head(d)
  })
  
  output$Caption1 <- renderUI({
    str.0 <- paste("&emsp;")
    str.1 <- paste("Table 1 : Example dataset for two gene groups (over and under-expressed) and two patient groups (Normal, Tumor).")
    HTML(paste(strong(str.1),str.0, str.0, str.0, sep = '<br/>'))
  })
  
  output$Eg2 <- renderTable({
    coln1 <- c("gene_id", "Groups","GSM9981", "GSM1870", "GSM4618", "GSM7689", "GSM8772", "GSM1121","GSM1250", "GSM3112", "GSM4987", "GSM1277")
    coln2 <- c(" ", " ", rep("MM", 5), rep("MUGS", 2), "NPC", rep("SM", 2))
    s.1 <- c("YWHAE|210996_s_at", "A", 1.47, 2.18, 5.87, 9.12, 7.34, 1.56, 3.0, 7.77, 3.40, 1.56 )
    s.2 <- c("YWHAE|201020_at", "A", 1.98, 7.93, 2.76, 9.11, 8.46, 0.98, 5.98, 8.19, 8.91, 5.98)
    s.3 <- c("YWHAH|33323_r_at", "A", 8.02, 8.00, 2.17, 10.12, 8.76, 9.76, 3.76, 0.02, 3.67, 7.94)
    s.4 <- c("YWHAB|208743_s_at", "A", 2.75, 5.99, 3.19, 11.86, 6.54, 8.17, 2.00, 0.99, 2.00, 1.17)
    s.5 <- c("YWHAQ|213699_s_at", "A", 9.35, 8.96, 6.67, 8.33, 3.98, 7.11, 1.67, 1.01, 5.18, 8.17)
    
    d.f <- rbind.data.frame(coln2, s.1, s.2, s.3, s.4, s.5)
    colnames(d.f) <- coln1
    head(d.f)
  })
  
  
  output$Caption2 <- renderUI({
    strg.0 <- paste("&emsp;")
    strg.1 <- paste("Table 2: Example dataset for one gene group (marked A) and four patient groups (MM, MUGS, NPC and, SM).")
    strg.2 <- paste(" Where applicable, both row and column dendrograms can be extracted in their specific tabs. Using the options on the right panel, dendrograms can be cut into desired no. of clusters (default at 2) and a pvalue of significance between the clusters can be determined using bootstrap method. ")
    HTML(paste(strong(strg.1),strg.0, strg.0,strg.2, strg.0, strg.0, sep = '<br/>'))
  })
  
  output$download_GW_Ex1 <- downloadHandler(
    filename= function() {paste('CoMMpassIA9.GW.Expression.data_truncated.csv')}, 
    content = function(file) {
      d <- readRDS("data/CoMMpassIA9.GW.Expression.data_truncated.rds")
      write.csv(d, file, row.names = FALSE) }
  )
  
  output$download_GW_Ex2 <- downloadHandler(
    filename= function() {paste('CoMMpassIA9_Expression_demo_data.csv')}, 
    content = function(file) {
      d <- readRDS("data/CoMMpassIA9.GW.Expression.data_truncated.rds")
      write.csv(d, file, row.names = FALSE) }
  )
  
  
  input_gw_data <- reactive({
    if(input$gw_file1 == 'GW_Example1'){
      d <- readRDS("data/CoMMpassIA9.GW.Expression.data_truncated.rds")
    }
    else if(input$gw_file1 == 'load_my_own_gw'){
      inFile <- input$gw_file2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { d = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
      else if(grepl(".txt", inFile[1])) { d = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
      else if(grepl(".rds", inFile[1])) { d = read.table(as.character(inFile$datapath)) }
      
    }
    else 
      return(NULL)
    # dim(data)
    Dataset <- data.frame(d)
    return(Dataset)
  })
  
  gw_data <- reactive ({
    if(!is.null(input_gw_data()))
    {
      data <- input_gw_data()
      data2 <- data[-1,]
      
      rownames(data2) = paste(data2$gene_id, data2$Groups, sep = "|")
      data2 <- data2[, c(-1, -2)]
      colnames(data2) = paste(names(data2), data[1,c(-1, -2)], sep = "||")
      data2 <- data.frame(data.matrix(data2))
      n= ncol(data)-2
     
      data2$var <- apply(data2[, c(1:n)], 1, var)
      data2$mad <- apply(data2[, c(1:n)], 1, mad)
      data2$IQR <- apply(data2[, c(1:n)], 1, IQR)
      data2$Rank.var <- rank(data2$var,ties.method= "min")
      data2$Rank.mad <- rank(data2$mad,ties.method= "min")
      data2$Rank.IQR <- rank(data2$IQR,ties.method= "min")
      data2$sumofranks <- data2$Rank.var + data2$Rank.mad + data2$Rank.IQR
     
      
      if(is.null(input$gw_subset)) {
        return()
      } else if(length(input$gw_subset) == 1) {
         if(input$gw_subset == 'VAR') {
          data2 <- data2[order(data2$var),]
          data2$var_percen <- ifelse(input$var_PercenChoice == 'Percentile Slider', quantile(data2$var, as.numeric(input$var_pslider)/100),quantile(data2$var, as.numeric(input$var_pInput)/100))
         } else if(input$gw_subset == 'MAD') {
          data2 <- data2[order(data2$mad),]
          data2$mad_percen <- ifelse(input$mad_PercenChoice == 'Percentile Slider', quantile(data2$mad, input$mad_pslider/100),quantile(data2$mad, input$mad_pInput/100))
         } else if(input$gw_subset == 'IQR') {
          data2 <- data2[order(data2$IQR),]
          data2$iqr_percen <- ifelse(input$iqr_PercenChoice == 'Percentile Slider', quantile(data2$IQR, input$iqr_pslider/100),quantile(data2$IQR, input$iqr_pInput/100))
         } 
      }
     else 
       if(length(input$gw_subset) > 1) {
         
        if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset) {
           data2 <- data2[order(data2$sumofranks),]
           data2$IMVA_percen <- ifelse(input$IMVA_PercenChoice == 'Percentile Slider', quantile(data2$sumofranks, as.integer(input$IMVA_pslider)/100),quantile(data2$sumofranks, as.integer(input$IMVA_pInput)/100))
        } else if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & !("IQR" %in% input$gw_subset) ) {
           data2$sumofranks_VM <- data2$Rank.var + data2$Rank.mad 
           data2 <- data2[order(data2$sumofranks_VM),]
           data2$var_mad_percen <- ifelse(input$var_mad_PercenChoice == 'Percentile Slider', quantile(data2$sumofranks_VM, input$var_mad_pslider/100),quantile(data2$sumofranks_VM, input$var_mad_pInput/100))
        } else if("VAR" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("MAD" %in% input$gw_subset)) {
           data2$sumofranks_VI <- data2$Rank.var + data2$Rank.IQR 
           data2 <- data2[order(data2$sumofranks_VI),]
           data2$var_iqr_percen <- ifelse(input$var_iqr_PercenChoice == 'Percentile Slider', quantile(data2$sumofranks_VI, input$var_iqr_pslider/100),quantile(data2$sumofranks_VI, input$var_iqr_pInput/100))
        } else  if("MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("VAR" %in% input$gw_subset)) {
           data2$sumofranks_MI <- data2$Rank.mad + data2$Rank.IQR 
           data2 <- data2[order(data2$sumofranks_MI),]
           data2$mad_iqr_percen <- ifelse(input$mad_iqr_PercenChoice == 'Percentile Slider', quantile(data2$sumofranks_MI, input$mad_iqr_pslider/100),quantile(data2$sumofranks_MI, input$mad_iqr_pInput/100))
         }  
           
         
     }
      return(data.frame(data2))
    }
  })
  
  output$geneSelector <- renderUI({
    selectizeInput(inputId = "Genes", "Choose Option:", as.list(getOSgenes()),options=list(maxOptions=getOSgenes())) 
  })
  
  output$dropdowngene <- renderText({ 
    paste("You have selected gene", input$Genes)
  })
  
  extracted_data <- reactive ({
    data2 <- gw_data()
    data <- input_gw_data()
    n= ncol(data)-2
    
    
    if(is.null(input$gw_subset)) {
      return()
    } else if(length(input$gw_subset) == 1) {
      if(input$gw_subset == 'VAR') {
        data3 <- data2[data2$var > data2$var_percen[1], 1:n]
      } else if(input$gw_subset == 'MAD') {
        data3 <- data2[data2$mad > data2$mad_percen[1], 1:n]
      } else if(input$gw_subset == 'IQR') {
        data3 <- data2[data2$IQR > data2$iqr_percen[1], 1:n]
      } 
    }
    else 
      if(length(input$gw_subset) > 1) {
        if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset) {
          data3 <- data2[data2$sumofranks > data2$IMVA_percen[1], 1:n]
        }else if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & !("IQR" %in% input$gw_subset)) {
          data3 <- data2[data2$sumofranks_VM > data2$var_mad_percen[1], 1:n]
        } else if("VAR" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("MAD" %in% input$gw_subset)) {
          data3 <- data2[data2$sumofranks_VI > data2$var_iqr_percen[1], 1:n]
        } else  if("MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("VAR" %in% input$gw_subset)) {
          data3 <- data2[data2$sumofranks_MI > data2$mad_iqr_percen[1], 1:n]
        }  
      }
    return(data.frame(data3))
  })
  
  getOSgenes <- reactive({
    if(!is.null(input_gw_data())) 
    {
      data <- input_gw_data()
      return(as.character(data[, 1]))
    }
    else 
      return(NULL)
  })
  
   output$n_selected <- renderUI({ 
      data3 <- extracted_data()
      st1 <- paste("The number of genes selected : ")
      st2 <- paste(nrow(data3))
      
      HTML(paste(st1, strong(st2)), sep = ' ')
     })
    
 
  
  output$Boxplot <- renderPlotly({
    data <- gw_data()
    n = ncol(input_gw_data())-2
    data$var <- apply(data[, c(1:n)], 1, var)
    data$mad <- apply(data[, c(1:n)], 1, mad)
    data$IQR <- apply(data[, c(1:n)], 1, IQR)
    y <- list(
      title = " ")
    pp <- plot_ly(data, y =~var, type = 'box', name = 'Var') %>%
      add_trace(y = ~mad, name = 'MAD')  %>%
      add_trace(y = ~IQR, name = 'IQR') %>%
      layout(yaxis = y)
    pp
  })
  
  output$GW_Scatter_LH <- renderPlotly({
      data <- gw_data()
      
      goi <- ifelse(input$Genes== "", NA, input$Genes)
       
      if(is.null(input$gw_subset)) {
        return()
      }
      
      if(length(input$gw_subset) == 1) {
        if(input$gw_subset == 'VAR') {   
           p <- plot_ly(data, y=~var, type = "scatter", mode = "markers", name = "Ordered variances") %>%
             add_trace(p, y = ~var_percen, line=list(dash=3, width= 1, color = "green"), mode = "lines", name = 'Percentile cut-off'  )
           if (!is.na(goi[1])) {
            p <- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$var, text = input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
           }
        p %>% layout(showlegend = TRUE)
        p
        }
        else if(input$gw_subset == 'MAD') {
          p<- plot_ly(data, y=~mad, type = "scatter", mode = "markers", name = "Ordered MAD") %>%
            add_trace(p, y = data$mad_percen, line=list(dash=3, width= 1, color = "green"), mode = "lines", name = 'Percentile cut-off' )
          if (!is.na(goi[1])) {
            p<- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$mad, text = input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
          }
          p %>% layout(showlegend = TRUE)
          p
        } else if(input$gw_subset == 'IQR') {
          p <- plot_ly(data, y=~IQR, type = "scatter", mode = "markers", name = "Ordered IQR") %>%
            add_trace(p, y = ~iqr_percen, line=list(dash=3, width= 1, color = "green"), mode = "lines", name ='Percentile cut-off' )
          if (!is.na(goi[1])) {
            p<- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$IQR, text = input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
          }
          p %>% layout(showlegend = TRUE)
          p
        } 
        
      } else 
        if(length(input$gw_subset) > 1) {
          
           if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset) {
            
            p <- plot_ly(data, y=~sumofranks, type = "scatter", mode = "markers", name = "Ordered sum of VAR, MAD and IQR") %>%
              add_trace(p, y = ~IMVA_percen, line=list(dash=3, width= 1, color = "green"), mode = "lines", name = 'Percentile cut-off' ) 
            if (!is.na(goi[1])) {
              p <- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$sumofranks, text= input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
            }
            p %>% layout(showlegend = TRUE)
            p
          }
          
          else if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & !("IQR" %in% input$gw_subset) ) {
            p <- plot_ly(data, y=~sumofranks_VM, type = "scatter", mode = "markers", name = "Ordered sum of VAR and MAD") %>%
              add_trace(p, y = ~var_mad_percen, line=list(dash=3, width= 1, color = "green" ), mode = "lines", name = 'Percentile cut-off'  )
            if (!is.na(goi[1])) {
              p <- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$sumofranks_VM, text= input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
            }
            p %>% layout(showlegend = TRUE)
            p
          }
          else if ("VAR" %in% input$gw_subset & "IQR" %in% input$gw_subset& !("MAD" %in% input$gw_subset)) {
         
            p <- plot_ly(data, y=~sumofranks_VI, type = "scatter", mode = "markers", name = "Ordered sum of VAR and IQR") %>%
              add_trace(p, y = ~var_iqr_percen, line=list(dash=3, width= 1, color = "green"), mode = "lines", name = 'Percentile cut-off'  )
            if (!is.na(goi[1])) {
              p <- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$sumofranks_VI, text= input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
            }
            p %>% layout(showlegend = TRUE)
            p
          }
          else if ("MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("VAR" %in% input$gw_subset)) {
            p <- plot_ly(data, y=~sumofranks_MI, type = "scatter", mode = "markers", name = "Ordered sum of MAD and IQR") %>%
              add_trace(p, y = ~mad_iqr_percen, line=list(dash=3, width= 1, color = "green"), mode = "lines", name = 'Percentile cut-off'  )
            if (!is.na(goi[1])) {
              p <- add_trace(p, x= which(grepl(input$Genes, rownames(data))), y= data[which(grepl(input$Genes, rownames(data))),]$sumofranks_MI, text= input$Genes, showlegend = TRUE, type = "scatter", mode = "markers", name = input$Genes, marker = list(color = "orange"))
            }
            p %>% layout(showlegend = TRUE)
            p
          } 
        }
      
   }) 
  
  extracted_data2 <- reactive({
    if(!is.null(extracted_data()))
    {
    data3 <- extracted_data()
    colheaders <- sub("^.*\\.","", colnames(data3))
    names(data3) <- gsub("\\.{2}.*", "", colnames(data3))
    data4 <- rbind(colheaders, data3)
    Groups <- sub("^.*\\|","", rownames(data4))
    Groups[1] <- " "
    gene_id <- sub("\\|.*$", "", rownames(data4)) 
    gene_id[1] <- " "
    
    #names(data4) <- gsub("\\.{2}.*", "", colnames(data4))
    data5 <- cbind(gene_id, Groups, data4)
    rownames(data5)[1] <- ""
    
    return(data.frame(data5))
     
    } 
    else
      return(NULL)
    
  })
  
  output$downloadSubset <- downloadHandler(
    filename= function() {paste(input$fname_subset, Sys.time(),'.csv', sep='')}, 
    content = function(file) {
      d <- extracted_data2()
      write.csv(d, file, row.names = FALSE) }
  )
    
  output$GW_subset_heatmap <- renderPlot({
    if(!is.null(extracted_data2()))
    {
      data5 <- extracted_data2()
      
      # sort columns based on colnames
      data <- data5[,order(data5[1, ])]
      data <- data[order(data[,1]),]
      
      ### gene names, column name as gene
      gene <- as.character(data$gene_id)
      gene <- gene[-1]
      row.groups <- as.character(as.vector(data[,2]))
      row.groups <- row.groups[-1]
      row.groups.name <- names(table(row.groups))
      number.row.groups <- length(row.groups.name)
      
      ### column groups
      col.groups <- as.character(as.vector(data[1,]))
      col.groups <- col.groups[c(-1, -2)] # calculate no. of column groups
      col.groups.name <- names(table(col.groups))
      number.col.groups <- length(col.groups.name)
      
      data <- data[-1, c(-1, -2)]
      rownames(data) <- gene
      data<-data[complete.cases(data),]
      data <- data.matrix(data)
      
      ## Set color palette
      col1 = colorRampPalette(c(input$low,input$mid,input$high))(299)
      colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
      
      ### Color vector for columns
      if(number.col.groups==1) { 
        cell <- c(rep(col.groups.name, number.col.groups))
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
      } else if(number.col.groups==2) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'  
      } else if(number.col.groups==3) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
      } else if(number.col.groups==4) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
      } else if(number.col.groups==5) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
      } else if(number.col.groups==6) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
      } else if(number.col.groups==7) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
      } else if(number.col.groups==8) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
      } else if(number.col.groups==9) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'darkorchid2'
      } else if(number.col.groups==10) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]]), rep(col.groups.name[10], table(col.groups)[[10]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'darkorchid2'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+table(col.groups)[[9]]+1:table(col.groups)[[10]]] <- 'maroon'
      }
      
      ### Color vector for rows
      
      if(number.row.groups==1) { 
        cell2 <- c(rep(row.groups.name, number.row.groups))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
      } else if(number.row.groups==2) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'  
      } else if(number.row.groups==3) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
      } else if(number.row.groups==4) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
      } else if(number.row.groups==5) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
      } else if(number.row.groups==6) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]]), rep(row.groups.name[6],table(row.groups)[[6]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+table(row.groups)[[5]]+1:table(row.groups)[[6]]] <- 'maroon'
      }
      
      ############# HEATPLOT2 EQUIVALENT HEATMAP2 CLUSTERING ###############
      z <- list()
      
      #data <- as.numeric(data)
      if(input$norm == "Z-Score") {
        z <- zClust(data, scale =input$norm2, zlim=c(input$inSlider[1],input$inSlider[2]))
      } else if (input$norm == "Modified Z-Score") { 
        z <- modzClust(data, scale =input$norm2, zlim=c(input$inSlider[1],input$inSlider[2]))
        check_z_mod  <<- z[[1]]
      } else if(input$norm == "none") {
        z[[1]] <- as.matrix(data)
      }
      
      if(input$dist == "pearson correlation") {
        if(input$dispRow == "No" & input$dispCol=='No') {
          hm <<- heatmap.2(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        } else if(input$dispRow == "No" & input$dispCol=='Yes' ) {
          hm <<- heatmap.2(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        } else if(input$dispRow == "Yes" & input$dispCol=='No') {
          hm <<- heatmap.2(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
          hm <<- heatmap.2(z[[1]], Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        }
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } 
      }  else {
        if(input$dispRow == "No" & input$dispCol=='No') {
          hm <<- heatmap.2(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        } else if(input$dispRow == "No" & input$dispCol=='Yes') {
          hm <<- heatmap.2(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        } else if(input$dispRow == "Yes" & input$dispCol=='No') {
          hm <<- heatmap.2(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
          hm <<- heatmap.2(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
        }
        
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        }
      }
      
      coldendo <- reactive({
        par(cex = input$sizeClable)
        dend1 <- as.dendrogram(hm$colDendrogram)
        d <- data.frame(v1 =hm$colInd, v2=1:length(hm$colInd))
        m <- data.frame(v3 = 1:length(cc1), v4 = cc1)
        
        colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
        colbar <- colbar[,2]
        labels_colors(dend1) <- as.character(colbar)
        plot(dend1)
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        }
      })
      
      output$plot1 <- renderPlot({
        coldendo()
      })
      
      
      colDen <- reactive({
        if(input$cutcolden == 'TRUE') {
          cuttable <- as.data.frame(cutree(as.hclust(hm$colDendrogram), k=as.numeric(input$cuttree))[as.hclust(hm$colDendrogram)$order])
          cuttable <- cbind.data.frame(rownames(cuttable), cuttable)
          names(cuttable)[1] <- "Sample"
          names(cuttable)[2] <- "Cluster"
          data_l1_l2 <- extracted_data2()
          data_l1_l2 <- data_l1_l2[1,c(-1, -2)]
          t_data_l1_l2 <- t(data_l1_l2)
          t_data_l1_l2 <- cbind.data.frame(rownames(t_data_l1_l2), t_data_l1_l2)
          names(t_data_l1_l2)[1] <- "Sample"
          names(t_data_l1_l2)[2] <- "Group"
          m.cut.data <- merge(cuttable, t_data_l1_l2, by = "Sample", sort= F)
          m.cut.data <- m.cut.data[, c(1, 3, 2)]
        }
        else {
          return(NULL)
        } 
        
      })
      
      output$display <- renderUI({
        if(input$cutcolden != 'TRUE') {
          return(br(strong(em("Please select Cut Col dendrogram?: = 'Yes' to display column clusters. Also select value at which you would like to cut the col dendogram (default is at k= 2)"))))
        }
      })
      
      output$df <- DT::renderDataTable({
        if(input$cutcolden == 'TRUE') {
          DT::datatable(colDen(), options = list(
            lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
            pageLength = 5))
        }
      })
      
      output$downloadCuttree <- downloadHandler(
        filename = function() {
          paste(paste(input$fname_HM, input$hclust, "clustering", input$dist, "distance", sep="_"), '_Col_Dendrogram_cutree_', 'k=', input$cuttree, '.csv', sep='') 
        },
        content = function(con) {
          write.csv(colDen(), con, quote=F, row.names = F)
        })
      
      rowdendo <- reactive({
        par(cex = input$sizeRlable)
        dend2 <- as.dendrogram(hm$rowDendrogram)
        dd <- data.frame(v1 =rev(hm$rowInd), v2=1:length(hm$rowInd))
        mm <- data.frame(v3 = 1:length(cc2), v4 = cc2)
        
        colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
        colbar2 <- colbar2[,2]
        labels_colors(dend2) <- rev(as.character(colbar2))
        plot(dend2, horiz = T)
        if(number.row.groups==1) {
          legend("topright", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==2) {
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==3) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==4) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==5) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==6) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } 
      })
      
      output$plot2 <- renderPlot({
        par(cex= input$sizeRlable)
        rowdendo()
      })
      
      rowDen <- reactive({
        if(input$cutrowden == 'TRUE') {
          cuttable2 <- as.data.frame(cutree(as.hclust(hm$rowDendrogram), k=as.integer(input$cuttree2))[as.hclust(hm$rowDendrogram)$order])
          cuttable2 <- cbind.data.frame(rownames(cuttable2), cuttable2)
          names(cuttable2)[1] <- "gene_id"
          names(cuttable2)[2] <- "Cluster"
          data_l1_l2_2 <- extracted_data2()
          data_l1_l2_2 <- data_l1_l2_2[-1, c(1,2)]
          m2 <- merge(cuttable2, data_l1_l2_2, by = "gene_id", sort= F)
          m2 <- m2[, c(1, 3, 2)]
        }
        else {
          return(NULL)
        }
        
      })
      
      output$display2 <- renderUI({
        if(input$cutrowden != 'TRUE') {
          return(br(strong(em("Please select Cut Row dendrogram?: = 'Yes' to display row clusters. Also select value at which you would like to cut the row dendrogram (default is at k= 2)"))))
        }
      })
      
      output$df2 <- DT::renderDataTable({
        if(input$cutrowden == 'TRUE') {
          DT::datatable(rowDen(), options = list(
            lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
            pageLength = 5))
        }
      })
      
      output$downloadCuttree2 <- downloadHandler(
        filename = function() {
          paste(paste(input$fname_HM, input$hclust, "clustering", input$dist, "distance", sep="_"), '_Row_Dendrogram_cutree_', 'k=', input$cuttree2, '.csv', sep='') 
        },
        content = function(con) {
          write.csv(rowDen(), con, quote=F, row.names = F)
        })
    
  
  output$downloadHM <- downloadHandler(
    
    filename <- function() {
      pdf_file <<- paste(input$fname_HM, input$hclust, "clustering", input$dist, "distance", sep="_")
      paste('NOJAH_', pdf_file, Sys.time(),'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file=paste(pdf_file,".pdf",sep="") , height= 10, width=10)
      plot.new()
      title("NOJAH: Clustering Analysis",cex.main=1.2, sub = "Genome Wide analysis" ,cex.main=1.2, col = "blue", font=3)
      df <- rbind.data.frame(c("Data Normalization Type", input$norm),
                             c("Normalized by (row/column/both)", input$norm2),
                             c("Distance Method", input$dist),
                             c("Clustering Method", input$hclust),
                             c("Scale", ifelse(input$norm == "none", paste(as.integer(min(data)), as.integer(max(data)), sep = ":"), paste(input$inSlider[1], input$inSlider[2], sep=":"))),
                             c("HeatMap colors", paste(input$low, input$mid, input$high, sep="-")))
      names(df)[1] <- "Parameters"
      names(df)[2] <- "Value Selected"
      grid.table(df, rows= NULL)
      
      eval(hm$call) # call the heatmap here
      if(number.col.groups==1) {
        legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==2) {
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==3) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==4) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==5) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==6) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==7) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==8) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7], col.groups.name[8] )), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==9) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7], col.groups.name[8], col.groups.name[9] )), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==10) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      }
      
      if(number.row.groups==1) {
        legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==2) {
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==3) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==4) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==5) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==6) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } 
      
      par(cex = 0.6*input$sizeClable)
      #par(mar=c(5,7,4,2))
      dend1 <- as.dendrogram(hm$colDendrogram)
      d <- data.frame(v1 =hm$colInd, v2=1:length(hm$colInd))
      m <- data.frame(v3 = 1:length(cc1), v4 = cc1)
      
      colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
      colbar <- colbar[,2]
      labels_colors(dend1) <- as.character(colbar)
      plot(dend1, main="Column Dendrogram")
      if(number.col.groups==1) {
        legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==2) {
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==3) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==4) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==5) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==6) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==7) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==8) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7], col.groups.name[8] )), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==9) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7], col.groups.name[8], col.groups.name[9] )), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.col.groups==10) {  
        legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6],  col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      }
      
      par(cex = input$sizeRlable)
      dend2 <- as.dendrogram(hm$rowDendrogram)
      dd <- data.frame(v1 =rev(hm$rowInd), v2=1:length(hm$rowInd))
      mm <- data.frame(v3 = 1:length(cc2), v4 = cc2)
      
      colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
      colbar2 <- colbar2[,2]
      labels_colors(dend2) <- rev(as.character(colbar2))
      plot(dend2, horiz = T, main="Row Dendrogram")
      if(number.row.groups==1) {
        legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==2) {
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==3) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==4) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==5) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.row.groups==6) {  
        legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } 
      dev.off()
      file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
     }
     )
    }
    else
      return(NULL)
  })
  
  Exp_input <- reactive({
    if(input$Exp_file == 'Exp_example'){
      d <- read.csv("data/Most_variable_extracted_Expression.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
    }
    else if(input$Exp_file == 'Exp_load_my_own'){
      inFile <- input$Exp_file2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { d = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
      else if(grepl(".txt", inFile[1])) { d = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
    }
    else 
      return(NULL)
    Dataset <- data.frame(d)
    return(as.data.frame(Dataset))
  })
  
  output$Exp_download <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_Expression_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds2 <- Exp_input()
      write.csv(ds2, file, row.names = FALSE)
    }
  )
  
  Variant_input <- reactive({
    if(input$Variant_file == 'Variant_example'){
      d2 <- read.csv("data/Most_variable_extracted_Variant_transposed.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
    }
    else if(input$Variant_file == 'Variant_load_my_own'){
      inFile <- input$Variant_file2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
    }
    else 
      return(NULL)
    Dataset2 <- data.frame(d2)
    return(as.data.frame(Dataset2))
  })
  
  output$Variant_download <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_Variant_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds3 <- Variant_input()
      write.csv(ds3, file, row.names = FALSE)
    }
  )
  
  CNV_input <- reactive({
    if(input$CNV_file == 'CNV_example'){
      d3 <- read.csv("data/Most_variable_extracted_CNV.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
    }
    else if(input$CNV_file == 'CNV_load_my_own'){
      inFile <- input$CNV_file2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { d3 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
      else if(grepl(".txt", inFile[1])) { d3 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d3)
    return(as.data.frame(Dataset3))
  })
  
  output$CNV_download <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_CNV_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds2 <- CNV_input()
      write.csv(ds2, file, row.names = FALSE)
    }
  )
  
  indiv <- reactive({
    if(!is.null(Exp_input()))
    {
     exp_data <- Exp_input()
     exp_data2 <- as.matrix(exp_data)
     cc1 <- consensus_clustering(dinput=exp_data2, mK=10, rep=as.integer(input$Exp_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$Exp_dist, iL=input$Exp_hclust, fL=input$Exp_hclust)
    }
     return(list(output= cc1[["output"]][[as.integer(input$Exp_opt_k)]]$consensusClass, data= cc1[["data"]], distance=cc1[["distance"]], order= cc1[["output"]][[as.integer(input$Exp_opt_k)]]$consensusTree$order))
   })
  
  output$Exp_cc <- renderPlot({
    exp_data <- Exp_input()
    exp_data2 <- as.matrix(exp_data)
    par(mfrow= c(1, 3))
    indiv()
  })
 
  output$Exp_sil <- renderPlot({
    silhouette_plot(data_use= indiv()[["data"]], opt_k=as.integer(input$Exp_opt_k), res=indiv()[["output"]], dist = indiv()[["distance"]] )
  })
  
  output$Exp_cc_dl <- downloadHandler(
      filename <- function(){
        pdf_file1 <<- paste("Expression", input$Exp_dist, input$Exp_hclust, sep = "_")
        paste("ConsensusClustring_Results_", pdf_file1,'.pdf', sep='')
        },
      content <- function(file) {
        pdf(file=paste(pdf_file1,".pdf",sep=""))
         exp_data <- Exp_input()
         exp_data2 <- as.matrix(exp_data)
         cc1 <- consensus_clustering(dinput=exp_data2, mK=10, rep=as.integer(input$Exp_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$Exp_dist, iL=input$Exp_hclust, fL=input$Exp_hclust)
         silhouette_plot(data_use= indiv()[["data"]], opt_k=as.integer(input$Exp_opt_k), res=indiv()[["output"]], dist = indiv()[["distance"]] )
        dev.off()
        file.copy(paste(pdf_file1,'.pdf', sep='') ,file, overwrite=TRUE)
        }
    )
  
  indiv2 <- reactive({
    if(!is.null(Variant_input()))
    {
      variant_data <- Variant_input()
      variant_data2 <- as.matrix(variant_data)
      cc2 <- consensus_clustering(dinput=variant_data2, mK=10, rep=as.integer(input$Variant_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$Variant_dist, iL=input$Variant_hclust, fL=input$Variant_hclust)
    }
    return(list(output= cc2[["output"]][[as.integer(input$Variant_opt_k)]]$consensusClass, data= cc2[["data"]], distance=cc2[["distance"]],order= cc2[["output"]][[as.integer(input$Variant_opt_k)]]$consensusTree$order))
  })
  
  output$Variant_cc <- renderPlot({
    variant_data <- Variant_input() ### check this!
    variant_data2 <- as.matrix(variant_data)
    par(mfrow= c(1, 3))
    indiv2()
  })
  
  output$Variant_sil <- renderPlot({
    silhouette_plot(data_use= indiv2()[["data"]], opt_k=as.integer(input$Variant_opt_k), res=indiv2()[["output"]], dist = indiv2()[["distance"]])
  })
  
  output$Variant_cc_dl <- downloadHandler(
    filename <- function(){
      pdf_file2 <<- paste("Variant", input$Variant_dist, input$Variant_hclust, sep = "_")
      paste("ConsensusClustring_Results_", pdf_file2,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file=paste(pdf_file2,".pdf",sep=""))
      variant_data <- Exp_input()
      variant_data2 <- as.matrix(variant_data)
      cc2 <- consensus_clustering(dinput=variant_data2, mK=10, rep=as.integer(input$Variant_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$Variant_dist, iL=input$Variant_hclust, fL=input$Variant_hclust)
      silhouette_plot(data_use= indiv2()[["data"]], opt_k=as.integer(input$Variant_opt_k), res=indiv2()[["output"]], dist = indiv2()[["distance"]] )
      dev.off()
      file.copy(paste(pdf_file2,'.pdf', sep='') ,file, overwrite=TRUE)
    }
  )
  
  indiv3 <- reactive({
    if(!is.null(CNV_input()))
    {
      CNV_data <- CNV_input()
      CNV_data2 <- as.matrix(CNV_data)
      cc3 <- consensus_clustering(dinput=CNV_data2, mK=10, rep=as.integer(input$CNV_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$CNV_dist, iL=input$CNV_hclust, fL=input$CNV_hclust)
    }
    return(list(output= cc3[["output"]][[as.integer(input$CNV_opt_k)]]$consensusClass, data= cc3[["data"]], distance=cc3[["distance"]], order= cc3[["output"]][[as.integer(input$CNV_opt_k)]]$consensusTree$order))
  })
  
  output$CNV_cc <- renderPlot({
    CNV_data <- CNV_input()
    CNV_data2 <- as.matrix(CNV_data)
    par(mfrow= c(1, 3))
    indiv3()
  })
  
  output$CNV_sil <- renderPlot({
    silhouette_plot(data_use= indiv3()[["data"]], opt_k=as.integer(input$CNV_opt_k), res=indiv3()[["output"]], dist = indiv3()[["distance"]] )
  })
  
  output$CNV_cc_dl <- downloadHandler(
    filename <- function(){
      pdf_file3 <<- paste("CNV", input$CNV_dist, input$CNV_hclust, sep = "_")
      paste("ConsensusClustring_Results_", pdf_file3,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file=paste(pdf_file3,".pdf",sep=""))
      CNV_data <- CNV_input()
      CNV_data2 <- as.matrix(CNV_data)
      cc3 <- consensus_clustering(dinput=CNV_data2, mK=10, rep=as.integer(input$CNV_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$CNV_dist, iL=input$CNV_hclust, fL=input$CNV_hclust)
      silhouette_plot(data_use= indiv3()[["data"]], opt_k=as.integer(input$CNV_opt_k), res=indiv3()[["output"]], dist = indiv3()[["distance"]] )
      dev.off()
      file.copy(paste(pdf_file3,'.pdf', sep='') ,file, overwrite=TRUE)
    }
  )
  
  
  coca_input <- reactive({
    if(input$coca_file == 'coca_example'){
      cc1 <- indiv()
      cc2 <- indiv2()
      cc3 <- indiv3()
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) { 
        d4 <- cbind.data.frame(cc1[["output"]], cc2[["output"]], cc3[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c("Expression", "Variant", "CNV")
      }  else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) { 
        d4 <- cbind.data.frame(cc1[["output"]], cc2[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c("Expression", "Variant")
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
        d4 <- cbind.data.frame(cc1[["output"]], cc3[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c("Expression", "CNV")
        
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
        d4 <- cbind.data.frame(cc2[["output"]], cc3[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c("Variant", "CNV")
        
      }
     }
    else if(input$coca_file == 'coca_load_my_own'){
      inFile <- input$coca_file2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { d4 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
      else if(grepl(".txt", inFile[1])) { d4 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T, row.names = 1) }
    }
    else 
      return(NULL)
    Dataset4 <- data.frame(d4)
    return(as.data.frame(Dataset4))
  })
  
  output$coca_download <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_coca_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds4 <- coca_input()
      write.csv(ds4, file, row.names = T)
    }
  )
  
  combined <- reactive({
    validate(
      need(!(input$coca_platform == "EXP" & input$coca_platform != "PROP" & input$coca_platform != "CNV"),"Please select atleast two platforms e.g. Expression and Variant"),
      need(!(input$coca_platform != "EXP" & input$coca_platform == "PROP" & input$coca_platform != "CNV"), "Please select atleast two platforms e.g. Expression and Variant"),
      need(!(input$coca_platform != "EXP" & input$coca_platform != "PROP" & input$coca_platform == "CNV"), "Please select atleast two platforms e.g. Expression and Variant")
    )
    if(!is.null(coca_input()))
    {
      cc = coca_input()
      cc1 <- list(output = cc[1,])
      cc2 <- list(output = cc[2,])
      cc3 <- list(output = cc[3,])
     
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        cc.final <- coca(cc= list(cc1, cc2, cc3), type = c("E", "V", "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", "V"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("V", "CNV"), opt_k= c(as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      }
    }
    dev.off()
    
    return(list(output= cc.final[["output"]], data= cc.final[["data"]], distance = cc.final[["distance"]]))
    
  })
  
  
  
  output$coca_cc <- renderPlot({
    #coca_data <- coca_input()
    #coca_data2 <- as.matrix(coca_data)
    par(mfrow= c(1, 3))
    combined()
  })
  
  output$coca_cc_dl <- downloadHandler(
    filename <- function(){
      pdf_file4 <<- paste("CoCA", input$coca_dist, input$coca_hclust, sep = "_")
      paste("ConsensusClustring_Results_", pdf_file4,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file=paste(pdf_file4,".pdf",sep=""))
      cc = coca_input()
      cc1 = list(output = cc[1,])
      cc2 = list(output = cc[2,])
      cc3 = list(output = cc[3,])
      
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        cc.final <- coca(cc= list(cc1, cc2, cc3), type = c("E", "V", "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", "V"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("V", "CNV"), opt_k= c(as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      }
      silhouette_plot(data_use= combined()[["data"]], opt_k=as.integer(input$coca_opt_k), res=combined()[["output"]], dist = combined()[["distance"]] )
      dev.off()
      file.copy(paste(pdf_file4,'.pdf', sep='') ,file, overwrite=TRUE)
    })
  
  combined2 <- reactive({
    
    if(!is.null(combined()))
    {
      cc.final <- combined()
      df = as.data.frame(cc.final[["output"]])
      dat = as.data.frame(cc.final[["data"]])
      tdf = t(df)
      df.dat = rbind(tdf, dat)
      rownames(df.dat)[1] <- "CC_Cluster"
    }
    return(as.data.frame(df.dat))
  })
  
  output$coca_sil <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_coca_clusters', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds5 <- combined2()
      ds5.t <- t(ds5)
      ds5.t2 <- ds5.t[, 1] 
      write.csv(ds5.t2, file, row.names = T)
    }
  )
  
  colbardata <- reactive({
    if(input$clinical == 'available'){
      d5 <- read.csv("data/Clinical_file_example.csv", header =T, sep =",", stringsAsFactors = T)
      
    }
    else if(input$clinical == 'load_my_own_c'){
      inFile <- input$clinical_file
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { d5 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = T, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d5 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = T, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset5 <- data.frame(d5)
    return(as.data.frame(Dataset5))
  })
  
  output$coca_heatmap <- renderPlot({
    if(!is.null(combined2())) {
      clust_data <- combined2()
      colbar_data <- colbardata()
      
      if(!is.null(colbardata()))
      {
      colbar_data <- colbar_data[order(colbar_data[,1]),]
      
      
      t.clust_data = as.data.frame(t(clust_data))
      t.clust_data$Sample = rownames(t.clust_data)
      m <- merge(t.clust_data, colbar_data, by = "Sample" )
      #clust_data = cbind(rownames(clust_data), clust_data)
      
      # sort columns based on colnames
      if(nrow(m) == ncol(clust_data))
      {
      data <- clust_data[,order(clust_data[1, ])]
      data <- data[order(rownames(data)),]
      
      data = cbind(rownames(data), data)
      
      ### gene names, column name as gene
      gene <- as.character(data[,1])
      gene <- gene[-1]
      row.groups <- as.character(as.vector(gsub("\\_.*","",data[,1])))
      row.groups <- row.groups[-1]
      row.groups.name <- names(table(row.groups))
      number.row.groups <- length(row.groups.name)
      
      ### column groups
      col.groups <- as.character(as.vector(data[1,-1]))
      #col.groups <- col.groups[c(-1)] # calculate no. of column groups
      col.groups.name <- names(table(col.groups))
      number.col.groups <- length(col.groups.name)
      
      data <- data[-1, c(-1)]
      rownames(data) <- gene
      
      data<-data[complete.cases(data),]
      data <- data.matrix(data)
      
      ## Set color palette 
      col1 = colorRampPalette(c("khaki","black","deepskyblue"))(299)
      colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
      
      ### Color vector for columns
      if(number.col.groups==1) { 
        cell <- c(rep(col.groups.name, number.col.groups))
        
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        colbars2 <- if(n.colbars ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==2) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'  
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==3) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        #colbars2 <- cbind(cc1, colbars(df2 = colbar_data))
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==4) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==5) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==6) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        colbars2 <- as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==7) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==8) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
      } else if(number.col.groups==9) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'darkorchid2'
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==10) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]]), rep(col.groups.name[10], table(col.groups)[[10]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'darkorchid2'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+table(col.groups)[[9]]+1:table(col.groups)[[10]]] <- 'maroon'
        colbars2 <- if(is.null(input$cli_feature)) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- colbars(df2 = colbar_data)[[3]]
        colors_used <- colbars(df2 = colbar_data)[[4]]
      }
      
      
      ### Color vector for rows
      
      if(number.row.groups==1) { 
        cell2 <- c(rep(row.groups.name, number.row.groups))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        rowbars2 <- as.matrix(t(cc2))
      } else if(number.row.groups==2) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray48'  
        rowbars2 <- as.matrix(t(cc2))
      } else if(number.row.groups==3) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray48'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        rowbars2 <- as.matrix(t(cc2))
      } else if(number.row.groups==4) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray48'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        rowbars2 <- as.matrix(t(cc2))      
      } else if(number.row.groups==5) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray48'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        rowbars2 <- as.matrix(t(cc2))
     } else if(number.row.groups==6) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]]), rep(row.groups.name[6],table(row.groups)[[6]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray48'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+table(row.groups)[[5]]+1:table(row.groups)[[6]]] <- 'maroon'
        rowbars2 <- as.matrix(t(cc2))
        
      }
      
      #heatmap clustering
     
      if(input$dist_4 == "pearson correlation") {
        if(input$dispRow_4 == "No" & input$dispCol_4=='No') {
          hm2 <- heatmap.3(data, labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(x) hclust(x,method=input$hclust_4), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1_4,cexCol =input$size2_4,  key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2, RowSideColors = rowbars2, RowSideColorsSize = 1) 
        } else if(input$dispRow_4 == "No" & input$dispCol_4=='Yes' ) {
          hm2 <- heatmap.3(data, labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(x) hclust(x,method=input$hclust_4), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1_4,cexCol =input$size2_4,key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2, RowSideColors = rowbars2, RowSideColorsSize = 1) 
        } else if(input$dispRow_4 == "Yes" & input$dispCol_4=='No') {
          hm2 <- heatmap.3(data, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(x) hclust(x,method=input$hclust_4), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1_4,cexCol =input$size2_4,  key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2,  ColSideColorsSize = 2, RowSideColors = rowbars2, RowSideColorsSize = 1) 
        } else if(input$dispRow_4 == "Yes" & input$dispCol_4=='Yes') {
          hm2 <- heatmap.3(data, Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), scale="none", hclust=function(x) hclust(x,method=input$hclust_4), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1_4,cexCol =input$size2_4,key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2, RowSideColors = rowbars2, RowSideColorsSize = 1) 
        }
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n" )
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "red", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "red", "orange", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "red", "orange", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n", inset = c(0, .1))
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray48"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray48", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray48", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray48", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray48", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } 
        
        if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class[1]), col = paste(colors_used[1]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = paste(c(colors_used[1], colors_used[2])), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }else if(number.colbar.class==3) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = paste(c(colors_used[1], colors_used[2], colors_used[3])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
        }
        
      }  else {
        if(input$dispRow_4 == "No" & input$dispCol_4=="No") {
          hm2 <- heatmap.3(data, labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(c) {hclust(c,method=input$hclust_4)}, distfun=function(c) {dist(c,method=input$dist_4)},cexRow=input$size1_4,cexCol =input$size2_4,key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2, RowSideColors = rowbars2, RowSideColorsSize = 1) 
        } else if(input$dispRow_4 == "No" & input$dispCol_4=='Yes') {
          hm2 <- heatmap.3(data, labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(c) {hclust(c,method=input$hclust_4)}, distfun=function(c) {dist(c,method=input$dist_4)},cexRow=input$size1_4,cexCol =input$size2_4,key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2, RowSideColors = rowbars2, RowSideColorsSize = 1) 
        } else if(input$dispRow_4 == "Yes" & input$dispCol_4=='No') {
          hm2 <- heatmap.3(data, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(c) {hclust(c,method=input$hclust_4)}, distfun=function(c) {dist(c,method=input$dist_4)},cexRow=input$size1_4,cexCol =input$size2_4,key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2,RowSideColors = rowbars2, RowSideColorsSize = 1) 
        } else if(input$dispRow_4 == "Yes" & input$dispCol_4=='Yes') {
          hm2 <- heatmap.3(data, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow_4))), Colv=eval(parse(text=paste(input$clust_bycol_4))), hclust=function(c) {hclust(c,method=input$hclust_4)}, distfun=function(c) {dist(c,method=input$dist_4)},cexRow=input$size1_4,cexCol =input$size2_4,key=TRUE,keysize=1.0, margin = c(input$inSlider2_4[1],input$inSlider2_4[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, ColSideColorsSize = 2,RowSideColors = rowbars2, RowSideColorsSize = 1) 
        }
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n" )
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "red"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "red", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "red", "orange", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "red", "orange", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n", inset = c(0, .1))
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray48"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray48", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray48", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n" )
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray48", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray48", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } 
        
        if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class[1]), col = paste(colors_used[1]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = paste(c(colors_used[1], colors_used[2])), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }else if(number.colbar.class==3) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = paste(c(colors_used[1], colors_used[2], colors_used[3])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
        }
      } 
      
    output$download_coca_HM <- downloadHandler(
      
      filename <- function() {
        pdf_file5 <<- paste("CoCA", input$hclust_4, "clustering", input$dist_4, "distance", sep="_")
        paste('CoCA_', pdf_file5, Sys.time(),'.pdf', sep='')
      },
      content <- function(file) {
        pdf(file=paste(pdf_file5,".pdf",sep="") , height= 10, width=10)
        plot.new()
        title("NOJAH: Clustering Analysis", sub= "CoC Analysis HeatMap",cex.main=1.2, col = "blue", font=3)
        df <- rbind.data.frame(c("Distance Method", input$dist_4),
                               c("Clustering Method", input$hclust_4),
                              # c("Scale", paste(input$inSlider_4[1], input$inSlider_4[2], sep=":")),
                               c("HeatMap colors", paste("Khaki","Black", "Blue", sep="-")))
        names(df)[1] <- "Parameters"
        names(df)[2] <- "Value Selected"
        grid.table(df, rows= NULL)
        
        eval(hm2$call) # call the heatmap here
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n" )
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "red"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "red", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "red", "orange", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "red", "orange", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "red", "orange", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n", inset = c(0, .1))
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray48"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray48", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray48", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n" )
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray48", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray48", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable, bty = "n")
        } 
        
        if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class[1]), col = paste(colors_used[1]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = paste(c(colors_used[1], colors_used[2])), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
        }else if(number.colbar.class==3) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = paste(c(colors_used[1], colors_used[2], colors_used[3])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
        }
        
        dev.off()
        
        file.copy(paste(pdf_file5,'.pdf', sep='') ,file, overwrite=TRUE)
      }
     )
    } else
      return(NULL)
      }
    }
 })
  
  output$com_text1 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Determination of number of clusters by Consensus Clustering")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
 output$com_text12 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Determination of number of clusters by Consensus Clustering")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
 output$com_text13 <- renderUI({
   hs1 <- paste("&emsp;")
   hs2 <- paste("Determination of number of clusters by Consensus Clustering")
   HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
 })
 output$com_text14 <- renderUI({
   hs1 <- paste("&emsp;")
   hs2 <- paste("Determination of number of clusters by Consensus Clustering")
   HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
 })
 
 
  
  output$com_text2 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Silhouette Plot")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text22 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Silhouette Plot")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text23 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Silhouette Plot")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text3 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("CoC HeatMap")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text31 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Cluster Interpretation")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$varmean_bxplot <- renderPlot({
   if(!is.null(combined())){
    if("EXP" %in% input$coca_platform) {
      cc_clust1 <- as.matrix(indiv()[["output"]])
      colnames(cc_clust1)[1] <- "x"
    
      expression <- Exp_input()
      ID1 <- 1:ncol(expression)
      expression2 <- as.matrix(rbind(ID1, t(cc_clust1), expression))
      expression_order <- indiv()[["order"]]
      exp_dist <- indiv()[["distance"]]
    }

    if("PROP" %in% input$coca_platform) {
      cc_clust2 <- as.matrix(indiv2()[["output"]])
      colnames(cc_clust2)[1] <- "x"
    
      variant <- Variant_input()
      ID2 <- 1:ncol(variant)
      variant2 <- as.matrix(rbind(ID2, t(cc_clust2), variant))
      variant_order <- indiv2()[["order"]]
      variant_dist <- indiv2()[["distance"]]
    }
    
    if("CNV" %in% input$coca_platform) { 
      cc_clust3 <- as.matrix(indiv3()[["output"]])
      colnames(cc_clust3)[1] <- "x"
    
      cnv <- CNV_input()
      ID3 <- 1:ncol(cnv)
      cnv2 <- as.matrix(rbind(ID3, t(cc_clust3), cnv))
      cnv_order <- indiv3()[["order"]]
      cnv_dist <- indiv3()[["distance"]]
    }
     
    #plots
    if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        plotMeans(data= list(expression2, variant2, cnv2), data_order=list(expression_order, variant_order, cnv_order), 
                  type= c("Expression", "Variant", "CNV"), dist= c(exp_dist, variant_dist, cnv_dist) )
    } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) {
        plotMeans(data= list(expression2, variant2), data_order=list(expression_order, variant_order) , 
                  type= c("Expression", "Variant"), dist= c(exp_dist, variant_dist))   
    } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
        plotMeans(data= list(expression2, cnv2), data_order=list(expression_order, cnv_order) , 
                type= c("Expression", "CNV"), dist= c(exp_dist, cnv_dist))     
    } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
        plotMeans(data= list(variant2, cnv2), data_order=list(variant_order, cnv_order) , 
                type= c("Variant", "CNV"), dist= c(variant_dist, cnv_dist)) 
    }
  } else {
      return(NULL)
  }
  })
  
  output$dl_coca_inter<- downloadHandler(  
    filename <- function() {
      pdf_file6 <<- "CoCA clustering Interpretation"
      paste('CoCA_', pdf_file6, Sys.time(),'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file=paste(pdf_file6,".pdf",sep="") , height= 8, width=10)
      if("EXP" %in% input$coca_platform) {
        cc_clust1 <- as.matrix(indiv()[["output"]])
        colnames(cc_clust1)[1] <- "x"
        
        expression <- Exp_input()
        ID1 <- 1:ncol(expression)
        expression2 <- as.matrix(rbind(ID1, t(cc_clust1), expression))
        expression_order <- indiv()[["order"]]
        exp_dist <- indiv()[["distance"]]
      }
      
      if("PROP" %in% input$coca_platform) {
        cc_clust2 <- as.matrix(indiv2()[["output"]])
        colnames(cc_clust2)[1] <- "x"
        
        variant <- Variant_input()
        ID2 <- 1:ncol(variant)
        variant2 <- as.matrix(rbind(ID2, t(cc_clust2), variant))
        variant_order <- indiv2()[["order"]]
        variant_dist <- indiv2()[["distance"]]
      }
      
      if("CNV" %in% input$coca_platform) { 
        cc_clust3 <- as.matrix(indiv3()[["output"]])
        colnames(cc_clust3)[1] <- "x"
        
        cnv <- CNV_input()
        ID3 <- 1:ncol(cnv)
        cnv2 <- as.matrix(rbind(ID3, t(cc_clust3), cnv))
        cnv_order <- indiv3()[["order"]]
        cnv_dist <- indiv3()[["distance"]]
      }
      
      #plots
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        plotMeans(data= list(expression2, variant2, cnv2), data_order=list(expression_order, variant_order, cnv_order), 
                  type= c("Expression", "Variant", "CNV"), dist= c(exp_dist, variant_dist, cnv_dist) )
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) {
        plotMeans(data= list(expression2, variant2), data_order=list(expression_order, variant_order) , 
                  type= c("Expression", "Variant"), dist= c(exp_dist, variant_dist))   
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
        plotMeans(data= list(expression2, cnv2), data_order=list(expression_order, cnv_order) , 
                  type= c("Expression", "CNV"), dist= c(exp_dist, cnv_dist))     
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
        plotMeans(data= list(variant2, cnv2), data_order=list(variant_order, cnv_order) , 
                  type= c("Variant", "CNV"), dist= c(variant_dist, cnv_dist)) 
      }
      dev.off()
      file.copy(paste(pdf_file6,'.pdf', sep='') ,file, overwrite=TRUE)
    }
  )
  
  
  #### Significance tesing of clusters
  output$download_Sig_Ex1 <- downloadHandler(
    filename= function() {paste('Example data set_CoMMpass_IA9_filtered_expression_ds.csv')}, 
    content = function(file) {
    	d <- read.csv("data/Most_variable_extracted_Expression_withHRgroups.csv", header = T, sep  = ",", stringsAsFactors = F)
        write.csv(d, file, row.names = FALSE) }
  )

  output$download_Sig_Ex2 <- downloadHandler(
    filename= function() {paste('Example data set_TCGA_BRCA_Meth.csv')}, 
    content = function(file) {
      d <- read.csv("data/BRCA.Example.data_made_up.csv", header = T, sep  = ",", stringsAsFactors = F)
      write.csv(d, file, row.names = FALSE) }
  )

  
  Sig_data_input <- reactive({
    if(input$Sig_file1 == 'Sig_Example1'){
	d <- read.csv("data/Most_variable_extracted_Expression_withHRgroups.csv", header = T, sep  = ",", stringsAsFactors = F)
      } else if(input$Sig_file1 == 'Sig_Example2') {
	d <- read.csv("data/BRCA.Example.data_made_up.csv", header = T, sep  = ",", stringsAsFactors = F)
    } else if(input$Sig_file1 == 'load_my_own_Sig'){
      inFile <- input$Sig_file2
      if (is.null(inFile))
        return(NULL)
      #else if(grepl(".xlsx", inFile[1])) { d = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
      else if(grepl(".txt", inFile[1])) { d = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
      else if(grepl(".rds", inFile[1])) { d = read.table(as.character(inFile$datapath)) }
      
    }
    else 
      return(NULL)
    # dim(data)
    Dataset <- data.frame(d)
    return(Dataset)
  })
  
  # plot HM for Significance testing
  
 
  
  output$Sig_plot <- renderPlot({
    
    if(!is.null(Sig_data_input()))
    {
      data <- Sig_data_input()
      
      # sort columns based on colnames
      data <- data[,order(data[1, ])]
      data <- data[order(data[,2]),]
      
      ### gene names, column name as gene
      gene <- as.character(data$gene_id)
      gene <- gene[(as.integer(input$DataR)-1):nrow(data)] #gene[(6-1):nrow(data)] 
      row.groups <- as.character(as.vector(data[(as.integer(input$DataR)-1): nrow(data),2])) #as.character(as.vector(data[(6-1):nrow(data),2])) 
      row.groups.name <- names(table(row.groups))
      number.row.groups <- length(row.groups.name)
      
      ### column groups
      col.groups <- as.character(as.vector(data[1,]))
      #col.groups <- col.groups[c(-1, -2)] # calculate no. of column groups
      col.groups <-  col.groups[as.numeric(input$DataC):ncol(data)] #col.groups[4:ncol(data)] 
      col.groups.name <- names(table(col.groups))
      number.col.groups <- length(col.groups.name)
      
      ## additional info on columns and rows
      colbar_data <-  as.matrix(t(data[1:(as.numeric(input$DataR)-2), c(1, (as.integer(input$DataC)-1):ncol(data))])) #t(data[1:(6-2), c(1, 3:ncol(data))]) 
      colnames(colbar_data) <- as.character(unlist(colbar_data[1,]))
      colbar_data <- as.matrix(colbar_data[c(-1, -2),])
      colnames(colbar_data)[1] <- "Groups"
      n.colbar_data <- ncol(colbar_data)
      
      rowbar_data <- as.matrix(data[(as.integer(input$DataR)-1):nrow(data), 1:(as.numeric(input$DataC)-1)]) #as.matrix(data[5:nrow(data), 1:3])
      rownames(rowbar_data) <- rowbar_data[,1]
      rowbar_data <- as.matrix(rowbar_data[,-1])
      #colnames(rowbar_data)[1] <- "Groups"
      n.rowbar_data <- ncol(rowbar_data)
      
      data <- data[(as.numeric(input$DataR)-1):nrow(data), as.numeric(input$DataC):ncol(data)] #data[5:nrow(data), 4:ncol(data)]
      rownames(data) <- gene
      data<-data[complete.cases(data),]
      data <- data.matrix(data)
   
    
      
      ## Set color palette
      col1 = colorRampPalette(c(input$Sig_low,input$Sig_mid,input$Sig_high))(299)
      colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
      
      ### Color vector for columns
      if(number.col.groups==1) { 
        cell <- c(rep(col.groups.name, number.col.groups))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==2) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue' 
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
       
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==3) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        #colbars2 <- cbind(cc1, colbars(df2 = colbar_data))
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==4) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==5) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==6) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==7) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==8) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==9) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'darkorchid2'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[4]]
      } else if(number.col.groups==10) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]]), rep(col.groups.name[10], table(col.groups)[[10]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'firebrick4'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'dodgerblue'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'khaki1'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'gray'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'brown'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'darkorchid2'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+table(col.groups)[[9]]+1:table(col.groups)[[10]]] <- 'maroon'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "Group"
        number.colbar.class <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else  colbars(df2 = colbar_data)[[4]]
      }
      
      ### Color vector for rows
      
      if(number.row.groups==1) { 
        cell2 <- c(rep(row.groups.name, number.row.groups))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[4]]
      } else if(number.row.groups==2) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[4]]
        
      } else if(number.row.groups==3) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[4]]
        
      } else if(number.row.groups==4) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[4]]
      } else if(number.row.groups==5) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[4]]
        
      } else if(number.row.groups==6) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]]), rep(row.groups.name[6],table(row.groups)[[6]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+table(row.groups)[[5]]+1:table(row.groups)[[6]]] <- 'maroon'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars(df2 = rowbar_data)[[4]]
        
      }
      
      ############# HEATPLOT2 EQUIVALENT HEATMAP2 CLUSTERING ###############
      z <- list()
      
      #data <- as.numeric(data)
      if(input$Sig_norm == "Z-Score") {
        z <- zClust(data, scale =input$Sig_norm2, zlim=c(input$Sig_inSlider[1],input$Sig_inSlider[2]))
      } else if (input$Sig_norm == "Modified Z-Score") { 
        z <- modzClust(data, scale =input$Sig_norm2, zlim=c(input$Sig_inSlider[1],input$Sig_inSlider[2]))
      } else if(input$Sig_norm == "none") {
        z[[1]] <- as.matrix(data)
      }
      
     
      if(input$Sig_dist == "pearson correlation") {
        if(input$Sig_dispRow == "No" & input$Sig_dispCol=='No') {
          Sig_hm <<- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,  key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "No" & input$Sig_dispCol=='Yes' ) {
          Sig_hm <<- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='No') {
          Sig_hm <<- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,  key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='Yes') {
          Sig_hm <<- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        }
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } 
      }  else {
        if(input$Sig_dispRow == "No" & input$Sig_dispCol=='No') {
          Sig_hm <<- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "No" & input$Sig_dispCol=='Yes') {
          Sig_hm <<- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='No') {
          Sig_hm <<- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='Yes') {
          Sig_hm <<- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        }
        
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        }
        
        if(number.row.groups==1) {
          legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(number.row.groups==2) {
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==3) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==4) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==5) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==6) {  
          legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        }
      }
      
      Sig_coldendo <- reactive({
        par(cex = input$Sig_sizeClable)
        dend1 <- as.dendrogram(Sig_hm$colDendrogram)
        d <- data.frame(v1 =Sig_hm$colInd, v2=1:length(Sig_hm$colInd))
        m <- data.frame(v3 = 1:length(cc1), v4 = cc1)
        
        colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
        colbar <- colbar[,2]
        labels_colors(dend1) <- as.character(colbar)
        plot(dend1)
        
        if(number.col.groups==1) {
          legend("topright", legend = paste(col.groups.name), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(number.col.groups==2) {
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(number.col.groups==3) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(number.col.groups==4) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(number.col.groups==5) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(number.col.groups==6) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==7) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        }
      })
      
      output$Sig_plot1 <- renderPlot({
        Sig_coldendo()
      })
      
      
      Sig_colDen <- reactive({
        if(input$Sig_cutcolden == 'TRUE') {
          cuttable <- as.data.frame(cutree(as.hclust(Sig_hm$colDendrogram), k=as.numeric(input$Sig_cuttree))[as.hclust(Sig_hm$colDendrogram)$order])
          cuttable <- cbind.data.frame(rownames(cuttable), cuttable)
          names(cuttable)[1] <- "Sample"
          names(cuttable)[2] <- "Cluster"
          data_l1_l2 <- Sig_data_input()
          data_l1_l2 <- data_l1_l2[1,c(-1, -2)]
          t_data_l1_l2 <- t(data_l1_l2)
          t_data_l1_l2 <- cbind.data.frame(rownames(t_data_l1_l2), t_data_l1_l2)
          names(t_data_l1_l2)[1] <- "Sample"
          names(t_data_l1_l2)[2] <- "Group"
          m.cut.data <- merge(cuttable, t_data_l1_l2, by = "Sample", sort= F)
          m.cut.data <- m.cut.data[, c(1, 3, 2)]
        }
        else {
          return(NULL)
        } 
        
      })
      
      output$Sig_display <- renderUI({
        if(input$Sig_cutcolden != 'TRUE') {
          return(br(strong(em("Please select Cut Col dendrogram?: = 'Yes' to display column clusters. Also select value at which you would like to cut the col dendogram (default is at k= 2)"))))
        }
      })
      
      output$Sig_df <- DT::renderDataTable({
        if(input$Sig_cutcolden == 'TRUE') {
          DT::datatable(Sig_colDen(), options = list(
            lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
            pageLength = 5))
        }
      })
      
      output$Sig_downloadCuttree <- downloadHandler(
        filename = function() {
          paste(paste(input$Sig_fname, input$Sig_hclust, "clustering", input$Sig_dist, "distance", sep="_"), '_Col_Dendrogram_cutree_', 'k=', input$Sig_cuttree, '.csv', sep='') 
        },
        content = function(con) {
          write.csv(Sig_colDen(), con, quote=F, row.names = F)
        })
      
      output$Sig_pv <- renderUI({
        if(input$Sig_cutcolden == 'TRUE' & number.col.groups >= 2 ){
          HTML(paste("<br/>", br(strong(em(paste("Would you want to assess gene set significance in the separation of specimens into two clusters? (Yes/No)")))), sep = "<br>")) 
        }
        else 
          return(NULL)
      })
      
      
      output$Sig_pvalue <- renderUI({
        input$Sig_goButton 
        
        isolate(
          
          if(input$Sig_cutcolden == 'TRUE') {
            pobs.col <- numeric()
            perms.col <- numeric()
            hc.cols <- as.hclust(Sig_hm$colDendrogram)
            cut <- as.data.frame(cutree(hc.cols, k=as.numeric(input$Sig_cuttree))[hc.cols$order])
            #cuttable <- cut_table()
            cut <- cbind.data.frame(rownames(cut), cut)
            names(cut)[1] <- "Sample"
            names(cut)[2] <- "Cluster"
            data_l1_l2 <- Sig_data_input()
            data_l1_l2 <- data_l1_l2[1,c(-1, -2)]
            t_data_l1_l2 <- t(data_l1_l2)
            t_data_l1_l2 <- cbind.data.frame(rownames(t_data_l1_l2), t_data_l1_l2)
            names(t_data_l1_l2)[1] <- "Sample"
            names(t_data_l1_l2)[2] <- "Group"
            mer <- merge(cut,t_data_l1_l2, by = "Sample")
            mer <- mer[, c(1, 3, 2)]
            
            if(input$Sig_pvalue_cal == TRUE) 
            {
              if(input$Sig_file3 == 'Sig_Exp.Example'){
                s_data <- readRDS("data/CoMMpassIA9_GW_Expression_data.rds")
              } else if(input$Sig_file3 == 'Sig_Meth.Example') {
		 s_data <- readRDS("data/Meth27K.GW.BRCA.Example.data.rds")
	      } else {
                inFile2 <- input$Sig_file4
                if (is.null(inFile2))
                  return(NULL)
                else if(grepl(".xlsx", inFile2[1])) { s_data = read.xlsx(as.character(inFile2$datapath), colNames = TRUE, rowNames = F) }
                else if(grepl(".csv", inFile2[1])) { s_data = read.csv(as.character(inFile2$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
                else if(grepl(".txt", inFile2[1])) { s_data = read.table(as.character(inFile2$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
                else if(grepl(".rds", inFile2[1])) { s_data = readRDS(as.character(inFile2$datapath)) }
              }
              
              # Create a Progress object
              progress <- shiny::Progress$new()
              progress$set(message = "Computing data", value = 0)
              # Close the progress when this reactive exits (even if there's an error)
              on.exit(progress$close())
              
              updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                  value <- progress$getValue()
                  value <- value + (progress$getMax() - value) / 10
                }
                progress$set(value = value, detail = detail)
              }
              
              
              
              # Bootstrap data, and pass in the updateProgress function so that it can update the progress indicator.
              Sig_b1 <<- bootstrapfun(obsdata=mer, samplingdata=s_data, distmethod = input$Sig_dist, clustmethod= input$Sig_hclust, norm= input$Sig_norm, scale=input$Sig_norm2, n=as.numeric(input$Sig_n), k=as.numeric(input$Sig_cuttree), n.iter=input$Sig_n_iter, zlim=c(input$Sig_inSlider[1],input$Sig_inSlider[2]), sampler = "Column", updateProgress )
              
              hstring1 <- paste("&emsp;")
              hstring2 <- paste("The p-value to test the gene set significance in the separation of specimens into 2 clusters is =", Sig_b1$p.value, sep = " ")
              if(Sig_b1$p.value <= 0.05) {
                hstring3  <- paste("The gene set cluster is statistically significant, i.e., a random sample of CpG probes/gene sets of the same number is Not able to separate the specimens when compared to the CpG probes/gene sets of interest of the same class")
              } else {
                hstring3  <- paste("The gene set cluster is NOT statistically significant, i.e., a random sample of CpG probes/gene sets of the same number is able to separate the specimens when compared to the CpG probes/gene sets of interest of the same class")
              }
              
              HTML(paste(hstring1, h5(strong(hstring2)), h5(em(hstring3)), hstring1, hstring1,  sep = '<br/>'))
              
              
            }
          }
          
        )
        
        
      })
      
      Sig_rowdendo <- reactive({
        par(cex = input$Sig_sizeRlable)
        dend2 <- as.dendrogram(Sig_hm$rowDendrogram)
        dd <- data.frame(v1 =rev(Sig_hm$rowInd), v2=1:length(Sig_hm$rowInd))
        mm <- data.frame(v3 = 1:length(cc2), v4 = cc2)
        
        colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
        colbar2 <- colbar2[,2]
        labels_colors(dend2) <- rev(as.character(colbar2))
        plot(dend2, horiz = T)
        if(number.row.groups==1) {
          legend("topright", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==2) {
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==3) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==4) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==5) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(number.row.groups==6) {  
          legend("topright", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } 
      })
      
      output$Sig_plot2 <- renderPlot({
        par(cex= input$Sig_sizeRlable)
        Sig_rowdendo()
      })
      
      Sig_rowDen <- reactive({
        if(input$Sig_cutrowden == 'TRUE') {
          cuttable2 <- as.data.frame(cutree(as.hclust(Sig_hm$rowDendrogram), k=as.integer(input$Sig_cuttree2))[as.hclust(Sig_hm$rowDendrogram)$order])
          cuttable2 <- cbind.data.frame(rownames(cuttable2), cuttable2)
          names(cuttable2)[1] <- "gene_id"
          names(cuttable2)[2] <- "Cluster"
          data_l1_l2_2 <- Sig_data_input()
          data_l1_l2_2 <- data_l1_l2_2[-1, c(1,2)]
          m2 <- merge(cuttable2, data_l1_l2_2, by = "gene_id", sort= F)
          m2 <- m2[, c(1, 3, 2)]
        }
        else {
          return(NULL)
        }
        
      })
      
      output$Sig_display2 <- renderUI({
        if(input$Sig_cutrowden != 'TRUE') {
          return(br(strong(em("Please select Cut Row dendrogram?: = 'Yes' to display row clusters. Also select value at which you would like to cut the row dendrogram (default is at k= 2)"))))
        }
      })
      
      output$Sig_df2 <- DT::renderDataTable({
        if(input$Sig_cutrowden == 'TRUE') {
          DT::datatable(Sig_rowDen(), options = list(
            lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
            pageLength = 5))
        }
      })
      
      output$Sig_downloadCuttree2 <- downloadHandler(
        filename = function() {
          paste(paste(input$Sig_fname, input$Sig_hclust, "clustering", input$Sig_dist, "distance", sep="_"), '_Row_Dendrogram_cutree_', 'k=', input$Sig_cuttree2, '.csv', sep='') 
        },
        content = function(con) {
          write.csv(Sig_rowDen(), con, quote=F, row.names = F)
        })
      
      
      output$Sig_pv2 <- renderUI({
        if(input$Sig_cutrowden == 'TRUE' & number.row.groups >= 2 ){
          HTML(paste("<br/>", paste("Would you want to assess significance of patients in the separation of genes into two clusters? (Yes/No)"), sep = "<br>")) 
        }
        else 
          return(NULL)
      })
      
      
      output$Sig_pvalue2 <- renderUI ({
        input$Sig_goButton2 
        
        isolate(
          
          if(input$Sig_cutrowden== TRUE){
            m2 <- Sig_rowDen()
            
            if(input$Sig_pvalue_cal2 == 'TRUE') 
            {
              if(input$Sig_file5 == 'Sig_Exp.Example2') {
                s_data <- readRDS("data/CoMMpassIA9_GW_Expression_data.rds")
             } else if(input$Sig_file5 == 'Sig_Meth.Example2') {
                s_data <- readRDS("data/Meth27K.GW.BRCA.Example.data.rds")
             } else {
                inFile3 <- input$Sig_file6
                if (is.null(inFile3))
                  return(NULL)
                #else if(grepl(".xlsx", inFile3[1])) { s_data = read.xlsx(as.character(inFile3$datapath), colNames = TRUE, rowNames = F) }
                else if(grepl(".csv", inFile3[1])) { s_data = read.csv(as.character(inFile3$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
                else if(grepl(".txt", inFile3[1])) { s_data = read.table(as.character(inFile3$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
                else if(grepl(".rds", inFile3[1])) { s_data = readRDS(as.character(inFile3$datapath)) }
                
              }  
              
              # Create a Progress object
              progress <- shiny::Progress$new()
              progress$set(message = "Computing data", value = 0)
              # Close the progress when this reactive exits (even if there's an error)
              on.exit(progress$close())
              
              updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                  value <- progress$getValue()
                  value <- value + (progress$getMax() - value) / 10
                }
                progress$set(value = value, detail = detail)
              }
             
              
              # Bootstrap data, and pass in the updateProgress function so that it can update the progress indicator.
              Sig_b2 <<- bootstrapfun(obsdata=m2, samplingdata=s_data, distmethod = input$Sig_dist, clustmethod= input$Sig_hclust, norm= input$Sig_norm, scale=input$Sig_norm2, n=as.numeric(input$Sig_n2), k=as.numeric(input$Sig_cuttree2), n.iter=input$Sig_n_iter2, zlim=c(input$Sig_inSlider[1],input$Sig_inSlider[2]), sampler = "Row", updateProgress )
              
              rhstring1 <- paste("&emsp;")
              rhstring2 <- paste("The p-value to test the sample significance in the separation of sample into 2 clusters is = =", Sig_b2$p.value, sep = " ")
              if(Sig_b2$Sig_p.value <= 0.05) {
                hstring3  <- paste("The cluster is statistically significant, i.e., a random sample of Sample sets of the same number is Not able to separate the gene sets when compared to the samples of interest of the same class")
              } else {
                hstring3  <- paste("The cluster is NOT statistically significant, i.e., a random sample of Sample sets of the same number is able to separate the gene sets when compared to the samples of interest of the same class")
              }
              
              HTML(paste(rhstring1, h5(strong(rhstring2)), h5(em(rhstring3)), rhstring1, rhstring1,  sep = '<br/>'))
              
              
              
            }
            
            
          }
          
        )
      })
      
      
      ############################
      # Download plots  #
      ############################
      output$Sig_downloadPlots <- downloadHandler(
        
        filename <- function() {
          Sig_pdf_file <<- paste(input$Sig_fname, input$Sig_hclust, "clustering", input$Sig_dist, "distance", sep="_")
          paste('NOJAH_', Sig_pdf_file, Sys.time(),'.pdf', sep='')
        },
        content <- function(file) {
          pdf(file=paste(Sig_pdf_file,".pdf",sep="") , height= 10, width=10)
          plot.new()
          title("NOJAH: Clustering Analysis",cex.main=1.2, sub = "Clustering Analysis with Significance", col = "blue", font=3)
          df <- rbind.data.frame(c("Data Normalization Type", input$Sig_norm),
                                 c("Normalized by (row/column/both)", input$Sig_norm2),
                                 c("Distance Method", input$Sig_dist),
                                 c("Clustering Method", input$Sig_hclust),
                                 c("Scale", ifelse(input$Sig_norm == "none", paste(as.integer(min(data)), as.integer(max(data)), sep = ":"), paste(input$Sig_inSlider[1], input$Sig_inSlider[2], sep=":"))),
                                 c("HeatMap colors", paste(input$Sig_low, input$Sig_mid, input$Sig_high, sep="-")))
          names(df)[1] <- "Parameters"
          names(df)[2] <- "Value Selected"
          grid.table(df, rows= NULL)
          
          eval(Sig_hm$call) # call the heatmap here
          if(number.col.groups==1) {
            legend("topright", legend = paste(col.groups.name), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==2) {
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==3) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==4) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==5) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==6) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==7) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==8) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==9) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==10) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          }
          
          if(number.row.groups==1) {
            legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==2) {
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==3) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==4) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==5) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==6) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
          
          par(cex = 0.6*input$Sig_sizeClable)
          #par(mar=c(5,7,4,2))
          dend1 <- as.dendrogram(Sig_hm$colDendrogram)
          d <- data.frame(v1 =Sig_hm$colInd, v2=1:length(Sig_hm$colInd))
          m <- data.frame(v3 = 1:length(cc1), v4 = cc1)
          
          colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
          colbar <- colbar[,2]
          labels_colors(dend1) <- as.character(colbar)
          plot(dend1, main="Column Dendrogram")
          if(number.col.groups==1) {
            legend("topright", legend = paste(col.groups.name), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==2) {
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==3) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==4) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==5) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==6) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==7) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==8) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==9) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==10) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          }
          
          par(cex = input$Sig_sizeRlable)
          dend2 <- as.dendrogram(Sig_hm$rowDendrogram)
          dd <- data.frame(v1 =Sig_hm$rowInd, v2=1:length(Sig_hm$rowInd))
          mm <- data.frame(v3 = 1:length(cc2), v4 = cc2)
          
          colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
          colbar2 <- colbar2[,2]
          labels_colors(dend2) <- as.character(colbar2)
          plot(dend2, horiz = T, main="Row Dendrogram")
          if(number.row.groups==1) {
            legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==2) {
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==3) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==4) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==5) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.row.groups==6) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
          
          if(input$Sig_goButton) {
            plot.new()
            
            title("Significance of Column Cluster", cex.main=2)
            df2 <- rbind.data.frame(c("Observed Data Fisher's exact p-value ", Sig_b1$p.obs),
                                    c("Sample Iterations", input$Sig_n_iter),
                                    c("Number of bootstrap Samples with replacement", input$Sig_n),
                                    c("Monte Carlo p-value", Sig_b1$p.value ), 
                                    c("Interpretation", ifelse(Sig_b1$p.value <= 0.05, paste("The CpG island/ gene set cluster is statistically significant", "i.e., a random sample of CpG probes/gene sets of the same number", " is Not able to separate the Sample groups when compared", "to the CpG islands/gene sets of interest", sep = '\n'), paste("The CpG island/ gene set cluster is NOT statistically significant", "i.e., a random sample of CpG probes/gene sets of the same number", " is able to separate the Sample groups when compared", "to the CpG islands/gene sets of interest", sep = '\n') ))
            )
            names(df2)[1] <- "Bootstrap approach parameters"
            names(df2)[2] <- "Value"
            grid.table(df2, rows= NULL)
          }
          
          if(input$Sig_goButton2) {
            plot.new()
            title("Significance of Row Cluster", cex.main=2)
            df3 <- rbind.data.frame(c("Observed Data Fisher's exact p-value", Sig_b2$p.obs),
                                    c("Sample Iterations", input$Sig_n_iter2),
                                    c("Number of bootstrap Samples with replacement", input$Sig_n2),
                                    c("Monte Carlo p-value", Sig_b2$p.value ),
                                    c("Interpretation", ifelse(Sig_b2$p.value <= 0.05, paste("The Sample cluster is statistically significant", "i.e., a random sample of Sample groups of the same number", " is Not able to separate the CpG islands/ gene sets when compared", "to the Sample groups of interest", sep = '\n'), paste("The Sample cluster is NOT statistically significant", "i.e., a random sample of Sample groups of the same number", " is able to separate the CpG islands/ gene sets when compared", "to the Sample groups of interest", sep = '\n') ))
            )
            names(df3)[1] <- "Bootstrap approach parameters"
            names(df3)[2] <- "Value"
            grid.table(df3, rows= NULL)
          }
          dev.off()
          
          file.copy(paste(Sig_pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
        }
      )
      
    } else {
      return(NULL)
    }
  })
  

}


shinyApp(ui = ui, server = server)