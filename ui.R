library(shiny)
library(shinythemes)
library(shinycssloaders)
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
library(matrixStats)
library(ggplot2)

fluidPage(
   tags$head(includeScript("google-analytics.js")),
   theme = shinytheme("spacelab"),
                
 
   navbarPage(h4(strong("NOJAH : NOt Just Another Heatmap")), 
              
  
    tabPanel("HomePage", 
             fluidRow(
               column(8, offset = 2, 
                      p("This interactive web application: NOt Just Another Heatmap (NOJAH) is developed in R with Shiny to"),
                      br(),
                      em("1) Perform Genome-Wide Heatmap (GWH) Analysis on any cancer genomic data set"),
                      br(),
                      em("2) Perform Combined results Clustering (CrC) Analysis for up to three different data types."),
                      br(),
                      em("3) Perform Significance of Cluster (SoC) Analysis using a robust bootstrap approach"), 
                      br(),
                      br(),
                      "The goal of this tool is to provide a one stop shop to perform genomic analyses.",
                      htmlOutput("ReadMe"), tableOutput("Eg"), htmlOutput("Caption1"), tableOutput("Eg2"), htmlOutput("Caption2")
                      ))),    
    tabPanel("Genome-Wide Heatmap (GWH) Analysis",
             fluidRow(
               column(2,
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(), 
                      br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(), 
                      br(), br(), br(), br(), br(),
              conditionalPanel("input.gw_conditionedPanels == 1",            
                wellPanel(h3("Input GW Data"),
                 selectInput("gw_file1",label= "Select an example dataset or upload your own with 'Load my own GW data.'", 
                             choices = c("Example TCGA BRCA Exp Data"="GW_Example1", 
                                         "Load my own data" = "load_my_own_gw")),
                 conditionalPanel("input.gw_file1 == 'load_my_own_gw'",
                                  fileInput('gw_file2', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                 conditionalPanel("input.gw_file1 == 'GW_Example1'", downloadButton('download_GW_Ex1', 'Download GW TCGA BRCA Expression DataSet'))
                 
                 ),
                #br(), br(), br(), br(), br(), br(), br(), br(), 
                br(), br(), br(), br(), br(),
                wellPanel(h3("Genome-Wide Dendrogram Display"),
                   radioButtons("PlotGW", "Plot GW Dendrogram (May take up to a few minutes to load)", c("No" = FALSE, "Yes" = TRUE))
                  ),
                 wellPanel(h3("Data Subsetting"),
                 checkboxGroupInput("gw_subset","Subset GW data by:",choices = list("Variance"= "VAR","Median Absoute Deviation" = "MAD", "Inter Quartile Range" = "IQR"), selected = c("IQR")),
                 conditionalPanel(condition = "$.inArray('VAR', input.gw_subset) > -1 & $.inArray('MAD', input.gw_subset) > -1 & $.inArray('IQR', input.gw_subset) > -1",
                                  radioButtons("IMVA_PercenChoice", "Percentile", c("Percentile Slider" = "Percentile Slider", "Manually Enter Percentile" = "Manually Enter Percentile")), ###
                                  conditionalPanel(condition = "input.IMVA_PercenChoice == 'Percentile Slider'",
                                                   sliderInput("IMVA_pslider", "Percentile Value:", 
                                                               min=0, max=100, value=75)),
                                  conditionalPanel(condition = "input.IMVA_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("IMVA_pInput", label = "Percentile value", min = 0, max = 100, value = 95, step = 5))),
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
                                                               min=0, max=100, value=99)),
                                  conditionalPanel(condition = "input.iqr_PercenChoice == 'Manually Enter Percentile'",
                                                   numericInput("iqr_pInput", label = "Percentile value", min = 0, max = 100, value = 75, step = 5))
                   ),
                  
                 actionButton("button1", "Run Analysis"),
                 h5("Hit button to update results after each change in input parameter(s)")
                 ),
                
               wellPanel(h3("Download subset data"),
                         textInput("fname_subset", "Type the file name you would like to save subset data as :", value = "Subset_data"),
                         downloadButton('downloadSubset', 'Download Subset data'))),
              conditionalPanel("input.gw_conditionedPanels == 2", 
               wellPanel(h3("Subset HeatMap Input"),
                         selectInput("gw_file_1",label= "Select an pre-existing dataset or upload your own with 'Load my own subset data.'", 
                                     choices = c("Pre-existing data subset"= "GW_Example_1", 
                                                 "Load my own data" = "load_my_own_gw_subset_1")),
                         conditionalPanel("input.gw_file_1 == 'load_my_own_gw_subset_1'",
                                          fileInput('gw_file_2', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                         conditionalPanel("input.gw_file_1 == 'GW_Example_1'", downloadButton('download_GW_Ex_1', 'Download subset Pre-existing DS')),
                         br(),
                         br(),
                         actionButton("button2", "Run Analysis"),
                         h5("Hit button to update results after each change in input parameter(s)")
               ),
               
               conditionalPanel("input.gw_file_1 == 'load_my_own_gw_subset_1'",
                                wellPanel(h3("Input Settings"),
                                          h4("Select first row and column where numeric data starts"),
                                          numericInput("DataR2", label = "Data row", min = 1, max = 10, value = 3),
                                          numericInput("DataC2", label = "Data column", min = 1, max = 10, value = 3)
                                )),
               
                wellPanel(h3("Download subset Heatmap"),
                  textInput("fname_HM", "Save Heatmap as :", value = "GW_subset_HM"),
                  downloadButton('downloadHM', 'Download HM'),
                  conditionalPanel("input.conditionedPanels == 2",   
                                 h5(downloadLink('downloadCuttree', 'Download Column clusters after cut-tree'))),
                  conditionalPanel("input.conditionedPanels == 3", 
                                 h5(downloadLink('downloadCuttree2', 'Download Row clusters after cut-tree')))) ),
              conditionalPanel("input.gw_conditionedPanels == 3",
                wellPanel (h3("Consensus Clustering Input"),
                selectInput("gw_file3",label= "Select an pre-existing dataset or upload your own with 'Load my own GW data.'", 
                            choices = c("Pre-existing data subset"= "GW_Example3", 
                                        "Load my own data" = "load_my_own_gw_subset")),
                conditionalPanel("input.gw_file3 == 'load_my_own_gw_subset'",
                                 fileInput('gw_file4', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                conditionalPanel("input.gw_file3 == 'GW_Example3'", downloadButton('download_GW_Ex3', 'Download subset Pre-existing DS')),
                br(),
                br(),
                actionButton("button3", "Run Analysis"),
                h5("Hit button to update results after each change in input parameter(s)")
               ),
              
                conditionalPanel("input.conditionedPanels == 1", 
                wellPanel(h3("Download Consensus Cluster Results"),
                          downloadButton("con_dl", label = "Download Consensus Cluster Results"),
                          h6("Download may take a while. Once complete, the result pdf file will automatically open.")) )),
              conditionalPanel("input.gw_conditionedPanels == 4",
                wellPanel (h3("Silhouette Input"),
                           selectInput("gw_file5",label= "Select an pre-existing dataset or upload your own with 'Load my own subset data.'", 
                                       choices = c("Pre-existing data subset"= "GW_Example5", 
                                                   "Load my own data" = "load_my_own_gw_subset2")),
                           conditionalPanel("input.gw_file5 == 'load_my_own_gw_subset2'",
                                            fileInput('gw_file6', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                           conditionalPanel("input.gw_file5 == 'GW_Example5'", downloadButton('download_GW_Ex5', 'Download subset Pre-existing silhouette DS')),
                           br(),
                           selectInput("gw_file7",label= "Select an pre-existing cluster or upload your own with 'Load my own clusters.'", 
                                       choices = c("Pre-existing data subset"= "GW_Example7", 
                                                   "Load my own clusters" = "load_my_own_gw_subset3")),
                           conditionalPanel("input.gw_file7 == 'load_my_own_gw_subset3'",
                                            fileInput('gw_file8', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                           conditionalPanel("input.gw_file7 == 'GW_Example7'", downloadButton('download_GW_Ex7', 'Download subset Pre-existing silhouette clusters')),
                           br(),
                           br(),
                           actionButton("button4", "Run Analysis"),
                           h5("Hit button to update results after each change in input parameter(s)")
                            
                ),
                wellPanel(h3("Download Silhouette Core Samples"),
                          downloadButton("download_GW_Ex6", label = "Download Silhouette Core Samples")
                          
                )
                
               
                ),
                
              conditionalPanel("input.gw_conditionedPanels == 5",
               wellPanel(h3("Input Core Data"),
                         selectInput("gw_file9",label= "Select an example dataset or upload your own with 'Load my own Core data.'", 
                                      choices = c("Pre-existing data subset"= "GW_Example9", 
                                     "Load my own data" = "load_my_own_core")),
                         conditionalPanel("input.gw_file9 == 'load_my_own_core'",
                                fileInput('gw_file10', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                         br(),
                        
                         #conditionalPanel("input.gw_file9 == 'Core_data'", downloadButton('download_cc_Ex1', 'Download DS'))
                         actionButton("button5", "Run Analysis"),
                         h5("Hit button to update results after each change in input parameter(s)")
                               ),
               conditionalPanel("input.gw_file9 == 'load_my_own_core'",
               wellPanel(h3("Input Settings"),
                         h4("Select first row and column where numeric data starts"),
                         numericInput("DataR1", label = "Data row", min = 1, max = 10, value = 5),
                         numericInput("DataC1", label = "Data column", min = 1, max = 10, value = 3)
               )),
                         
                         
              
                conditionalPanel("input.cc_conditionedPanels == 1 | input.cc_conditionedPanels == 2 | input.cc_conditionedPanels == 3",  
                #wellPanel(
                #         actionButton("button5", "Run Analysis")),
                
                wellPanel(h3("Core Sample based Downloads"),
                          br(), 
                          textInput("cc_fname_subset", "Type the file name you would like to save subset data as :", value = "Core_Subset_data"),
                          downloadButton('dlCoreSubset', 'Download Subset data'),
                          br(), br(),
                          textInput("cc_fname_HM", "Save HeatMap as :", value = "Core_subset_HM"),
                          downloadButton('cc_downloadHM', 'Download HM'),
                          conditionalPanel("input.cc_conditionedPanels == 2",   
                                           h5(downloadLink('cc_downloadCuttree', 'Download Column clusters after cut-tree'))),
                          conditionalPanel("input.cc_conditionedPanels == 3", 
                                           h5(downloadLink('cc_downloadCuttree2', 'Download Row clusters after cut-tree')))
                         
                          
                          ))
              )),               
                column(8, p(HTML("<p> <h3> A Genome-Wide Heatmap can be very dense. Given the limitation with the computational power required to construct a genome wide heatmap, NOJAH showcases a <font color = '#338AFF'> <i> Genome-Wide Dendrogram. </i> </font> </h3> <br> 
                                 <p> <h3> <strong> <u> Genome-Wide Heatmap Analysis workflow is divided into four main subparts: </strong> </u> </h3> 
                                 <p> <font color = '#338AFF'> <h3> 1. <i> Define Core Features with Most Variable Approach </i> <br>
                                 <p> 2. <i> Heatmap of Core Features </i> <br>
                                 <p> 3. <i> Define Cluster Number </i> <br>
                                 <p> 4. <i> Define Core Samples </i> </font> <br>
                                 <p> Heatmap is <i> <font color = '#338AFF'> updated</i> </font> based on the Core Features with Core Samples. <br> <br>
                                 <p> When using the analysis workflow, each step of the workflow is intended to be used sequentially i.e. the output of step 1 is fed into step 2 as input and so on. However each of these components can also be used independently.For example, if only consensus clustering needs to be performed then the 'Cluster Number' tab can be used. </h3> <br> <br> <br>"
                                 )),
                       
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
                tabPanel("Core Features",
                     conditionalPanel(condition= "input.PlotGW == 'TRUE'" , 
                        h3(strong("Genome-Wide Dendrogram")),
                        plotOutput("gw_dend", height = 500) %>% withSpinner(color = "#0275D8") #color="#0dc5c1")
                        ),
                     h3(strong("Measures of Spread")),
                     #plotlyOutput("Boxplot") %>% withSpinner(color="#0275D8"), 
                     plotOutput("Boxplot") %>% withSpinner(color="#0275D8"),
                     h3(strong("Define Core Features with Most Variable Approach")),
                     h5(em("To see the position of your 'gene of interest', enter gene id in the text box to the right")),
                     #plotlyOutput("GW_Scatter_LH") %>% withSpinner(color="#0275D8"),
                     plotOutput("GW_Scatter_LH") %>% withSpinner(color="#0275D8"),
                     htmlOutput("n_selected"), 
                     conditionalPanel(condition="input.button1 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=1),
                tabPanel("Heatmap",
                     h3(strong("Heatmap of Core Features")),
                     tabsetPanel(type = "tabs", 
                       tabPanel("Heatmap", withSpinner(plotOutput("GW_subset_heatmap", height = 1200)), 
                                conditionalPanel(condition="input.button2 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=1), 
                       tabPanel("Column Dendrogram",withSpinner(plotOutput("plot1", height = 800)), htmlOutput("display"), br(), DT::dataTableOutput("df"), 
                                conditionalPanel(condition="(input.button2 & input.button2) & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=2), 
                       tabPanel("Row Dendrogram",withSpinner(plotOutput("plot2", height = 800)), htmlOutput("display2"), br(), DT::dataTableOutput("df2"),  
                                conditionalPanel(condition="input.button2 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value =3), 
                       id = "conditionedPanels"), value =2),
                tabPanel("Cluster Number",
                     h3(strong("Define Number of Clusters")),
                     plotOutput("gw_cc", height = 600, width = 2400) %>% withSpinner(color = "#0275D8"), 
                     conditionalPanel(condition="input.button3 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value = 3),
                tabPanel("Core Samples",
                     h3(strong("Define Core Samples")),
                     plotOutput("gw_core_sam", height = 1800) %>% withSpinner(color = "#0275D8"), 
                     conditionalPanel(condition="input.button4 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value= 4),
                tabPanel("Updated Heatmap",
                     h3(strong("Heatmap of Core Features with Core Samples")),
                     tabsetPanel(type = "tabs", 
                                 tabPanel("Heatmap", withSpinner(plotOutput("cc_GW_subset_heatmap", height = 1200)), 
                                          conditionalPanel(condition="input.button5 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=1), 
                                 tabPanel("Column Dendrogram",withSpinner(plotOutput("cc_plot1", height = 800)), htmlOutput("cc_display"), br(), DT::dataTableOutput("cc_df"), 
                                          conditionalPanel(condition="input.button5 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=2), 
                                 tabPanel("Row Dendrogram",withSpinner(plotOutput("cc_plot2", height = 800)), htmlOutput("cc_display2"), br(), DT::dataTableOutput("cc_df2"), 
                                          conditionalPanel(condition="input.button5 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value =3),
                                 id = "cc_conditionedPanels"), value= 5),
                tabPanel("Workflow",
                      h3(strong("Workflow for Genome-Wide Heatmap (GWH) Analysis")),
                      DiagrammeROutput("diagram", height = 1000), 
                      conditionalPanel(condition="input.button3 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value = 6),
                
                id= "gw_conditionedPanels")
      
                     ),
               column(2, 
                
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                      br(), br(), br(), br(), br(), br(), br(),
                   
                    conditionalPanel("input.gw_conditionedPanels == 1",
                      conditionalPanel(condition= "input.PlotGW == 'TRUE'" ,
                      ########## GW dendrogram options ##########
                      wellPanel(  
                      h4("Genome-Wide Dendrogram Options"),
                      selectInput("GW_norm", "Normalization Type",
                                  c("Z-Score", "Modified Z-Score", "none"), selected = "Z-Score"),
                      selectInput("GW_norm2", "Normalize by:",
                                  c("row", "col", "both"), selected = "row"),
                      sliderInput("GW_inSlider", "Scale Range",
                                  min = -10, max = 20, value = c(-2, 2)),
                      selectInput("GW_dist", "Distance Method",
                                  c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson correlation"),
                      selectInput("GW_hclust", "Agglomerative Linkage Method", 
                                  c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "ward.D"),
                      sliderInput("gw_dend_size", "Dendrogram Label size:", min=0, max=1, value=0.4)
                      )),
                     br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                     br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                     
                     br(), br(), br(),  br(), br(), br(), br(), br(), br(),
                     #wellPanel(h3("Gene Display"),
                    #        radioButtons("Sgene", "Would you like to annotate genes on the scatter plot ?", c("No" = FALSE, "Yes" = TRUE))
                    # ),
                     #conditionalPanel(condition = "input.Sgene == 'TRUE'",
                                   #uiOutput(outputId="geneSelector", width=NULL)
                                   textInput("Genes", "Enter the 'gene id' you wish to annotate. Please note that input is case sensitive. Entering only gene name instead of gene id will throw error. For example, when using preloaded data enter: SCGB1D2>ENSG00000124935.3", value = ""),
                                   textOutput("text1")
                     #)
                  ),
                    
                 conditionalPanel("input.gw_conditionedPanels == 2",
                    conditionalPanel("input.conditionedPanels == 1", 
                       wellPanel(  
                        ########## HeatMap Clustering options ##########
                        h4("Heat Map Options"),
                        selectInput("norm", "Normalization Type",
                             c("Z-Score", "Modified Z-Score", "none")),
                        selectInput("norm2", "Normalize by:",
                             c("row", "col", "both"), selected = "both"),
                        sliderInput("inSlider", "Scale Range",
                             min = -10, max = 20, value = c(-2, 2)),
                        sliderInput("inSlider2", "Plot Margin dimensions",
                             min = 0, max = 20, value = c(14, 16)) ),
                        wellPanel(  
                          h4("Clustering Measures"),
                          selectInput("dist", "Distance Method",
                             c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson correlation"),
                          selectInput("hclust", "Agglomerative Linkage Method", 
                             c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "ward.D"),
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
                            ) ),
                    ##br(), br(),
                 conditionalPanel("input.gw_conditionedPanels == 3",
                    conditionalPanel("input.conditionedPanels == 1", 
                    wellPanel(h3("Consensus Clustering Options"),
                                               selectInput("con_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson"),
                                               selectInput("con_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                                               numericInput("con_maxK", "Maximum number of clusters:", 10),
                                               radioButtons("con_pItems", "No. of iterations", choices = c(10, 100, 500, 1000), selected = 100)),
                                     
                    wellPanel(h3("Choose Optimal Number of clusters"),
                            numericInput("con_opt_k", "Optimal k is ", value = 2))
                        )),
                 conditionalPanel("input.gw_conditionedPanels == 4",
                    wellPanel(h3("Silhouette Options"),
                              conditionalPanel("input.gw_file5 == 'load_my_own_gw_subset2'",
                                   selectInput("sil_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean")
                               ),
                              radioButtons("sil_choice", "Remove samples based on:", choices= c("Fixed value" = "Fixed value", "Percentile"= "Percentile", "Change point" = "Change point")),
                              conditionalPanel("input.sil_choice == 'Fixed value'",
                                    conditionalPanel("input.con_opt_k >= 1", sliderInput("upto_slider1", "Remove Samples within Cluster1 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.15)),
                                    conditionalPanel("input.con_opt_k >= 2", sliderInput("upto_slider2", "Remove Samples within Cluster2 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.34)),
                                    conditionalPanel("input.con_opt_k >= 3", sliderInput("upto_slider3", "Remove Samples within Cluster3 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0)),
                                    conditionalPanel("input.con_opt_k >= 4", sliderInput("upto_slider4", "Remove Samples within Cluster4 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0)),
                                    conditionalPanel("input.con_opt_k >= 5", sliderInput("upto_slider5", "Remove Samples within Cluster5 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0)),
                                    conditionalPanel("input.con_opt_k >= 6", sliderInput("upto_slider6", "Remove Samples within Cluster6 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0)),
                                    conditionalPanel("input.con_opt_k >= 7", sliderInput("upto_slider7", "Remove Samples within Cluster7 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0)),
                                    conditionalPanel("input.con_opt_k >= 8", sliderInput("upto_slider8", "Remove Samples within Cluster8 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0)),
                                    conditionalPanel("input.con_opt_k >= 9", sliderInput("upto_slider9", "Remove Samples within Cluster9 with Silhoutte width less than:", min = 0, max = 1.0, value = 0.0))
                             ),
                              conditionalPanel("input.sil_choice == 'Percentile'",
                                    #numericInput("sil_percen", "Remove Samples with Silhoutte width less than what percentile:", value = 10, min = 0, max = 50, step = 5)
                                    #conditionalPanel("input.sil_choice == 'Fixed value'",
                                    conditionalPanel("input.con_opt_k >= 1", sliderInput("sil_percen1", "Remove Samples within Cluster1 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 2", sliderInput("sil_percen2", "Remove Samples within Cluster2 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 3", sliderInput("sil_percen3", "Remove Samples within Cluster3 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 4", sliderInput("sil_percen4", "Remove Samples within Cluster4 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 5", sliderInput("sil_percen5", "Remove Samples within Cluster5 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 6", sliderInput("sil_percen6", "Remove Samples within Cluster6 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 7", sliderInput("sil_percen7", "Remove Samples within Cluster7 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 8", sliderInput("sil_percen8", "Remove Samples within Cluster8 with percentile silhoutte width less than:", min = 0, max = 50, value = 10)),
                                    conditionalPanel("input.con_opt_k >= 9", sliderInput("sil_percen9", "Remove Samples within Cluster9 with percentile silhoutte width less than:", min = 0, max = 50, value = 10))
                                    ),
                              conditionalPanel("input.sil_choice == 'Change point'",
                                    selectInput("changes","Changes Using:", c('mean','var','meanvar'),selected="meanvar", selectize = FALSE),
                                    selectInput("method","Method:", c('AMOC','PELT','BinSeg','SeqNeigh'),selected="BinSeg", selectize = FALSE),
                                    numericInput("max","Max Q Allowed:", min=1, max=60, step=5, value=60),
                                    conditionalPanel("input.con_opt_k >= 1", numericInput("sil_cp1", "Remove Samples within Cluster1 with silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 2", numericInput("sil_cp2", "Remove Samples within Cluster2 with silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 3", numericInput("sil_cp3", "Remove Samples within Cluster3 with silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 4", numericInput("sil_cp4", "Remove Samples within Cluster4 with Silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 5", numericInput("sil_cp5", "Remove Samples within Cluster5 with Silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 6", numericInput("sil_cp6", "Remove Samples within Cluster6 with Silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 7", numericInput("sil_cp7", "Remove Samples within Cluster7 with Silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 8", numericInput("sil_cp8", "Remove Samples within Cluster8 with Silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1)),
                                    conditionalPanel("input.con_opt_k >= 9", numericInput("sil_cp9", "Remove Samples within Cluster9 with Silhoutte width greater than change point:", value = 1, min = 0, max = 10, step = 1))
                                    
                              )
                              
                    ) ),
                 conditionalPanel("input.gw_conditionedPanels == 5",
                    conditionalPanel("input.cc_conditionedPanels == 1", 
                                     wellPanel(  
                                       ########## HeatMap Clustering options ##########
                                       h4("Heat Map Options"),
                                       selectInput("cc_norm", "Normalization Type",
                                                   c("Z-Score", "Modified Z-Score", "none")),
                                       selectInput("cc_norm2", "Normalize by:",
                                                   c("row", "col", "both"), selected = "both"),
                                       sliderInput("cc_inSlider", "Scale Range",
                                                   min = -10, max = 20, value = c(-2, 2)),
                                       sliderInput("cc_inSlider2", "Plot Margin dimensions",
                                                   min = 0, max = 20, value = c(14, 16)) ),
                                     wellPanel(  
                                       h4("Clustering Measures"),
                                       selectInput("cc_dist", "Distance Method",
                                                   c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson correlation"),
                                       selectInput("cc_hclust", "Agglomerative Linkage Method", 
                                                   c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "ward.D"),
                                       conditionalPanel("input.conditionedPanels==1 ", 
                                                        radioButtons("cc_clust_byrow", "Row dendrogram", inline = TRUE, c("TRUE", "FALSE")) ,
                                                        radioButtons("cc_clust_bycol", "Col dendrogram", inline = TRUE, c("TRUE", "FALSE")) ,
                                                        radioButtons("cc_dispRow", "Display Row labels?:", inline = TRUE,c("No", "Yes")),
                                                        conditionalPanel("input.cc_dispRow == 'Yes'", sliderInput("cc_size1", "If yes, Row Label font size", min = 0.01, max = 3, value = 0.5)),
                                                        radioButtons("cc_dispCol", "Display Col labels?:", inline = TRUE, c("No", "Yes")),
                                                        conditionalPanel("input.cc_dispCol == 'Yes'", sliderInput("cc_size2", "If yes, Col Label font size", min = 0.01, max = 3, value = 0.5) ))
                                     ), 
                                     wellPanel(
                                       h4("Heatmap colors"),
                                       colourInput("cc_low", "low", "green", returnName = TRUE, palette = "limited", showColour = "background"),
                                       colourInput("cc_mid", "mid", "black", returnName = TRUE, palette = "limited", showColour = "background"),
                                       colourInput("cc_high", "high", "red", returnName = TRUE, palette = "limited", showColour = "background") )
                    ),
                    conditionalPanel("input.cc_conditionedPanels ==2 | input.cc_conditionedPanels == 3",
                                      wellPanel(
                                        conditionalPanel("input.cc_conditionedPanels==2", 
                                                         sliderInput("cc_sizeClable", "Adjust Column Label font size", min = 0.01, max = 3, value = 0.8) ),
                                        conditionalPanel("input.cc_conditionedPanels==3", 
                                                         sliderInput("cc_sizeRlable", "Adjust Row Label font size", min = 0.01, max = 3, value = 0.42) ),
                                        
                                        # Column Dendrogram tab
                                        conditionalPanel("input.cc_conditionedPanels==2", 
                                                         radioButtons("cc_cutcolden", "Cut Col dendrogram?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE)) ),
                                        conditionalPanel("input.cc_conditionedPanels==2 & input.cc_cutcolden == 'TRUE'", 
                                                         numericInput("cc_cuttree", "Cut Col Dendrogram at:", 2)),
                                        
                                        # Row Dendrogram tab
                                        conditionalPanel("input.cc_conditionedPanels==3", 
                                                         radioButtons("cc_cutrowden", "Cut Row dendrogram?:", inline = TRUE, c("No" = FALSE, "Yes" = TRUE)) ),
                                        conditionalPanel("input.cc_conditionedPanels==3 & input.cc_cutrowden == 'TRUE'", 
                                                         numericInput("cc_cuttree2", "Cut Row Dendrogram at:", 2))
                                        
                                      )
                    )))
             )
              
             
    ),
    tabPanel("Combined Results Clustering (CrC) Analysis", 
             fluidRow(
               column(2,
                      wellPanel(
                        h3("Input file"),
                        conditionalPanel("input.cPanels1 == 1",
                                         selectInput("Exp_file",label= "Select the example Expression ds or upload your own with 'Load my own'", 
                                                     choices = c("Example TCGA BRCA Exp ds"="Exp_example", #"Example CoMMpass Exp ds"="Exp_example1", 
                                                                 "Load my own data" = "Exp_load_my_own")),
                                         conditionalPanel("input.Exp_file == 'Exp_load_my_own'",
                                                          fileInput('Exp_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.Exp_file == 'Exp_example'",
                                                          downloadButton('Exp_download', 'Download TCGA BRCA Expression data'))#,
                                         #conditionalPanel("input.Exp_file == 'Exp_example1'",
                                         #                  downloadButton('Exp_download1', 'Download CoMMpass Expression data'))
                                         ),
                        conditionalPanel("input.cPanels1 == 2",
                                         selectInput("Variant_file",label= "Select the example Meth/Variant ds or upload your own with 'Load my own'", 
                                                     choices = c("Example TCGA BRCA Meth ds"="Variant_example", #"Example CoMMpass Variant ds"="Variant_example1", 
                                                                 "Load my own data" = "Variant_load_my_own")),
                                         conditionalPanel("input.Variant_file == 'Variant_load_my_own'",
                                                          fileInput('Variant_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.Variant_file == 'Variant_example'",
                                                          downloadButton('Variant_download', 'Download TCGA BRCA Meth data'))#,
                                         #conditionalPanel("input.Variant_file == 'Variant_example1'",
                                         #                   downloadButton('Variant_download1', 'Download CoMMpass Variant data'))
                                         ),
                        conditionalPanel("input.cPanels1 == 3",
                                         selectInput("CNV_file",label= "Select the example CNV ds or upload your own with 'Load my own'", 
                                                     choices = c("Example TCGA BRCA CNV ds"= "CNV_example", #"Example CoMMpass CNV ds"="CNV_example1", 
                                                                 "Load my own data" = "CNV_load_my_own")),
                                         conditionalPanel("input.CNV_file == 'CNV_load_my_own'",
                                                          fileInput('CNV_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.CNV_file == 'CNV_example'",
                                                          downloadButton('CNV_download', 'Download TCGA BRCA CNV data'))#,
                                         #conditionalPanel("input.CNV_file == 'CNV_example1'",
                                         #                  downloadButton('CNV_download1', 'Download coMMpass CNV data'))
                                         ),
                        conditionalPanel("input.cPanels1 == 4",
                                         selectInput("coca_file",label= "Select the example CrC analysis ds or upload your own with 'Load my own'", 
                                                     choices = c("Example ds File"="coca_example", "Load my own data" = "coca_load_my_own")),
                                         conditionalPanel("input.coca_file == 'coca_load_my_own'",
                                                          fileInput('coca_file2', 'Choose file to upload (maximum size 500 MB).', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                                         conditionalPanel("input.coca_file == 'coca_example'",
                                                          downloadButton('coca_download', 'Download CrC analysis data')))
                        ),
                      conditionalPanel("input.cPanels1 == 2",
                       wellPanel(h3("Data type"),
                                 radioButtons("dt_choice", "Choose data used", choices = c("Methylation" = "M", "Variant" = "V"), selected = "M")
                     )),
                       wellPanel(h3("Clustering Options"),
                          conditionalPanel("input.cPanels1 == 1",      
                            selectInput("Exp_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson"),
                            selectInput("Exp_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"),
                            radioButtons("Exp_pItems", "No. of iterations", choices = c(50,100, 500, 1000), selected = 100),
                            br(),
                            actionButton("button6", "Run Analysis"),
                            h5("Hit button to update results after each change in input parameter(s)")),
                          conditionalPanel("input.cPanels1 == 2",      
                            selectInput("Variant_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson"), #"manhattan"
                            selectInput("Variant_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "average"), #"average"
                            radioButtons("Variant_pItems", "No. of iterations", choices = c(50,100, 500, 1000), selected = 100),
                            br(),
                            actionButton("button7", "Run Analysis"),
                            h5("Hit button to update results after each change in input parameter(s)")),
                          conditionalPanel("input.cPanels1 == 3",      
                            selectInput("CNV_dist","Distance Measure",choices = c("pearson", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "canberra"),
                            selectInput("CNV_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "mcquitty"),
                            radioButtons("CNV_pItems", "No. of iterations", choices = c(50, 100, 500, 1000), selected = 100),
                            br(),
                            actionButton("button8", "Run Analysis"),
                            h5("Hit button to update results after each change in input parameter(s)")),
                          conditionalPanel("input.cPanels1 == 4",
                            checkboxGroupInput("coca_platform","Perform CrC analysis using", choices = list("Expression"= "EXP","Methylation/Variant" = "PROP", "Copy Number" = "CNV"), selected = c("EXP", "PROP", "CNV")),
                            h6("Select atleast two platforms to run CrC analysis"),
                            selectInput("coca_dist","Distance Measure",choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean"),
                            selectInput("coca_hclust", "Clustering Method", c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "ward.D"),
                            radioButtons("coca_pItems", "No. of iterations", choices = c(50, 100, 500, 1000), selected = 100),
                            br(),
                            actionButton("button9", "Run Analysis"),
                            h5("Hit button to update results after each change in input parameter(s)"))
                          
                      ),
                      conditionalPanel("input.cPanels1 == 4",
                            wellPanel(h4("Clinical markers"),
                               selectInput("clinical",label= "Select a clinical feature from below or add your own using 'Load your own'.", 
                                 choices = c("Features available from TCGA BRCA"="available", #"Features available from coMMpass ds"="available1", 
                                             "Load my own clinical data" = "load_my_own_c")),
                               conditionalPanel("input.clinical == 'available'",
                                 downloadButton('download_clinical', 'Download TCGA BRCA clinical ds Example')),
                               #conditionalPanel("input.clinical == 'available1'",
                              #                    downloadButton('download_clinical1', 'Download coMMpass clinical ds Example')),
                               br(), br(),
                               checkboxInput("cli_feature","Add clinical features (Yes/No)", value = FALSE),
                              # conditionalPanel("input.cli_feature == 'TRUE'",
                               numericInput("cli_feature_no", "Number of features", min = 1, max = 5, step =1, value = 1), #), 
                               conditionalPanel("input.clinical == 'load_my_own_c'",
                                    fileInput('clinical_file', 'Choose file to upload (maximum size 500 MB).', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')))
                      )),
                     wellPanel(h3("Download Results"),
                          h4("Consensus Clustering"),
                          conditionalPanel("input.cPanels1 == 1", 
                                  downloadButton("Exp_cc_dl", label = "Download Expression clusters"),
                                  h6("Download may take a while. Once complete, the result pdf file will automatically open.")),
                          conditionalPanel("input.cPanels1 == 2", 
                                  downloadButton("Variant_cc_dl", label = "Download Meth/Variant clusters"), 
                                  h6("Download may take a while. Once complete, the result pdf file will automatically open.")),
                          conditionalPanel("input.cPanels1 == 3", 
                                  downloadButton("CNV_cc_dl", label = "Download CNV clusters"),
                                  h6("Download may take a while. Once complete, the result pdf file will automatically open.")),
                          conditionalPanel("input.cPanels1 == 4", 
                                 downloadButton("coca_cc_dl", label = "Download CrC analysis"),
                                 h6("Download may take a while. Once complete, the result pdf file will automatically open."),
                                 h4("Sample Clusters"),
                                 downloadButton("coca_sil", label = "Download CrC Sample clusters"),
                                 br(),br(),
                                 h4("CoC HeatMap"),
                                 downloadButton('download_coca_HM', 'Download CrC Analysis HM'),
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
                                  tabPanel("Expresssion", htmlOutput("com_text1"), withSpinner(plotOutput("Exp_cc", height = 500)), br(), htmlOutput("com_text2"), withSpinner(plotOutput("Exp_sil", height = 800)), 
                                      conditionalPanel(condition="input.button6 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=1), 
                                  tabPanel("Methylation/Variant", htmlOutput("com_text12"), withSpinner(plotOutput("Variant_cc", height = 500)), br(), htmlOutput("com_text22"), withSpinner(plotOutput("Variant_sil", height = 800)), 
                                           conditionalPanel(condition="input.button7 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=2), 
                                  tabPanel("Copy Number", htmlOutput("com_text13"), withSpinner(plotOutput("CNV_cc", height = 500)), br(), htmlOutput("com_text23"), withSpinner(plotOutput("CNV_sil", height = 800)),  
                                           conditionalPanel(condition="input.button8 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value =3),
                                  
                                  
                                  tabPanel("Cluster of Clusters (CoC) Analysis", htmlOutput("com_text14"), withSpinner(plotOutput("coca_cc", height = 500)), br(),htmlOutput("com_text3"), withSpinner(plotOutput("coca_heatmap", height = 1200)), 
                                           conditionalPanel("input.coca_file != 'coca_load_my_own'" ,
                                           htmlOutput("com_text31"), withSpinner(plotOutput("varmean_bxplot")), htmlOutput("com_text32"), htmlOutput("com_text33"), rHandsontableOutput("foo")), 
                                           conditionalPanel(condition="input.button9 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value = 4),
                                  id = "cPanels1")
                      
                ),
               column(2,
                      wellPanel(h3("Choose Optimal Number of clusters"),
                                conditionalPanel("input.cPanels1 == 1", numericInput("Exp_opt_k", "Optimal k is ",value = 2)),
                                conditionalPanel("input.cPanels1 == 2", numericInput("Variant_opt_k", "Optimal k is ",value = 2)),
                                conditionalPanel("input.cPanels1 == 3", numericInput("CNV_opt_k", "Optimal k is ",value = 2)),
                                conditionalPanel("input.cPanels1 == 4", numericInput("coca_opt_k", "Optimal k is ",value = 2))
                      ),
                      conditionalPanel("input.cPanels1 == 4", 
                                       br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       br(), br(), br(), br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(),
                                       #br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                      
                                       wellPanel(  
                                         h4("Clustering Measures"),
                                         selectInput("dist_4", "Distance Method",
                                                     c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean"),
                                         selectInput("hclust_4", "Agglomerative Linkage Method", 
                                                     c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "ward.D"),
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
                                                     min = 0, max = 20, value = c(13, 15)) ),
                                       
                                       br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                       br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                       br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                       br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                       br(),br(),br(),br(),br(),br(),br(),br(),
                                       
                                       conditionalPanel("input.coca_file != 'coca_load_my_own'" ,
                                       wellPanel(
                                         h4("Contingency Table(s) options:"),
                                         radioButtons("stratify_by", "Starify by:", inline = FALSE, c("Expression", "Methylation/Variant", "Copy Number"), selected = "Copy Number"))
                                       )
                                       
                                       #wellPanel(
                                      #   h4("Heat Map colors"), 
                                      #   colourInput("low_4", "low", "cyan", returnName = TRUE, palette = "limited", showColour = "background"),
                                      #   colourInput("mid_4", "mid", "khaki1", returnName = TRUE, palette = "limited", showColour = "background"),
                                      #   colourInput("high_4", "high", "chartreuse", returnName = TRUE, palette = "limited", showColour = "background") )
                      )))                    
                      
             ),
    tabPanel("Significance of Clusters (SoC) Analysis",
             #tags$iframe(src= "http://shinygispa.winship.emory.edu/CASH", width = 1900, height = 1000)
             fluidRow(
               column(2,
                      wellPanel(h3("Input Data to test significance of clusters "),
                                selectInput("Sig_file1",label= "Select an example dataset or upload your own with 'Load my own GW data.'", 
                                            choices = c("Example TCGA BRCA Expression Data"="Sig_Example1", #"Example CoMMpass  RNASeq gene Expression Data"="Sig_Example2", 
                                                        "Load my own data" = "load_my_own_Sig")),
                                conditionalPanel("input.Sig_file1 == 'load_my_own_Sig'",
                                                 fileInput('Sig_file2', 'Choose file to upload (maximum size 50 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))), 
                                conditionalPanel("input.Sig_file1 == 'Sig_Example1'", downloadButton('download_Sig_Ex1', 'Download Example DataSet')),
				                        conditionalPanel("input.Sig_file1 == 'Sig_Example2'", downloadButton('download_Sig_Ex2', 'Download Example DataSet'))

                      ),
                      wellPanel(h3("Input Settings"),
                                h4("Select first row and column where numeric data starts"),
                                numericInput("DataR", label = "Data row", min = 1, max = 15, value = 5),
                                numericInput("DataC", label = "Data column", min = 1, max = 10, value = 3),
                                #br(),
                                actionButton("button10", "Run Analysis"),
                                h5("Hit button to update results after each change in input parameter(s)")
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
                                  tabPanel("HeatMap", withSpinner(plotOutput("Sig_plot", width = 1300, height = 1300 )), 
                                           conditionalPanel(condition="input.button10 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=2), 
                                  tabPanel("Column Dendrogram", withSpinner(plotOutput("Sig_plot1", height= 600, width = 1500)), htmlOutput("Sig_display"), br(), DT::dataTableOutput("Sig_df"), htmlOutput("Sig_pv"), htmlOutput("Sig_pvalue"),  
                                           conditionalPanel(condition="input.button10 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value=3), 
                                  tabPanel("Row Dendrogram", withSpinner(plotOutput("Sig_plot2", height = 800, width = 1500)), htmlOutput("Sig_display2"), br(), DT::dataTableOutput("Sig_df2"), htmlOutput("Sig_pv2"), htmlOutput("Sig_pvalue2"),  
                                           conditionalPanel(condition="input.button10 & $('html').hasClass('shiny-busy')", tags$div("Loading...",id="loadmessage")), value =4),
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
                                                     c("pearson correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "pearson correlation"),
                                         selectInput("Sig_hclust", "Agglomerative Linkage Method", 
                                                     c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), selected = "ward.D"),
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
                                                          selectInput("Sig_file3", label= "Select a dataset or upload your own with 'Load my own data.'", choices = c("Example TCGA BRCA Exp Sampling Data" ="Sig_Exp.Example", #"Example CoMMPass Exp Sampling Data" ="Sig_Exp.Example1",  
                                                                                                                                                                      "Load my own sampling data" = "Sig_load_my_own_s_data"))),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_file3 == 'Sig_load_my_own_s_data'",
                                                          fileInput('Sig_file4', 'Choose file to upload to sample from to estimate significance of separation (Maximum size 100 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))) ,
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                                          numericInput("Sig_n", "Gene-Set size:", 605)),
                                         #conditionalPanel("input.cPanel2s==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                        #                  numericInput("Sig_n_iter", "No. of iterations for bootstrap:", 1000)),
                                         conditionalPanel("input.cPanels2==3 & input.Sig_cutcolden == 'TRUE' & input.Sig_pvalue_cal == 'TRUE'", 
                                                          numericInput("Sig_n_iter", "No. of bootstrap samples:", 1000)),
                                         
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
                                                          selectInput("Sig_file5", label= "Select a dataset or upload your own with 'Load my own data.'", choices = c("Example TCGA BRCA Exp Sampling Data" ="Sig_Exp.Example2", #"Example CoMMpass Exp Sampling Data" ="Sig_Exp.Example21", 
                                                                                                                                                                      "Load my own sampling data" = "load_my_own_s_data2"))),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_file5 == 'load_my_own_s_data2'",
                                                          fileInput('Sig_file6', 'Choose file to upload to sample from to estimate significance of separation (Maximum size 75 MB)', accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv'))) ,
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE'", 
                                                          numericInput("Sig_n2", "Sample-set size:", 1000)),
                                         conditionalPanel("input.cPanels2==4 & input.Sig_cutrowden == 'TRUE' & input.Sig_pvalue_cal2 == 'TRUE'", 
                                                          numericInput("Sig_n_iter2", "No. of bootstrap samples:", 1000)),
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

