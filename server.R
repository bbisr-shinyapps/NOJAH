options(shiny.maxRequestSize=1000*1024^2)
options(shiny.sanitize.errors = FALSE) 
#options(bitmapType='cairo')

source("functions/zClust.R")
source("functions/modzClust.R")
source("functions/consensus_clustering.R")
source("functions/silhouette_plot.R")
source("functions/coca.R")
source("functions/helper.R")
source("functions/colbars.R")
source("functions/plotMeans.R")
source("functions/bxplot.R")
source("functions/bootstrapfun.R")
source("functions/contingencyfun.R")
source("functions/pvaluefunc.R")
source("functions/col_color.R")
source("functions/row_color.R")
source("functions/silhouette_plot2.R")
source("functions/silhouette_plot3.R")
source("functions/colbars_gw.R")

function(input, output) {
  
  output$ReadMe <- renderUI({
    str0 <- paste("&emsp;")
    str00 <- paste("DATA INPUT")
    str1 <- paste("Data should be input as a .txt or .csv file. The first two rows of the data file have information about the patients/specimens and their response/subtype; all remaining rows have gene expression data, one row per gene.  In the case of Microarray gene expression data in which there are several probes corresponding to a single gene, a unique identifier would need to be created to separately identify each probe such as, 'Gene 1_p1', 'Gene1_p2' indicating Gene 1 has two probes.  The columns represent the different experimental samples. A maximum of up to 10 different sample groups and 6 different gene groups may be used with this tool.")
    str2 <- paste("Data Format")
    str3 <- paste("&emsp; 	1.	The first line of the file contains the gene identifier 'gene_id' (column 1), gene group classification 'Groups' (column 2) followed by the patient IDs e.g. TCGA.01.1A2B, one per column, starting at column 3. Column 1 gene identifier has to be labelled 'gene_id' and column 2 header should be labelled 'Groups' for using this tool. Other titles may cause the program to display errors. ")
    str4 <- paste("&emsp; 	2.	The second line of the file contains the patient response classification e.g. Fav/Unf for favorable outcome group vs the unfavorable outcome group or Normal/Tumor, etc., in alphabetical order, starting at column 3. The first two columns for this row should be blank.")
    str5 <- paste("&emsp; 	3.	The remaining lines contain gene expression measurements one line per gene, described in the format below.")
    str6 <- paste("&emsp;&emsp;  a) Column_1. This should contain the gene name, for the user's reference. Each row should be unique. For microarray data, for example, gene and probes can be combined into a single identifier using any delimitor except the '|'. Delimitors such as >,;:#%&(!)_+ are acceptable.")
    str7 <- paste("&emsp;&emsp;  b)   Column_ 2.  This should contain the gene group classification e.g. O/U for Over-expressed/Under-expressed or Hyper/Hypo for hypermethylated/hypomethylated in alphabetical order. If only one gene group, use any alphabet e.g. A or even na for each row instead.")
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
    s.1 <- c("YWHAE>210996_s_at", "na", 1.47, 2.18, 5.87, 9.12, 7.34, 1.56, 3.0, 7.77, 3.40, 1.56 )
    s.2 <- c("YWHAE>201020_at", "na", 1.98, 7.93, 2.76, 9.11, 8.46, 0.98, 5.98, 8.19, 8.91, 5.98)
    s.3 <- c("YWHAH>33323_r_at", "na", 8.02, 8.00, 2.17, 10.12, 8.76, 9.76, 3.76, 0.02, 3.67, 7.94)
    s.4 <- c("YWHAB>208743_s_at", "na", 2.75, 5.99, 3.19, 11.86, 6.54, 8.17, 2.00, 0.99, 2.00, 1.17)
    s.5 <- c("YWHAQ>213699_s_at", "na", 9.35, 8.96, 6.67, 8.33, 3.98, 7.11, 1.67, 1.01, 5.18, 8.17)
    
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
    filename= function() {paste('TCGA.BRCA.Expression.csv')}, 
    content = function(file) {
      d <- read.csv("data/BRCA_Expression_gene.normalized_log2T_RNAseq_Tumor_Normal_150samples.csv")
      write.csv(d, file, row.names = FALSE) }
  )
  
  
  input_gw_data <- reactive({
    if(input$gw_file1 == 'GW_Example1'){
      d <- readRDS("data/CoMMpassIA9.GW.Expression.data_truncated.rds")
    } else if(input$gw_file1 == 'GW_Example2'){
      d <- read.csv("data/BRCA_Expression_gene.normalized_log2T_RNAseq_Tumor_Normal_150samples.csv", header = TRUE, sep = ",", stringsAsFactors = F)
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
   
    Dataset <- data.frame(d)
    return(Dataset)
  })
  
  output$gw_dend <- renderPlot({
    data <- input_gw_data()
    data2 <- data.frame(data[-1,], stringsAsFactors = F)
    
    rownames(data2) = paste(data2$gene_id, data2$Group, sep = "|")
    data2 <- data.frame(data2[, c(-1, -2)], stringsAsFactors = F)
    colnames(data2) = paste(names(data2), as.vector(unlist(data[1,c(-1, -2)])), sep = "||")
    data2 <- data.frame(as.matrix(data2), stringsAsFactors = F)
    data2 <- data.frame(apply(data2, 2, function(x) as.numeric(as.character(x))))
    rownames(data2) = paste(data$gene_id[-1], data$Group[-1], sep = "|")
    
    data3 <- data.matrix(data2)
    getVar <- apply(data3, 1, var)
    data4 <- data3[getVar != 0, ]
    
    col.groups <- as.character(as.vector(data[1,]))
    col.groups <- col.groups[c(-1, -2)] # calculate no. of column groups
    col.groups.name <- names(table(col.groups))
    number.col.groups <- length(col.groups.name)
    
    cc1 = vector()
    
    ## Set color palette
    col1 <- colorRampPalette(c(input$low,input$mid,input$high))(299)
    colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
    
    cc1 <- col_color(col1 = col1, col.groups= col.groups, number.col.groups= number.col.groups, col.groups.name= col.groups.name)
    
    z <- list()
    
    if(input$GW_norm == "Z-Score") {
      z <- zClust(data4, scale =input$GW_norm2, zlim=c(input$GW_inSlider[1],input$GW_inSlider[2]))
    } else if (input$GW_norm == "Modified Z-Score") { 
      z <- modzClust(data4, scale =input$GW_norm2, zlim=c(input$GW_inSlider[1],input$GW_inSlider[2]))
    } else if(input$GW_norm == "none") {
      z[[1]] <- as.matrix(data4)
    }
    
    # scaling and normalization
    if(input$dist == "pearson correlation") {
     dist.func =  function(x) as.dist((1-cor(t(x))))
     hc = hclust(dist.func(t(z[[1]])), method = input$hclust)
    } else {
     hc = hclust(dist(t(z[[1]]), method = input$dist), method = input$hclust)
    }
   
    par(cex=input$gw_dend_size)
    dend <- as.dendrogram(hc)
    d <- data.frame(v1 = hc$order, v2=1:length(hc$order))
    m <- data.frame(v3 = 1:length(cc1), v4 = cc1)
    
    colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
    colbar <- colbar[,2]
    labels_colors(dend) <- as.character(colbar)
    plot(dend)
    
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
  
  gw_data <- reactive ({
    if(!is.null(input_gw_data()))
    {
      data <- input_gw_data()
      data2 <- data.frame(data[-1,], stringsAsFactors = F)
      
      rownames(data2) = paste(data2$gene_id, data2$Group, sep = "|")
      data2 <- data.frame(data2[, c(-1, -2)], stringsAsFactors = F)
      colnames(data2) = paste(names(data2), as.vector(unlist(data[1,c(-1, -2)])), sep = "||")
      data2 <- data.frame(as.matrix(data2), stringsAsFactors = F)
      data2 <- data.frame(apply(data2, 2, function(x) as.numeric(as.character(x))))
      rownames(data2) = paste(data$gene_id[-1], data$Group[-1], sep = "|")
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
    } else
      return(NULL)
  })
  
  output$geneSelector <- renderUI({
    selectizeInput(inputId = "Genes", "Choose Option:", as.list(getOSgenes()),options=list(maxOptions=getOSgenes())) 
  })
  
  output$dropdowngene <- renderText({ 
    paste("You have selected gene", input$Genes)
  })
  
  extracted_data <- reactive ({
    
  if(!is.null(input_gw_data())) {
    data2 <- gw_data()
    data <- input_gw_data()
    n= ncol(data)-2
   
    if(is.null(input$gw_subset)) {
      return(NULL)
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
  } else 
    return(NULL)
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
    plot_ly(data, y = data$var, type = 'box', name = 'Var') %>%
      add_trace(y = data$mad, name = 'MAD')  %>%
      add_trace(y = data$IQR, name = 'IQR') %>%
      layout(yaxis = y)
    #pp
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
    if (input$gw_file_1 == 'GW_Example_1') {
        data3 <- extracted_data()
        colheaders <- sub("^.*\\.","", colnames(data3))
        names(data3) <- gsub("\\.{2}.*", "", colnames(data3))
        data4 <- rbind(colheaders, data3)
        Groups <- sub("^.*\\|","", rownames(data4))
        Groups[1] <- " "
        gene_id <- sub("\\|.*$", "", rownames(data4)) 
        gene_id[1] <- " "
      
        data5 <- cbind(gene_id, Groups, data4)
        rownames(data5)[1] <- ""
      } else if(input$gw_file_1 == "load_my_own_gw_subset_1") {
        inFile_1 <- input$gw_file_2
        if (is.null(inFile_1))
          return(NULL)
        else if(grepl(".csv", inFile_1[1])) { data5 = read.csv(as.character(inFile_1$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
        else if(grepl(".txt", inFile_1[1])) { data5 = read.table(as.character(inFile_1$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
        # else if(grepl(".rds", inFile_1[1])) { data5 = read.rds(as.character(inFile_1$datapath)) }
      } 
      return(data.frame(data5))
  
   })
  
  output$downloadSubset <- downloadHandler(
    filename= function() {paste(input$fname_subset, Sys.time(),'.csv', sep='')}, 
    content = function(file) {
      d <- extracted_data2()
      write.csv(d, file, row.names = FALSE) }
  )
  
  
  
  mv_hm_data = reactive({
    if(is.null(extracted_data2()))
    {
    data <- NULL
    } else {
      data5 <- extracted_data2()
  
      # sort columns based on colnames
      data <- data5[,order(data5[1, ])]
      data <- data[order(data[,2]),]
      
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
      
      cc1 = vector()
      cc2 = vector()
      
      ## Set color palette
      col1 <- colorRampPalette(c(input$low,input$mid,input$high))(299)
      colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
      
      ### Color vector for columns
      cc1 <- col_color(col1 = col1, col.groups= col.groups, number.col.groups= number.col.groups, col.groups.name= col.groups.name)
     
      ### Color vector for rows
      cc2 <- row_color(col1 = col1, row.groups= row.groups, number.row.groups= number.row.groups, row.groups.name= row.groups.name)
      
      ############# HEATPLOT2 EQUIVALENT HEATMAP2 CLUSTERING ###############
      z <- list()
      
      #data <- as.numeric(data)
      if(input$norm == "Z-Score") {
        z <- zClust(data, scale =input$norm2, zlim=c(input$inSlider[1],input$inSlider[2]))
      } else if (input$norm == "Modified Z-Score") { 
        z <- modzClust(data, scale =input$norm2, zlim=c(input$inSlider[1],input$inSlider[2]))
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
      
     hmdata <- list(data = data, z = list(z), hm = list(hm), col1=col1, cc1=cc1, cc2 = cc2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name)
    }
  })
  
  tgplot <- function(z, col1, cc1, cc2, number.col.groups, number.row.groups, row.groups.name , col.groups.name ) {  
                    
    if(input$dist == "pearson correlation") {
      if(input$dispRow == "No" & input$dispCol=='No') {
        hm <- heatmap.2(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
      } else if(input$dispRow == "No" & input$dispCol=='Yes' ) {
        hm <- heatmap.2(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
      } else if(input$dispRow == "Yes" & input$dispCol=='No') {
        hm <- heatmap.2(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
      } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
        hm <- heatmap.2(z[[1]], Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
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
        hm <- heatmap.2(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
      } else if(input$dispRow == "No" & input$dispCol=='Yes') {
        hm <- heatmap.2(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
      } else if(input$dispRow == "Yes" & input$dispCol=='No') {
        hm <- heatmap.2(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
      } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
        hm <- heatmap.2(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2)
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
  }
  

  
  output$GW_subset_heatmap <- renderPlot({
      if(is.null(mv_hm_data())) {
        hm <- NULL
      } else {    
        tgplot(z= mv_hm_data()$z[[1]], col1 = mv_hm_data()$col1, cc1 = mv_hm_data()$cc1, cc2 = mv_hm_data()$cc2, 
               number.col.groups =  mv_hm_data()$number.col.groups, number.row.groups = mv_hm_data()$number.row.groups, 
               row.groups.name = mv_hm_data()$row.groups.name, col.groups.name = mv_hm_data()$col.groups.name)
      }
  })
  
  coldendo <- reactive({
        par(cex = input$sizeClable)
        dend1 <- as.dendrogram(mv_hm_data()$hm[[1]]$colDendrogram)
        d <- data.frame(v1 =mv_hm_data()$hm[[1]]$colInd, v2=1:length(mv_hm_data()$hm[[1]]$colInd))
        m <- data.frame(v3 = 1:length(mv_hm_data()$cc1), v4 = mv_hm_data()$cc1)
        
        colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
        colbar2 <- colbar[,2]
        labels_colors(dend1) <- as.character(colbar2)
        plot(dend1)
        
        if(mv_hm_data()$number.col.groups==1) {
          legend("topright", legend = paste(mv_hm_data()$col.groups.name), col = "darkblue", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==2) {
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==3) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==4) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==5) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4], mv_hm_data()$col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==6) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4], mv_hm_data()$col.groups.name[5], mv_hm_data()$col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==7) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4], mv_hm_data()$col.groups.name[5], mv_hm_data()$col.groups.name[6], mv_hm_data()$col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==8) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4], mv_hm_data()$col.groups.name[5], mv_hm_data()$col.groups.name[6], mv_hm_data()$col.groups.name[7], mv_hm_data()$col.groups.name[8])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(mv_hm_data()$number.col.groups==9) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4], mv_hm_data()$col.groups.name[5], mv_hm_data()$col.groups.name[6], mv_hm_data()$col.groups.name[7], mv_hm_data()$col.groups.name[8], mv_hm_data()$col.groups.name[9])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable9)
        } else if(mv_hm_data()$number.col.groups==10) {  
          legend("topright", legend = paste(c(mv_hm_data()$col.groups.name[1], mv_hm_data()$col.groups.name[2], mv_hm_data()$col.groups.name[3], mv_hm_data()$col.groups.name[4], mv_hm_data()$col.groups.name[5], mv_hm_data()$col.groups.name[6], mv_hm_data()$col.groups.name[7], mv_hm_data()$col.groups.name[8], mv_hm_data()$col.groups.name[9], mv_hm_data()$col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        }
      })
      
      output$plot1 <- renderPlot({
        coldendo()
      })
      
      
      colDen <- reactive({
        if(input$cutcolden == 'TRUE') {
          cuttable <- as.data.frame(cutree(as.hclust(mv_hm_data()$hm[[1]]$colDendrogram), k=as.numeric(input$cuttree))[as.hclust(mv_hm_data()$hm[[1]]$colDendrogram)$order])
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
        dend2 <- as.dendrogram(mv_hm_data()$hm[[1]]$rowDendrogram)
        dd <- data.frame(v1 =rev(mv_hm_data()$hm[[1]]$rowInd), v2=1:length(mv_hm_data()$hm[[1]]$rowInd))
        mm <- data.frame(v3 = 1:length(mv_hm_data()$cc2), v4 = mv_hm_data()$cc2)
        
        colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
        colbar2 <- colbar2[,2]
        labels_colors(dend2) <- rev(as.character(colbar2))
        plot(dend2, horiz = T)
        if(mv_hm_data()$number.row.groups==1) {
          legend("topright", legend = paste(mv_hm_data()$row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(mv_hm_data()$number.row.groups==2) {
          legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(mv_hm_data()$number.row.groups==3) {  
          legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(mv_hm_data()$number.row.groups==4) {  
          legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3], mv_hm_data()$row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(mv_hm_data()$number.row.groups==5) {  
          legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3], mv_hm_data()$row.groups.name[4], mv_hm_data()$row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } else if(mv_hm_data()$number.row.groups==6) {  
          legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3], mv_hm_data()$row.groups.name[4], mv_hm_data()$row.groups.name[5], mv_hm_data()$row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
        } 
      })
      
      output$plot2 <- renderPlot({
        par(cex= input$sizeRlable)
        rowdendo()
      })
      
      rowDen <- reactive({
        if(input$cutrowden == 'TRUE') {
          cuttable2 <- as.data.frame(cutree(as.hclust(mv_hm_data()$hm[[1]]$rowDendrogram), k=as.integer(input$cuttree2))[as.hclust(mv_hm_data()$hm[[1]]$rowDendrogram)$order])
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
                                 c("Scale", ifelse(input$norm == "none", paste(as.integer(min(mv_hm_data()$data)), as.integer(max(mv_hm_data()$data)), sep = ":"), paste(input$inSlider[1], input$inSlider[2], sep=":"))),
                                 c("HeatMap colors", paste(input$low, input$mid, input$high, sep="-")))
          names(df)[1] <- "Parameters"
          names(df)[2] <- "Value Selected"
          grid.table(df, rows= NULL)
          
          #eval(hm$call)
          #eval(mv_hm_data()$hm[[1]]$call) # call the heatmap here
          #eval(hm$call)
          #tgplot()
          tgplot(z= mv_hm_data()$z[[1]], col1 = mv_hm_data()$col1, cc1 = mv_hm_data()$cc1, cc2 = mv_hm_data()$cc2, 
                 number.col.groups =  mv_hm_data()$number.col.groups, number.row.groups = mv_hm_data()$number.row.groups, 
                 row.groups.name = mv_hm_data()$row.groups.name, col.groups.name = mv_hm_data()$col.groups.name)
          
          par(cex = 0.6*input$sizeClable)
          coldendo()
          
          #par(cex = input$sizeRlable)
          rowdendo()
          
          dev.off()
          file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
        })
      
      
    
    consen <- reactive({
      if (input$gw_file3 == "GW_Example3")
        {
          mv_data <- mv_hm_data()$data
         
      } else if(input$gw_file3 == "load_my_own_gw_subset")
      {
        inFile4 <- input$gw_file4
        if (is.null(inFile4))
          return(NULL)
        else if(grepl(".csv", inFile4[1])) { mv_data = read.csv(as.character(inFile4$datapath), header = TRUE, sep = ",", stringsAsFactors = F, row.names = 1) }
        else if(grepl(".txt", inFile4[1])) { mv_data = read.table(as.character(inFile4$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, row.names = 1) }
       # else if(grepl(".rds", inFile4[1])) { mv_data = read.rds(as.character(inFile4$datapath)) }
        
      } else
        return(NULL)
      
      mv_data2 <- as.matrix(mv_data)
        
      con <- consensus_clustering(dinput=mv_data2, mK=10, rep=as.integer(input$con_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$con_dist, iL=input$con_hclust, fL=input$con_hclust)
      
      
      return(list(output= con[["output"]][[as.integer(input$con_opt_k)]]$consensusClass, data= con[["data"]], distance=con[["distance"]]))
     #  return(as.data.frame(mv_data2))
    })
    
    output$download_GW_Ex3 <- downloadHandler(
      filename= function() {paste('Pre_existing subset for consensus.csv')}, 
      content = function(file) {
        d <-  mv_hm_data()$data
        write.csv(d, file, row.names = TRUE) }
    )
      
   output$gw_cc <- renderPlot ({
     par(mfrow= c(1, 3))
     consen()
     
  })
   
   sil_data <- reactive({
     # if(!is.null(consen())){
     if(input$gw_file5 == "GW_Example5") {
    
     data = consen()$data
     sample = colnames(consen()$data)
     
     colors = mv_hm_data()$cc1
     df = cbind.data.frame(sample, colors, stringsAsFactors = F)
     colnames(df) <- c("Sample","colors")
     
     sil = silhouette_plot2(data_use= consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], upto_width= as.numeric(input$upto_slider), cols = df)
     
     colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
     par(mfrow = c(1, 2))
     plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
     plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
     
     return(list(core = sil$core.samples, clust = sil$sk3, order = consen()[["output"]]))
     
     } else if(input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3"){
       inFile6 <- input$gw_file6
       if (is.null(inFile6))
         return(NULL)
       else if(grepl(".csv", inFile6[1])) { mv_data = read.csv(as.character(inFile6$datapath), header = TRUE, sep = ",", stringsAsFactors = F, row.names = 1) }
       else if(grepl(".txt", inFile6[1])) { mv_data = read.table(as.character(inFile6$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, row.names = 1) }
       # else if(grepl(".rds", inFile6[1])) { mv_data = read.rds(as.character(inFile6$datapath)) }
       
       inFile8 <- input$gw_file8
       if (is.null(inFile8))
         return(NULL)
       else if(grepl(".csv", inFile8[1])) { c_order = read.csv(as.character(inFile8$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
       else if(grepl(".txt", inFile8[1])) { c_order = read.table(as.character(inFile8$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
       # else if(grepl(".rds", inFile8[1])) { c_order = read.rds(as.character(inFile8$datapath)) }
       
       mv_data2 <- as.data.frame(mv_data)
       c_order2 <- as.data.frame(c_order)
       
       k= unique(c_order2$Cluster)[length(unique(c_order2$Cluster))]
       
       sil = silhouette_plot3(data_use= mv_data2, opt_k= k, res= c_order2, dist = input$sil_dist, upto_width= as.numeric(input$upto_slider), cols = c("orange", "darkblue", "black", "maroon", "violet", "plum2")[1:k])
       
       colors = c("cyan", "khaki1", "pink", "plum3", "purple", "darkgreen", "hotpink","brown", "darkorchid2",  "maroon")
       par(mfrow = c(1, 2))
       plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = c("orange", "darkblue", "black", "maroon", "violet", "plum2")[1:k])
       plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = c("orange", "darkblue", "black", "maroon", "violet", "plum2")[1:k]) #sil[["sk3.col"]]
       
       return(list(core = sil$core.samples, clust = sil$sk3, data = mv_data2, c_order = c_order2))
     }
       
     #colors = c("cyan", "khaki1", "pink", "plum3", "purple", "darkgreen", "hotpink","brown", "darkorchid2",  "maroon")
     #par(mfrow = c(1, 2))
     #plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
     #plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
     
     #return(list(core = sil$core.samples, clust = sil$sk3, order = consen()[["output"]]))
   })
   
   output$download_GW_Ex5 <- downloadHandler(
     filename= function() {paste('Core samples from Silhouette.csv')}, 
     content = function(file) {
       d <-  mv_hm_data()$data 
       write.csv(d, file, row.names = TRUE) }
   )
   
   output$download_GW_Ex6 <- downloadHandler(
     filename= function() {paste('Core samples from Silhouette.csv')}, 
     content = function(file) {
       d <-  core_mv_hm_data()$data_use
       d2 = t(d[1, c(-1, -2)])
       write.csv(d2, file, row.names = TRUE) }
   )
   
   output$download_GW_Ex7 <- downloadHandler(
     filename= function() {paste('Sample clusters.csv')}, 
     content = function(file) {
       d <-  sil_data()$order
       write.csv(d, file, row.names = TRUE) }
   )
    
    output$gw_core_sam <- renderPlot ({
    #  if(!is.null(sil_data()))
     # {
        #sil_data()
        if(input$gw_file5 == "GW_Example5") {
          
          check_consen <<- consen()
          
          data = consen()$data
          sample = colnames(consen()$data)
          check_sample <<- sample
          
          colors = mv_hm_data()$cc1
          check_colors <<- colors
          df = cbind.data.frame(sample, colors, stringsAsFactors = F)
          colnames(df) <- c("Sample","colors")
          
          sil = silhouette_plot2(data_use= consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], upto_width= as.numeric(input$upto_slider), cols = df)
          
          colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "lightpink","steelblue", "darkorchid2",  "yellowgreen", "violetred")
          par(mfrow = c(1, 2))
          plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
          plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
          
        } else if(input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3"){
          inFile6 <- input$gw_file6
          if (is.null(inFile6))
            return(NULL)
          else if(grepl(".csv", inFile6[1])) { mv_data = read.csv(as.character(inFile6$datapath), header = TRUE, sep = ",", stringsAsFactors = F, row.names = 1) }
          else if(grepl(".txt", inFile6[1])) { mv_data = read.table(as.character(inFile6$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, row.names = 1) }
          # else if(grepl(".rds", inFile6[1])) { mv_data = read.rds(as.character(inFile6$datapath)) }
          
          inFile8 <- input$gw_file8
          if (is.null(inFile8))
            return(NULL)
          else if(grepl(".csv", inFile8[1])) { c_order = read.csv(as.character(inFile8$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
          else if(grepl(".txt", inFile8[1])) { c_order = read.table(as.character(inFile8$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
          # else if(grepl(".rds", inFile8[1])) { c_order = read.rds(as.character(inFile8$datapath)) }
          
          mv_data2 <- as.data.frame(mv_data)
          c_order2 <- as.data.frame(c_order)
         
          k= unique(c_order2$Cluster)[length(unique(c_order2$Cluster))]
         
          sil = silhouette_plot3(data_use= mv_data2, opt_k= k, res= c_order2, dist = input$sil_dist, upto_width= as.numeric(input$upto_slider), cols = c("orange", "darkblue", "black", "maroon", "violet", "plum2")[1:k])
          
          colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "lightpink","steelblue", "darkorchid2",  "yellowgreen", "violetred")
          par(mfrow = c(1, 2))
          plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = c("orange", "darkblue", "black", "maroon", "violet", "plum2")[1:k])
          plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = c("orange", "darkblue", "black", "maroon", "violet", "plum2")[1:k]) #sil[["sk3.col"]]
          
         
        }
        
      #} else
      #  return(NULL)
    })
    
    output$con_dl <- downloadHandler(
      filename <- function(){
        pdf_file_con <<- paste(input$con_dist, input$con_hclust, sep = "_")
        paste("ConsensusClustring_Results_", pdf_file_con,'.pdf', sep='')
      },
      content <- function(file) {
        pdf(file=paste(pdf_file_con,".pdf",sep=""))
        mv_data <- mv_hm_data()$data
        mv_data2 <- as.matrix(mv_data)
        con <- consensus_clustering(dinput=mv_data2, mK=10, rep=as.integer(input$con_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$con_dist, iL=input$con_hclust, fL=input$con_hclust)
        #silhouette_plot2(data_use= consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], upto_width= as.numeric(input$upto_slider))
        dev.off()
        file.copy(paste(pdf_file_con,'.pdf', sep='') ,file, overwrite=TRUE)
      }
    )
    
    core_mv_hm_data = reactive({
      if(is.null(consen()))
      {
        data <- NULL
      } else {
        
        if(input$gw_file5 == "GW_Example5") {
        data5 <- extracted_data2()
       
        samples_to_include <- sil_data()$core
        check_samples_to_include <<- samples_to_include

        cc_clust <- as.data.frame(consen()$output)
        colnames(cc_clust) <- "Cluster"
        cc_clust$Sample <- rownames(cc_clust)
        check_cc_clust <<- cc_clust
       
        # sort columns based on colnames
        data1 <- data5[,order(data5[1, ])]
        data1 <- data1[order(data1[,2]),]
        data_use = data1[, colnames(data1) %in% samples_to_include]
        rownames(data_use) <- paste(data1[,1], data1[,2], sep = "|")
        
        cc_clust2 <- cc_clust[cc_clust$Sample %in% samples_to_include, ]
        cc_clust3 <- cc_clust2[match(colnames(data_use), cc_clust2$Sample),] ######
        df2 <- t(cc_clust3)
        
        data2 <- rbind(df2[-2,], data_use)
        rownames(data2)[2] <- "Groups"
        
        data3 <- data2[, order(data2[1,])]
       
        gene_id <- sub('\\|.*', '', rownames(data2))
        Groups <- sub('.*\\|', '', rownames(data2))
        
        data3 <- data.frame(gene_id, Groups, data3)
        data3[1:2,c(1:2)] <-  " "
        
        #data3 <- as.data.frame(data3[,order(data3[1,3:ncol(data3)])])
        check_data3 <<- data3
        
        ### gene names, column name as gene
        gene <- as.character(data3$gene_id)
        gene <- gene[(4-1):nrow(data3)] #gene[(6-1):nrow(data3)] 
        row.groups <- as.character(as.vector(data3[(4-1): nrow(data3),2])) #as.character(as.vector(data3[(6-1):nrow(data3),2])) 
        row.groups.name <- names(table(row.groups))
        number.row.groups <- length(row.groups.name)
        
        
        ### column groups
        col.groups <- as.character(as.vector(data3[1,]))
        #col.groups <- col.groups[c(-1, -2)] # calculate no. of column groups
        col.groups <-  col.groups[3:ncol(data3)] #col.groups[4:ncol(data3)] 
        col.groups.name <- names(table(col.groups))
        number.col.groups <- length(col.groups.name)
        
        ## additional info on columns and rows
        colbar_data <-  as.matrix(t(data3[1:2, 3:ncol(data3)])) #t(data[1:(6-2), c(1, 3:ncol(data))]) 
        #colnames(colbar_data) <- as.character(unlist(colbar_data[1,]))
        #colbar_data <- as.matrix(colbar_data[c(-1, -2),])
        colnames(colbar_data)[1] <- "Groups"
        n.colbar_data <- ncol(colbar_data)
        
        # rowbar_data <- as.matrix(data[(4-1):nrow(data), 1:(4-1)]) #as.matrix(data[5:nrow(data), 1:3])
        #  rownames(rowbar_data) <- rowbar_data[,1]
        rowbar_data <- as.matrix(data3[c(-1,-2),2])
        colnames(rowbar_data)[1] <- "Groups"
        n.rowbar_data <- ncol(rowbar_data)
        
        data <- data3[(4-1):nrow(data3), 3:ncol(data3)] #data[5:nrow(data), 4:ncol(data)]
        rownames(data) <- gene
        data<-data[complete.cases(data),]
        data <- data.matrix(data)
        
        
        
        } else if(input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3") {
          data5 = sil_data()$data
          samples_to_include <- sil_data()$core
          cc_clust = sil_data()$c_order
          
        # sort columns based on colnames
        data1 <- data5[,order(data5[1, ])]
        data1 <- data1[order(data1[,2]),]
        data_use = data1[, colnames(data1) %in% samples_to_include]
        
        
        # sort columns based on colnames
        data1 <- data5[,order(data5[1, ])]
        data1 <- data1[order(data1[,2]),]
        data_use = data1[, colnames(data1) %in% samples_to_include]
        
        cc_clust2 <- cc_clust[cc_clust$Sample %in% samples_to_include, ]
        cc_clust3 <- cc_clust2[match(colnames(data_use), cc_clust2$Sample),] ######
        df2 <- t(cc_clust3)
        
        data2 <- rbind(df2[-1,], data_use)
        #rownames(data2)[2] <- "Groups"
        
        
        gene_id <- sub('\\.|*', '', rownames(data2))
        Groups <- rep("A", length(gene_id))
        
        data3 <- data.frame(gene_id, Groups, data2)
        data3[1,c(1:2)] <-  " "
       
        ### gene names, column name as gene
        gene <- as.character(data3$gene_id)
        gene <- gene[(3-1):nrow(data3)] #gene[(6-1):nrow(data3)] 
        row.groups <- as.character(as.vector(data3[(3-1): nrow(data3),2])) #as.character(as.vector(data3[(6-1):nrow(data3),2])) 
        row.groups.name <- names(table(row.groups))
        number.row.groups <- length(row.groups.name)
      
        
        ### column groups
        col.groups <- as.character(as.vector(data3[1,]))
        #col.groups <- col.groups[c(-1, -2)] # calculate no. of column groups
        col.groups <-  col.groups[3:ncol(data3)] #col.groups[4:ncol(data3)] 
        col.groups.name <- names(table(col.groups))
        number.col.groups <- length(col.groups.name)
        
        ## additional info on columns and rows
        colbar_data <-  as.matrix(t(data3[1, 3:ncol(data3)])) #t(data[1:(6-2), c(1, 3:ncol(data))]) 
        #colnames(colbar_data) <- as.character(unlist(colbar_data[1,]))
        #colbar_data <- as.matrix(colbar_data[c(-1, -2),])
        colnames(colbar_data)[1] <- "Groups"
        n.colbar_data <- ncol(colbar_data)
        
       # rowbar_data <- as.matrix(data[(4-1):nrow(data), 1:(4-1)]) #as.matrix(data[5:nrow(data), 1:3])
      #  rownames(rowbar_data) <- rowbar_data[,1]
        rowbar_data <- as.matrix(data3[c(-1),2])
        colnames(rowbar_data)[1] <- "Groups"
        n.rowbar_data <- ncol(rowbar_data)
        
        data <- data3[(3-1):nrow(data3), 3:ncol(data3)] #data[5:nrow(data), 4:ncol(data)]
        rownames(data) <- gene
        data<-data[complete.cases(data),]
        data <- data.matrix(data)
      }
        
        cc1 = vector()
        cc2 = vector()
        
        ## Set color palette
        col1 <- colorRampPalette(c(input$cc_low,input$cc_mid,input$cc_high))(299)
        colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
        
        check_number.row.groups = number.row.groups
        
        ### Color vector for columns
        if(number.col.groups==1) { 
          cell <- c(rep(col.groups.name, number.col.groups))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==2) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'                                                                                                                                   
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1' 
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==3) {
          #colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine")
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          #colbars2 <- cbind(cc1, colbars_gw(df2 = colbar_data))
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==4) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==5) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'aquamarine'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==6) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'aquamarine'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'coral'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==7) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'aquamarine'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'coral'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'lightpink'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==8) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'aquamarine'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'coral'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'lightpink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'steelblue'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==9) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'aquamarine'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'coral'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'lightpink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'steelblue'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'yellowgreen'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[4]]
        } else if(number.col.groups==10) {
          cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]]), rep(col.groups.name[10], table(col.groups)[[10]])) 
          cc1 <- rep(col1[50], length(cell))
          cc1[1:table(col.groups)[[1]]] <- 'cyan'
          cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'khaki1'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'pink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'plum3'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'aquamarine'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'coral'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'lightpink'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- 'steelblue'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- 'yellowgreen'
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+table(col.groups)[[9]]+1:table(col.groups)[[10]]] <- 'violetred'
          colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw(df2 = colbar_data)[[1]]))
          colnames(colbars2)[1] <- "CC"
          number.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[2]]
          names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[3]]
          colors_used <- if(n.colbar_data ==1) NULL else  colbars_gw(df2 = colbar_data)[[4]]
        }
        
        ### Color vector for rows
        
        if(number.row.groups==1) { 
          cell2 <- c(rep(row.groups.name, number.row.groups))
          cc2 <- rep(col1[50], length(cell2))
          cc2[1:table(row.groups)[[1]]] <- 'black'
          rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw(df2 = rowbar_data, colscheme = 1)[[1]])))
          number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[2]]
          names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[3]]
          colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[4]]
          check_rowbars2 <<- rowbars2
        } else if(number.row.groups==2) {
          cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]))
          cc2 <- rep(col1[50], length(cell2))
          cc2[1:table(row.groups)[[1]]] <- 'black'
          cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
          rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw(df2 = rowbar_data, colscheme = 1)[[1]])))
          number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[2]]
          names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[3]]
          colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[4]]
          
        } else if(number.row.groups==3) {
          cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]])) 
          cc2 <- rep(col1[50], length(cell2))
          cc2[1:table(row.groups)[[1]]] <- 'black'
          cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
          rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw(df2 = rowbar_data, colscheme = 1)[[1]])))
          number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[2]]
          names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[3]]
          colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[4]]
          
        } else if(number.row.groups==4) {
          cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]])) 
          cc2 <- rep(col1[50], length(cell2))
          cc2[1:table(row.groups)[[1]]] <- 'black'
          cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
          rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw(df2 = rowbar_data, colscheme = 1)[[1]])))
          number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[2]]
          names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[3]]
          colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[4]]
        } else if(number.row.groups==5) {
          cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]])) 
          cc2 <- rep(col1[50], length(cell2))
          cc2[1:table(row.groups)[[1]]] <- 'black'
          cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
          rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw(df2 = rowbar_data, colscheme = 1)[[1]])))
          number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[2]]
          names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[3]]
          colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[4]]
          
        } else if(number.row.groups==6) {
          cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]]), rep(row.groups.name[6],table(row.groups)[[6]])) 
          cc2 <- rep(col1[50], length(cell2))
          cc2[1:table(row.groups)[[1]]] <- 'black'
          cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
          cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+table(row.groups)[[5]]+1:table(row.groups)[[6]]] <- 'maroon'
          rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw(df2 = rowbar_data, colscheme = 1)[[1]])))
          number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[2]]
          names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[3]]
          colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw(df2 = rowbar_data)[[4]]
          
        }
        
        ############# HEATPLOT2 EQUIVALENT HEATMAP2 CLUSTERING ###############
        z <- list()
        
        #data <- as.numeric(data)
        if(input$cc_norm == "Z-Score") {
          z <- zClust(data, scale =input$cc_norm2, zlim=c(input$cc_inSlider[1],input$cc_inSlider[2]))
        } else if (input$cc_norm == "Modified Z-Score") { 
          z <- modzClust(data, scale =input$cc_norm2, zlim=c(input$cc_inSlider[1],input$cc_inSlider[2]))
        } else if(input$cc_norm == "none") {
          z[[1]] <- as.matrix(data)
        }
        
        
        if(input$cc_dist == "pearson correlation") {
          if(input$cc_dispRow == "No" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,  key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "No" & input$cc_dispCol=='Yes' ) {
            cc_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,  key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='Yes') {
            cc_hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          }
          
          if(number.col.groups==1) {
            legend("topright", legend = paste(col.groups.name), col = "cyan", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==2) {
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("cyan", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==3) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("cyan", "khaki1", "pink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==4) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("cyan", "khaki1", "pink", "plum3"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==5) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==6) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==7) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==8) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==9) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==10) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          }
          
          if(number.row.groups==1) {
            legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==2) {
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==3) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "lightpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==4) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "lightpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==5) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "lightpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==6) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "lightpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } 
        }  else {
          if(input$cc_dispRow == "No" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "No" & input$cc_dispCol=='Yes') {
            cc_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='Yes') {
            cc_hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
          }
          
          
          if(number.col.groups==1) {
            legend("topright", legend = paste(col.groups.name), col = "cyan", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==2) {
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("cyan", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==3) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("cyan", "khaki1", "pink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==4) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("cyan", "khaki1", "pink", "plum3"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==5) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==6) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==7) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==8) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==9) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==10) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          }
          
          if(number.row.groups==1) {
            legend("bottomleft", legend = paste(row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
          } else if(number.row.groups==2) {
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==3) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==4) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4])), col = c("black", "gray", "hotpink", "steelblue1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==5) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          } else if(number.row.groups==6) {  
            legend("bottomleft", legend = paste(c(row.groups.name[1], row.groups.name[2], row.groups.name[3], row.groups.name[4], row.groups.name[5], row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1","cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeRlable)
          }
        }
        
        hmdata <- list(data = data,data_use= data3, z = list(z), hm = list(cc_hm), col1=col1, cc1=cc1, cc2 = cc2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name)
        
      }
    })
    
    output$cc_GW_subset_heatmap <- renderPlot({
      if(is.null(core_mv_hm_data())) {
        hm <- NULL
      } else {    
        core_mv_hm_data()
        #tgplot(z= core_mv_hm_data()$z[[1]], col1 = core_mv_hm_data()$col1, cc1 = core_mv_hm_data()$cc1, cc2 = core_mv_hm_data()$cc2, 
        #       number.col.groups =  core_mv_hm_data()$number.col.groups, number.row.groups = core_mv_hm_data()$number.row.groups, 
        #       row.groups.name = core_mv_hm_data()$row.groups.name, col.groups.name = core_mv_hm_data()$col.groups.name)
      }
    })
    
    cc_coldendo <- reactive({
      par(cex = input$cc_sizeClable)
      dend1 <- as.dendrogram(core_mv_hm_data()$hm[[1]]$colDendrogram)
      d <- data.frame(v1 =core_mv_hm_data()$hm[[1]]$colInd, v2=1:length(core_mv_hm_data()$hm[[1]]$colInd))
      m <- data.frame(v3 = 1:length(core_mv_hm_data()$cc1), v4 = core_mv_hm_data()$cc1)
      
      colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
      colbar2 <- colbar[,2]
      labels_colors(dend1) <- as.character(colbar2)
      plot(dend1)
      
      if(core_mv_hm_data()$number.col.groups==1) {
        legend("topright", legend = paste(core_mv_hm_data()$col.groups.name), col = "cyan", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==2) {
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2])), col = c("cyan", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==3) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3])), col = c("cyan", "khaki1", "pink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==4) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4])), col = c("cyan", "khaki1", "pink", "plum3"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==5) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==6) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==7) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==8) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==9) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8], core_mv_hm_data()$col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==10) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8], core_mv_hm_data()$col.groups.name[9], core_mv_hm_data()$col.groups.name[10])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "lightpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      }
    })
    
    output$cc_plot1 <- renderPlot({
      cc_coldendo()
    })
    
    
    cc_colDen <- reactive({
      if(input$cc_cutcolden == 'TRUE') {
        cuttable <- as.data.frame(cutree(as.hclust(core_mv_hm_data()$hm[[1]]$colDendrogram), k=as.numeric(input$cc_cuttree))[as.hclust(core_mv_hm_data()$hm[[1]]$colDendrogram)$order])
        cuttable <- cbind.data.frame(rownames(cuttable), cuttable)
        
        names(cuttable)[1] <- "Sample"
        names(cuttable)[2] <- "Cluster"
        data_l1_l2 <- core_mv_hm_data()$data_use #extracted_data2()
       
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
    
    output$cc_display <- renderUI({
      if(input$cc_cutcolden != 'TRUE') {
        return(br(strong(em("Please select Cut Col dendrogram?: = 'Yes' to display column clusters. Also select value at which you would like to cut the col dendogram (default is at k= 2)"))))
      }
    })
    
    output$cc_df <- DT::renderDataTable({
      if(input$cc_cutcolden == 'TRUE') {
        DT::datatable(cc_colDen(), options = list(
          lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
          pageLength = 5))
      }
    })
    
    output$cc_downloadCuttree <- downloadHandler(
      filename = function() {
        paste(paste(input$cc_fname_HM, input$cc_hclust, "clustering", input$cc_dist, "distance", sep="_"), '_Col_Dendrogram_cutree_', 'k=', input$cc_cuttree, '.csv', sep='') 
      },
      content = function(con) {
        write.csv(cc_colDen(), con, quote=F, row.names = F)
      })
    
    cc_rowdendo <- reactive({
      par(cex = input$cc_sizeRlable)
      dend2 <- as.dendrogram(core_mv_hm_data()$hm[[1]]$rowDendrogram)
      dd <- data.frame(v1 =rev(core_mv_hm_data()$hm[[1]]$rowInd), v2=1:length(core_mv_hm_data()$hm[[1]]$rowInd))
      mm <- data.frame(v3 = 1:length(core_mv_hm_data()$cc2), v4 = core_mv_hm_data()$cc2)
      
      colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
      colbar2 <- colbar2[,2]
      labels_colors(dend2) <- rev(as.character(colbar2))
      plot(dend2, horiz = T)
      if(core_mv_hm_data()$number.row.groups==1) {
        legend("topright", legend = paste(mv_hm_data()$row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
      } else if(mv_hm_data()$number.row.groups==2) {
        legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
      } else if(mv_hm_data()$number.row.groups==3) {  
        legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
      } else if(mv_hm_data()$number.row.groups==4) {  
        legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3], mv_hm_data()$row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
      } else if(mv_hm_data()$number.row.groups==5) {  
        legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3], mv_hm_data()$row.groups.name[4], mv_hm_data()$row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
      } else if(mv_hm_data()$number.row.groups==6) {  
        legend("topright", legend = paste(c(mv_hm_data()$row.groups.name[1], mv_hm_data()$row.groups.name[2], mv_hm_data()$row.groups.name[3], mv_hm_data()$row.groups.name[4], mv_hm_data()$row.groups.name[5], mv_hm_data()$row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeRlable)
      } 
    })
    
    output$cc_plot2 <- renderPlot({
      par(cex= input$cc_sizeRlable)
      cc_rowdendo()
    })
    
    cc_rowDen <- reactive({
      if(input$cc_cutrowden == 'TRUE') {
        cuttable2 <- as.data.frame(cutree(as.hclust(core_mv_hm_data()$hm[[1]]$rowDendrogram), k=as.integer(input$core_cuttree2))[as.hclust(mv_hm_data()$hm[[1]]$rowDendrogram)$order])
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
    
    output$cc_display2 <- renderUI({
      if(input$cc_cutrowden != 'TRUE') {
        return(br(strong(em("Please select Cut Row dendrogram?: = 'Yes' to display row clusters. Also select value at which you would like to cut the row dendrogram (default is at k= 2)"))))
      }
    })
    
    output$cc_df2 <- DT::renderDataTable({
      if(input$cc_cutrowden == 'TRUE') {
        DT::datatable(cc_rowDen(), options = list(
          lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
          pageLength = 5))
      }
    })
    
    output$cc_downloadCuttree2 <- downloadHandler(
      filename = function() {
        paste(paste(input$cc_fname_HM, input$cc_hclust, "clustering", input$cc_dist, "distance", sep="_"), '_Row_Dendrogram_cutree_', 'k=', input$cc_cuttree2, '.csv', sep='') 
      },
      content = function(con) {
        write.csv(cc_rowDen(), con, quote=F, row.names = F)
      })
    
    
    output$cc_downloadHM <- downloadHandler(
      
      filename <- function() {
        core_pdf_file <<- paste(input$cc_fname_HM, input$cc_hclust, "clustering", input$cc_dist, "distance", sep="_")
        paste('NOJAH_Core_', core_pdf_file, Sys.time(),'.pdf', sep='')
      },
      content <- function(file) {
        pdf(file=paste(core_pdf_file,".pdf",sep="") , height= 10, width=10)
        plot.new()
        title("NOJAH: Clustering Analysis",cex.main=1.2, sub = "Genome Wide analysis" ,cex.main=1.2, col = "blue", font=3)
        df <- rbind.data.frame(c("Data Normalization Type", input$cc_norm),
                               c("Normalized by (row/column/both)", input$cc_norm2),
                               c("Distance Method", input$cc_dist),
                               c("Clustering Method", input$cc_hclust),
                               c("Scale", ifelse(input$cc_norm == "none", paste(as.integer(min(core_mv_hm_data()$data)), as.integer(max(core_mv_hm_data()$data)), sep = ":"), paste(input$cc_inSlider[1], input$cc_inSlider[2], sep=":"))),
                               c("HeatMap colors", paste(input$cc_low, input$cc_mid, input$cc_high, sep="-")))
        names(df)[1] <- "Parameters"
        names(df)[2] <- "Value Selected"
        grid.table(df, rows= NULL)
        
        #eval(hm$call)
        #eval(mv_hm_data()$hm[[1]]$call) # call the heatmap here
        #eval(hm$call)
        #tgplot()
        core_mv_hm_data()
        
        par(cex = 0.6*input$cc_sizeClable)
        coldendo()
        
        par(cex = input$cc_sizeRlable)
        rowdendo()
        
        dev.off()
        file.copy(paste(core_pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
      })
    
    
   
    output$dlCoreSubset <- downloadHandler(
      filename = function() {
        paste(paste(input$cc_fname_subset, input$cc_hclust, "clustering", input$cc_dist, "distance", sep="_"), '.csv', sep='') 
      },
      content = function(con) {
        
        if(input$gw_file5 == "GW_Example5") {
          
          data5 <- extracted_data2()
          samples_to_include <- sil_data()$core
          
          cc_clust <- as.data.frame(consen()$output)
          colnames(cc_clust) <- "Cluster"
          cc_clust$Sample <- rownames(cc_clust)
          
          # sort columns based on colnames
          data1 <- data5[,order(data5[1, ])]
          data1 <- data1[order(data1[,2]),]
          data_use = data1[, colnames(data1) %in% samples_to_include]
          rownames(data_use) <- paste(data1[,1], data1[,2], sep = "|")
          
          cc_clust2 <- cc_clust[cc_clust$Sample %in% samples_to_include, ]
          cc_clust3 <- cc_clust2[match(colnames(data_use), cc_clust2$Sample),] ######
          df2 <- t(cc_clust3)
          
          data2 <- rbind(df2[-2,], data_use)
          rownames(data2)[2] <- "Groups"
          
          gene_id <- sub('\\|.*', '', rownames(data2))
          Groups <- sub('.*\\|', '', rownames(data2))
          
          data3 <- data.frame(gene_id, Groups, data2)
          data3[1:2,c(1:2)] <-  " "
        
        } else if(input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3") {
          data5 = sil_data()$data
          samples_to_include <- sil_data()$core
          cc_clust = sil_data()$c_order
          
          # sort columns based on colnames
          data1 <- data5[,order(data5[1, ])]
          data1 <- data1[order(data1[,2]),]
          data_use = data1[, colnames(data1) %in% samples_to_include]
          
          
          # sort columns based on colnames
          data1 <- data5[,order(data5[1, ])]
          data1 <- data1[order(data1[,2]),]
          data_use = data1[, colnames(data1) %in% samples_to_include]
          
          cc_clust2 <- cc_clust[cc_clust$Sample %in% samples_to_include, ]
          cc_clust3 <- cc_clust2[match(colnames(data_use), cc_clust2$Sample),] ######
          df2 <- t(cc_clust3)
          
          data2 <- rbind(df2[-1,], data_use)
          #rownames(data2)[2] <- "Groups"
          
          
          gene_id <- sub('\\.|*', '', rownames(data2))
          Groups <- rep("A", length(gene_id))
          
          data3 <- data.frame(gene_id, Groups, data2)
          data3[1,c(1:2)] <-  " "
        }
        
      
        
      write.csv(data3, con, quote=F, row.names = F)    
    })
    

###########################################################################
  
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
          
          cc1 = vector()
          cc2 = vector()
          
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
      
      cc1 = vector()
      cc2 = vector()
      
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
          Sig_hm <<- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
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


