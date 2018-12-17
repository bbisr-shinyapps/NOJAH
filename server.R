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
source("functions/colbars_gw2.R")
source("functions/sil_function.R")
source("functions/sil_function2.R")
source("functions/sil_cpt.R")
source("functions/sil_cpt2.R")
source("functions/cptGenes.R")

function(input, output) {
  
  output$ReadMe <- renderUI({
    str0 <- paste("&emsp;")
    str00 <- paste("DATA INPUT")
    str1 <- paste("Data should be input as a TXT or a CSV file. The first two rows of the data file have information about the patients/specimens and their response/subtype; all remaining rows have gene expression data, one row per gene.  In the case of Microarray gene expression data in which there are several probes corresponding to a single gene, a unique identifier would need to be created to separately identify each probe such as, 'Gene 1_p1', 'Gene1_p2' indicating Gene 1 has two probes.  The columns represent the different experimental samples. A maximum of up to 10 different sample groups and 6 different gene groups may be used with this tool.")
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
    str.1 <- paste("Table 1 : Example dataset for two gene groups (over and under-expressed) and two patient groups (Normal, Tumor). Numeric data starts at row three and column three.")
    HTML(paste(strong(str.1),str.0, str.0, str.0, sep = '<br/>'))
  })
  
  output$Eg2 <- renderTable({
    coln1 <- c("gene_id", "Groups","GSM9981", "GSM1870", "GSM4618", "GSM7689", "GSM8772", "GSM1121","GSM1250", "GSM3112", "GSM4987", "GSM1277")
    coln2 <- c(" ", " ", rep("MM", 5), rep("MUGS", 2), "NPC", rep("SM", 2))
    coln3 <- c("Subtype", "", rep("Classical", 2), rep("Neural", 1), rep("Pro-neural", 2), "Classical", rep("Mesenchymal", 1), "Classical", "Neural", "Neural")
    s.1 <- c("YWHAE>210996_s_at", "na", 1.47, 2.18, 5.87, 9.12, 7.34, 1.56, 3.0, 7.77, 3.40, 1.56 )
    s.2 <- c("YWHAE>201020_at", "na", 1.98, 7.93, 2.76, 9.11, 8.46, 0.98, 5.98, 8.19, 8.91, 5.98)
    s.3 <- c("YWHAH>33323_r_at", "na", 8.02, 8.00, 2.17, 10.12, 8.76, 9.76, 3.76, 0.02, 3.67, 7.94)
    s.4 <- c("YWHAB>208743_s_at", "na", 2.75, 5.99, 3.19, 11.86, 6.54, 8.17, 2.00, 0.99, 2.00, 1.17)
    s.5 <- c("YWHAQ>213699_s_at", "na", 9.35, 8.96, 6.67, 8.33, 3.98, 7.11, 1.67, 1.01, 5.18, 8.17)
    
    d.f <- rbind.data.frame(coln2, coln3, s.1, s.2, s.3, s.4, s.5)
    colnames(d.f) <- coln1
    head(d.f)
  })
  
  
  output$Caption2 <- renderUI({
    strg.0 <- paste("&emsp;")
    strg.1 <- paste("Table 2: Example dataset for one gene group (marked na) and four patient groups (MM, MUGS, NPC and, SM). Additonal subtype information is appended above the numeric data. Numeric data starts at row four and column three. Any missing entries in subtype should be coded as none. NA or blank may display errors.")
    strg.2 <- paste(" Where applicable, both row and column dendrograms can be extracted in their specific tabs. Using the options on the right panel, dendrograms can be cut into desired no. of clusters (default at 2) and a pvalue of significance between the clusters can be determined using bootstrap method. ")
    HTML(paste(strong(strg.1),strg.0, strg.0,strg.2, strg.0, strg.0, sep = '<br/>'))
  })
  
  output$download_GW_Ex1 <- downloadHandler(
    filename= function() {paste('TCGA.BRCA.Expression.csv')}, 
    content = function(file) {
      d <- read.csv("data/TCGA_BRCA_Early_Late_OS_6.9years.csv", header = TRUE, sep = ",", stringsAsFactors = F)
      write.csv(d, file, row.names = FALSE) }
  )
  
  
  
  input_gw_data <- reactive({
    if(input$gw_file1 == 'GW_Example1'){
      data <- read.csv("data/TCGA_BRCA_Early_Late_OS_6.9years.csv", header = TRUE, sep = ",", stringsAsFactors = F)
    }
    else if(input$gw_file1 == 'load_my_own_gw'){
      inFile <- input$gw_file2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".csv", inFile[1])) { data = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
      else if(grepl(".txt", inFile[1])) { data = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
      else if(grepl(".rds", inFile[1])) { data = readRDS(as.character(inFile$datapath)) }
      
    }
   
    n = ncol(data)-2
    
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
    
    # box plot data
    data2 <- cbind.data.frame(data2, var=as.numeric(getVar)) #apply(data2[, c(1:n)], 1, var)
    data2$mad <- apply(data2[, c(1:n)], 1, mad)
    data2$IQR <- apply(data2[, c(1:n)], 1, IQR)
    
    
    #ngenes_extracted = nrow(data5)
   
    dat = list(data=data, n= n, dendo_data = data4, colgroups = col.groups, colgroupsname= col.groups.name, numbercolgroups = number.col.groups, bp_data = data2,  
               genes = as.character(data[, 1]),  n_genes = nrow(data4),n_samples= ncol(data4)
               # sp_data = data2.2, extracted_data = data5, n_genes = nrow(data4),n_samples= ncol(data4), percen = percen, ngenes_extracted = nrow(data5)
               )
     
    })
  
  
  scatter_plot_data <- reactive({
    
    data2.2 <- input_gw_data()$bp_data
    n  <- input_gw_data()$n
    
    check_data2.2 <<- data2.2
    #scatter plot data
   # data2.2 <- data2
    data2.2$Rank.var <- rank(data2.2$var,ties.method= "min")
    data2.2$Rank.mad <- rank(data2.2$mad,ties.method= "min")
    data2.2$Rank.IQR <- rank(data2.2$IQR,ties.method= "min")
    data2.2$sumofranks <- data2.2$Rank.var + data2.2$Rank.mad + data2.2$Rank.IQR
    
    if(is.null(input$gw_subset)) {
      return()
    } else if(length(input$gw_subset) == 1) {
      if(input$gw_subset == 'VAR') {
        data2.2 <- data2.2[order(data2.2$var),]
        data2.2$var_percen <- ifelse(input$var_PercenChoice == 'Percentile Slider', quantile(data2.2$var, as.numeric(input$var_pslider)/100),quantile(data2.2$var, as.numeric(input$var_pInput)/100))
        percen <- ifelse(input$var_PercenChoice == 'Percentile Slider', input$var_pslider,  input$var_pInput)
      } else if(input$gw_subset == 'MAD') {
        data2.2 <- data2.2[order(data2.2$mad),]
        data2.2$mad_percen <- ifelse(input$mad_PercenChoice == 'Percentile Slider', quantile(data2.2$mad, input$mad_pslider/100),quantile(data2.2$mad, input$mad_pInput/100))
        percen <- ifelse(input$mad_PercenChoice == 'Percentile Slider', input$mad_pslider,  input$mad_pInput)
      } else if(input$gw_subset == 'IQR') {
        data2.2 <- data2.2[order(data2.2$IQR),]
        data2.2$iqr_percen <- ifelse(input$iqr_PercenChoice == 'Percentile Slider', quantile(data2.2$IQR, input$iqr_pslider/100),quantile(data2.2$IQR, input$iqr_pInput/100))
        percen <- ifelse(input$iqr_PercenChoice == 'Percentile Slider', input$iqr_pslider,  input$iqr_pInput)
      } 
      
    }
    else 
      if(length(input$gw_subset) > 1) {
        
        if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset) {
          data2.2 <- data2.2[order(data2.2$sumofranks),]
          data2.2$IMVA_percen <- ifelse(input$IMVA_PercenChoice == 'Percentile Slider', quantile(data2.2$sumofranks, as.integer(input$IMVA_pslider)/100),quantile(data2.2$sumofranks, as.integer(input$IMVA_pInput)/100))
          percen <- ifelse(input$IMVA_PercenChoice == 'Percentile Slider', input$IMVA_pslider,  input$IMVA_pInput)
          
        } else if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & !("IQR" %in% input$gw_subset) ) {
          data2.2$sumofranks_VM <- data2.2$Rank.var + data2.2$Rank.mad 
          data2.2 <- data2.2[order(data2.2$sumofranks_VM),]
          data2.2$var_mad_percen <- ifelse(input$var_mad_PercenChoice == 'Percentile Slider', quantile(data2.2$sumofranks_VM, input$var_mad_pslider/100),quantile(data2.2$sumofranks_VM, input$var_mad_pInput/100))
          percen <- ifelse(input$var_mad_PercenChoice == 'Percentile Slider', input$var_mad_pslider, input$var_mad_pInput)
          
        } else if("VAR" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("MAD" %in% input$gw_subset)) {
          data2.2$sumofranks_VI <- data2.2$Rank.var + data2.2$Rank.IQR 
          data2.2 <- data2.2[order(data2.2$sumofranks_VI),]
          data2.2$var_iqr_percen <- ifelse(input$var_iqr_PercenChoice == 'Percentile Slider', quantile(data2.2$sumofranks_VI, input$var_iqr_pslider/100),quantile(data2.2$sumofranks_VI, input$var_iqr_pInput/100))
          percen <- ifelse(input$var_iqr_PercenChoice == 'Percentile Slider', input$var_iqr_pslider,  input$var_iqr_pInput)
          
        } else  if("MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("VAR" %in% input$gw_subset)) {
          data2.2$sumofranks_MI <- data2.2$Rank.mad + data2.2$Rank.IQR 
          data2.2 <- data2.2[order(data2.2$sumofranks_MI),]
          data2.2$mad_iqr_percen <- ifelse(input$mad_iqr_PercenChoice == 'Percentile Slider', quantile(data2.2$sumofranks_MI, input$mad_iqr_pslider/100),quantile(data2.2$sumofranks_MI, input$mad_iqr_pInput/100))
          percen <- ifelse(input$mad_iqr_PercenChoice == 'Percentile Slider', input$mad_iqr_pslider,  input$mad_iqr_pInput)
          
        }  
        
        
      }
    
    if(is.null(input$gw_subset)) {
      return(NULL)
    } else if(length(input$gw_subset) == 1) {
      if(input$gw_subset == 'VAR') {
        data5 <- data2.2[data2.2$var > data2.2$var_percen[1], 1:n]
      } else if(input$gw_subset == 'MAD') {
        data5 <- data2.2[data2.2$mad > data2.2$mad_percen[1], 1:n]
      } else if(input$gw_subset == 'IQR') {
        data5 <- data2.2[data2.2$IQR > data2.2$iqr_percen[1], 1:n]
      } 
    }
    else 
      if(length(input$gw_subset) > 1) {
        if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset) {
          data5 <- data2.2[data2.2$sumofranks > data2.2$IMVA_percen[1], 1:n]
        }else if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & !("IQR" %in% input$gw_subset)) {
          data5 <- data2.2[data2.2$sumofranks_VM > data2.2$var_mad_percen[1], 1:n]
        } else if("VAR" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("MAD" %in% input$gw_subset)) {
          data5 <- data2.2[data2.2$sumofranks_VI > data2.2$var_iqr_percen[1], 1:n]
        } else  if("MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("VAR" %in% input$gw_subset)) {
          data5 <- data2.2[data2.2$sumofranks_MI > data2.2$mad_iqr_percen[1], 1:n]
        }  
      }
    
    dat  = list(sp_data = data2.2, extracted_data = data5, percen = percen, ngenes_extracted = nrow(data5))
  
  })
  
  dend1 <- eventReactive(input$button1, {
    data4 <- input_gw_data()$dendo_data
    
    col.groups <- input_gw_data()$colgroups
    col.groups.name <- input_gw_data()$colgroupsname
    number.col.groups <- input_gw_data()$numbercolgroups
    
    cc1 = vector()
    
    ## Set color palette
    col1 <- colorRampPalette(c("red","black","green"))(299)
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
    if(input$GW_dist == "pearson correlation") {
      dist.func =  function(x) as.dist((1-cor(t(x))))
      hc = hclust(dist.func(t(z[[1]])), method = input$GW_hclust)
    } else {
      hc = hclust(dist(t(z[[1]]), method = input$GW_dist), method = input$GW_hclust)
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
    } else if(number.col.groups==2) {
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2])), col = c("darkblue", "grey"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
    } else if(number.col.groups==3) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3])), col = c("darkblue", "grey", "orange"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==4) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4])), col = c("darkblue", "grey", "orange", "yellow"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==5) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5])), col = c("darkblue", "grey", "orange", "yellow", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==6) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==7) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==8) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==9) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    } else if(number.col.groups==10) {  
      legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$gw_dend_size)
    }
    
  })
  
  output$gw_dend <- renderPlot({ #renderCachedPlot({ #
     dend1()
  }#,
 # cacheKeyExpr = { 
  #  list(input$gw_file1, input$GW_norm) #, input$GW_norm2, input$GW_inSlider, input$GW_dist, input$GW_hclust, input$gw_dend_size)
  #}, cache = memoryCache()
  )
  
  
  
  #output$geneSelector <- renderUI({
  #  #selectizeInput(inputId = "Genes", "Choose Option:", as.list(getOSgenes()),options=list(maxOptions=getOSgenes())) 
  #  selectizeInput(inputId = "Genes", "Choose Option:", as.list(input_gw_data()$genes[1:605]),options=list(maxOptions=input_gw_data()$n_genes[1:605])) 
  #})
  
  #output$dropdowngene <- renderText({ 
  #  paste("You have selected gene", input$Genes)
  #})
  
  
  #getOSgenes <- eventReactive(input$button1, {
  #  if(!is.null(input_gw_data()$genes)) 
  #  {
  #    return(input_gw_data()$genes)
  #  }
  #  else 
  #    return(NULL)
  #})
  
  gs <- eventReactive(input$button1, {
    
    n_sel <<- scatter_plot_data()$ngenes_extracted
    isolate(
      HTML(paste("The number of genes selected : ", strong(n_sel)), sep = ' ')
    )
  })
  
  
  output$n_selected <- renderUI({ 
    gs()
    
  })
  
  bplot <- eventReactive(input$button1, {
    data <- input_gw_data()$bp_data
    
    boxplot(data$var, data$mad, data$IQR, names = c("VAR", "MAD", "IQR"), border = c("blue", "orange", "green"),col = c("blue", "orange", "green"))
   #beanplot(data$var, data$mad, data$IQR, col = c("orange", "blue", "green"), border = "black")
  })
  
 # output$Boxplot <- renderPlotly({
  output$Boxplot <- renderPlot({
       bplot()
    
  })
  
 rt <- eventReactive(input$button1, {
   data <- scatter_plot_data()$sp_data #input_gw_data()$sp_data
   gene_names <- sub("\\|.*", "", rownames(data))
   
   goi <- ifelse(input$Genes== "", NA, input$Genes)
   
   if(!is.na(goi) & (any(goi == gene_names)== FALSE))
   {
     paste("Please select valid gene id")
   }
 })  
  
 output$text1 <- renderText({
   rt()
  
   })
 
 
 splot <- eventReactive(input$button1, {
 # splot <- reactive ({ 
   data <- scatter_plot_data()$sp_data #input_gw_data()$sp_data
   
   gene_names <- sub("\\|.*", "", rownames(data))
   
   goi <- ifelse(input$Genes== "", NA, input$Genes)
   #goi <- ifelse(input$Sgene == "", NA, input$Sgene)
    
    ding <<- goi
    dong <<- data
    
    if(is.null(input$gw_subset)) {
      return()
    }
    
    if(length(input$gw_subset) == 1) {
      if(input$gw_subset == 'VAR') {   
       
        plot(data$var, col = "blue", ylab= "VAR", xlab = "", pch = 20)
        if(!is.na(goi[1])) {
          poi <- which(goi == gene_names)
          points(x= poi, y = data[poi,]$var, col = "black", pch =20, cex= 1.5)  
        }
        abline(h= data$var_percen[1], col = "red", lty =2)
  
      }
      else if(input$gw_subset == 'MAD') {
        plot(data$mad, col = "darkgoldenrod1", ylab= "MAD", xlab = "", pch = 20)
        if(!is.na(goi[1])) {
          poi <- which(goi == gene_names) 
          points(x= poi, y = data[poi,]$mad, col = "black", pch =20, cex= 1.5)  
        }
        abline(h= data$mad_percen[1], col = "red", lty =2)
        
      } else if(input$gw_subset == 'IQR') {
        plot(data$IQR, col = "green", ylab= "IQR", xlab = "", pch = 20)
        if(!is.na(goi[1])) {
          poi <- which(goi == gene_names)
          points(x= poi, y = data[poi,]$IQR, col = "black", pch =20, cex= 1.5) 
          }
        abline(h= data$iqr_percen[1], col = "red", lty =2)
      } 
      
    } else 
      if(length(input$gw_subset) > 1) {
        
        if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset) {
          
          plot(data$sumofranks, col = "grey", ylab= "Ordered sum of ranks of VAR, MAD and IQR", xlab = "", pch = 20)
          if(!is.na(goi[1])) {
            poi <- which(grepl(input$Genes, rownames(data))) 
            points(x= poi, y = data[poi,]$sumofranks, col = "black", pch =20, cex= 1.5)  
          }
          abline(h= data$IMVA_percen[1], col = "red", lty =2)
        }
        
        else if("VAR" %in% input$gw_subset & "MAD" %in% input$gw_subset & !("IQR" %in% input$gw_subset) ) {
          
          plot(data$sumofranks_VM, col = "grey", ylab= "Ordered sum of ranks of VAR and MAD", xlab = "", pch = 20)
          if(!is.na(goi[1])) {
            poi <- which(goi == gene_names)
            points(x= poi, y = data[poi,]$sumofranks_VM, col = "black", pch =20, cex= 1.5)  
          }
          abline(h= data$var_mad_percen[1], col = "red", lty =2)
        }
        else if ("VAR" %in% input$gw_subset & "IQR" %in% input$gw_subset& !("MAD" %in% input$gw_subset)) {
          
          plot(data$sumofranks_VI, col = "grey", ylab= "Ordered sum of ranks of VAR and IQR", xlab = "", pch = 20)
          if(!is.na(goi[1])) {
            poi <- which(goi == gene_names)
            points(x= poi, y = data[poi,]$sumofranks_VI, col = "black", pch =20, cex= 1.5)  
          }
          abline(h= data$var_iqr_percen[1], col = "red", lty =2)
          
          }
        else if ("MAD" %in% input$gw_subset & "IQR" %in% input$gw_subset & !("VAR" %in% input$gw_subset)) {
          
          plot(data$sumofranks_MI, col = "grey", ylab= "Ordered sum of ranks of MAD and IQR", xlab = "", pch = 20)
          if(!is.na(goi[1])) {
            poi <- which(goi == gene_names)
            points(x= poi, y = data[poi,]$sumofranks_MI, col = "black", pch =20, cex= 1.5)  
          }
          abline(h= data$mad_iqr_percen[1], col = "red", lty =2)
        } 
      }
    
  })
  
  
  output$GW_Scatter_LH <- renderPlot({
     splot()
  }) 
  
  
  extracted_data2 <- reactive({
    if (input$gw_file_1 == 'GW_Example_1') {
       
        data3 <- scatter_plot_data()$extracted_data #input_gw_data()$extracted_data
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
        else if(grepl(".rds", inFile_1[1])) { data5 = readRDS(as.character(inFile_1$datapath)) }
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
    if(is.null(extracted_data2())) {
        return(NULL)
    } else {
      
      if(input$gw_file_1 == "GW_Example_1") {
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
      
      colbars2 = matrix()
      rowbars2 = matrix()
      
      ## Set color palette
      col1 <- colorRampPalette(c(input$low,input$mid,input$high))(299)
      colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
      
      ### Color vector for columns
      colbars2 <- col_color(col1 = col1, col.groups= col.groups, number.col.groups= number.col.groups, col.groups.name= col.groups.name)
      colbars2 <- as.matrix(colbars2, nrow =ncol(data))
     
      ### Color vector for rows
      rowbars2 <- row_color(col1 = col1, row.groups= row.groups, number.row.groups= number.row.groups, row.groups.name= row.groups.name)
      rowbars2 <- as.matrix(rowbars2, ncol = nrow(data))
      rowbars2 <- t(rowbars2)
      
      number.colbar.class <- NULL
      names.colbar.class <- NULL
      colors_used <- NULL
      
      number.rowbar.class <- NULL
      names.rowbar.class <- NULL
      colors_used2 <- NULL
      
    } else if(input$gw_file_1 == "load_my_own_gw_subset_1") {
      inFile_2 <- input$gw_file_2
      if (is.null(inFile_2))
        return(NULL)
      else if(grepl(".csv", inFile_2[1])) { data1 = read.csv(as.character(inFile_2$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
      else if(grepl(".txt", inFile_2[1])) { data1 = read.table(as.character(inFile_2$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
      else if(grepl(".rds", inFile_2[1])) { data1 = readRDS(as.character(inFile_2$datapath)) }
      
      #dding <<- data1
      
      data1 <- data1[,order(data1[1, ])]
      data2 <- data1[order(data1[,2]),]
      
      # remove rows with zero variance
      #data2 <- data2[apply(data2[(as.integer(input$DataR2)-1):nrow(data2),as.numeric(input$DataC2):ncol(data2)], 1, var) != 0, ]
      
      
      ### row groups
      gene <- as.character(data2$gene_id)
      gene <- gene[(as.integer(input$DataR2)-1):nrow(data2)] #gene[(6-1):nrow(data)] 
      row.groups <- as.character(as.vector(data2[(as.integer(input$DataR2)-1): nrow(data2),2])) #as.character(as.vector(data[(6-1):nrow(data),2])) 
      row.groups.name <- names(table(row.groups))
      number.row.groups <- length(row.groups.name)
      
      ### column groups
      col.groups <- as.character(as.vector(data2[1,]))
      col.groups <-  col.groups[as.numeric(input$DataC2):ncol(data2)] #col.groups[4:ncol(data)] 
      col.groups.name <- names(table(col.groups))
      number.col.groups <- length(col.groups.name)
      check_number.col.groups <<- number.col.groups
      
      ## additional info on columns and rows
      colbar_data <-  as.matrix(t(data2[1:(as.numeric(input$DataR2)-2), c(1, (as.integer(input$DataC2)-1):ncol(data2))])) #t(data[1:(6-2), c(1, 3:ncol(data))]) 
      colnames(colbar_data) <- as.character(unlist(colbar_data[1,]))
      colbar_data <- as.matrix(colbar_data[c(-1, -2),])
      colnames(colbar_data)[1] <- "Groups"
      n.colbar_data <- ncol(colbar_data)
      
      check_colbar_data <<- colbar_data
      
      rowbar_data <- as.matrix(data2[(as.integer(input$DataR2)-1):nrow(data2), 1:(as.numeric(input$DataC2)-1)]) #as.matrix(data[5:nrow(data), 1:3])
      rownames(rowbar_data) <- rowbar_data[,1]
      rowbar_data <- as.matrix(rowbar_data[,-1])
      #colnames(rowbar_data)[1] <- "Groups"
      n.rowbar_data <- ncol(rowbar_data)
      
      data3 <- data2[(as.numeric(input$DataR2)-1):nrow(data2), as.numeric(input$DataC2):ncol(data2)] #data[5:nrow(data), 4:ncol(data)]
      rownames(data3) <- gene
      data4 <-data3[complete.cases(data3),]
      data <- data.matrix(data4)
      
      cc1 = vector()
      cc2 = vector()
      
      ## Set color palette
      col1 <- colorRampPalette(c(input$low,input$mid,input$high))(299)
      colors <- c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)) # check slider
      
      ### Color vector for columns
      if(number.col.groups==1) { 
        cell <- c(rep(col.groups.name, number.col.groups))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
      } else if(number.col.groups==2) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'                                                                                                                                   
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey' 
        check_colbar_data <<- colbar_data
        check_cc1 <<- cc1
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        check_colbars <<- colbars2
        
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
      } else if(number.col.groups==3) {
        #colors = c("darkblue", "grey", "orange", "plum3", "mediumaquamarine")
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        #colbars2 <- cbind(cc1, colbars_gw(df2 = colbar_data))
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
      } else if(number.col.groups==4) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
      } else if(number.col.groups==5) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
      } else if(number.col.groups==6) {
        cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]])) 
        cc1 <- rep(col1[50], length(cell))
        cc1[1:table(col.groups)[[1]]] <- 'darkblue'
        cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'grey'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- 'orange'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- 'yellow'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- 'purple'
        cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- 'darkgreen'
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
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
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
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
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else colbars_gw2(df2 = colbar_data)[[4]]
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
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[4]]
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
        colbars2 <- if(n.colbar_data ==1) as.matrix(cc1) else as.matrix(cbind(cc1, colbars_gw2(df2 = colbar_data)[[1]]))
        colnames(colbars2)[1] <- "CC"
        number.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) NULL else  colbars_gw2(df2 = colbar_data)[[4]]
      }
      
      ### Color vector for rows
      
      if(number.row.groups==1) { 
        cell2 <- c(rep(row.groups.name, number.row.groups))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw2(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[4]]
        check_rowbars2 <<- rowbars2
      } else if(number.row.groups==2) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]))
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw2(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[4]]
        
      } else if(number.row.groups==3) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw2(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[4]]
        
      } else if(number.row.groups==4) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw2(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[4]]
      } else if(number.row.groups==5) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw2(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[4]]
        
      } else if(number.row.groups==6) {
        cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]]), rep(row.groups.name[6],table(row.groups)[[6]])) 
        cc2 <- rep(col1[50], length(cell2))
        cc2[1:table(row.groups)[[1]]] <- 'black'
        cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- 'gray'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- 'hotpink'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- 'brown1'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- 'cyan'
        cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+table(row.groups)[[5]]+1:table(row.groups)[[6]]] <- 'maroon'
        rowbars2 <- if(n.rowbar_data ==1) as.matrix(t(cc2)) else as.matrix(t(cbind(cc2, colbars_gw2(df2 = rowbar_data, colscheme = 1)[[1]])))
        number.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[2]]
        names.rowbar.class <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[3]]
        colors_used2 <- if(n.rowbar_data ==1) NULL else colbars_gw2(df2 = rowbar_data)[[4]]
        
      }
    }
      
      
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
          hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$dispRow == "No" & input$dispCol=='Yes' ) {
          hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$dispRow == "Yes" & input$dispCol=='No') {
          hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
          hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
       
       if(!is.null(number.colbar.class)) {
         if(number.colbar.class==1) {
           legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.colbar.class==2) {
           legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.colbar.class==3) {  
           legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.colbar.class==4) {  
           legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.colbar.class==5) {  
           legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.colbar.class==6) {  
           legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } 
       }
       
       if(!is.null(number.rowbar.class)) {
         if(number.rowbar.class== 1) {
           legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.rowbar.class==2) {
           legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.rowbar.class==3) {  
           legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.rowbar.class==4) {  
           legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.rowbar.class==5) {  
           legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } else if(number.rowbar.class==6) {  
           legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
         } 
       }
      }  else {
        if(input$dispRow == "No" & input$dispCol=='No') {
          hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$dispRow == "No" & input$dispCol=='Yes') {
          hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$dispRow == "Yes" & input$dispCol=='No') {
          hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
          hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
        
        if(!is.null(number.colbar.class)) {
          if(number.colbar.class==1) {
            legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==2) {
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==3) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==4) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==5) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==6) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
        }
        
        if(!is.null(number.rowbar.class)) {
          if(number.rowbar.class== 1) {
            legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==2) {
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==3) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==4) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==5) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==6) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
        }
      }
    
   
     hmdata <- list(data = data, z = list(z), hm = list(hm), col1=col1, cc1= colbars2, cc2 = rowbars2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name, 
                    number.colbar.class = number.colbar.class, names.colbar.class = names.colbar.class, colors_used = colors_used, number.rowbar.class= number.rowbar.class, names.rowbar.class = names.rowbar.class,  colors_used2 = colors_used2)
    }
  })
  
  tgplot <- function(z, col1, cc1, cc2, number.col.groups, number.row.groups, row.groups.name , col.groups.name, number.colbar.class, names.colbar.class, colors_used, number.rowbar.class, names.rowbar.class, colors_used2) {  
                
    if(input$dist == "pearson correlation") {
      if(input$dispRow == "No" & input$dispCol=='No') {
        hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
      } else if(input$dispRow == "No" & input$dispCol=='Yes' ) {
        hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
      } else if(input$dispRow == "Yes" & input$dispCol=='No') {
        hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,  key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
      } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
        hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
      
      if(!is.null(number.colbar.class)) {
        if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==3) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==4) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==5) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==6) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } 
      }
      
      if(!is.null(number.rowbar.class)) {
        if(number.rowbar.class== 1) {
          legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==2) {
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==3) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==4) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==5) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==6) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } 
      }
    }  else {
      if(input$dispRow == "No" & input$dispCol=='No') {
        hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
      } else if(input$dispRow == "No" & input$dispCol=='Yes') {
        hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
      } else if(input$dispRow == "Yes" & input$dispCol=='No') {
        hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
      } else if(input$dispRow == "Yes" & input$dispCol=='Yes') {
        hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$clust_byrow))), Colv=eval(parse(text=paste(input$clust_bycol))), hclust=function(c) {hclust(c,method=input$hclust)}, distfun=function(c) {dist(c,method=input$dist)},cexRow=input$size1,cexCol =input$size2,key=TRUE,keysize=1.0, margin = c(input$inSlider2[1],input$inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=cc1, RowSideColors = cc2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
      
      if(!is.null(number.colbar.class)) {
        if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==3) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==4) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==5) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==6) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } 
      }
      
      if(!is.null(number.rowbar.class)) {
        if(number.rowbar.class== 1) {
          legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==2) {
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==3) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==4) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==5) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.rowbar.class==6) {  
          legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } 
      }
      
    }
  }
  
  hmplot <- eventReactive(input$button2, {
  if(is.null(mv_hm_data())) {
    hm <- NULL
  } else {    
    tgplot(z= mv_hm_data()$z[[1]], col1 = mv_hm_data()$col1, cc1 = mv_hm_data()$cc1, cc2 = mv_hm_data()$cc2, 
           number.col.groups =  mv_hm_data()$number.col.groups, number.row.groups = mv_hm_data()$number.row.groups, 
           row.groups.name = mv_hm_data()$row.groups.name, col.groups.name = mv_hm_data()$col.groups.name,
           number.colbar.class = mv_hm_data()$number.colbar.class, names.colbar.class = mv_hm_data()$names.colbar.class, 
           colors_used = mv_hm_data()$colors_used, 
           number.rowbar.class= mv_hm_data()$number.rowbar.class, names.rowbar.class = mv_hm_data()$names.rowbar.class,  
           colors_used2 = mv_hm_data()$colors_used2 )
  }
  })
  
  output$GW_subset_heatmap <- renderPlot({
    hmplot()
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
          m2$order <- 1:nrow(m2)
          m3 <- m2[order(-m2$order),]
          rownames(m3) <- 1:nrow(m3)
          m4 <- m3[, c(1,2, 3)]
          
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
                 row.groups.name = mv_hm_data()$row.groups.name, col.groups.name = mv_hm_data()$col.groups.name,
                 number.colbar.class = mv_hm_data()$number.colbar.class, names.colbar.class = mv_hm_data()$names.colbar.class, 
                 colors_used = mv_hm_data()$colors_used, 
                 number.rowbar.class= mv_hm_data()$number.rowbar.class, names.rowbar.class = mv_hm_data()$names.rowbar.class,  
                 colors_used2 = mv_hm_data()$colors_used2)
          
          
          if(input$clust_bycol == TRUE) {
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
            
          }
          
          if(input$clust_byrow == TRUE){
          #par(cex = input$sizeRlable)
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
          }  
          
          dev.off()
          file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
        })
      
      
      output$download_GW_Ex_1 <- downloadHandler(
        filename= function() {paste('Example_subset_data.csv')}, 
        content = function(file) {
          d <- extracted_data2()
          
          write.csv(d, file, row.names = FALSE) }
      )
      
      
    
    consen <- eventReactive(input$button3, {
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
        else if(grepl(".rds", inFile4[1])) { mv_data = readRDS(as.character(inFile4$datapath)) }
        
      } else
        return(NULL)
      
      mv_data2 <- as.matrix(mv_data)
        
      con <- consensus_clustering(dinput=mv_data2, mK=as.integer(input$con_maxK), rep=as.integer(input$con_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$con_dist, iL=input$con_hclust, fL=input$con_hclust)
      
      
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
     
     upto_slider <- list()
     
     
     if(input$sil_choice == 'Fixed value'){
       upto_slider = list(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, input$upto_slider9)
       check_upto_slider <<- upto_slider
       
       sil = silhouette_plot2(data_use= consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], upto_width= upto_slider, cols = df)
     
       colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
       par(mfrow = c(1, 2))
       plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
       plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
    
     } else if(input$sil_choice == 'Percentile') {
       
       upto_percen = list(input$sil_percen1, input$sil_percen2, input$sil_percen3, input$sil_percen4, input$sil_percen5, input$sil_percen6, input$upto_slider7, input$sil_percen8, input$sil_percen9)
       
       sil = sil_function(data_use = consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], percen = upto_percen, cols = df)
       
       check_percen <<- as.numeric(input$sil_percen)
       check_sil <<- sil
       
       colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
       par(mfrow = c(1, 2))
       plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
       plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
       
       
       
     } else if(input$sil_choice == 'Change point'){
      
       upto_sil_cp = list(input$sil_cp1, input$sil_cp2, input$sil_cp3, input$sil_cp4, input$sil_cp5, input$sil_cp6, input$sil_cp7, input$sil_cp8, input$sil_cp9)
       
       sil = sil_cpt(data_use = consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], sil_cp = upto_sil_cp, cols = df, cpt_data = input$changes, cpt_method = input$method, cpt_max = input$max)
       
       colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
       par(mfrow = c(1, 2))
       plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
       plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
       
       
     }
     
     return(list(data = consen()[["data"]], core = sil$core.samples, clust = sil$sk3, c_order = consen()[["output"]] ))
     
     } else if(input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3"){
       inFile6 <- input$gw_file6
       if (is.null(inFile6))
         return(NULL)
       else if(grepl(".csv", inFile6[1])) { mv_data = read.csv(as.character(inFile6$datapath), header = TRUE, sep = ",", stringsAsFactors = F, row.names = 1) }
       else if(grepl(".txt", inFile6[1])) { mv_data = read.table(as.character(inFile6$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, row.names = 1) }
       else if(grepl(".rds", inFile6[1])) { mv_data = readRDS(as.character(inFile6$datapath)) }
       
       inFile8 <- input$gw_file8
       if (is.null(inFile8))
         return(NULL)
       else if(grepl(".csv", inFile8[1])) { c_order = read.csv(as.character(inFile8$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
       else if(grepl(".txt", inFile8[1])) { c_order = read.table(as.character(inFile8$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
       else if(grepl(".rds", inFile8[1])) { c_order = readRDS(as.character(inFile8$datapath)) }
       
       mv_data2 <- as.data.frame(mv_data)
       c_order2 <- as.data.frame(c_order)
       
       colnames(c_order2)[2] <- "Cluster"
       
        k= unique(c_order2$Cluster)[length(unique(c_order2$Cluster))]
       
       if(input$sil_choice == 'Fixed value'){
         
       upto_slider = list(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, input$upto_slider9)
       check_upto_slider <<- upto_slider
         
       sil = silhouette_plot3(data_use= mv_data2, opt_k= k, res= c_order2, dist = input$sil_dist, upto_width= upto_slider, cols = c("orange", "darkblue", "black", "maroon", "violet", "plum2", "lightyellow", "lightblue", "green")[1:k])
         
       colors = c("cyan", "khaki1", "pink", "plum3", "purple", "darkgreen", "hotpink","brown", "darkorchid2",  "maroon")
       par(mfrow = c(1, 2))
       plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:k])
       plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:k]) #sil[["sk3.col"]]
       
       } else if(input$sil_choice == 'Percentile') {
         
         upto_percen = list(input$sil_percen1, input$sil_percen2, input$sil_percen3, input$sil_percen4, input$sil_percen5, input$sil_percen6, input$upto_slider7, input$sil_percen8, input$sil_percen9)
         check_upto_percen <<- upto_percen
         
         sil = sil_function2(data_use = mv_data2, opt_k=k, res=c_order2, dist = input$sil_dist, percen = upto_percen, cols = colors[1:k])
         check_sil <<- sil
         
         colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:k])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[sil[['k']]]) #sil[["sk3.col"]]
         
       } else if(input$sil_choice == 'Change point'){
         
         upto_sil_cp = list(input$sil_cp1, input$sil_cp2, input$sil_cp3, input$sil_cp4, input$sil_cp5, input$sil_cp6, input$sil_cp7, input$sil_cp8, input$sil_cp9)
         
         sil = sil_cpt2(data_use = mv_data2, opt_k=k, res=c_order2, dist = input$sil_dist, sil_cp = upto_sil_cp, cols = colors[1:k], cpt_data = input$changes, cpt_method = input$method, cpt_max = input$max)
         
         colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:k])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[sil[['k']]]) #sil[["sk3.col"]]
         
         
       }
   
       return(list(data= mv_data2, core = sil$core.samples, clust = sil$sk3, data = mv_data2, c_order = c_order2))
     }
       
     #colors = c("cyan", "khaki1", "pink", "plum3", "purple", "darkgreen", "hotpink","brown", "darkorchid2",  "maroon")
     #par(mfrow = c(1, 2))
     #plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
     #plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
     
     #return(list(core = sil$core.samples, clust = sil$sk3, order = consen()[["output"]]))
   })
   
   output$download_GW_Ex5 <- downloadHandler(
     filename= function() {paste('Core Sample Data set from Silhouette.csv')}, 
     content = function(file) {                           
       d <-  mv_hm_data()$data 
       write.csv(d, file, row.names = TRUE) }
   )
   
   output$download_GW_Ex6 <- downloadHandler(
     filename= function() {paste('Core Sample with clusters.csv')}, 
     content = function(file) {
       #d <-  core_mv_hm_data()$data
       #d2 = t(d[1, c(-1, -2)])
       d <-  data.frame(Cluster= sil_data()$c_order)
       d$Sample <- rownames(d)
       
       d2 <- d[,c(2, 1)]
       
       
       write.csv(d, file, row.names = FALSE) }
   )
   
   output$download_GW_Ex7 <- downloadHandler(
     filename= function() {paste('Core Sample with clusters.csv')}, 
     content = function(file) {
       
       d <-  data.frame(Cluster = sil_data()$c_order)
       d$Sample <- rownames(d)
       d2 <- d[,c(2, 1)]
       
       write.csv(d, file, row.names = FALSE) }
   )
   
   sil_plot <- eventReactive(input$button4, {
     if(input$gw_file5 == "GW_Example5") {
       
       data = consen()$data
       sample = colnames(consen()$data)
       check_sample <<- sample
       
       colors = mv_hm_data()$cc1
       check_colors <<- colors
       df = cbind.data.frame(sample, colors, stringsAsFactors = F)
       colnames(df) <- c("Sample","colors")
       
       upto_slider <- list()
       upto_percen <- list()
       upto_sil_cp <- list()
       
       check_df <<- df
       check_consen <<- consen()
       check_sil_choice <<- input$sil_choice
     
       if(input$sil_choice == 'Fixed value'){
         upto_slider = list(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, input$upto_slider9)
         check_upto_slider <<- upto_slider
         
         sil = silhouette_plot2(data_use= consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], upto_width= upto_slider, cols = df)
         
         colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
         
       } else if(input$sil_choice == 'Percentile') {
         
         upto_percen = list(input$sil_percen1, input$sil_percen2, input$sil_percen3, input$sil_percen4, input$sil_percen5, input$sil_percen6, input$upto_slider7, input$sil_percen8, input$sil_percen9)
         check_upto_percen <<- upto_percen
         
         sil = sil_function(data_use = consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], percen = upto_percen, cols = df)
         
         
         check_sil <<- sil
         
         colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
       
       }  else if(input$sil_choice == 'Change point'){
         
         upto_sil_cp = list(input$sil_cp1, input$sil_cp2, input$sil_cp3, input$sil_cp4, input$sil_cp5, input$sil_cp6, input$sil_cp7, input$sil_cp8, input$sil_cp9)
         
         sil = sil_cpt(data_use = consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], sil_cp = upto_sil_cp, cols = df, cpt_data = input$changes, cpt_method = input$method, cpt_max = input$max)
         
         colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = sil[["sk2.col"]])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = colors[1:sil[["k"]]]) #sil[["sk3.col"]]
         
         
       }
     } else if(input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3"){
       inFile6 <- input$gw_file6
       if (is.null(inFile6))
         return(NULL)
       else if(grepl(".csv", inFile6[1])) { mv_data = read.csv(as.character(inFile6$datapath), header = TRUE, sep = ",", stringsAsFactors = F, row.names = 1) }
       else if(grepl(".txt", inFile6[1])) { mv_data = read.table(as.character(inFile6$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, row.names = 1) }
       else if(grepl(".rds", inFile6[1])) { mv_data = readRDS(as.character(inFile6$datapath)) }
       
       inFile8 <- input$gw_file8
       if (is.null(inFile8))
         return(NULL)
       else if(grepl(".csv", inFile8[1])) { c_order = read.csv(as.character(inFile8$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
       else if(grepl(".txt", inFile8[1])) { c_order = read.table(as.character(inFile8$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
       else if(grepl(".rds", inFile8[1])) { c_order = readRDS(as.character(inFile8$datapath)) }
       
       mv_data2 <- as.data.frame(mv_data)
       c_order2 <- as.data.frame(c_order)
       names(c_order2)[2] <- "Sample" 
       
       check_mvdata <<- mv_data2
       check_corder2 <<- c_order2
       
       k = length(unique(c_order2$Cluster)) # unique(c_order2$Cluster)[length(unique(c_order2$Cluster))]
       
       if(input$sil_choice == 'Fixed value'){
         
         upto_slider = list(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, input$upto_slider9)
         check_upto_slider <<- upto_slider
         
         sil = silhouette_plot3(data_use= mv_data2, opt_k= k, res= c_order2, dist = input$sil_dist, upto_width= upto_slider, cols = c("orange", "darkblue", "black", "maroon", "violet", "plum2", "lightyellow", "lightblue", "green")[1:k])
      
         #colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "lightpink","steelblue", "darkorchid2",  "yellowgreen", "violetred")
        par(mfrow = c(1, 2))
        plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[1:k])
        plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[sil[['k']]]) #sil[["sk3.col"]]
       
       } else if(input$sil_choice == 'Percentile') {
         
         upto_percen = list(input$sil_percen1, input$sil_percen2, input$sil_percen3, input$sil_percen4, input$sil_percen5, input$sil_percen6, input$upto_slider7, input$sil_percen8, input$sil_percen9)
         check_upto_percen <<- upto_percen
         
         sil = sil_function2(data_use = mv_data2, opt_k=k, res=c_order2, dist = input$sil_dist, percen = upto_percen, cols = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[1:k])
         check_sil <<- sil
         
        #colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[1:k])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col =  c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[sil[['k']]]) #sil[["sk3.col"]]
         
       } else if(input$sil_choice == 'Change point'){
         
         upto_sil_cp = list(input$sil_cp1, input$sil_cp2, input$sil_cp3, input$sil_cp4, input$sil_cp5, input$sil_cp6, input$sil_cp7, input$sil_cp8, input$sil_cp9)
         
         sil = sil_cpt2(data_use = mv_data2, opt_k=k, res=c_order2, dist = input$sil_dist, sil_cp = upto_sil_cp, cols = c("orange", "darkblue", "black", "maroon", "violet", "plum2", "lightyellow", "lightblue", "green")[1:k], cpt_data = input$changes, cpt_method = input$method, cpt_max = input$max)
         
         #colors = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")
         par(mfrow = c(1, 2))
         plot(sil[["sk2"]],  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[1:k])
         plot(sil[["sk3"]],  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = c("cyan", "khaki1", "pink", "plum3", "mediumaquamarine", "coral", "hotpink","brown", "darkorchid2",  "maroon")[sil[['k']]]) #sil[["sk3.col"]]
         
         
       }
     }
   })
    
    output$gw_core_sam <- renderPlot({
      sil_plot()
    #  if(!is.null(sil_data()))
     # {
        #sil_data()
        
        
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
        con <- consensus_clustering(dinput=mv_data2, mK=as.integer(input$con_maxK), rep=as.integer(input$con_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$con_dist, iL=input$con_hclust, fL=input$con_hclust)
        #silhouette_plot2(data_use= consen()[["data"]], opt_k=as.integer(input$con_opt_k), res=consen()[["output"]], dist = consen()[["distance"]], upto_width= as.numeric(input$upto_slider))
        dev.off()
        file.copy(paste(pdf_file_con,'.pdf', sep='') ,file, overwrite=TRUE)
      }
    )
    
    
    core_mv_hm_data = eventReactive(input$button5, {
      #if(is.null(consen()))
      if(input$gw_file9 == "load_my_own_core")
      {
        #data <- NULL
        inFile10 <- input$gw_file10
        if (is.null(inFile10))
          return(NULL)
        else if(grepl(".csv", inFile10[1])) { data1 = read.csv(as.character(inFile10$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
        else if(grepl(".txt", inFile10[1])) { data1 = read.table(as.character(inFile10$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
        else if(grepl(".rds", inFile10[1])) { data1 = readRDS(as.character(inFile10$datapath)) }
        
        data1 <- data1[,order(data1[1, ])]
        data2 <- data1[order(data1[,2]),]
        
        check_d2 <<- data2
        
        ### row groups
        gene <- as.character(data2$gene_id)
        gene <- gene[(as.integer(input$DataR1)-1):nrow(data2)] #gene[(6-1):nrow(data)] 
        row.groups <- as.character(as.vector(data2[(as.integer(input$DataR1)-1): nrow(data2),2])) #as.character(as.vector(data[(6-1):nrow(data),2])) 
        row.groups.name <- names(table(row.groups))
        number.row.groups <- length(row.groups.name)
        
        ### column groups
        col.groups <- as.character(as.vector(data2[1,]))
        #col.groups <- col.groups[c(-1, -2)] # calculate no. of column groups
        col.groups <-  col.groups[as.numeric(input$DataC1):ncol(data2)] #col.groups[4:ncol(data)] 
        col.groups.name <- names(table(col.groups))
        number.col.groups <- length(col.groups.name)
        
        name_clust <<- paste(col.groups.name, collapse= ", ")
        count_clust <<- paste(table(col.groups), collapse = ", ")
        
        ## additional info on columns and rows
        colbar_data <-  as.matrix(t(data2[1:(as.numeric(input$DataR1)-2), c(1, (as.integer(input$DataC1)-1):ncol(data2))])) #t(data[1:(6-2), c(1, 3:ncol(data))]) 
        colnames(colbar_data) <- as.character(unlist(colbar_data[1,]))
        colbar_data <- as.matrix(colbar_data[c(-1, -2),])
        colnames(colbar_data)[1] <- "Groups"
        n.colbar_data <- ncol(colbar_data)
        
        #check_colbar_data <<- colbar_data
        
        rowbar_data <- as.matrix(data2[(as.integer(input$DataR1)-1):nrow(data2), 1:(as.numeric(input$DataC1)-1)]) #as.matrix(data[5:nrow(data), 1:3])
        rownames(rowbar_data) <- rowbar_data[,1]
        rowbar_data <- as.matrix(rowbar_data[,-1])
        #colnames(rowbar_data)[1] <- "Groups"
        n.rowbar_data <- ncol(rowbar_data)
        
        data3 <- data2[(as.numeric(input$DataR1)-1):nrow(data2), as.numeric(input$DataC1):ncol(data2)] #data[5:nrow(data), 4:ncol(data)]
        rownames(data3) <- gene
        data4 <-data3[complete.cases(data3),]
        data <- data.matrix(data4)
        check_d <<- data
        
      } else {
        
        if((input$gw_file5 == "GW_Example5" & input$gw_file7 == 'GW_Example7') & input$gw_file9 == "GW_Example9") {
        data5 <- extracted_data2()
       
        samples_to_include <- sil_data()$core #################
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
        
        name_clust <<- paste(col.groups.name, collapse= ", ")
        count_clust <<- paste(table(col.groups), collapse = ", ")
        
        ## additional info on columns and rows
        colbar_data <-  as.matrix(t(data3[1:2, 3:ncol(data3)])) #t(data[1:(6-2), c(1, 3:ncol(data))]) 
        colnames(colbar_data)[1] <- "Groups"
        n.colbar_data <- ncol(colbar_data)
        
        rowbar_data <- as.matrix(data3[c(-1,-2),2])
        colnames(rowbar_data)[1] <- "Groups"
        n.rowbar_data <- ncol(rowbar_data)
        
        data <- data3[(4-1):nrow(data3), 3:ncol(data3)] #data[5:nrow(data), 4:ncol(data)]
        rownames(data) <- gene
        data <-data[complete.cases(data),]
        data <- data.matrix(data)
        
    } else if((input$gw_file5 == "load_my_own_gw_subset2" & input$gw_file7 == "load_my_own_gw_subset3") & input$gw_file9 == "GW_Example9") {
          data5 = sil_data()$data                                                                                                                                                                                                                                                            
          samples_to_include <- sil_data()$core
          cc_clust = sil_data()$c_order
          
          names(cc_clust)[1] <- "Sample"
          
          check_data5 <<- data5
          check_samples_to_include <<- samples_to_include 
          check_cc_clust <<- cc_clust
          
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
        
        data2.2 <- data2[,order(data2[1,])]
        
        data3 <- data.frame(gene_id, Groups, data2.2)
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
        
        name_clust <<- paste(col.groups.name, collapse= ", ")
        count_clust <<- paste(table(col.groups), collapse = ", ")
        
        
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
          check_ccolbar_data <<- colbar_data
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
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
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
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
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
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
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
          cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- 'hotpink'
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
        check_z <<- z
        
        if(input$cc_dist == "pearson correlation") {
          if(input$cc_dispRow == "No" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,  key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "No" & input$cc_dispCol=='Yes' ) {
            cc_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,  key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize =1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='Yes') {
            cc_hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==8) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==9) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==10) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
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
          
          if(!is.null(number.colbar.class)) {
            if(number.colbar.class==1) {
              legend("right", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==2) {
              legend("right", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==3) {  
              legend("right", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==4) {  
              legend("right", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==5) {  
              legend("right", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==6) {  
              legend("right", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } 
          }
          
          
          if(!is.null(number.rowbar.class)) {
            if(number.rowbar.class==1) {
              legend("bottomright", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==2) {
              legend("bottomright", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==3) {  
              legend("bottomright", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==4) {  
              legend("bottomright", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==5) {  
              legend("bottomright", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==6) {  
              legend("bottomright", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } 
          }
          
        }  else {
          if(input$cc_dispRow == "No" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "No" & input$cc_dispCol=='Yes') {
            cc_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='No') {
            cc_hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
          } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='Yes') {
            cc_hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==8) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==9) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.col.groups==10) {  
            legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
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
          
          if(!is.null(number.colbar.class)) {
            if(number.colbar.class==1) {
              legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==2) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==3) {  
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==4) {  
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==5) {  
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.colbar.class==6) {  
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } 
          }
          
          
          if(!is.null(number.rowbar.class)) {
            if(number.rowbar.class==1) {
              legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==2) {
              legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==3) {  
              legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==4) {  
              legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==5) {  
              legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } else if(number.rowbar.class==6) {  
              legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
            } 
          }
        }
        
        eval(cc_hm$call)
        
        #hmdata <- list(data = data,data_use= data3, z = list(z), hm = list(cc_hm), col1=col1, cc1=cc1, cc2 = cc2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name, colbars2= colbars2, rowbars2= rowbars2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name, number.colbar.class = number.colbar.class, names.colbar.class = names.colbar.class, colors_used = colors_used, number.rowbar.class= number.rowbar.class, names.rowbar.class = names.rowbar.class,  colors_used2 = colors_used2)
        hm_data <- return(list(data = data, data_use= data3, z = list(z), hm = list(cc_hm), col1=col1, colbars2 = colbars2, rowbars2 = rowbars2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name, number.colbar.class = number.colbar.class, names.colbar.class = names.colbar.class, colors_used = colors_used, number.rowbar.class= number.rowbar.class, names.rowbar.class = names.rowbar.class,  colors_used2 = colors_used2 ))
        
    })
    
    
    tgplot2 <- function(z, col1, colbars2, rowbars2, number.col.groups, number.row.groups, row.groups.name , col.groups.name,  number.colbar.class, names.colbar.class, colors_used, number.rowbar.class, names.rowbar.class, colors_used2){
      if(input$cc_dist == "pearson correlation") {
        if(input$cc_dispRow == "No" & input$cc_dispCol=='No') {
          cc_hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,  key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$cc_dispRow == "No" & input$cc_dispCol=='Yes' ) {
          cc_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='No') {
          cc_hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,  key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='Yes') {
          cc_hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$cc_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
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
        
        if(!is.null(number.colbar.class)) {
          if(number.colbar.class==1) {
            legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==2) {
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==3) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==4) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==5) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==6) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
        }
        
        
        if(!is.null(number.rowbar.class)) {
          if(number.rowbar.class==1) {
            legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==2) {
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==3) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==4) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==5) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==6) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
        }
        
      }  else {
        
        check_number.col.groups <- number.col.groups
        check_number.row.groups <- number.row.groups
        
        check_number.colbar.class <- number.colbar.class
        check_number.rowbar.class <- number.rowbar.class
        
        if(input$cc_dispRow == "No" & input$cc_dispCol=='No') {
          cc_hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$cc_dispRow == "No" & input$cc_dispCol=='Yes') {
          cc_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='No') {
          cc_hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
        } else if(input$cc_dispRow == "Yes" & input$cc_dispCol=='Yes') {
          cc_hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$cc_clust_byrow))), Colv=eval(parse(text=paste(input$cc_clust_bycol))), hclust=function(c) {hclust(c,method=input$cc_hclust)}, distfun=function(c) {dist(c,method=input$cc_dist)},cexRow=input$cc_size1,cexCol =input$cc_size2,key=TRUE,keysize=1.0, margin = c(input$cc_inSlider2[1],input$cc_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 2, RowSideColorsSize = 1)
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
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==8) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==9) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.col.groups==10) {  
          legend("topright", legend = paste(c(col.groups.name[1], col.groups.name[2], col.groups.name[3], col.groups.name[4], col.groups.name[5], col.groups.name[6], col.groups.name[7], col.groups.name[8], col.groups.name[9], col.groups.name[10])), col = c("darkblue", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
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
        
        if(!is.null(number.colbar.class)) {
          if(number.colbar.class==1) {
            legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==2) {
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==3) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==4) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==5) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.colbar.class==6) {  
            legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
        }
        
        
        if(!is.null(number.rowbar.class)) {
          if(number.rowbar.class==1) {
            legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==2) {
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==3) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==4) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==5) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(number.rowbar.class==6) {  
            legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
        }
        
      }
    }
    
    output$cc_GW_subset_heatmap <- renderPlot({
      if(is.null(core_mv_hm_data())) {
        cc_hm <- NULL
      } else {    
        core_mv_hm_data()
        #tgplot2(z= core_mv_hm_data()$z[[1]], col1 = core_mv_hm_data()$col1, colbars2 = core_mv_hm_data()$colbars2, rowbars2  = core_mv_hm_data()$rowbars2, 
        #       number.col.groups =  core_mv_hm_data()$number.col.groups, number.row.groups = core_mv_hm_data()$number.row.groups, 
        #       row.groups.name = core_mv_hm_data()$row.groups.name, col.groups.name = core_mv_hm_data()$col.groups.name,
        #       number.colbar.class = core_mv_hm_data()$number.colbar.class, names.colbar.class = core_mv_hm_data()$names.colbar.class, 
        #       colors_used = core_mv_hm_data()$colors_used, 
        #       number.rowbar.class= core_mv_hm_data()$number.rowbar.class, names.rowbar.class = core_mv_hm_data()$names.rowbar.class,  
        #       colors_used2 = core_mv_hm_data()$colors_used2)
      }
    })
    
    cc_coldendo <- reactive({
      par(cex = input$cc_sizeClable)
      dend1 <- as.dendrogram(core_mv_hm_data()$hm[[1]]$colDendrogram)
      d <- data.frame(v1 =core_mv_hm_data()$hm[[1]]$colInd, v2=1:length(core_mv_hm_data()$hm[[1]]$colInd))
      
      m <- data.frame(v3 = 1:length(core_mv_hm_data()$colbars2[, 1]), v4 = core_mv_hm_data()$colbars2[, 1])
      
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
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==8) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==9) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8], core_mv_hm_data()$col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
      } else if(core_mv_hm_data()$number.col.groups==10) {  
        legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8], core_mv_hm_data()$col.groups.name[9], core_mv_hm_data()$col.groups.name[10])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
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
        data_l1_l2 <-  extracted_data2()
        check_d1 <<- data_l1_l2 #extracted_data2() #core_mv_hm_data()$data_use 
       
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
      mm <- data.frame(v3 = 1:length(core_mv_hm_data()$rowbars2[, 1]), v4 = core_mv_hm_data()$rowbars2[,1])
      
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
        cuttable2 <- as.data.frame(cutree(as.hclust(core_mv_hm_data()$hm[[1]]$rowDendrogram), k=as.integer(input$cc_cuttree2))[as.hclust(core_mv_hm_data()$hm[[1]]$rowDendrogram)$order])
        cuttable2 <- cbind.data.frame(rownames(cuttable2), cuttable2)
        names(cuttable2)[1] <- "gene_id"
        names(cuttable2)[2] <- "Cluster"
        check1 <<- cuttable2
        data_l1_l2_2 <- extracted_data2() #core_mv_hm_data()$data_use #extracted_data2() #
        data_l1_l2_2 <- data_l1_l2_2[-1, c(1,2)]
        check2 <<-  data_l1_l2_2
        m2 <- merge(cuttable2, data_l1_l2_2, by = "gene_id", sort= F)
        m2 <- m2[, c(1, 3, 2)]
        m2$order <- 1:nrow(m2)
        m3 <- m2[order(-m2$order),]
        rownames(m3) <- 1:nrow(m3)
        m4 <- m3[, c(1,2, 3)]
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
        
        #eval(cc_hm$call)
        #eval(mv_hm_data()$hm[[1]]$call) # call the heatmap here
        #eval(hm$call)
        tgplot2(z= core_mv_hm_data()$z[[1]], col1 = core_mv_hm_data()$col1, colbars2 = core_mv_hm_data()$colbars2, rowbars2  = core_mv_hm_data()$rowbars2, 
                number.col.groups =  core_mv_hm_data()$number.col.groups, number.row.groups = core_mv_hm_data()$number.row.groups, 
                row.groups.name = core_mv_hm_data()$row.groups.name, col.groups.name = core_mv_hm_data()$col.groups.name,
                number.colbar.class = core_mv_hm_data()$number.colbar.class, names.colbar.class = core_mv_hm_data()$names.colbar.class, 
                colors_used = core_mv_hm_data()$colors_used, 
                number.rowbar.class= core_mv_hm_data()$number.rowbar.class, names.rowbar.class = core_mv_hm_data()$names.rowbar.class,  
                colors_used2 = core_mv_hm_data()$colors_used2)
        
        if(input$cc_clust_byrow == TRUE){
        par(cex = 0.6*input$cc_sizeClable)
          par(cex = input$cc_sizeClable)
          dend1 <- as.dendrogram(core_mv_hm_data()$hm[[1]]$colDendrogram)
          d <- data.frame(v1 =core_mv_hm_data()$hm[[1]]$colInd, v2=1:length(core_mv_hm_data()$hm[[1]]$colInd))
          m <- data.frame(v3 = 1:length(core_mv_hm_data()$colbars2[,1]), v4 = core_mv_hm_data()$colbars2[,1])
          
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
            legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
          } else if(core_mv_hm_data()$number.col.groups==8) {  
            legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
          } else if(core_mv_hm_data()$number.col.groups==9) {  
            legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8], core_mv_hm_data()$col.groups.name[9])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
          } else if(core_mv_hm_data()$number.col.groups==10) {  
            legend("topright", legend = paste(c(core_mv_hm_data()$col.groups.name[1], core_mv_hm_data()$col.groups.name[2], core_mv_hm_data()$col.groups.name[3], core_mv_hm_data()$col.groups.name[4], core_mv_hm_data()$col.groups.name[5], core_mv_hm_data()$col.groups.name[6], core_mv_hm_data()$col.groups.name[7], core_mv_hm_data()$col.groups.name[8], core_mv_hm_data()$col.groups.name[9], core_mv_hm_data()$col.groups.name[10])), col = c("cyan", "khaki1", "pink", "plum3", "aquamarine", "coral", "hotpink", "steelblue", "yellowgreen", "violetred"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$cc_sizeClable)
          }
        }
        
        if(input$cc_clust_bycol == TRUE){
        par(cex = input$cc_sizeRlable)
          par(cex = input$cc_sizeRlable)
          dend2 <- as.dendrogram(core_mv_hm_data()$hm[[1]]$rowDendrogram)
          dd <- data.frame(v1 =rev(core_mv_hm_data()$hm[[1]]$rowInd), v2=1:length(core_mv_hm_data()$hm[[1]]$rowInd))
          mm <- data.frame(v3 = 1:length(core_mv_hm_data()$rowbars2[, 1]), v4 = core_mv_hm_data()$rowbars2[,1])
          
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
        }
        
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
    
    
    output$diagram <- renderDiagrammeR({
      mermaid(paste(
        paste0("graph ", "TB"),
        if(input$button1){
          paste(
            paste0("A[Input Data: <br/> <br/> No. of features/genes: ", input_gw_data()$n_genes, "<br/> No. of samples: ", input_gw_data()$n_samples, "<br/> No. of sample groups: ", input_gw_data()$numbercolgroups, "] --> B[Filtering based on: <br/> <br/>", as.character(paste(input$gw_subset, collapse = '<br/>')), "]"), 
            paste0("B --> C[Percentile cut-off: greater than P", as.character(scatter_plot_data()$percen), " <br/> Features selected: ", paste(n_sel, collapse= "AND <br/>"), "]"),
            sep = "\n", collapse = "" )
        }, if(input$button2){
          paste(
            paste0("C --> D[Heatmap Parameters: <br/> <br/> Normalized by: ", as.character(input$norm), "<br/> Scaling: ",as.character(input$norm2), "<br/> Distance: ", as.character(input$dist), "<br/> Clustering: ", as.character(input$hclust), "]"),
            sep = "\n", collapse = "" )
        }, if(input$button3){
          paste(
            paste0("D --> E[Consensus Clustering Parameters: <br/>  <br/> Distance: ", as.character(input$con_dist), "<br/> Clustering: ", as.character(input$con_hclust), "<br/> Cluster Algorithm: hierarchial clustering <br/> Number of Iterations: ",  as.character(input$con_pItems), "<br/> Proportion of items to sample = 0.8 <br/> Proportion of features to sample = 1]"),
            sep = "\n", collapse = "" )
        }, if(input$button4){
          paste(
            paste0("E --> F[Optimal Number of Clusters: <br/> <br/>", as.character(input$con_opt_k),"]"), 
            paste0("F --> G[Silhouette Plot: <br/> <br/> Clusters: ", as.character(name_clust), "<br/> Samples: ", as.character(count_clust), "<br/> Silhouette width cut-off: ", 
               if(input$sil_choice == 'Fixed value'){  
                   if(input$con_opt_k >= 1) {paste(input$upto_slider1, input$upto_slider2, sep = ",")} 
                   else if(input$con_opt_k >= 2) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3,  sep = ",")}
                   else if(input$con_opt_k >= 3) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, sep = ",")}
                   else if(input$con_opt_k >= 4) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, sep = ",")}
                   else if(input$con_opt_k >= 5) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, sep = ",")}
                   else if(input$con_opt_k >= 6) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, sep = ",")}
                   else if(input$con_opt_k >= 7) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, sep = ",")}
                   else if(input$con_opt_k >= 8) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, input$upto_slider9, sep = ",")}
                   else if(input$con_opt_k >= 9) {paste(input$upto_slider1, input$upto_slider2, input$upto_slider3, input$upto_slider4, input$upto_slider5, input$upto_slider6, input$upto_slider7, input$upto_slider8, input$upto_slider9, input$upto_slider10, sep = ",")}
              } else if(input$sil_choice == 'Percentile'){
                   if(input$con_opt_k >= 1) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), sep = ",")} 
                   else if(input$con_opt_k >= 2) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 3) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4,  sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 4) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4, sep = ""), paste("P",  input$sil_percen5, sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 5) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4, sep = ""), paste("P",  input$sil_percen5, sep = ""), paste("P", input$sil_percen6,  sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 6) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4, sep = ""), paste("P",  input$sil_percen5, sep = ""), paste("P", input$sil_percen6,  sep = ""), paste("P", input$sil_percen7, sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 7) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4, sep = ""), paste("P",  input$sil_percen5, sep = ""), paste("P", input$sil_percen6,  sep = ""), paste("P", input$sil_percen7, sep = ""), paste("P", input$sil_percen8, sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 8) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4, sep = ""), paste("P",  input$sil_percen5, sep = ""), paste("P", input$sil_percen6,  sep = ""), paste("P", input$sil_percen7, sep = ""), paste("P", input$sil_percen8, sep = ""), paste("P", input$sil_percen9,  sep = ""), sep = ",")}
                   else if(input$con_opt_k >= 9) {paste(paste("P", input$sil_percen1, sep = ""), paste("P", input$sil_percen2, sep = ""), paste("P", input$sil_percen3, sep = ""), paste("P", input$sil_percen4, sep = ""), paste("P",  input$sil_percen5, sep = ""), paste("P", input$sil_percen6,  sep = ""), paste("P", input$sil_percen7, sep = ""), paste("P", input$sil_percen8, sep = ""), paste("P", input$sil_percen9,  sep = ""), paste("P", input$sil_percen10,  sep = ""), sep = ",")}
              }  else if(input$sil_choice == 'Change point'){
                if(input$con_opt_k >= 1) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), sep = ",")} 
                else if(input$con_opt_k >= 2) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), sep = ",")}
                else if(input$con_opt_k >= 3) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP",  input$sil_cp4,  sep = ""), sep = ",")}
                else if(input$con_opt_k >= 4) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP", input$sil_cp4, sep = ""), paste("CP",  input$sil_cp5, sep = ""), sep = ",")}
                else if(input$con_opt_k >= 5) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP", input$sil_cp4, sep = ""), paste("CP",  input$sil_cp5, sep = ""), paste("CP", input$sil_cp6,  sep = ""), sep = ",")}
                else if(input$con_opt_k >= 6) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP", input$sil_cp4, sep = ""), paste("CP",  input$sil_cp5, sep = ""), paste("CP", input$sil_cp6,  sep = ""), paste("CP", input$sil_cp7, sep = ""), sep = ",")}
                else if(input$con_opt_k >= 7) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP", input$sil_cp4, sep = ""), paste("CP",  input$sil_cp5, sep = ""), paste("CP", input$sil_cp6,  sep = ""), paste("CP", input$sil_cp7, sep = ""), paste("CP", input$sil_cp8, sep = ""), sep = ",")}
                else if(input$con_opt_k >= 8) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP", input$sil_cp4, sep = ""), paste("CP",  input$sil_cp5, sep = ""), paste("CP", input$sil_cp6,  sep = ""), paste("CP", input$sil_cp7, sep = ""), paste("CP", input$sil_cp8, sep = ""), paste("CP", input$sil_cp9,  sep = ""), sep = ",")}
                else if(input$con_opt_k >= 9) {paste(paste("CP", input$sil_cp1, sep = ""), paste("CP", input$sil_cp2, sep = ""), paste("CP", input$sil_cp3, sep = ""), paste("CP", input$sil_cp4, sep = ""), paste("CP",  input$sil_cp5, sep = ""), paste("CP", input$sil_cp6,  sep = ""), paste("CP", input$sil_cp7, sep = ""), paste("CP", input$sil_cp8, sep = ""), paste("CP", input$sil_cp9,  sep = ""), paste("CP", input$sil_cp10,  sep = ""), sep = ",")}
              }
                   
                   , "]"),
            #if(input$con_opt_k >= 1) { paste(input$upto_slider1)} 
            sep = "\n", collapse = "" )
        }, if(input$button5){
          paste(
            paste0("G --> H[Updated Heatmap Parameters: <br/> <br/> Normalized by: ", as.character(input$cc_norm), "<br/> Scaling: ", input$cc_norm2, "<br/> Distance: ", as.character(input$cc_dist), "<br/> Clustering: ", as.character(input$cc_hclust), "]"),
            sep = "\n", collapse = "" )
        }, sep = "\n"))
     
      
    })
    
    
    
    
    

####################### Combine results cluster (CrC) Tab  #######################
  
  Exp_input <- reactive({
    if(input$Exp_file == 'Exp_example'){
      d <- read.csv("data/MV_TCGA_BRCA_Early_Late_OS_6.9years_IQR_99P_subset_data.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
    } else if(input$Exp_file == 'Exp_example1'){
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
      paste('TCGA_BRCA_Early_Late_Expression_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds2 <- read.csv("data/MV_TCGA_BRCA_Early_Late_OS_6.9years_IQR_99P_subset_data.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)

      write.csv(ds2, file, row.names = TRUE)
    }
  )
  
  output$Exp_download1 <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_Expression_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds2 <- read.csv("data/Most_variable_extracted_Expression.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
      
      write.csv(ds2, file, row.names = TRUE)
    }
  )
  
  Variant_input <- reactive({
    if(input$Variant_file == 'Variant_example'){
     # d2 <- read.csv("data/MV_TCGA_BRCA_Variant_proportion_25T_Early_Late_IMVA_95p.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
      d2 <- read.csv("data/MV_TCGA_BRCA_450K_methylation_tranf_betaValues_25_Tumors_Early_Late_NAs_rmvd_MAD_IQR_99p8_percentile_use.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)      
    } else if(input$Variant_file == 'Variant_example1') {
        d2 <- read.csv("data/Most_variable_extracted_Variant_transposed.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
    } else if(input$Variant_file == 'Variant_load_my_own'){
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
      paste('TCGA_BRCA_Methylation_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
     # ds3 <- Variant_input()
       #ds3 <- read.csv("data/MV_TCGA_BRCA_Variant_proportion_25T_Early_Late_IMVA_95p.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
       ds3 <- read.csv("data/MV_TCGA_BRCA_450K_methylation_tranf_betaValues_25_Tumors_Early_Late_NAs_rmvd_MAD_IQR_99p8_percentile_use.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)      
        
        write.csv(ds3, file, row.names = TRUE)
    }
  )
  
  output$Variant_download1 <- downloadHandler(
    filename <- function() {
      paste('CoMMpass_IA9b_560pt_Variant_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
       ds3 <- read.csv("data/Most_variable_extracted_Variant_transposed.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
       write.csv(ds3, file, row.names = TRUE)
    }
  )
  
  CNV_input <- reactive({
    if(input$CNV_file == 'CNV_example'){
    # d3 <- read.csv("data/MV_TCGA_BRCA_450K_methylation_tranf_betaValues_25_Tumors_Early_Late_NAs_rmvd_MAD_IQR_99p8_percentile.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
      d3 <- read.csv("data/MV_TCGA_BRCA_CNV_data_Var_96P.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
      } else if(input$CNV_file == 'CNV_example1'){
      d3 <- read.csv("data/Most_variable_extracted_CNV_updated.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
    } else if(input$CNV_file == 'CNV_load_my_own'){
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
      paste('TCGA_BRCA_CNV_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds2 <- read.csv("data/MV_TCGA_BRCA_CNV_data_Var_96P.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)
      
      write.csv(ds2, file, row.names = TRUE)
    }
  )
  
  output$CNV_download1 <- downloadHandler(
   filename <- function() {
      paste('CoMMpass_IA9b_560pt_CNV_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds2 <- read.csv("data/Most_variable_extracted_CNV_updated.csv", header =T, sep =",", stringsAsFactors = F, row.names = 1)

      write.csv(ds2, file, row.names = TRUE)
    }
  )
  
  indiv <- eventReactive(input$button6, {
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
  
  exp1 <- eventReactive(input$button6, {
    silhouette_plot(data_use= indiv()[["data"]], opt_k=as.integer(input$Exp_opt_k), res=indiv()[["output"]], dist = indiv()[["distance"]] )
  })
  
  output$Exp_sil <- renderPlot({
    exp1()
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
  
  indiv2 <- eventReactive(input$button7, {
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
  
  var1 <- eventReactive(input$button7, {
    silhouette_plot(data_use= indiv2()[["data"]], opt_k=as.integer(input$Variant_opt_k), res=indiv2()[["output"]], dist = indiv2()[["distance"]])
  })
  
  output$Variant_sil <- renderPlot({
    var1()
  })
  
  output$Variant_cc_dl <- downloadHandler(
    filename <- function(){
      pdf_file2 <<- paste(ifelse(input$dt_choice == "V", "Variant", "Methylation"), input$Variant_dist, input$Variant_hclust, sep = "_")
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
  
  indiv3 <- eventReactive(input$button8, {
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
  
  cnv1 <- eventReactive(input$button8, {
    silhouette_plot(data_use= indiv3()[["data"]], opt_k=as.integer(input$CNV_opt_k), res=indiv3()[["output"]], dist = indiv3()[["distance"]] )
  })
  
  output$CNV_sil <- renderPlot({
    cnv1()
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
  
  
  output$download_clinical <- downloadHandler(
    
    filename <- function() {
      paste('TCGA_BRCA_Clinical_file_example', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds4 <- read.csv("data/TCGA_BRCA_Clinical_file_example.csv", header =T, sep =",", stringsAsFactors = F)
      write.csv(ds4, file, row.names = F)
    }
  )
  
  output$download_clinical1 <- downloadHandler(
    
    filename <- function() {
      paste('CoMMpass_Clinical_file_example', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds4 <-  read.csv("data/Clinical_file_example.csv", header =T, sep =",", stringsAsFactors = F)
      write.csv(ds4, file, row.names = F)
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
        rownames(d4) <- c("Expression", ifelse(input$dt_choice== "V", "Variant", "Meth"), "CNV")
      }  else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) { 
        d4 <- cbind.data.frame(cc1[["output"]], cc2[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c("Expression", ifelse(input$dt_choice== "V", "Variant", "Meth"))
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
        d4 <- cbind.data.frame(cc1[["output"]], cc3[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c("Expression", "CNV")
        
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
        d4 <- cbind.data.frame(cc2[["output"]], cc3[["output"]])
        d4 <- t(d4)
        rownames(d4) <- c(ifelse(input$dt_choice== "V", "Variant", "Meth"), "CNV")
        
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
      paste('CrC_ds', Sys.time(),'.csv', sep='')
    },
    content <- function(file) {
      ds4 <- coca_input()
      write.csv(ds4, file, row.names = T)
    }
  )
  
  combined <- eventReactive(input$button9, {
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
      
      check_cc1 <<- cc1
      check_cc2 <<- cc2
      check_cc3 <<- cc3
      
      cc.final <- NULL
      
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        cc.final <- coca(cc= list(cc1, cc2, cc3), type = c("E", ifelse(input$dt_choice== "V", "V", "M"), "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", ifelse(input$dt_choice== "V", "V", "M")), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c(ifelse(input$dt_choice== "V", "V", "M"), "CNV"), opt_k= c(as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      }
    }
    
    consensus_clustering(dinput=cc.final[["data_in"]], mK=10, rep=as.integer(input$coca_pItems), pI=0.8, pF= 1, cAlg="hc", dist=input$coca_dist, iL=input$coca_hclust, fL=input$coca_hclust)
    
    return(list(output= cc.final[["output"]], data= cc.final[["data"]], distance = cc.final[["distance"]]))
    
  }, ignoreNULL = FALSE)
  
  
  
  output$coca_cc <- renderPlot({
    #coca_data <- coca_input()
    #coca_data2 <- as.matrix(coca_data)
    par(mfrow= c(1, 3))
    combined()
  })
  
  output$coca_cc_dl <- downloadHandler(
    filename <- function(){
      pdf_file4 <<- paste("CrC", input$coca_dist, input$coca_hclust, sep = "_")
      paste("ConsensusClustring_Results_", pdf_file4,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file=paste(pdf_file4,".pdf",sep=""))
      cc = coca_input()
      cc1 = list(output = cc[1,])
      cc2 = list(output = cc[2,])
      cc3 = list(output = cc[3,])
      
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        cc.final <- coca(cc= list(cc1, cc2, cc3), type = c("E", ifelse(input$dt_choice== "V", "V", "M"), "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", ifelse(input$dt_choice== "V", "V", "M")), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$Variant_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c("E", "CNV"), opt_k= c(as.integer(input$Exp_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in% input$coca_platform)) {
        cc.final <- coca(cc= list(cc1, cc2), type = c(ifelse(input$dt_choice== "V", "V", "M"), "CNV"), opt_k= c(as.integer(input$Variant_opt_k), as.integer(input$CNV_opt_k)), coca_reps = as.integer(input$coca_pItems), clust= input$coca_hclust, dist= input$coca_dist, coca_opt_k= as.integer(input$coca_opt_k))
      }
      silhouette_plot(data_use= combined()[["data"]], opt_k=as.integer(input$coca_opt_k), res=combined()[["output"]], dist = combined()[["distance"]] )
      dev.off()
      file.copy(paste(pdf_file4,'.pdf', sep='') ,file, overwrite=TRUE)
    })
  
  combined2 <- reactive({
    
    if(!is.null(combined()))
    {
      cc.final <- combined()
      check_cc.final <<- cc.final
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
      paste('CrC_clusters_', Sys.time(),'.csv', sep='')
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
      d5 <- read.csv("data/TCGA_BRCA_Clinical_file_example.csv", header =T, sep =",", stringsAsFactors = F)
      
    }else if(input$clinical == 'available1'){
      d5 <- read.csv("data/Clinical_file_example.csv", header =T, sep =",", stringsAsFactors = F)
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
        
          m2 <- m[match(colnames(data), m$Sample),]
          colbar_data <- m2[,c(1, (ncol(m2)+1-as.integer(input$cli_feature_no)):ncol(m2) )]
        
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
            number.colbar.class <- colbars(df2 = colbar_data)[[2]]
            names.colbar.class <- colbars(df2 = colbar_data)[[3]]
            colors_used <- colbars(df2 = colbar_data)[[4]]
          } else if(number.col.groups==2) {
            cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
            cc1 <- rep(col1[50], length(cell))
            cc1[1:table(col.groups)[[1]]] <- 'darkblue'
            cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- 'red'  
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            colbars2 <- if(input$cli_feature == FALSE) as.matrix(cc1) else as.matrix(cbind(cc1, colbars(df2 = colbar_data)[[1]]))
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
            
            if(input$cli_feature == TRUE){
            if(number.colbar.class==1) {
              legend("bottomright", legend = paste(names.colbar.class[1]), col = paste(colors_used[1]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
            }else if(number.colbar.class==2) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = paste(c(colors_used[1], colors_used[2])), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
            }else if(number.colbar.class==3) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = paste(c(colors_used[1], colors_used[2], colors_used[3])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
            }else if(number.colbar.class==4) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = paste(c(colors_used[1], colors_used[2], colors_used[3], colors_used[4])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
            } else if(number.colbar.class==5) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = paste(c(colors_used[1], colors_used[2], colors_used[3], colors_used[4], colors_used[5])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
            }
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
            
            if(input$cli_feature == TRUE){
            if(number.colbar.class==1) {
              legend("bottomright", legend = paste(names.colbar.class[1]), col = paste(colors_used[1]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
            }else if(number.colbar.class==2) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = paste(c(colors_used[1], colors_used[2])), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
            }else if(number.colbar.class==3) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = paste(c(colors_used[1], colors_used[2], colors_used[3])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
            }else if(number.colbar.class==4) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = paste(c(colors_used[1], colors_used[2], colors_used[3], colors_used[4])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
            } else if(number.colbar.class==5) {
              legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = paste(c(colors_used[1], colors_used[2], colors_used[3], colors_used[4], colors_used[5])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
            }
            }
          } 
          
          output$download_coca_HM <- downloadHandler(
            
            filename <- function() {
              pdf_file5 <<- paste("CrC", input$hclust_4, "clustering", input$dist_4, "distance", sep="_")
              paste('CrC_', pdf_file5, Sys.time(),'.pdf', sep='')
            },
            content <- function(file) {
              pdf(file=paste(pdf_file5,".pdf",sep="") , height= 10, width=10)
              plot.new()
              title("NOJAH: Clustering Analysis", sub= "CrC Analysis HeatMap",cex.main=1.2, col = "blue", font=3)
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
              
              if(input$cli_feature == TRUE){
                if(number.colbar.class==1) {
                  legend("bottomright", legend = paste(names.colbar.class[1]), col = paste(colors_used[1]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
                }else if(number.colbar.class==2) {
                  legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = paste(c(colors_used[1], colors_used[2])), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9, bty = "n")
                }else if(number.colbar.class==3) {
                  legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = paste(c(colors_used[1], colors_used[2], colors_used[3])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
                }else if(number.colbar.class==4) {
                  legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = paste(c(colors_used[1], colors_used[2], colors_used[3], colors_used[4])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
                } else if(number.colbar.class==5) {
                  legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = paste(c(colors_used[1], colors_used[2], colors_used[3], colors_used[4], colors_used[5])), lty= 1, lwd = 10, pt.cex = 1, cex= 0.9, bty = "n")
                }
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
    hs1.2 <- paste("Individual data types (such as, Expression, Methylation/Variant, CNV) must be clustered using the corresponding tabs or outside NOJAH (input data as load my own dataset) prior to performing Cluster of Cluster (CoC) analysis on this tab.")
    hs2 <- paste("Determination of number of clusters by Consensus Clustering")
    HTML(paste(hs1, hs1.2, h3(strong(hs2)), hs1, sep = '<br/>'))
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
    hs2 <- paste("CrC HeatMap")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text31 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Cluster Interpretation")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text32 <- renderUI({
    hs1 <- paste("&emsp;")
    hs2 <- paste("Contingency Table(s)")
    HTML(paste(h3(strong(hs2)), hs1, sep = '<br/>'))
  })
  
  output$com_text33 <- renderUI({
    hs1 <- paste("&emsp;")
    if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
    if(input$stratify_by == "Copy Number") { hs2 <- paste("Stratified by Copy Number") }
    else if(input$stratify_by == "Expression") { hs2 <- paste("Stratified by Expression") }
    else if(input$stratify_by == "Methylation/Variant") { hs2 <- paste("Stratified by Methylation/Variant") }
    HTML(paste(h4(hs2), hs1, sep = '<br/>'))
    }
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
                  type= c("Expression", ifelse(input$dt_choice== "V", "Variant", "Methylation"), "CNV"), dist= c(exp_dist, variant_dist, cnv_dist) )
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) {
        plotMeans(data= list(expression2, variant2), data_order=list(expression_order, variant_order) , 
                  type= c("Expression", ifelse(input$dt_choice== "V", "Variant", "Methylation")), dist= c(exp_dist, variant_dist))   
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
        plotMeans(data= list(expression2, cnv2), data_order=list(expression_order, cnv_order) , 
                  type= c("Expression", "CNV"), dist= c(exp_dist, cnv_dist))     
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
        plotMeans(data= list(variant2, cnv2), data_order=list(variant_order, cnv_order) , 
                  type= c(ifelse(input$dt_choice== "V", "Variant", "Methylation"), "CNV"), dist= c(variant_dist, cnv_dist)) 
      }
    } else {
      return(NULL)
    }
  })
  
  output$foo <- renderRHandsontable({
    if(!is.null(combined())){
     
      if("EXP" %in% input$coca_platform) {
        cc_clust1 <- as.matrix(indiv()[["output"]])
        cc_clust1 <- cbind(Sample = rownames(cc_clust1), cc_clust1)
        colnames(cc_clust1)[2] <- "Exp.Clust"
      }
      
      if("PROP" %in% input$coca_platform) {
        cc_clust2 <- as.matrix(indiv2()[["output"]])
        cc_clust2 <- cbind(Sample =  rownames(cc_clust2), cc_clust2)
        colnames(cc_clust2)[2] <- "Variant.Clust"
       }
      
      if("CNV" %in% input$coca_platform) { 
        cc_clust3 <- as.matrix(indiv3()[["output"]])
        cc_clust3 <- cbind(Sample = rownames(cc_clust3), cc_clust3)
        colnames(cc_clust3)[2] <- "CNV.Clust"
      }
      
      check_cc_clust1 <<- cc_clust1
      check_cc_clust2 <<- cc_clust2
      check_cc_clust3 <<- cc_clust3
      
      m0 <- merge(cc_clust1, cc_clust2, by = "Sample")
      m1 <- merge(m0, cc_clust3, by = "Sample")
      
      observed <- list()
      
      #plots
      if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform) {
        
        if(input$stratify_by == "Copy Number") {
        
        for(i in 1:length(unique(m1$CNV.Clust))) {
           observed[[i]] <- as.data.frame.matrix(table(m1$Exp.Clust, m1$Variant.Clust, m1$CNV.Clust, dnn = c("EXP", "VAR"))[, ,i])
          
       }
        check_obsCNV <<- observed
        obs <- cbind.data.frame(observed)
        rhandsontable(
          data.frame(obs),
          rowHeaders = paste("E", 1:length(unique(m1$Exp.Clust)), sep = ""),
          colHeaders = sort(apply(expand.grid(paste("CNV", 1:length(unique(m1$CNV.Clust)), sep = ""), paste(ifelse(input$dt_choice == "V", "V", "M"), 1:length(unique(m1$Variant.Clust)), sep = "")), 1, paste, collapse=""))
            #c(rep(paste("V", 1:length(unique(m1$Variant.Clust)), sep = ""), length(unique(m1$Variant.Clust))))
          ) %>%
          hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
      
       } else if(input$stratify_by == "Expression") {
         check_in_E <<- TRUE
        
        for(i in 1:length(unique(m1$Exp.Clust))) {
          observed[[i]] <- as.data.frame.matrix(table(m1$Exp.Clust, m1$Variant.Clust, m1$CNV.Clust, dnn = c("VAR", "CNV"))[i, , ])
         }
        check_obsE <<- observed
        obs <- cbind.data.frame(observed)
        rhandsontable(
          data.frame(obs),
          rowHeaders = paste(ifelse(input$dt_choice == "V", "V", "M"), 1:length(unique(m1$Variant.Clust)), sep = ""),
          colHeaders = sort(apply(expand.grid(paste("E", 1:length(unique(m1$Exp.Clust)), sep = ""), paste("CNV", 1:length(unique(m1$CNV.Clust)), sep = "")), 1, paste, collapse=""))
          #c(rep(paste("V", 1:length(unique(m1$Variant.Clust)), sep = ""), length(unique(m1$Variant.Clust))))
          
        ) %>%
          hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
        
       } else if(input$stratify_by == "Methylation/Variant") {
         
         for(i in 1:length(unique(m1$Variant.Clust))) {
           observed[[i]] <- as.data.frame.matrix(table(m1$Exp.Clust, m1$Variant.Clust, m1$CNV.Clust, dnn = c("E", "CNV"))[, i, ])
         }
         check_obsV <<- observed
         obs <- cbind.data.frame(observed)
         rhandsontable(
           data.frame(obs),
           rowHeaders = paste("E", 1:length(unique(m1$Exp.Clust)), sep = ""),
           colHeaders = sort(apply(expand.grid(paste(ifelse(input$dt_choice == "V", "V", "M"), 1:length(unique(m1$Variant.Clust)), sep = ""), paste("CNV", 1:length(unique(m1$CNV.Clust)), sep = "")), 1, paste, collapse=""))
           #c(rep(paste("V", 1:length(unique(m1$Variant.Clust)), sep = ""), length(unique(m1$Variant.Clust))))
           
         ) %>%
           hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
       }
        
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) {
            observed <- as.data.frame.matrix(table(m1$Exp.Clust, m1$Variant.Clust))
            
            rhandsontable(
              observed,
              rowHeaders = paste("E", 1:length(unique(m1$Exp.Clust)), sep = ""),
              colHeaders = paste(ifelse(input$dt_choice == "V", "V", "M"), 1:length(unique(m1$Variant.Clust)), sep = ""),
              readOnly = TRUE
            ) %>%
              hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE, width = 800, 
                               allowedTags = "<em><b><span><strong><a><big>")
            
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
           observed <- as.data.frame.matrix(table(m1$Exp.Clust, m1$CNV.Clust))
           
           rhandsontable(
             observed,
             rowHeaders = paste("E", 1:length(unique(m1$Exp.Clust)), sep = ""),
             colHeaders = paste("CNV", 1:length(unique(m1$CNV.Clust)), sep = ""),
             readOnly = TRUE
           ) %>%
             hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
  
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
           observed <- as.data.frame.matrix(table(m1$Variant.Clust, m1$CNV.Clust))
        
        rhandsontable(
          observed,
          rowHeaders = paste(ifelse(input$dt_choice == "V", "V", "M"), 1:length(unique(m1$Variant.Clust)), sep = ""),
          colHeaders = paste("CNV", 1:length(unique(m1$CNV.Clust)), sep = ""),
          readOnly = TRUE
        ) %>%
          hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
              
      }
      
      
    } else {
      return(NULL)
    }
  })
  
 
  output$dl_coca_inter<- downloadHandler(  
    filename <- function() {
      pdf_file6 <<- "CrC clustering Interpretation"
      paste('CrC_', pdf_file6, Sys.time(),'.pdf', sep='')
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
                  type= c("Expression", ifelse(input$dt_choice== "V", "Variant", "Methylation"), "CNV"), dist= c(exp_dist, variant_dist, cnv_dist) )
      } else if("EXP" %in% input$coca_platform & "PROP" %in% input$coca_platform & !("CNV" %in%input$coca_platform ) ) {
        plotMeans(data= list(expression2, variant2), data_order=list(expression_order, variant_order) , 
                  type= c("Expression", ifelse(input$dt_choice== "V", "Variant", "Methylation")), dist= c(exp_dist, variant_dist))   
      } else if("EXP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("PROP" %in%input$coca_platform )) { 
        plotMeans(data= list(expression2, cnv2), data_order=list(expression_order, cnv_order) , 
                  type= c("Expression", "CNV"), dist= c(exp_dist, cnv_dist))     
      } else if("PROP" %in% input$coca_platform & "CNV" %in% input$coca_platform & !("EXP" %in%input$coca_platform )) { 
        plotMeans(data= list(variant2, cnv2), data_order=list(variant_order, cnv_order) , 
                  type= c(ifelse(input$dt_choice== "V", "Variant", "Methylation"), "CNV"), dist= c(variant_dist, cnv_dist)) 
      }
      dev.off()
      file.copy(paste(pdf_file6,'.pdf', sep='') ,file, overwrite=TRUE)
    }
  )
  
  
  ################################ Significance tesing of clusters Tab ##############################3
  output$download_Sig_Ex1 <- downloadHandler( 
    filename= function() {paste('Example data set_TCGA_BRCA_subset_data.csv')}, 
    content = function(file) {
      d <- read.csv("data/TCGA_BRCA_Early_Late_OS_6.9years_IQR_99P_subset_data_SOC_use.csv", header = T, sep  = ",", stringsAsFactors = F)
      write.csv(d, file, row.names = FALSE) }
  )
  
  output$download_Sig_Ex2 <- downloadHandler( 
    filename= function() {paste('Example data set_CoMMpass_IA9_filtered_expression_ds.csv')}, 
    content = function(file) {
      d <- read.csv("data/Most_variable_extracted_Expression_withHRgroups_updated.csv", header = T, sep  = ",", stringsAsFactors = F)
      write.csv(d, file, row.names = FALSE) }
  )
  
  
  Sig_data_input <- reactive({
    if(input$Sig_file1 == 'Sig_Example2'){
      d <- read.csv("data/Most_variable_extracted_Expression_withHRgroups_updated.csv", header = T, sep  = ",", stringsAsFactors = F) 
    } else if(input$Sig_file1 == 'Sig_Example1') {
      d <- read.csv("data/TCGA_BRCA_Early_Late_OS_6.9years_IQR_99P_subset_data_SOC_use.csv", header = T, sep  = ",", stringsAsFactors = F)
    } else if(input$Sig_file1 == 'load_my_own_Sig'){
      inFile <- input$Sig_file2
      if (is.null(inFile))
        return(NULL)
      #else if(grepl(".xlsx", inFile[1])) { d = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F) }
      else if(grepl(".txt", inFile[1])) { d = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F) }
      else if(grepl(".rds", inFile[1])) { d = readRDS(as.character(inFile$datapath)) }
      
    }
    else 
      return(NULL)
    # dim(data)
    Dataset <- data.frame(d)
    return(Dataset)
  })
  
  # plot HM for Significance testing
  
  
  
  sig_hm_plot <- reactive({
    if(!is.null(Sig_data_input()))
    {
      data <- Sig_data_input()
      
      check <<- data
       
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
      ncg <<- number.col.groups
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
        colors_used <- if(n.colbar_data ==1) as.matrix(cc1) else colbars(df2 = colbar_data)[[4]]
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
        check_colbars2 <<- colbars2
        check_cc1 <<- cc1
        number.colbar.class <- if(n.colbar_data ==1) 1 else colbars(df2 = colbar_data)[[2]]
        names.colbar.class <- if(n.colbar_data ==1) 1 else colbars(df2 = colbar_data)[[3]]
        colors_used <- if(n.colbar_data ==1) cc1 else colbars(df2 = colbar_data)[[4]]
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
          Sig_hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,  key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "No" & input$Sig_dispCol=='Yes' ) {
          Sig_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='No') {
          Sig_hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,  key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='Yes') {
          Sig_hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
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
        
        if(!is.null(number.colbar.class)) {
         if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==3) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==4) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==5) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==6) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } 
        }
        
        
      }  else {
        if(input$Sig_dispRow == "No" & input$Sig_dispCol=='No') {
          Sig_hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "No" & input$Sig_dispCol=='Yes') {
          Sig_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='No') {
          Sig_hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
        } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='Yes') {
          Sig_hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 1, RowSideColorsSize = 1)
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
         k_ncc <<- number.colbar.class
        
         if(!is.null(number.colbar.class)) { 
        if(number.colbar.class==1) {
          legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==2) {
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==3) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==4) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==5) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } else if(number.colbar.class==6) {  
          legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
        } 
         }
      }
      
      hm_data <- return(list(data = data, z = list(z), hm = list(Sig_hm ), col1=col1, colbars2 = colbars2, rowbars2 = rowbars2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name, number.colbar.class = number.colbar.class, names.colbar.class = names.colbar.class, colors_used = colors_used, number.rowbar.class= number.rowbar.class, names.rowbar.class = names.rowbar.class,  colors_used2 = colors_used2 ))
    } 
    else {
      return(NULL)
    }
    
  })
  
  tgplot3 <- function(z, col1, colbars2, rowbars2, number.col.groups, number.row.groups, row.groups.name , col.groups.name, number.colbar.class, names.colbar.class, colors_used, number.rowbar.class, names.rowbar.class,  colors_used2 ) { 
    if(input$Sig_dist == "pearson correlation") {
      if(input$Sig_dispRow == "No" & input$Sig_dispCol=='No') {
        Sig_hm <- heatmap.3(z[[1]], labRow = NA, labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,  key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
      } else if(input$Sig_dispRow == "No" & input$Sig_dispCol=='Yes' ) {
        Sig_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
      } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='No') {
        Sig_hm <- heatmap.3(z[[1]], labCol= NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,  key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
      } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='Yes') {
        Sig_hm <- heatmap.3(z[[1]], Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), scale="none", hclust=function(x) hclust(x,method=input$Sig_hclust), distfun=function(x) as.dist((1-cor(t(x)))),cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
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
      
      if(!is.null(number.colbar.class)) {
      if(number.colbar.class==1) {
        legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==2) {
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==3) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==4) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==5) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==6) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } 
      }
      
       
      if(!is.null(number.rowbar.class)) {
      if(number.rowbar.class==1) {
        legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==2) {
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==3) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==4) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==5) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==6) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } 
      }
      
    }  else {
      if(input$Sig_dispRow == "No" & input$Sig_dispCol=='No') {
        Sig_hm <- heatmap.3(z[[1]], labRow=NA, labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
      } else if(input$Sig_dispRow == "No" & input$Sig_dispCol=='Yes') {
        Sig_hm <- heatmap.3(z[[1]], labRow=NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
      } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='No') {
        Sig_hm <- heatmap.3(z[[1]], labCol = NA, scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
      } else if(input$Sig_dispRow == "Yes" & input$Sig_dispCol=='Yes') {
        Sig_hm <- heatmap.3(z[[1]], scale="none", Rowv=eval(parse(text=paste(input$Sig_clust_byrow))), Colv=eval(parse(text=paste(input$Sig_clust_bycol))), hclust=function(c) {hclust(c,method=input$Sig_hclust)}, distfun=function(c) {dist(c,method=input$Sig_dist)},cexRow=input$Sig_size1,cexCol =input$Sig_size2,key=TRUE,keysize=1.0, margin = c(input$Sig_inSlider2[1],input$Sig_inSlider2[2]), density.info=c("none"),trace=c("none"),col=col1,ColSideColors=colbars2, RowSideColors = rowbars2, ColSideColorsSize = 3, RowSideColorsSize = 1)
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
      
     
      if(!is.null(number.colbar.class)) {
      if(number.colbar.class==1) {
        legend("bottomright", legend = paste(names.colbar.class), col = colors_used[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==2) {
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2])), col = c(colors_used[1], colors_used[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==3) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3])), col = c(colors_used[1], colors_used[2], colors_used[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==4) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==5) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.colbar.class==6) {  
        legend("bottomright", legend = paste(c(names.colbar.class[1], names.colbar.class[2], names.colbar.class[3], names.colbar.class[4], names.colbar.class[5], names.colbar.class[6])), col = c(colors_used[1], colors_used[2], colors_used[3],colors_used[4], colors_used[5], colors_used[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } 
      }
      
      if(!is.null(number.rowbar.class)) {
      if(number.rowbar.class== 1) {
        legend("right", legend = paste(names.rowbar.class), col = colors_used2[1], lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==2) {
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2])), col = c(colors_used2[1], colors_used2[2]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==3) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3])), col = c(colors_used2[1], colors_used2[2], colors_used2[3]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==4) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==5) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } else if(number.rowbar.class==6) {  
        legend("right", legend = paste(c(names.rowbar.class[1], names.rowbar.class[2], names.rowbar.class[3], names.rowbar.class[4], names.rowbar.class[5], names.rowbar.class[6])), col = c(colors_used2[1], colors_used2[2], colors_used2[3],colors_used2[4], colors_used2[5], colors_used2[6]), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
      } 
      }
      
    }
    
    #hm_data <- return(list(data = data, z = list(z), hm = list(Sig_hm ), col1=col1, cc1= colbars2, cc2 = rowbars2, number.col.groups = number.col.groups,number.row.groups=number.row.groups, row.groups.name=row.groups.name, col.groups.name= col.groups.name))
  } 
  

  sig_hmplot <- eventReactive(input$button10, {
     if(is.null(sig_hm_plot())) {
        Sig_hm <- NULL
     } else {    #sig_hm_plot
       tgplot3(z= sig_hm_plot()$z[[1]], col1 = sig_hm_plot()$col1, colbars2 = sig_hm_plot()$colbars2, rowbars2 = sig_hm_plot()$rowbars2, 
           number.col.groups = sig_hm_plot()$number.col.groups, number.row.groups = sig_hm_plot()$number.row.groups, 
           row.groups.name = sig_hm_plot()$row.groups.name, col.groups.name = sig_hm_plot()$col.groups.name, 
           number.colbar.class = sig_hm_plot()$number.colbar.class, names.colbar.class = sig_hm_plot()$names.colbar.class, 
           colors_used = sig_hm_plot()$colors_used, 
           number.rowbar.class= sig_hm_plot()$number.rowbar.class, names.rowbar.class = sig_hm_plot()$names.rowbar.class,  
           colors_used2 = sig_hm_plot()$colors_used2 )
     }
   })

  
  output$Sig_plot <- renderPlot({
    #input$button10
    sig_hmplot()
  })
  
  
  Sig_coldendo <- eventReactive(input$button10, {
        par(cex = input$Sig_sizeClable)
        dend1 <- as.dendrogram(sig_hm_plot()$hm[[1]]$colDendrogram)
        d <- data.frame(v1 =sig_hm_plot()$hm[[1]]$colInd, v2=1:length(sig_hm_plot()$hm[[1]]$colInd))
        m <- data.frame(v3 = 1:length(sig_hm_plot()$colbars2), v4 = sig_hm_plot()$colbars2)
        
        colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
        colbar <- colbar[,2]
        labels_colors(dend1) <- as.character(colbar)
        plot(dend1)
        
        if(sig_hm_plot()$number.col.groups==1) {
          legend("topright", legend = paste(sig_hm_plot()$col.groups.name), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(sig_hm_plot()$number.col.groups==2) {
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(sig_hm_plot()$number.col.groups==3) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(sig_hm_plot()$number.col.groups==4) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(sig_hm_plot()$number.col.groups==5) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeClable)
        } else if(sig_hm_plot()$number.col.groups==6) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(sig_hm_plot()$number.col.groups==7) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(sig_hm_plot()$number.col.groups==8) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7], sig_hm_plot()$col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        } else if(sig_hm_plot()$number.col.groups==9) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7], sig_hm_plot()$col.groups.name[8], sig_hm_plot()$col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable9)
        } else if(sig_hm_plot()$number.col.groups==10) {  
          legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7], sig_hm_plot()$col.groups.name[8], sig_hm_plot()$col.groups.name[9], sig_hm_plot()$col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$sizeClable)
        }
      })
      
      output$Sig_plot1 <- renderPlot({
        Sig_coldendo()
      })
      
      
    Sig_colDen <- eventReactive(input$button10, {
        if(input$Sig_cutcolden == 'TRUE') {
          cuttable <- as.data.frame(cutree(as.hclust(sig_hm_plot()$hm[[1]]$colDendrogram), k=as.numeric(input$Sig_cuttree))[as.hclust(sig_hm_plot()$hm[[1]]$colDendrogram)$order])
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
        if(input$Sig_cutcolden == 'TRUE' & sig_hm_plot()$number.col.groups >= 2 ){
          HTML(paste("<br/>", br(strong(em(paste("Would you want to assess gene set significance in the separation of specimens into two clusters? (Yes/No)")))), sep = "<br>")) 
        }
        else 
          return(NULL)
      })
      
      
      output$Sig_pvalue <- renderUI({
        input$Sig_goButton 
        
        Sig_b1 <<- NULL
        
        isolate(
          
          if(input$Sig_cutcolden == 'TRUE') {
            pobs.col <- numeric()
            perms.col <- numeric()
            hc.cols <- as.hclust(sig_hm_plot()$hm[[1]]$colDendrogram)
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
                s_data <- read.csv("data/GW_TCGA_BRCA_Early_Late_OS_6.9years_Core_samples_for_SoC_use.csv", header = T, sep = ",", stringsAsFactors = F)
              } else if(input$Sig_file3 == 'Sig_Exp.Example1') {
                s_data <- readRDS("data/CoMMpassIA9_GW_Expression_data_SoC.rds")
                #s_data <- readRDS("data/Meth27K.GW.BRCA.Example.data.rds")
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
              
              
              check_obsdata <<- mer
              check_samplingdata <<- s_data
              check_distmethod <<- input$Sig_dist
              check_clustmethod <<- input$Sig_hclust
              check_norm <<- input$Sig_norm 
              check_scale <<- input$Sig_norm2 
              check_n <<- as.numeric(input$Sig_n) 
              check_k <<- as.numeric(input$Sig_cuttree)
              check_n.iter <<- input$Sig_n_iter 
              check_zlim <<- c(input$Sig_inSlider[1],input$Sig_inSlider[2]) 
              check_sampler <<-  "Column"
              
              
              Sig_b1 <<- bootstrapfun(obsdata=mer, samplingdata=s_data, distmethod = input$Sig_dist, clustmethod= input$Sig_hclust, norm= input$Sig_norm, scale=input$Sig_norm2, n=as.numeric(input$Sig_n), k=as.numeric(input$Sig_cuttree), n.iter=input$Sig_n_iter, zlim=c(input$Sig_inSlider[1],input$Sig_inSlider[2]), sampler = "Column", updateProgress )
              
              hstring1 <- paste("&emsp;")
              hstring2 <- paste("The p-value to test the gene set significance in the separation of specimens into 2 clusters is =", Sig_b1$p.value, sep = " ")
              if(Sig_b1$p.value <= 0.05) {
                hstring3  <- paste("The gene set of interest is able to separate sample groups, outperforming ", input$Sig_n_iter, " random samples of genes of the same number in this regard.", sep = "")
              } else {
                hstring3  <- paste(input$Sig_n_iter, " random samples of genes of the same number are able to separate sample groups, outperforming the gene set of interest in this regard.", sep = "")
              }
              
              HTML(paste(hstring1, h5(strong(hstring2)), h5(em(hstring3)), hstring1, hstring1,  sep = '<br/>'))
              
              
            }
          }
          
        )
        
        
      })
      
      Sig_rowdendo <- eventReactive(input$button10, {
        par(cex = input$Sig_sizeRlable) 
        dend2 <- as.dendrogram(sig_hm_plot()$hm[[1]]$rowDendrogram)
        dd <- data.frame(v1 =rev(sig_hm_plot()$hm[[1]]$rowInd), v2=1:length(sig_hm_plot()$hm[[1]]$rowInd))
        mm <- data.frame(v3 = 1:length(sig_hm_plot()$rowbars2), v4 = sig_hm_plot()$rowbars2)
        
        colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
        colbar2 <- colbar2[,2]
        labels_colors(dend2) <- rev(as.character(colbar2))
        plot(dend2, horiz = T)
        if(sig_hm_plot()$number.row.groups==1) {
          legend("topright", legend = paste(sig_hm_plot()$row.groups.name[1]), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(sig_hm_plot()$number.row.groups==2) {
          legend("topright", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(sig_hm_plot()$number.row.groups==3) {  
          legend("topright", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(sig_hm_plot()$number.row.groups==4) {  
          legend("topright", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3], sig_hm_plot()$row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(sig_hm_plot()$number.row.groups==5) {  
          legend("topright", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3], sig_hm_plot()$row.groups.name[4], sig_hm_plot()$row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } else if(sig_hm_plot()$number.row.groups==6) {  
          legend("topright", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3], sig_hm_plot()$row.groups.name[4], sig_hm_plot()$row.groups.name[5], sig_hm_plot()$row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 2*input$Sig_sizeRlable)
        } 
      })
      
      output$Sig_plot2 <- renderPlot({
        par(cex= input$Sig_sizeRlable)
        Sig_rowdendo()
      })
      
      Sig_rowDen <- reactive({
        if(input$Sig_cutrowden == 'TRUE') {
          cuttable2 <- as.data.frame(cutree(as.hclust(sig_hm_plot()$hm[[1]]$rowDendrogram), k=as.integer(input$Sig_cuttree2))[as.hclust(sig_hm_plot()$hm[[1]]$rowDendrogram)$order])
          cuttable2 <- cbind.data.frame(rownames(cuttable2), cuttable2)
          names(cuttable2)[1] <- "gene_id"
          names(cuttable2)[2] <- "Cluster"
          data_l1_l2_2 <- Sig_data_input()
          data_l1_l2_2 <- data_l1_l2_2[-1, c(1,2)]
          m2 <- merge(cuttable2, data_l1_l2_2, by = "gene_id", sort= F)
          m2 <- m2[, c(1, 3, 2)]
          m2$order <- 1:nrow(m2)
          m3 <- m2[order(-m2$order),]
          rownames(m3) <- 1:nrow(m3)
          m4 <- m3[, c(1,2, 3)]
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
        if(input$Sig_cutrowden == 'TRUE' & sig_hm_plot()$number.row.groups >= 2 ){
          HTML(paste("<br/>", paste("Would you want to assess significance of patients in the separation of genes into two clusters? (Yes/No)"), sep = "<br>")) 
        }
        else 
          return(NULL)
      })
      
      
      
      output$Sig_pvalue2 <- renderUI ({
        input$Sig_goButton2 
        sig_b2 <<- NULL
        
        isolate(
          
          if(input$Sig_cutrowden== TRUE){
            m2 <- Sig_rowDen()
            
            if(input$Sig_pvalue_cal2 == 'TRUE') 
            {
              if(input$Sig_file5 == 'Sig_Exp.Example2') {
                s_data <- read.csv("data/GW_TCGA_BRCA_Early_Late_OS_6.9years_Core_samples_for_SoC_use.csv", header = T, sep = ",", stringsAsFactors = F)
              } else if(input$Sig_file5 == 'Sig_Exp.Example21') {
                s_data <- readRDS("data/CoMMpassIA9_GW_Expression_data_SoC.rds")
                #s_data <- readRDS("data/Meth27K.GW.BRCA.Example.data.rds")
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
                hstring3  <- paste("The sample set of interest is able to separate gene groups, outperforming ", input$Sig_n_iter2, " random samples of the same number in this regard.", sep = "")
              } else {
                hstring3  <- paste(input$Sig_n_iter2, " random samples of genes of the same number are able to separate sample groups, outperforming the gene set of interest in this regard.", sep = "")
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
                                 c("Scale", ifelse(input$Sig_norm == "none", paste(as.integer(min(sig_hm_plot()$data)), as.integer(max(sig_hm_plot()$data)), sep = ":"), paste(input$Sig_inSlider[1], input$Sig_inSlider[2], sep=":"))),
                                 c("HeatMap colors", paste(input$Sig_low, input$Sig_mid, input$Sig_high, sep="-")))
          names(df)[1] <- "Parameters"
          names(df)[2] <- "Value Selected"
          grid.table(df, rows= NULL)
          
          #eval(Sig_hm$call) # call the heatmap here
          #eval(sig_hm_plot()$hm$call)
          
          tgplot3(z= sig_hm_plot()$z[[1]], col1 = sig_hm_plot()$col1, colbars2 = sig_hm_plot()$colbars2, rowbars2 = sig_hm_plot()$rowbars2, 
                  number.col.groups = sig_hm_plot()$number.col.groups, number.row.groups = sig_hm_plot()$number.row.groups, 
                  row.groups.name = sig_hm_plot()$row.groups.name, col.groups.name = sig_hm_plot()$col.groups.name, 
                  number.colbar.class = sig_hm_plot()$number.colbar.class, names.colbar.class = sig_hm_plot()$names.colbar.class, 
                  colors_used = sig_hm_plot()$colors_used, 
                  number.rowbar.class= sig_hm_plot()$number.rowbar.class, names.rowbar.class = sig_hm_plot()$names.rowbar.class,  
                  colors_used2 = sig_hm_plot()$colors_used2 )
          
          if(input$Sig_clust_bycol == "TRUE") {
          par(cex = 0.6*input$Sig_sizeClable)
          #par(mar=c(5,7,4,2))
          dend1 <- as.dendrogram(sig_hm_plot()$hm[[1]]$colDendrogram)
          d <- data.frame(v1 =sig_hm_plot()$hm[[1]]$colInd, v2=1:length(sig_hm_plot()$hm[[1]]$colInd))
          m <- data.frame(v3 = 1:length(sig_hm_plot()$colbars2), v4 = sig_hm_plot()$colbars2)
          
          colbar <- data.frame(v1=d$v2, v4=m[match(d$v1, m$v3), 2])
          colbar <- colbar[,2]
          labels_colors(dend1) <- as.character(colbar)
          plot(dend1, main="Column Dendrogram")
          if(sig_hm_plot()$number.col.groups==1) {
            legend("topright", legend = paste(sig_hm_plot()$col.groups.name[1]), col = "firebrick4", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==2) {
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2])), col = c("firebrick4", "dodgerblue"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==3) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3])), col = c("firebrick4", "dodgerblue", "khaki1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==4) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4])), col = c("firebrick4", "dodgerblue", "khaki1", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==5) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==6) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==7) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==8) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7], sig_hm_plot()$col.groups.name[8])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==9) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7], sig_hm_plot()$col.groups.name[8], sig_hm_plot()$col.groups.name[9])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.col.groups==10) {  
            legend("topright", legend = paste(c(sig_hm_plot()$col.groups.name[1], sig_hm_plot()$col.groups.name[2], sig_hm_plot()$col.groups.name[3], sig_hm_plot()$col.groups.name[4], sig_hm_plot()$col.groups.name[5], sig_hm_plot()$col.groups.name[6], sig_hm_plot()$col.groups.name[7], sig_hm_plot()$col.groups.name[8], sig_hm_plot()$col.groups.name[9], sig_hm_plot()$col.groups.name[10])), col = c("firebrick4", "dodgerblue", "khaki1", "gray", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          }
          }
          
          if(input$Sig_clust_byrow == "TRUE") {
          par(cex = input$Sig_sizeRlable)
          dend2 <- as.dendrogram(sig_hm_plot()$hm[[1]]$rowDendrogram)
          dd <- data.frame(v1 =sig_hm_plot()$hm[[1]]$rowInd, v2=1:length(sig_hm_plot()$hm[[1]]$rowInd))
          mm <- data.frame(v3 = 1:length(sig_hm_plot()$rowbars2), v4 = sig_hm_plot()$rowbars2)
          
          colbar2 <- data.frame(v1=dd$v2, v4=mm[match(dd$v1, mm$v3), 2])
          colbar2 <- colbar2[,2]
          labels_colors(dend2) <- as.character(colbar2)
          plot(dend2, horiz = T, main="Row Dendrogram")
          if(sig_hm_plot()$number.row.groups==1) {
            legend("bottomleft", legend = paste(sig_hm_plot()$row.groups.name), col = "black", lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.row.groups==2) {
            legend("bottomleft", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2])), col = c("black", "gray"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.row.groups==3) {  
            legend("bottomleft", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3])), col = c("black", "gray", "hotpink"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.row.groups==4) {  
            legend("bottomleft", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3], sig_hm_plot()$row.groups.name[4])), col = c("black", "gray", "hotpink", "brown1"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.row.groups==5) {  
            legend("bottomleft", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3], sig_hm_plot()$row.groups.name[4], sig_hm_plot()$row.groups.name[5])), col = c("black", "gray", "hotpink", "brown1", "cyan"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } else if(sig_hm_plot()$number.row.groups==6) {  
            legend("bottomleft", legend = paste(c(sig_hm_plot()$row.groups.name[1], sig_hm_plot()$row.groups.name[2], sig_hm_plot()$row.groups.name[3], sig_hm_plot()$row.groups.name[4], sig_hm_plot()$row.groups.name[5], sig_hm_plot()$row.groups.name[6])), col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"), lty= 1, lwd = 10, pt.cex = 1, cex = 0.9)
          } 
          }
          
          if(input$Sig_goButton) {
            plot.new()
            
            title("Significance of Column Cluster", cex.main=2)
            df2 <- rbind.data.frame(c("Fisher's Exact Test p-value ", ifelse(Sig_b1$p.obs < 0.001, "<0.001", ifelse(Sig_b1$p.obs < 0.01, "<0.01", round(Sig_b1$p.obs, 3)))),
                                    c("Gene-Set size", input$Sig_n),
                                    c("Number of bootstrap samples", input$Sig_n_iter),
                                    c("Monte Carlo p-value", ifelse(Sig_b1$p.value < 0.001, "<0.001", ifelse(Sig_b1$p.value < 0.01, "<0.01", round(Sig_b1$p.value, 3)))),#Sig_b1$p.value ), 
                                    c("Interpretation", ifelse(Sig_b1$p.value <= 0.05, paste("The gene set of interest is able to", "separate sample groups,", paste("outperforming ", input$Sig_n_iter, " random samples of genes", sep = ""), "of the same number in this regard.", sep = '\n'), paste(paste(input$Sig_n_iter, " random samples of genes", sep = ""), "of the same number are able", "to separate sample groups outperforming", "the gene set of interest in this regard.", sep = '\n') ))
                                    
            )
            names(df2)[1] <- "Summary"
            names(df2)[2] <- "Value"
            grid.table(df2, rows= NULL)
          }
          
          if(input$Sig_goButton2) {
            plot.new()
            title("Significance of Row Cluster", cex.main=2)
            df3 <- rbind.data.frame(c("Fisher's Exact Test p-value", ifelse(Sig_b2$p.obs < 0.001, "<0.001", ifelse(Sig_b2$p.obs < 0.01, "<0.01", round(Sig_b2$p.obs, 3)))),
                                    c("Sample-set size", input$Sig_n2),
                                    c("Number of bootstrap samples", input$Sig_n_iter2),
                                    c("Monte Carlo p-value", ifelse(Sig_b2$p.value < 0.001, "<0.001", ifelse(Sig_b2$p.value < 0.01, "<0.01", round(Sig_b2$p.value, 3)))),
                                    c("Interpretation", ifelse(Sig_b2$p.value <= 0.05, paste("The sample set of interest is able to", "separate gene groups,", paste("outperforming", input$Sig_n_iter2, " random genes of the ", sep = ""), "same number in this regard.", sep = '\n'), paste(paste(input$Sig_n_iter2, " random samples of the same", sep = ""), "number are able ", "to separate gene groups outperforming", "the sample groups of interest in this regard.", sep = '\n') ))
                                   
            )
            names(df3)[1] <- "Summary"
            names(df3)[2] <- "Value"
            grid.table(df3, rows= NULL)
          }
          dev.off()
          
          file.copy(paste(Sig_pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
        })
  
  
}


