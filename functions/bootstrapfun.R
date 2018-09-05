bootstrapfun <- function(obsdata, samplingdata, distmethod, clustmethod, norm, scale, n, k, n.iter, zlim, sampler, updateProgress = NULL) #n is the sample size for the expected data
{
  classfn.count.obs <- integer()
  names(obsdata)[1] <- "Sample"
  obsdata$Sample <- gsub("[[:punct:]]", "-", obsdata$Sample)
  
  # estimate no. of col groups for obs
  classfn.obs <- as.character(unlist(obsdata$Group)) 
  classfn.count.obs <- length(table(unique(classfn.obs)))
  classfn.names.obs <- names(table(classfn.obs))
  
  if(classfn.count.obs < 2){
    stop("Need atleast two groups to calculate p-values using contigency table analysis")
  } else if(classfn.count.obs > 1) {
    obsdata$Sample <- ifelse(grepl(classfn.names.obs[1], obsdata$Group) , paste(obsdata$Sample, ".G1", sep=""), paste(obsdata$Sample, ".G2", sep=""))
    
    obstable <- contingencyfun(obsdata)
    pobs <- pvaluefunc(obstable)
    if(pobs <= 0.05)
    {
      b1 <- list()
      contab <- list()
      my.matrix <- matrix()
      perms.pvalue <- numeric()
      
      if(n <= nrow(samplingdata))
      {
        for(i in 1:n.iter) #change to no. of interations required
        {
          Sys.sleep(0.25)
          #Remove non-data columns
          if(sampler == "Row") {
            samplingdata2 <- setNames(data.frame(t(samplingdata[,-1])), samplingdata[,1]) 
            samplingdata2 <- cbind.data.frame(rownames(samplingdata2), samplingdata2)
            #colnames(samplingdata2) <- samplingdata2[1,]
          } else if(sampler == "Column") {
            samplingdata2 <- samplingdata
          }
          
          samplingdata2<- samplingdata2[,c(-1, -2)]
          # estimate no. of col groups
          classfn <- as.character(unlist(samplingdata2[1,])) 
          classfn.count <- table(classfn)
          classfn.names <- names(table(classfn))
          
          colnames(samplingdata2) <- ifelse(grepl(classfn.names[1], samplingdata2[1,]) , paste(colnames(samplingdata2), ".G1", sep=""), paste(colnames(samplingdata2), ".G2", sep=""))
          
          #Remove non-data rows
          samplingdata3<- samplingdata2[-1,]
          
          # randomly select user selected no. of rows from a dataset with replacement
          bootdata <- samplingdata3[sample(1:nrow(samplingdata3), n, replace=TRUE),]
          
          # Check for missing data values
          # complete.cases(bootdata2[[1]])
          bootdata2 <- bootdata[complete.cases(bootdata),]
          bootdata3 <- data.matrix(bootdata2)
          
          if(distmethod == "pearson correlation") {
          bootdata3 <- bootdata3[apply(bootdata3, 1, var) != 0, ]
          }
          
          # Call HEATPLOT2 EQUIVALENT HEATMAP2 CLUSTERING function
          if(norm == "Z-Score") {
            mat <- zClust(bootdata3, scale, zlim)
          } else if(norm == "Modified Z-Score") {
            mat <- modzClust(bootdata3, scale, zlim )
          } else if( norm == "none") {
            mat <- as.matrix(bootdata3)
          }
          
          
          
          if(distmethod == "pearson correlation") {
            # hm <- heatmap.2(as.matrix(mat[[1]]), Rowv=T, Colv=T, scale="none", hclust=function(x) hclust(x,method=clustmethod), distfun=function(x) as.dist((1-cor(t(x)))), margin = c(7,9), density.info=c("none"),trace=c("none"))
            # hc <- as.dendrogram(hm$colDendrogram)
            hc <- hclust(as.dist(1-cor(mat[[1]])))
            
          } else {
            # hm <- heatmap.2(as.matrix(mat[[1]]), scale="none", Rowv=T, Colv=T, hclust=function(c) {hclust(c,method=clustmethod)}, distfun=function(c) {dist(c,method=distmethod)}, margin = c(7,9), density.info=c("none"),trace=c("none"))
            #  hc <- as.dendrogram(hm$colDendrogram)
            hc <- hclust(dist(t(mat$data), method= distmethod), method=clustmethod)
            
          }
          
          ##    Hierarchical clustering of samples
          #hc <- hclust(hm$colDendrogram)
          
          #hclust(dist(t(z$data),method= distmethod), method = clustmethod) # will need to change depending on the clustring used
          # plot(hc)
          
          ### cut the tree into input number of clusters
          memb <- cutree(hc, k)
          # memb <- cutree(as.hclust(hm$colDendrogram), k)[as.hclust(hm$colDendrogram)$order]
          memb <- cbind(Row.Names = rownames(memb), memb, row.names = NULL)
          memb <- cbind.data.frame(rownames(memb),1:nrow(memb), memb) 
          memb
          
          contab[[i]] <- contingencyfun(memb)
          
          if (is.function(updateProgress)) {
            text <- paste0("Iteration #: ", i, " of ", n.iter)
            updateProgress(detail = text) }
          print(i)
        }
        
        # summarize each iteration row wise in a data frame
        my.matrix <- data.frame(matrix(unlist(contab), nrow=n.iter, byrow=T)) #no. of iterations
        
        perms.pvalue <- apply(my.matrix, 1, function(x) pvaluefunc(x))
        # my.matrix$Test <- length(which(my.matrix$p.value >= pobs))
        
        No <- ifelse(perms.pvalue > pobs, 1, 0)
        
        #calculate the pvalue to test significance of the CpG sites compared to  random sets of the same no. in separating T from N
        pvalue <- 1- (length(which(No==1)+1)/(n.iter+1))
        
        #hist(perms.pvalue, main = "Histogram")
        #abline(v= pobs)
        
        return(list(p.value= pvalue, p.obs = pobs, perms=  perms.pvalue) )
        
      }
      else
        stop("Sample size selected is larger than the dataset itself! Please select smaller sample size")
    }
    
    else
      stop("p-value from the obs data is not statistically significant, thus testing for random sites to assess their significance does not make sense")
  } 
  
}
