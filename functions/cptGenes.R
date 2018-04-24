cptGenes <- function (x,cpt_data,cpt_method,cpt_max,n){
  ## Identify the changepoints for the data
  # x : data within which you wish to identify the changepoints
  # i : numeric column index of the data frame to sort it by
  # b : sorting order, ascending (FALSE) or descending (TRUE)
  perDiffcpt = function(x,cpt_data,cpt_method,cpt_max){
    max <- length(x)/2
    if (cpt_method == "AMOC") {
      
      if(cpt_data == "mean"){
        changepoints <- cpt.mean(x,method=cpt_method,Q=1)
      }else if(cpt_data == "var"){
        changepoints <- cpt.var(x,method=cpt_method,Q=1)
      }else{
        changepoints <- cpt.meanvar(x,method=cpt_method,Q=1)
      }
      
    }else if (cpt_method == "PELT" || cpt_method == "SegNeigh" || cpt_method == "BinSeg") {
      
      if(max > cpt_max){
        if(cpt_data == "mean"){
          changepoints <- cpt.mean(x,method=cpt_method,Q=cpt_max)
        } else if (cpt_data == "var") {
          changepoints <- cpt.var(x,method=cpt_method,Q=cpt_max)
        }else {
          changepoints <- cpt.meanvar(x,method=cpt_method,Q=cpt_max)
        }
      }else if(max <= cpt_max){
        if(cpt_data == "mean"){
          changepoints <- cpt.mean(x,method=cpt_method,Q=max)
        } else if (cpt_data == "var") {
          changepoints <- cpt.var(x,method=cpt_method,Q=max)
        }else{
          changepoints <- cpt.meanvar(x,method=cpt_method,Q=max)
        }
      }
      
    }
    return ( cpts(changepoints) )
  }
  
  ## Add the estimated changepoint locations to the data frame
  # x : A data frame containing the data used for change point estimation
  # y : vector containing the changepoint locations for the data supplied
  cptAdd = function(x,y){
    level=character()
    l = 1
    for(i in 1:length(y) ) {
      if(i==1){
        for(j in 1:y[i]){
          level[j]=l
        }
      }else{
        a <- (y[i-1])+1
        for(j in a:y[i]){
          level[j]=l
        }
      }
      if(i==length(y)){
        l=1000
        a <- y[i]+1
        for(j in a:length(x)){
          level[j]=1000
        }
      }
      
      l=l+1
    }
    return(as.numeric(level))
  }
  
  ##########Sort the data
  sortData = function(x,i,b) {
    if(missing(x)){
      stop("input data is missing!")
    }
    if(missing(b)){
      b=FALSE
    }
    if(missing(i)){
      i=1
    }
    return( x[ order(x[,i], decreasing=b), ] )
  }
  
  ######### END: sub-functions used by the cptSamples main function #########
  
  #identify the changepoints
  cpts = perDiffcpt(x,cpt_data,cpt_method,cpt_max) ## subject to change
  #check for the number of changepoints identified
  if(length(cpts)<1){
    warning("No changepoints identified in the data set!")
    cpt_out=NULL
  }else{
    #add estimated changepoint locations to the data
    changepoints <- cptAdd(x,cpts)
    #append changepoints to the original data frame
    cpt_out <- data.frame(x,changepoints)
   
    #generate the data frame with all changepoints
    #format the data file for reading usability
    #any changepoint other than 1 is designated as '0'
    #plot changepoints within the data identified
    colnames(cpt_out)[ncol(cpt_out)-1] <- "Succ.Diff"
    #cpt_out [, ncol(cpt_out) ][ cpt_out [ ,ncol(cpt_out) ] == 1000 ] <- "NA"
    rownames(cpt_out) <- NULL
    #cpt_prop <- cpt_out
    
    #cpt_results <- cpt_out
    return(cpt_out)
  }
   
}


  