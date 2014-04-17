

cluster.rdist <- function(vec1, vec2){ 
    return(sqrt(sum((vec1-vec2)^2)))
}

# Construct a distance matrix. Since distance is symmetric, only the
# Upper-Right is calculated. Values on the diagonal will be 0.
cluster.dmatrix <- function(data){   
    nrecords <-  nrow(data)  
    results <- matrix( nrow=nrecords, ncol=nrecords) 
    for(i in 1:nrecords){ 
        vec1 <- data[i,]
        for(j in i:nrecords){
            vec2 <- data[(j),]
            results[i,j] <- cluster.rdist(vec1, vec2 )
        }  
    }  
    return(results)
}

# Takes a distance matrix and returns in sorted order the distance of the
# specified observation to all others in sorted order. The returned value
# is an (n-1) by 2 matrix where ''n' is the number of observations in
# the dmatrix. The two columns contain in sorted order the pointers and 
# the distances respectivsly. 
# 
# The dmatirx provided as input should specify the Upper-Right values.
# Lower-left values are ignored as symmetry  is assumed. Values on the
# diagonal should be 0.
cluster.neighbors <- function(dmatrix, i){
    # grab both col and row, including the 'NA' values in lower left
    foo <- c(dmatrix[,i],dmatrix[i,])
    # drop the NA and the diagonal
    distances <- subset(foo, foo > 0.0)
    # create the parallel 'pointer' vector. Special treatment required
    # if 'i' is the 1st or last observation.
    if(i == 1){
        pointers <- c( (i+1): (length(distances)+1) )
    }else if (i == ncol(dmatrix)){
        pointers <- c(1:(i-1) )
    }else{
        pointers <- c(1:(i-1),(i+1):(length(distances)+1))
    }	
    # now sort both by the distances (lowest to highest)
    unsorted <- cbind(pointers,distances)
    sorted <- unsorted[order(unsorted[,2],decreasing=FALSE),]
    return (sorted)
}

#' Takes a distance matrix and returns in sorted order the distance to
#' all other data points. The returned values are therefore two matrix:
#' one with the distance values and one with the point (i.e., observation) 
#'number.
#'
#' The dmatrix provided as input should specify the Upper-Right values.
#' Lower-left values are ignored as symmetry  is assumed. Values on the
#' diagonal should be 0.
cluster.nmatrix <- function(dmatrix){
    cnt <- nrow(dmatrix)-1 
    nmatrixP <- matrix(nrow=cnt,ncol=cnt+1) # pointers
    nmatrixD <- matrix(nrow=cnt ,ncol=cnt+1) # distancees
    for(i in 1: nrow(dmatrix)){
        foo <- cluster.neighbors(dmatrix,i)  
        nmatrixP[,i]<-foo[,1]
        nmatrixD[,i]<-foo[,2] 
    }
    return(list(dist=nmatrixD, ptr=nmatrixP))
}

bmrk.elbowPt <- function( dataSeq, min=FALSE, trim=0.01){
    margin <- as.integer(trim*nrow(dataSeq))    
    x1 <- dataSeq[(1+margin),1]
    y1 <- dataSeq[(1+margin),2]
    x2 <- dataSeq[(nrow(dataSeq)-margin),1]
    y2 <- dataSeq[(nrow(dataSeq)-margin),2]
    
    a <- (y2-y1)/(x2-x1)
    b <- -1
    c <- y1 - (x1 * a)
    denom <-   (a*a)+(b*b) 
    denomSqRt <- sqrt(denom)
    dist <-  rep(-1, nrow(dataSeq))
    for(i in (1+margin): (nrow(dataSeq)-margin)){ 
        dist[i]<- abs((a*dataSeq[i,1]) +(b*dataSeq[i,2]) + c) / denomSqRt
    }
    if(min){
        ptIndex <- which.min(dist[ (1+margin): (nrow(dataSeq)-margin)])
        distance <- dist[ptIndex]
        
    }else{
        ptIndex <- which.max(dist)
        distance <- dist[ptIndex]
    }
    # now calculate point of intersect along the baseline
    x0 <- dataSeq[ptIndex,1]
    y0 <- dataSeq[ptIndex,2]
    intersectPt.X <- ( (b*((b*x0)-(a*y0)))-(a*c))/denom
    intersectPt.Y <- ( (a*((a*y0)-(b*x0) ))-(b*c))/denom
    intersectPt <- list(X=intersectPt.X, Y=intersectPt.Y)
    results <- list(dist=distance, intersect =intersectPt, index=ptIndex, base.dist=dist[ (1+margin): (nrow(dataSeq)-margin)])
    return(results)
}

bmrk.epsAnalysis <- function(dataSet, trim=0.0, kVec = c(5,10,20,30,40)){
    system.time(distanceMatrix <-cluster.dmatrix(dataSet  ))   
    sorted <- cluster.nmatrix(distanceMatrix)
    rd2 <- sorted$dist 
    xMatrix <- matrix(nrow=nrow(rd2)+1, ncol=length(kVec))
    for (i in 1:length(kVec)){
        k <- kVec[i]
        rowK <- rd2[k,] 
        x <- sort(rowK)
        xMatrix[,i]<-x
    }
    # For each k-curve, add 'knee of the curve' analysis using 'farthest point 
    # from line' approach...
    eVec <- vector(mode='list',length=length(kVec))
    for(i in 1:length(kVec)){
        kCurve <- xMatrix[,i]
        dataSeq <- cbind(seq(0, length(kCurve)-1), kCurve)
        eRes <- bmrk.elbowPt( dataSeq, trim=trim) 
        eVec[[i]] <- eRes
    }
    return(list(dmatrix=distanceMatrix, xmatrix=xMatrix, evector=eVec, kvalues=kVec) )
}



bmrk.cluster<- function(dataSet , epsRes ){ 
    xMatrix <- epsRes$xmatrix
    eVector <- epsRes$evector
    kvalues <- epsRes$kvalues
    analysis <- matrix(ncol=5, nrow=0)
    for(i in 1:length(eVector)){  
        # need the Y value
        testRes <- eVector[[i]]
        dPtIdx <- testRes$index
        kCurve <- xMatrix[,i]
        eps <- kCurve[dPtIdx ]
        k <- kvalues[i]
        ds <- dbscan(dataSet , eps, MinPts=k)  
        d2base <- testRes$base.dist
        # normalize data to range [0,1]
        dMax <- max(d2base)
        # since its distance metric it is a given that dMin = 0
        d2base <- d2base * (1.0 / dMax)
        analysis <- rbind( analysis, c(k, eps, testRes$dist, mean(d2base), sd(d2base)))
    }
    colnames(analysis) <- c("k", "eps", "max(d)", "mean(d)", "sDev(d)")
    return(analysis)
}

# Starts out the same as clusterGen03 as far as the test data goes but then:
#   1) runs epsAnalysis ONLY on the entire data set
#   2) 
bmrk.dbscan <- function(size=1, seed=1234){
    # Uncomment next line to clear the entire workspace and start with
    # a 100% clean slate....
    # rm(list=ls()) 
    library('fpc')
    source('./TestDataGenerator.R', echo=FALSE)  
    
    noiseSize = 500*size
    clusterSize = 300*size
    set.seed(seed)
    # noise is same for all so do that upfront....
    noise <-genUniformXY(xCtr=5, yCtr=1, ptCnt=noiseSize, scale=12.0)
    colnames(noise) <- c("x", "y")
    
    # Generate 3 NON-INTERSECTING clusters that vary only in the Y coordinate
    #  band and the density. Number of points and the X coordinate band are 
    # identical.... 
    x <- c(0, 10)
    y1 <- c(-3, -3)
    sd1 <- 0.6    
    cluster1 <- genLineSeq(x, y1, sd1, iCnt=clusterSize)  
    # cluster 2 is higher density BUT same point count 
    y2 <- c(3, 3)
    sd2 <- 0.4    
    cluster2 <- genLineSeq(x, y2, sd2, iCnt=clusterSize)  
    # cluster 3 is highest density BUT same point count 
    y3 <- c(0, 0)
    sd3 <- 0.2    
    cluster3 <- genLineSeq(x, y3, sd3, iCnt=clusterSize)   
    # last is the whole cheese ball..
    dataSet3 <- rbind(cluster1,cluster2,cluster3,  noise) 
    
    # default analysis is to use  kVec = c(5,10,20,30,40)
    kValues <-  c(3, 5, 7, 20 ,25)
    
    epsRes3 <- bmrk.epsAnalysis(dataSet3, trim=0.1, kVec=kValues)  
    analysis <- bmrk.cluster(dataSet3, epsRes3  ) 
   
}

benchmark <- function(size=1, seed=1234){
    system.time(bmrk.dbscan(size, seed))
}