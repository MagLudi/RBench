# Functions used to support testing, demonstration, and validation of other
# components and libraries. Note that these functions ARE NOT test functions but
# utilities that support the creation tests, primarily by generating test data.
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generates a sequence of X-Y data points in which the X values are specified as
# an input and the Y values are randomly generated.
# 
# Algorithm: 
#
# Step 1) The X values are input as a list of 'inflection points'.
# Intermediate X coordinates are then generated via linear interpolation. [NOTE:
# currently the nummber of intermediate points (iCnt) is defaults to 5].
#
# Step 2) The Y values are randomly generated for the inflextion points. 
# Intermediate values are then generated via linear interpolation, after which
# noise is added. The noise has a random normal deviation with mean 0 
# and a standard deviation specified as an input parameter.
#............................
genXYSeq <-function(x, noise, iCnt=5){ 
    ptCnt <- 1+(length(x)-1)*iCnt
    y <- rnorm(length(x))
    i1 <- approx(x,y, n=ptCnt)
    yNoisy <- i1$y + rnorm(length(i1$y), sd=noise)  
    i1$y <- yNoisy
    return (i1)
}

# Generates a sequence of X-Y data points in which a sequence of X-Y coordinates are 
# specify a sequence of connected line segements and then data points are randomly 
# generated along the line segements. The numer of points per segment is specified by 
# iCnt. The randomness is specified by the 'jitter' parameter which is used as the
# standard deviation in a rnorm distribution.
genLineSeq <-function(x, y, jitter, iCnt=5){ 
    ptCnt <- 1+(length(x)-1)*iCnt 
    i1 <- approx(x,y, n=ptCnt)
    yNoisy <- i1$y + rnorm(length(i1$y), sd=jitter)  
    i1$y <- yNoisy
    xyPts <- cbind(i1$x,i1$y)
    return (xyPts)
}

genUniformXY <-function(xCtr=0, yCtr=0, ptCnt=10, scale=1.0){
    xpts <-((runif(ptCnt) - 0.5)*scale) + xCtr
    ypts <-((runif(ptCnt) - 0.5)*scale) + yCtr
    xyPts <- cbind(xpts,ypts)
    return(xyPts)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generates a sequence of X data points in which the X values are specified as
# an input and multiple Y values are randomly generated using noise values 
# specified as an input parameter.
#
# Example:
#   x <- c(0:10)
#   sd <- c(0.0, 0.2, 0.5, 0.7)
#   data <- genMVSeq(x, sd, iCnt=5)
#   write.csv(data, file="./PseudoData/TestSeq01.csv")
#
#   will generate 55 'X' points, each with 4 Y values with increasing levels
#   of random noise, and then save it in CSV format.
#............................
genMVSeq <-function(x, noise, iCnt=5){ 
  base <- genXYSeq(x,noise[1],iCnt)
  message(c("Base length=", length(base$x) )) 
  vNames <- c("x",  paste("y{" ,noise[1], "}", sep=""))
  for(i in 2:length(noise)){
    nextBase <-  genXYSeq(x,noise[i],iCnt)
    base <- merge(base, nextBase, by="x")
    nextCol <- paste("y{",noise[i],"}",sep="")
    vNames <- c(vNames, nextCol)
  }
  colnames(base)<- vNames
  return(base)
}