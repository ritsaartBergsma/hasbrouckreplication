library(xts)
library(lubridate)
library(timeDate)
library(tsDyn)
library(vars)
library(urca)
library(quantmod)
library(pracma)
library(systemfit)
#library(irtoys)
#set standard timezone
Sys.setenv(TZ='GMT')

cleanedDataPath <- "C:\\Users\\Ritsaart\\Documents\\R\\Hasbrouck SAS\\Cleaned Data\\"
outputPath <- "C:\\Users\\Ritsaart\\Documents\\R\\Hasbrouck SAS\\Output\\"

#This function creates a list of equations which form a VECM system. The z variables
#represent the cointegration relationship which are defined as:z_{x} = p1 - p_{x + 1} where is
#one of the other variables.The lag coefficients are represented as A_{i}_{j}_{k} where k is
#the lag number i and j represent the row and column coordinates in the coefficient matrix.

#input
#lagNumbers: column vector containing all the numbers of the lags to include
#so for examples 1:3 will result in a model with 3 lags.
#pricesnames: names of the variables which need to be estimated
vecmSystem <- function(pricenames,lagNumbers) {
  systemFormulas <- vector("list", length(pricenames))
  for (i in 1:length(pricenames)){
    a <- strcat(c(pricenames[i],"d", " ~ ",0), collapse = "");
    for (p in  1:(length(pricenames) -1)) {
      a <- strcat(c(a," + ","z",p),collapse = "");
    }
    for (k in lagNumbers){
      for (j in 1:length(pricenames)) {
        a <- strcat(c(a," + ",pricenames[j],"dL",k),collapse = "");
      }
    }
    systemFormulas[[i]] <- as.formula(a);
  }
  return(systemFormulas)
}

# Transforms the data which needed for the estimation of the VECM model. Input is
# the data and the lags which to include.
systemData <- function(pVector,lagNumbers) {
  pricenames <- colnames(pVector);
  
  pricenamesL <- paste(pricenames,"L",sep = "",collapse = NULL);
  
  dP <- diff(pVector,1);
  colnames(dP) <- paste(pricenames,"d",sep = "",collapse = NULL);
  dPL <- timeVector;
  for (i in lagNumbers){
    dPL1 <- lag(dP, i);
    colnames(dPL1) <- paste(pricenames,strcat(c("dL",i),collapse = ""),sep = "", collapse = NULL);
    dPL <- merge(dPL1,dPL);
  }
  
  
  pL = lag(pVector,1);
  colnames(pL) <- paste(pricenames,"L",sep = "",collapse = NULL);
  z <- timeVector;
  for (j in 2:length(pricenamesL)){
    zj <- pL[,pricenamesL[1]] - pL[,pricenamesL[j]] - (mean(pL[,pricenamesL[1]],na.rm = TRUE) - mean(pL[,pricenamesL[j]],na.rm = TRUE));
    z <- merge(zj,z);
  }
  colnames(z) <- paste("z",seq((length(pricenamesL)-1),1),sep = "" , collapse = NULL);
  
  modelData <- merge.xts(dP,merge.xts(dPL,merge.xts(pL,z)));
  modelData <- modelData[rowSums(is.na(modelData)) < 1,];
  
  return(modelData)
}

informationShare <- function(coefficientMatrix,lagNumbers,iterations,covMatrix) {
  numDims <- nrow(coefficientMatrix);
  alpha <- coefficientMatrix[,1:numDims-1];
  beta <- cbind(rep(1,numDims-1),-1 * diag(numDims-1));
  coefficient <- vector("list", length(lagNumbers))
  
  for (i in 1:length(lagNumbers)) {
    coefficient[[i]] <- coefficientMatrix[,(numDims + (i-1) * numDims) : ((numDims - 1) + i*numDims) ];
  }
  
  C <- matrix(0,numDims,numDims);
  for (shockToVariable in 1:numDims){
    irf <- matrix(0,numDims,max(lagNumbers) + 1);
    irf[shockToVariable,ncol(irf)] <- 1;
    p <- irf[,ncol(irf)];
    for (j in 1:iterations) {
      p <- p + alpha %*% (beta %*% p);
      
      a <- matrix(0,numDims,1);
      for (k in 1:length(lagNumbers)) {
        a <- a +  coefficient[[k]] %*% (irf[,ncol(irf) - (lagNumbers[k] -1)] - irf[,ncol(irf) - lagNumbers[k]]); 
      }
      p <- p + a;
      irf <- cbind(irf,p);
    }
    C[,shockToVariable] <- p;
  }
  
  c <- C[1,];
  sigmaW <- c %*% covMatrix %*% c;
  infoShare <- matrix(0,1,numDims);
  for (p in 1:numDims){
    infoShare[p] <- ((c[p] ^ 2) * covMatrix[p,p]) / sigmaW;
  }
  
  return(infoShare)
}

MVArepresenation <- function(coefficientMatrix,lagNumbers,shockToVariable){
  numDims <- nrow(coefficientMatrix);
  alpha <- coefficientMatrix[,1:numDims-1];
  beta <- cbind(rep(1,numDims-1),-1 * diag(numDims-1));
  coefficient <- vector("list", length(lagNumbers))
  
  for (i in 1:length(lagNumbers)) {
    coefficient[[i]] <- coefficientMatrix[,(numDims + (i-1) * numDims) : ((numDims - 1) + i*numDims) ];
  }
  
  irf <- matrix(0,numDims,max(lagNumbers) + 1);
  irf[shockToVariable,ncol(irf)] <- 1;
  p <- irf[,ncol(irf)];
  for (j in 1:600) {
    
    p <- p + alpha %*% (beta %*% p);
    
    a <- matrix(0,numDims,1);
    for (k in 1:length(lagNumbers)) {
      a <- a +  coefficient[[k]] %*% (irf[,ncol(irf) - (lagNumbers[k] -1)] - irf[,ncol(irf) - lagNumbers[k]]); 
    }
    
    p <- p + a;
    irf <- cbind(irf,p);
    
  }
  irf <- irf[,(max(lagNumbers)+1):ncol(irf)];
  return(irf)
}

plotIRF <- function(matrix,titleString) {
  plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen",rgb(r=0.5,g=0.5,b=0.0))
  
  plot(matrix[1,], type="l", col= plot_colors[1], xlab="seconds",
       ylab="Price impact", cex.lab=0.8, lwd=2, ylim = range(c(0,1)))
  for (i in 1:nrow(matrix)){ 
    lines(matrix[i,],type = "l",col = plot_colors[i])
  }
  
  #legend("topleft", colnames(matrix), cex=0.8, col=plot_colors, 
  #       lty=1:3, lwd=2, bty="n" );
  title(main=titleString)
}

Data <- read.csv("sasPriceData.csv", header=TRUE, stringsAsFactors=FALSE)
timeVector <- as.POSIXct(Data[,1],format = "%H:%M:%S")
pVector <- as.xts(Data[,2:4], order.by = timeVector);
pVector <- pVector[,c("p1","p2","p3")];
pricenames <- colnames(pVector);


nMA <- 2;
lagNumbers <- c(1:nMA);

systemFormulas <- vecmSystem(pricenames,lagNumbers)
modelData <- systemData(pVector,lagNumbers)
fitsur <- systemfit(systemFormulas,data = modelData,method = "SUR");

coefficientMatrix <- matrix(coef(fitsur), nrow = length(pricenames), byrow = TRUE);

irf <- MVArepresenation(coefficientMatrix,lagNumbers,3);
plotIRF(irf,"Shock to p1")
informationShare(coefficientMatrix,lagNumbers,2000,fitsur$residCovEst)


