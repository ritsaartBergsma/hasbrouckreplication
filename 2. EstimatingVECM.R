library(xts)
library(lubridate)
library(timeDate)
library(tsDyn)
library(quantmod)
library(systemfit)
#library(irtoys)
#set standard timezone
Sys.setenv(TZ='GMT')

#User input
#lagnumbers (either fill in nMA or choose a custom vector of which lags to include)
nMA <- 10;
lagNumbers <- c(1:nMA);

#Paths
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

#This function returns a list of 2 vectors containing the minimum and maximum informationShare
#respectively.
#input 
#CoefficientMatrix: a matrix of coefficients of a VECM system
#lagNumbers: vector of lags in the system
#iterations: number of periods to forecast for estimating the C matrix
#systemFormulas: list of equations which together form the VECM model
#modelData: data used to estimate the VECM model
informationShare <- function(coefficientMatrix,lagNumbers,iterations,systemFormulas,modelData) {
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
  # Information shares
  #To calculate the max and minimum information shares we have to re-order the 
  #covariance matrix through reordering the euqations which need to be estimated
  #to the first equation, the largest amount of variance will be contributed to the last equation 
  #in systemfit (apperently).And thus this will be used to compute the max information share for 
  #that particular price. the i variable represents the equation which will come first in the re-
  #ordering and the j variable the last. 
  max.infShare <- matrix(0,1,numDims)
  min.infShare <- matrix(0,1,numDims)
  for(i in 1:numDims) {
    j <- numDims - (i - 1);
    #If there is an odd number of equations, i and j can become equal resulting in an
    #invalid VECM system.
    if(i != j) {
      #First create an index for the reordering
      nums <- Filter(function(x) x != j & x != i, 1:length(pricenames))
      index <- c(i,nums,j)
      systemFormulasRearranged <- systemFormulas[index]
      fitsur <- systemfit(systemFormulasRearranged,data = modelData,method = "SUR");
      omega <- as.matrix(fitsur$residCovEst);
      sigmaW <- c[index] %*% omega %*% c[index];
      F <- chol(omega)
      min.infShare[i] <- ((c[index] %*% F) ^ 2)[1] / sigmaW;
      max.infShare[j] <- ((c[index] %*% F) ^ 2)[numDims] / sigmaW
      #An odd number of equations means that for the middle variable, both the
      #max and minimum information share need to be estimated differently
    } else {
      #Create two indices, one with the middle index number in the beginning for the maximum
      #information share and one with the middle index number last for the minimum information
      #share.
      nums <- Filter(function(x) x != ceiling(length(pricenames) / 2), 1:numDims)
      indexMax <- c(ceiling(numDims / 2),nums);
      indexMin <- c(nums,ceiling(numDims / 2));
      listIndex <- list(indexMax,indexMin)
      
      #Both the min and max informationshare has to be calculated.
      for (x in 1:2){
        index <- listIndex[[x]];
        systemFormulasRearranged <- systemFormulas[index]
        fitsur <- systemfit(systemFormulasRearranged,data = modelData,method = "SUR");
        omega <- as.matrix(fitsur$residCovEst);
        sigmaW <- c[index] %*% omega %*% c[index];
        F <- chol(omega)
        if(x == 1){
          min.infShare[ceiling(numDims / 2)] <- ((c[index] %*% F) ^ 2)[1] / sigmaW;
        } else {
          max.infShare[ceiling(numDims / 2)] <- ((c[index] %*% F) ^ 2)[numDims] / sigmaW
        }
      }
    }
  }
  return(list(min.infShare,max.infShare))
}

#This function returns the irf of a VECM system in matrix form
#input
#coefficientMatrix: matrix of coefficients of a VECM model
#lagNumbers: lags used in the model
#shockToVariable: number of the variable who receives the shock in the
#first period.
impulseResponseFunction <- function(coefficientMatrix,lagNumbers,shockToVariable){
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

#Plots an impulse response function. The function should be supplied as a matrix.
#Also a title and legend contents are required.
plotIRF <- function(matrix,titleString,legendnames) {
  plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen",rgb(r=0.5,g=0.5,b=0.0))
  
  plot(matrix[1,], type="l", col= plot_colors[1], xlab="seconds",
       ylab="Price impact", cex.lab=0.8, lwd=2, ylim = range(c(0,1)))
  for (i in 1:nrow(matrix)){ 
    lines(matrix[i,],type = "l",col = plot_colors[i])
  }
  
  legend("topleft", legendnames, cex=0.8, col=plot_colors, 
         lty=1:3, lwd=2, bty="n" );
  title(main=titleString)
}

#Function created a dataframe with descriptive statistics also reported in
#Hasbrouck (2003).
informationShareSummary <- function(infList,names) {
  infShare <- as.data.frame(matrix(unlist(infList),ncol = length(names), byrow = TRUE))
  Median <- sapply(infShare,median);
  Mean <- sapply(infShare,mean);
  St.Dev <- sapply(infShare,sd);
  SEM <- sapply(infShare, function(x) sd(x) / sqrt(length(x)));
  descriptiveStatistics <- rbind(Median,Mean,St.Dev,SEM);
  colnames(descriptiveStatistics) <- names;
  return(descriptiveStatistics)
}

Data <- read.csv(paste(cleanedDataPath,"pVector.csv", sep=""), header=TRUE, stringsAsFactors=FALSE)

dateVector <- as.timeDate(strftime(Data[,"Index"],format = "%Y-%m-%d"));
dates <- unique(dateVector);
timeVector <- as.POSIXct(Data[,"Index"],tz = "GMT", format = "%Y-%m-%d %H:%M:%S");


lagNumbers <- c(1:nMA);
pricenames <- colnames(Data[,2:ncol(Data)])
systemFormulas <- vecmSystem(pricenames,lagNumbers)

Data.xts <- log(as.xts(Data[,2:ncol(Data)], order.by = timeVector)) * 10000;


irfLists <- list();
for(i in 1:length(pricenames)){
  mylist <- list();
  irfLists <- append(irfLists,list(mylist));
}

informationShareListMax <-list();
informationShareListMin <- list();

counter <- 1;
for(i in 1:length(dates)) {
  a <- dates[i];
  pVector <- na.locf(Data.xts[dateVector == a]);
  #Check if pVector is not empty (even after filtering out holidays, there still seems to be
  #a day in the sample which does not contain any trading data)
  if(!(sum(is.na(pVector)) == nrow(pVector) * ncol(pVector))){
    modelData <- systemData(pVector,lagNumbers)
    fitsur <- systemfit(systemFormulas,data = modelData,method = "SUR");
    coefficientMatrix <- matrix(coef(fitsur), nrow = length(pricenames), byrow = TRUE);
    for(j in 1:length(pricenames)){
      irf <- impulseResponseFunction(coefficientMatrix,lagNumbers,j);
      irfLists[[j]][[counter]] <- irf;
    }
    ifShare <- informationShare(coefficientMatrix,lagNumbers,600,systemFormulas,modelData)
    informationShareListMin[[counter]] <- ifShare[[1]];
    informationShareListMax[[counter]] <- ifShare[[2]];
    counter <- counter + 1;
  }
}

#Taking averages of the impulse response functions
for(i in 1:length(pricenames)){
  #Take average of all irfs in list
  impulseResponse <- Reduce("+",irfLists[[i]]) * (1 / (counter -1))
  plotIRF(impulseResponse,paste("Shock to ",pricenames[i],sep = ""),pricenames);
  plotName <- paste(paste("Shock to ",pricenames[i],sep = ""),".png",sep = "");
  dev.copy(png,paste(outputPath,plotName,sep = ""))
  dev.off()
}

infShareMin <- informationShareSummary(informationShareListMin,paste(pricenames, "Min", sep = " "));
infShareMax <- informationShareSummary(informationShareListMax,paste(pricenames, "Max", sep = " "));
summaryTable <- NULL;
colnamesTable <- NULL;
for(i in 1:length(pricenames)){
  summaryTable <- cbind(summaryTable,infShareMin[,i],infShareMax[,i])
  colnamesTable <- c(colnamesTable,colnames(infShareMin)[i],colnames(infShareMax)[i])
}
#Remove dots from column headers
colnames(summaryTable) <- gsub("\\."," ",colnamesTable)
summaryTable
xtable(summaryTable,caption = "Information Shares",digits = 3)