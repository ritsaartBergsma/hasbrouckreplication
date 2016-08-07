library(xts)
library(lubridate)
library(timeDate)
library(systemfit)
library(xtable)
library(pracma)
library(foreach)
library(doParallel)
library(iterators)
#set standard timezone
Sys.setenv(TZ='GMT')

#Paths
cleanedDataPath <- "C:\\Users\\Ritsaart\\Documents\\R\\Hasbrouck SAS\\Cleaned Data\\"
outputPath <- "C:\\Users\\Ritsaart\\Documents\\R\\Hasbrouck SAS\\Output\\"

excecuteParalell = TRUE;

if(excecuteParalell) {
  #set number of cores to use (default is total number of cores - 1)
  #is always at least 1.
  no_cores <- max(detectCores() - 1,1);
}else{
  no_cores <- 1;
}
  #initiate cluster
  cl <- makeCluster(no_cores);
  registerDoParallel(cl,no_cores);

#User input
#lagnumbers (either fill in nMA or choose a custom vector of which lags to include)
nMA <- 10;
lagNumbers <- c(1:nMA);

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
  timeVector <- as.POSIXct(pVector,tz = "GMT", format = "%Y-%m-%d %H:%M:%S");
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
    zj <- pL[,pricenamesL[1]] - pL[,pricenamesL[j]] - 
      (mean(pL[,pricenamesL[1]],na.rm = TRUE) - mean(pL[,pricenamesL[j]],na.rm = TRUE));
    z <- merge(zj,z);
  }
  colnames(z) <- paste("z",seq((length(pricenamesL)-1),1),sep = "" , collapse = NULL)
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
#date: If the time series does not converge, an error message is printed
#with this date.
impulseResponseFunction <- function(coefficientMatrix,lagNumbers,iterations,shockToVariable,date){
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
  
  for (j in 1:iterations) {
    p <- p + alpha %*% (beta %*% p);
    a <- matrix(0,numDims,1);
    
    for (k in 1:length(lagNumbers)) {
      a <- a +  coefficient[[k]] %*% (irf[,ncol(irf) - (lagNumbers[k] -1)] - irf[,ncol(irf) - lagNumbers[k]]); 
    }
    p <- p + a;
    irf <- cbind(irf,p);
  }
  failedDate <- NULL;
  irf <- irf[,(max(lagNumbers)+1):ncol(irf)];
  
  #Check for convergence
  if(abs(irf[1,(iterations + 1)]-mean(irf[2:numDims,(iterations + 1)]) > 0.01)){
    failedDate <- date;
  }
  return(list(irf,failedDate))
}

#Plots an impulse response function. The function should be supplied as a matrix.
#Also a title and legend contents are required.
plotIRF <- function(matrix,titleString,legendnames) {
  #If more than four time series need to be plottes, add extra colours
  plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen",rgb(r=0.5,g=0.5,b=0.0))
  
  plot(matrix[1,], type="l", col= plot_colors[1], xlab="seconds",
       ylab="Price Impact", cex.lab=0.8, lwd=2, ylim = range(c(0,1)))
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

#Read data
Data <- read.csv(paste(cleanedDataPath,"pVector.csv", sep=""), header=TRUE, stringsAsFactors=FALSE)
dateVector <- as.timeDate(strftime(Data[,"Index"],format = "%Y-%m-%d"));
timeVector <- as.POSIXct(Data[,"Index"],tz = "GMT", format = "%Y-%m-%d %H:%M:%S");

#names of the time series
pricenames <- colnames(Data[,2:ncol(Data)])
systemFormulas <- vecmSystem(pricenames,lagNumbers)

#tranform to log prices in basispoints
Data.xts <- log(as.xts(Data[,2:ncol(Data)], order.by = timeVector)) * 10000;

#dates in the sample to estimate the model on
#filter out dates with zero observations (could be adjusted to filter out days
#with fewer than a certain number of observations)
observationsPerDay <- apply.daily(as.xts((!is.na(Data.xts)),order.by = timeVector), FUN = sum);
dates <- as.timeDate(strftime(observationsPerDay[observationsPerDay != 0],format = "%Y-%m-%d"));

#Create a list of pVectors for each day
pVectorList <- lapply(dates,function (a) na.locf(Data.xts[dateVector == a]));

#uncomment for debugging (include only the first two pvectors)
#pVectorList <- lapply(1:2,function (d) pVectorList[[d]]);

#foreach is a parallel for loop which outputs a list of two 
estimates <- foreach(pVector = pVectorList, .packages = c("xts","lubridate","timeDate","systemfit","pracma")) %dopar% {
  modelData <- systemData(pVector,lagNumbers)
  fitsur <- systemfit(systemFormulas,data = modelData,method = "SUR");
  coefficientMatrix <- matrix(coef(fitsur), nrow = length(pricenames), byrow = TRUE);
  
  #Calculate the impulse response function for a shock in every price and also get the date
  #if the irf fails to converge
  irfPlusFailedToConvergeDate <- lapply(1:length(pricenames),function(j)
                                 impulseResponseFunction(coefficientMatrix,lagNumbers,600,j,
                                 as.timeDate(strftime(pVector,format = "%Y-%m-%d"))[1])) #parse date
  irf <- lapply(1:length(pricenames),function (d) irfPlusFailedToConvergeDate[[d]][[1]]);
  failedToConvergeDates <- lapply(1:length(pricenames),function (d)
                           irfPlusFailedToConvergeDate[[d]][[2]]);
  ifShare <- informationShare(coefficientMatrix,lagNumbers,600,systemFormulas,modelData)
  list(irf,ifShare,failedToConvergeDates)
}

irfList <- lapply(estimates, function(d) d[[1]]);

#This list contains the dates on which the impulse response function
#failes to converge
failedToConverge <- unlist(lapply(estimates, function(d) d[[3]]));
print("failed to converge on:");
failedToConverge
#Create an impulse response graph for every irf 
for(i in 1:length(pricenames)){
  #Take average of all irfs in list
  impulseResponse <- Reduce("+",lapply(1:length(pVectorList),function (d) irfList[[d]][[i]])) * 
   (1 / (length(pVectorList)));
  plotIRF(impulseResponse,paste("Shock to ",pricenames[i],sep = ""),pricenames);
  plotName <- paste(paste("Shock to ",pricenames[i],sep = ""),".png",sep = "");
  dev.copy(png,paste(outputPath,plotName,sep = ""))
  dev.off()
}

infShareMin <- informationShareSummary(lapply(1:length(pVectorList),function (d)
  estimates[[d]][[2]][[1]]),paste(pricenames, "Min", sep = " "));
infShareMax <- informationShareSummary(lapply(1:length(pVectorList),function (d)
  estimates[[d]][[2]][[2]]),paste(pricenames, "Max", sep = " "));

#The next section formates the table in the format Hasbrouck (2003) uses. However, the header are
#too large for an A4 format. Adjusting for multiple row headers, is very complicated and
#is something for later on.
summaryTable <- cbind(infShareMin,infShareMax)[,unlist(sapply(1:length(pricenames), 
                      function(i) c(i,i + length(pricenames)),simplify = FALSE))]

#Remove dots from column headers
colnames(summaryTable) <- gsub("\\.", " " , colnames(summaryTable))
xtable(summaryTable,caption = "Information Shares",digits = 3)

#end cluster, important because R might crash if not done
stopCluster(cl);
