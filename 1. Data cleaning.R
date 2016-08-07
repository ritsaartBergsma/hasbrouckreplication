library(xts)
library(lubridate)
library(timeDate)
#set standard timezone
Sys.setenv(TZ='GMT')

#Use double backward slashes when defining the paths
inputPath <- "C:\\Users\\Ritsaart\\Documents\\R\\Hasbrouck SAS\\Input\\";
cleanedDataPath <- "C:\\Users\\Ritsaart\\Documents\\R\\Hasbrouck SAS\\Cleaned Data\\"

# input data
ETF.SP500.March.Quote <- read.csv(paste(inputPath,"2000_quotes-2000-03-SPY.A.csv",sep = ""), header=TRUE, stringsAsFactors=FALSE)
ETF.SP500.April.Quote <- read.csv(paste(inputPath,"2000_quotes-2000-04-SPY.A.csv",sep = ""), header=TRUE, stringsAsFactors=FALSE)
ETF.SP500.May.Quote <- read.csv(paste(inputPath,"2000_quotes-2000-05-SPY.A.csv",sep= ""), header=TRUE, stringsAsFactors=FALSE)

ETF.SP500.March.trade <- read.csv(paste(inputPath,"2000_trades-2000-03-SPY.csv",sep= ""),header=TRUE, stringsAsFactors=FALSE)
ETF.SP500.April.trade <- read.csv(paste(inputPath,"2000_trades-2000-04-SPY.csv",sep = ""),header=TRUE, stringsAsFactors=FALSE)
ETF.SP500.May.trade <- read.csv(paste(inputPath,"2000_trades-2000-05-SPY.csv", sep = ""),header=TRUE, stringsAsFactors=FALSE)

pit_price.SP500.March.set1 <- read.csv(paste(inputPath,"2000_trades-2000-03-2SPH0.csv",sep = ""),header=TRUE, stringsAsFactors=FALSE)
pit_price.SP500.March.set2 <- read.csv(paste(inputPath,"2000_trades-2000-03-2SPM0.csv",sep = ""),header=TRUE, stringsAsFactors=FALSE)
pit_price.SP500.April <- read.csv(paste(inputPath,"2000_trades-2000-04-2SPM0.csv",sep = ""),header=TRUE, stringsAsFactors=FALSE)
pit_price.SP500.May <- read.csv(paste(inputPath,"2000_trades-2000-05-2SPM0.csv", sep = ""),header=TRUE, stringsAsFactors=FALSE)

E_mini.SP500.March.set1 <- read.csv(paste(inputPath,"2000_trades-2000-03-ESH0.csv", sep = ""),header=TRUE, stringsAsFactors=FALSE)
E_mini.SP500.March.set2 <- read.csv(paste(inputPath,"2000_trades-2000-03-ESM0.csv", sep = ""),header=TRUE, stringsAsFactors=FALSE)
E_mini.SP500.April <- read.csv(paste(inputPath,"2000_trades-2000-04-ESM0.csv", sep = ""),header=TRUE, stringsAsFactors=FALSE)
E_mini.SP500.May <- read.csv(paste(inputPath,"2000_trades-2000-05-ESM0.csv", sep = ""),header=TRUE, stringsAsFactors=FALSE)

#Simply moving average filter
#input:
#data: data frame for filtering
#filterVar: variable which needs to be filtered
#thresholdPercentages: percentage of the mean of the filterVar at which the
#filterVar can deviate from the moving average
#n: number of periods for the moving average
movingAverageFilterCenteredOud <- function(data,filterVar,thresholdPercentages,n){
  dates <- as.POSIXct(as.character(data[,"Date.L."]),format = "%Y%m%d",tz = "GMT");
  threshold <- mean(data[,filterVar],na.rm = TRUE) *  (thresholdPercentages/100);
  print(threshold);
  movingAverage <- as.numeric(filter(data[,filterVar],rep(1/n,n), sides=2))
  movingAverage[1:((n - 1) / 2 )] <- data[1:((n - 1) / 2 ),filterVar];
  data2 <- data[(data[,filterVar] - movingAverage < threshold) & (movingAverage - data[,filterVar] < threshold),];
  return(data2)
}

#moving average filter which filters per day (and thus avoids deleting the first few entries of each day).
#input: data.frame with price data, variable which needs to be filtered on, threshold above which values
# are deleted and the number of periods used for the moving average.
movingAverageFilterCentered <- function(data,filterVar,thresholdPercentages,n){
  dates <- as.POSIXct(as.character(data[,"Date.L."]),format = "%Y%m%d",tz = "GMT");
  threshold <- mean(data[,filterVar],na.rm = TRUE) *  (thresholdPercentages/100);
  a <- NULL;
  for (i in unique(dates)) {
    data2 <- data[dates == i,];
    movingAverage <- as.numeric(filter(data2[,filterVar],rep(1/n,n), sides=2))
    movingAverage[1:((n - 1) / 2 )] <- data2[1:((n - 1) / 2 ),filterVar];
    data3 <- data2[(data2[,filterVar] - movingAverage < threshold) & (movingAverage - data2[,filterVar] < threshold),];
    a <- rbind(a,data3);
  }
  return(a)
}

filterQuote <- function(quoteData,threshold,nPeriods) {
  #remove negative spreads and spreads greater than 1
  quoteData1 <- subset(quoteData, (Ask.Price - Bid.Price) > 0 & (Ask.Price - Bid.Price) <1);
  #remove prices 50 cents greater than moving avarage
  #
  quoteData2 <- movingAverageFilterCentered(quoteData1,"Ask.Price",threshold,nPeriods)
  quoteData3 <- movingAverageFilterCentered(quoteData2,"Bid.Price",threshold,nPeriods)
  #remove empty values
  quoteData4 <- quoteData3[quoteData3["Quote.Time"] != "" & !is.na(quoteData3[,"Ask.Price"]) & !is.na(quoteData3[,"Bid.Price"]) & !is.na(quoteData3[,"Date.L."]),];
  return(quoteData4)
}

filterTradeData <- function(tradeData,threshold,nPeriods) {
  #filter out values 50 cent larger than moving avarage of the past 10 days
  step1 <- movingAverageFilterCentered(tradeData,"Price",threshold,nPeriods)
  #filter out unusual trades
  step2 <- step1[step1["Type"] == "Trade",];
  step3 <- step2[step2["Exch.Time"] != "" & !is.na(step2[,"Price"]) & !is.na(step2[,"Date.L."]),];
  return(step3)
}

# function takes to data frames does a roll over on the rollOverDate specified as "%Y%m%d".
rollOver <- function(dataSet1,dataSet2,rollOverDate){
  data1 <- dataSet1[dataSet1[,"Date.L."] <= rollOverDate,]
  data2 <- dataSet2[dataSet2[,"Date.L."] > rollOverDate,]
  data3 <- rbind(data1,data2)
  return(data3)
}

# This function creates a POSIXct object which contains all time on which
# the market was pen between the startDate and endate
createSequence <- function(startDate,endDate,marketOpen,marketClose){
  start <- as.POSIXct(paste(as.character(startDate),marketOpen, sep = " ", collapse = NULL), tz = "GMT", format = "%Y%m%d %H:%M:%S")
  end <- as.POSIXct(paste(as.character(endDate),marketClose, sep = " ", collapse = NULL), tz = "GMT", format = "%Y%m%d %H:%M:%S")
  timeSequence <- seq(from=start,to=end,by = "sec")
  
  #remove weekends
  timeSequence <- timeSequence[wday(timeSequence) != 1 & wday(timeSequence) != 7]
  
  
  #remove holidays
  holidayDates <- as.timeDate(strftime(holidayNYSE(unique(year(timeSequence))),format = "%Y-%m-%d"))
  dateVector <- as.timeDate(strftime(timeSequence,format = "%Y-%m-%d"))
  timeSequence <- timeSequence[!is.element(dateVector@Data,holidayDates@Data)]
  
  #remove times before and after market close
  dateVector <- as.character(date(timeSequence))
  marketOpenVector <- as.POSIXct(paste(dateVector,marketOpen, sep = " ", collapse = NULL), tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
  marketCloseVector <-as.POSIXct(paste(dateVector,marketClose, sep = " ", collapse = NULL), tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
  timeSequence <- timeSequence[timeSequence >= marketOpenVector & timeSequence <= marketCloseVector]
  return(timeSequence)
}

#This function creates an xts object from the dataframe objects. It also filters out all 
#price data before the market open and after the market close.
createXTS <- function(data, dateVar,timeVar, priceVar,marketOpen,marketClose){

  timeVector <- paste(as.character(data[,dateVar]),data[,timeVar], sep = " ", collapse = NULL);
  options(digits.secs = 3);
  data$timeDate <- as.POSIXct(timeVector,tz = "GMT", format = "%Y%m%d %H:%M:%S");
  
  #filter out trades before opening and after close
  
  marketOpenVector <- as.POSIXct(paste(as.character(data[,dateVar]),marketOpen, sep = " ", collapse = NULL), tz = "GMT", format = "%Y%m%d %H:%M:%S");
  marketCloseVector <-as.POSIXct(paste(as.character(data[,dateVar]),marketClose, sep = " ", collapse = NULL), tz = "GMT", format = "%Y%m%d %H:%M:%S");
  
  dataFiltered <- data[data[,"timeDate"] >= marketOpenVector & data[,"timeDate"] <= marketCloseVector,];
  time_series <- xts(dataFiltered[,priceVar], order.by = dataFiltered$timeDate);
  
  #add market openings
  #time_series <- addMarketOpenings(time_series,marketOpen)
  return(time_series)
}

# Do the roll over
E_mini.SP500.March <- rollOver(E_mini.SP500.March.set1,E_mini.SP500.March.set2,20000308);
pit_price.SP500.March <- rollOver(pit_price.SP500.March.set1,pit_price.SP500.March.set2,20000308);

# Aggregate the monthly data.
ETF.SP500.Quote <- rbind(ETF.SP500.March.Quote, ETF.SP500.April.Quote, ETF.SP500.May.Quote);
ETF.SP500.trade <- rbind(ETF.SP500.March.trade,ETF.SP500.April.trade,ETF.SP500.May.trade);
pit_price.SP500 <- rbind(pit_price.SP500.March,pit_price.SP500.April,pit_price.SP500.May);
E_mini.SP500 <- rbind(E_mini.SP500.March,E_mini.SP500.April,E_mini.SP500.May);

# Filter the data
ETF.SP500.Quote <- filterQuote(ETF.SP500.Quote,0.2,11);
ETF.SP500.trade <- filterTradeData(ETF.SP500.trade,0.2,11);
pit_price.SP500 <- filterTradeData(pit_price.SP500,0.2,11);
E_mini.SP500.trade <- filterTradeData(E_mini.SP500,0.2,11);

# create a midpoint variable
ETF.SP500.Quote$midpoint <- (ETF.SP500.Quote$Ask.Price + ETF.SP500.Quote$Bid.Price) / 2;

# Create XTS objects
ETF.SP500.Quote.xts <- createXTS(ETF.SP500.Quote,"Date.L.","Quote.Time","midpoint", "14:30:00", "21:00:00");
ETF.SP500.trade.xts <- createXTS(ETF.SP500.trade,"Date.L.","Exch.Time","Price", "14:30:00", "21:00:00");
pit_price.SP500.xts <- createXTS(pit_price.SP500,"Date.L.","Exch.Time","Price", "14:30:00", "21:00:00");
E_mini.SP500.trade.xts <- createXTS(E_mini.SP500.trade,"Date.L.","Exch.Time","Price", "14:30:00", "21:00:00");

# Aggregate the xts objects on seconds
ETF.SP500.Quote.xts <- ETF.SP500.Quote.xts[endpoints(ETF.SP500.Quote.xts, on = "seconds", k = 1)];
ETF.SP500.trade.xts <- ETF.SP500.trade.xts[endpoints(ETF.SP500.trade.xts, on = "seconds", k = 1)];
pit_price.SP500.xts <- pit_price.SP500.xts[endpoints(pit_price.SP500.xts, on = "seconds", k = 1)];
E_mini.SP500.trade.xts <- E_mini.SP500.trade.xts[endpoints(E_mini.SP500.trade.xts, on = "seconds", k = 1)];

# Create a sequence of all time at which the market is open between the first and the last day
# of the sample data.
marketTimes <- createSequence(20000301,20000531,"14:30:00","21:00:00");

#Merge all the data with this timesequence.
pVector <- merge.xts(ETF.SP500.Quote.xts,merge.xts(ETF.SP500.trade.xts,
                                                   merge.xts(pit_price.SP500.xts,merge.xts(E_mini.SP500.trade.xts,marketTimes,join = "right"),
                                                             join = "right"),join = "right"),join = "right")

colnames(pVector) <- c("ETF Quote Midpoint","ETF Trade Price","Pit Contract Price","E-mini Contract Price")
write.zoo(pVector,file = paste(cleanedDataPath,"pVector.csv",sep = ""),sep = ",");


