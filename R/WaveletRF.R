

#'@title Wavelet Transform Using Maximal Overlap Discrete Wavelet Transform (MODWT) Algorithm
#' @description Transforms the time series data by using hybrid MODWT algorithm.
#' @param ts Univariate time series
#' @param Wvlevels The level of wavelet decomposition
#' @param WFilter  Wavelet filter use in the decomposition
#' @param bndry The boundary condition of wavelet decomposition:'periodic' or 'reflection'
#' @param FFlag The FastFlag condition of wavelet decomposition: True or False
#' @import fracdiff forecast stats wavelets
#'
#' @return
#' \itemize{
#'   \item WaveletSeries - The wavelet trasnform of the series
#' }
#' @export
#'
#' @examples
#' data<-rnorm(200,mean=20,sd=3)
#' Wavelet<-WaveletFitting(ts=data,Wvlevels=3,WFilter='haar',bndry='periodic',FFlag=TRUE)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2007. Forecasting time series using wavelets. Internationa Journal of Wavelets, Multiresolution and Inforamtion Processing, 5, 709 to 724

#' \item Percival D. B. and Walden A. T. 2000. Wavelet Methods for Time-Series Analysis. Cambridge Univ. Press, U.K.

#' \item Paul R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319 to 327.
#' }
WaveletFitting <- function(ts,Wvlevels,WFilter='haar',bndry='periodic',FFlag=TRUE)
{
  mraout <- wavelets::modwt(ts, filter=WFilter, n.levels=Wvlevels,boundary=bndry, fast=FFlag)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  return(list(WaveletSeries=WaveletSeries))
}


#' @title Wavelet-RF Hybrid Model for Forecasting
#' @description The Wavelet Decomposition followed by Random Forest Regression (RF) models have been applied for time series forecasting. The maximum overlap discrete wavelet transform (MODWT) algorithm was chosen as it works for any length of the series. The series is first divided into training and testing sets. In each of the wavelet decomposed series, the  supervised machine learning approach namely random forest was employed to train the model. This package also provides accuracy metrics in the form of Root Mean Square Error (RMSE) and Mean Absolute Prediction Error (MAPE).
#' @param ts Univariate time series
#' @param tlag Number of lags
#' @param Waveletlevels The level of wavelet decomposition
#' @param WaveletFilter Wavelet filter use in the decomposition
#' @param boundary The boundary condition of wavelet decomposition
#' @param FastFlag The FastFlag condition of wavelet decomposition: True or False
#' @param SplitRatio Training and testing data split
#'
#' @import fracdiff forecast stats wavelets randomForest tsutils
#' @return
#' \itemize{
#'   \item TrainFittedValue - Fitted value of train data
#'   \item TestPredictedValue - Predicted value of test data
#'   \item AccuracyTable - RMSE and MAPE of train and test data
#' }
#' @export
#'
#' @examples
#'data<-rnorm(200,mean=20,sd=3)
#'WRF<-WaveletFittingRF(ts=data,tlag=2,Waveletlevels=3)
#' @references
#' \itemize{
#'\item Ding, Y., Zhang, W., Zhao, X., Zhang, L. and Yan, F., 2020. A hybrid random forest method fusing wavelet transform and variable importance for the quantitative analysis of K in potassic salt ore using laser-induced breakdown spectroscopy. Journal of Analytical Atomic Spectrometry, 35(6), 1131-1138.

#' \item Rezaali, M., Fouladi-Fard, R., Mojarad, H., Sorooshian, A., Mahdinia, M. and Mirzaei, N., 2021. A wavelet-based random forest approach for indoor BTEX spatiotemporal modeling and health risk assessment. Environmental Science and Pollution Research, 28(18), 22522-22535.

#' \item Paul, R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319-327.

#' \item Paul, R. K. and Birthal, P.S. 2015. Investigating rainfall trend over India using wavelet technique. Journal of Water and Climate Change, 7, 365 to 378.

#' \item Paul, R. K. 2015. ARIMAX-GARCH-WAVELET Model for forecasting volatile data. Model Assisted Statistics and Application, 10, 243 to 252.

#' }

WaveletFittingRF<- function(ts,tlag=ACF,Waveletlevels,WaveletFilter='haar',boundary='periodic',FastFlag=TRUE,SplitRatio=0.8)
{
  WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,WFilter=WaveletFilter,bndry=boundary,FFlag=FastFlag)$WaveletSeries
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ANN model to the Wavelet Coef                #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    N<-length(ts)
    limit<-(exp(2*1.96/sqrt(N-3)-1))/(exp(2*1.96/sqrt(N-3)+1))
    ACF<-sum(abs(acf(WS[,1])$acf)>limit)
    if(tlag==ACF)
    {lag<-ACF
    } else
    {
      lag<-tlag
    }

    lag_x_a<-lagmatrix(ts,lag=c(0:lag))
    lag_x<-lag_x_a[-c(1:lag),]

    train <- lag_x[c(1:(nrow(lag_x)*SplitRatio)),]
    test<- lag_x[-c(1:(nrow(lag_x)*SplitRatio)),]

    train_y<-train[,1]
    test_y<-test[,1]

    WaveletRFFit <- randomForest::randomForest(train[,1]~.,data=train)

    WaveletRFPredict <- stats::predict(WaveletRFFit)
    WaveletRFForecast <- stats::predict(WaveletRFFit,test)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletRFPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,WaveletRFForecast)
  }

  TrainFittedValue <- rowSums(AllWaveletPrediction,na.rm = T)
  TestPredictedValue <- rowSums(AllWaveletForecast,na.rm = T)
  AccuracyTable<-matrix(nrow=2, ncol=2)
  AccuracyTable[1,1]<-round(sqrt(mean((train_y-TrainFittedValue)^2)), digits = 4)
  AccuracyTable[1,2]<-round(mean(abs((train_y-TrainFittedValue)/train_y)), digits = 4)
  AccuracyTable[2,1]<-round(sqrt(mean((test_y-TestPredictedValue)^2)), digits = 4)
  AccuracyTable[2,2]<-round(mean(abs((test_y-TestPredictedValue)/test_y)), digits = 4)
  row.names(AccuracyTable)<-c("Train", "Test")
  colnames(AccuracyTable)<-c("RMSE", "MAPE")
  return(list(TrainFittedValue=TrainFittedValue,TestPredictedValue=TestPredictedValue, AccuracyTable=AccuracyTable))
}

