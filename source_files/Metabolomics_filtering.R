## Code from Courtney Schiffman on github to filter metabolomics data

## 'calc.icc' is a function that calculates ICC values for each features using biological and QC samples
## filled - a filled dataframe (as from xcms) of logged feature abundances. The dataframe should have features in the rows and samples in the columns (including QC samples), with the names of features being the rownames of the matrix and the sample names being the column names of the matrix. 
## biosamples - chacater string of names corresponding to the biological samples in the filled dataframe.
## qcsamples - chacater string of names corresponding to the QC samples in the filled dataframe.

calc.icc <- function(filled,biosamples,qcsamples){
  vars1 <- foreach(i=(1:nrow(filled)),.packages='nlme',.combine='rbind')%dopar%{
    reps <- factor(c(1:length(biosamples),rep(length(biosamples)+1,length(qcsamples))))
    data <- data.frame(y=c(as.numeric(filled[i,biosamples]),as.numeric(filled[i,qcsamples])),reps=reps)
    mm <- lme(y~ 1, random = ~1|reps, data=data,na.action = na.omit)
    return(as.numeric(VarCorr(mm)[1:2]))
  }
  ICC <- apply(vars1,1,function(x) x[1]/sum(x))
  names(ICC) <- rownames(filled)
  return(ICC)
}


## 'boxplot.icc' is a function that creates boxplots and boxplot statistics for high and low quality features that can be used to determine filtering thresholds

boxplot.icc <- function(icc,high,low){
  bp.high <- boxplot(icc[names(icc)%in%high],main="ICC, High quality")
  bp.low <- boxplot(icc[names(icc)%in%low],main="ICC, Low quality")
  boxplot(icc[names(icc)%in%c(high,low)] ~ as.factor(names(icc[names(icc)%in%c(high,low)])%in%high),notch=T
          ,names=c('Low Quality','High Quality'),main='Box Plot of Intra-class Correlation Coefficient',
          cex.lab=1.4,cex.main=1.5,ylab='icc',cex.axis=1.4,varwidth=T)
  return(list(LowQuality=bp.low$stats,HighQuality=bp.high$stats))
  
}

## 'density.icc' is a function that creates density plots of ICC values for high and low quality features

density.icc <- function(icc,high,low){
  plot(density(icc[names(icc)%in%high]*100),main='Density Plot of ICC Values',
       cex.lab=1.4,cex.main=1.5,ylab='ICC',cex.axis=1.4,xlab='',col='red')
  lines(density(icc[names(icc)%in%low]*100))
  legend('topleft',col=c('red','black'),lwd=2,legend=c('High Quality','Low Quality'),cex=1.2)
}

## 'filter.icc' is a function that identifies which features pass the specified icc threshold.

filter.icc <- function(icc,threshold){
  names(icc[icc<threshold])
}

## Mean centering
center_scale <- function(x) {
  scale(x, scale = FALSE)
}
