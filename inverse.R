source("data_multiDec.R")
source("specPdgrm.R")
source("multiDec_algebra.R")

library(signal)
library(Matrix)

########################
### Inverse function ###
########################


########################################################################
inverse = function(fs=4096, wData, detectors=c("LHO","LLO","VIR"), psd,
                          skyPosition, actPlot=TRUE, verbose=TRUE){
  ######################################################################
  # Inputs :  fs: sampling frequency
  #           wData: matrix of time-domain whitened data, one column per detector
  #           detectors: vector of detectors used in the analysis
  #           psd: power spectral densities of the detectors noises, one column per detector
  #                 Must be sampled at frequencies fs*[0:1/2, ]
  #           skyPosition: vector (dec, ra, t)
  #
  # Output :  hplus, hcross
  ######################################################################
  
  ### Command line arguments ###
  nDet=length(detectors);
  if (length(wData) < nDet){
    stop("Number of detectors inconsistent with the data provided")
  }
  nDat=length(wData[,1]);
  
  ### Detectors information ###
  dec=skyPosition[1];
  ra=skyPosition[2];
  t=skyPosition[3];
  F=antenna_patterns(dec,ra,t,0,detectors);
  delays=time_delays(dec,ra,t,detectors);
  # Reset reference position to first detector
  delays=delays-delays[1]
  if (verbose){
    print("Antenna response matrix : F")
    print(F)
  }
  
  ### Data process ###
  # Time alignment
  for (k in 1:nDet){
    wData[,k]=wData[,k]+delays[k]
  }
  
  # Fourier transform
  dataf=mvfft(wData);
  
  ### Resolution : standard likelihood ###
  fhplus = rep(0,T);
  fhcross = rep(0,T);
  
  diag = rep(0,nDet);
  
  for (f in 1:T){
    
    # Noise-spectrum-weighted antenna responses
    for (k in 1:nDet){
      diag[k] = 1/sqrt(psd[f,k])
    }
    W = diag(diag);
    wF = W %*% F
    
    # Dominant Polarization Frame
    wFp=wF[,1];
    wFc=wF[,2];
    DPF=convertToDPF(wFp,wFc)
    
    # Inverse
    d=dataf[f,]
    res = solve(Conj(t(DPF)) %*% DPF, Conj(t(DPF)) %*% d)
    fhplus[f] = res[1]
    fhcross[f] = res[2]
  }
  
  hplus = Re(fft(fhplus, inverse = TRUE))/T;
  hcross = Re(fft(fhcross, inverse = TRUE))/T;
  
  ### Plots ###
  if (actPlot){
    gw_filename=paste("Waveforms/",signal_name,sep="")
    sXX = read.table(gw_filename) # V1 time, V2 hplus, V3 hcross
    colnames(sXX) = c ("time","hplus","hcross")
    
    plot(sXX$time,sXX$hplus,type='l',xlab="Time [s]",ylab="hplus")
    points(hplus,type='l',col='red',pch=2)
    leg = c("True", "Estimated")
    col = c("black","red")
    legend(x=sXX$time[1]*1.1,y=max(sXX$hplus)*.9,legend=leg,col=col,pch=c(1,2))
    
    plot(sXX$time,sXX$hcross,type='l',xlab="Time [s]",ylab="hcross")
    points(hcross,type='l',col='red',pch=2)
    leg = c("True", "Estimated")
    col = c("black","red")
    legend(x=sXX$time[1]*1.1,y=max(sXX$hcross)*.9,legend=leg,col=col,pch=c(1,2))
    
    spec=specPdgrm(hplus,sXX$time,l=100,p=90,fs,main="Estimated hplus");
    #out=findGmodes(spec, l=100, p=90, gmode="right", actPlot = TRUE);
    spec=specPdgrm(hcross,sXX$time,l=100,p=90,fs,main="Estimated hcross");
    #out=findGmodes(spec, l=100, p=90, gmode="right", actPlot = TRUE);
  }
  
  return(list(hplus=hplus,hcross=hcross))
}