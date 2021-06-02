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
  if (dim(wData)[2] != nDet){
    stop("Number of detectors inconsistent with the data provided")
  }
  
  ### Detectors information ###
  dec=skyPosition[1];
  ra=skyPosition[2];
  t=skyPosition[3];
  F=antenna_patterns(dec,ra,t,0,detectors);
  delays=time_delays(dec,ra,t,detectors);
  # Reset reference position to first detector
  delays=delays-delays[1]
  delayLengths=round(delays*fs);
  if (verbose){
    print("Antenna response matrix : F")
    print(F)
  }
  
  ### Data process ###
  # Samples to ignore at the start and at the end (20ms)
  transient = floor(0.02*fs);
  start=transient;
  end=length(wData[,1])-transient-1;
  L=end-start+1;
  auxdata=zeros(L,nDet);
  # Time alignment
  for (k in 1:nDet){
    auxdata[,k]=wData[(start+delayLengths[k]):(end+delayLengths[k]),k];
  }
  
  # Fourier transform
  dataf=mvfft(auxdata);
  
  ### Resolution : standard likelihood ###
  fhplus = rep(0,L);
  fhcross = rep(0,L);
  
  diag = rep(0,nDet);
  
  for (f in 1:L){
    # Noise-spectrum-weighted antenna responses
    for (k in 1:nDet){
      diag[k] = 1/sqrt(psd[f,k])
    }
    W = diag(diag);
    wF = W %*% F
    
    # Dominant Polarization Frame
    wFp=wF[,1];
    wFc=wF[,2];
    DPF=convertToDPF(wFp,wFc);
    
    # Inverse
    d=dataf[f,];
    res = solve(Conj(t(DPF)) %*% DPF, Conj(t(DPF)) %*% d);
    fhplus[f] = res[1];
    fhcross[f] = res[2];
  }
  
  hplus = Re(fft(fhplus, inverse = TRUE))/L;
  hcross = Re(fft(fhcross, inverse = TRUE))/L;
  
  ### Plots ###
  if (actPlot){
    gw_filename=paste("Waveforms/",signal_name,sep="")
    sXX = read.table(gw_filename) # V1 time, V2 hplus, V3 hcross
    colnames(sXX) = c ("time","hplus","hcross")
    fs_orig = round(1/(sXX$time[2]-sXX$time[1]));
    n = length(sXX$time);
    ind0 = floor((L-n)/2)+1;
    indn = L-floor((L-n)/2);
    
    plot(sXX$time,sXX$hplus,type='l',xlab="Time [s]",ylab="hplus")
    points(sXX$time,hplus[ind0:indn],type='l',col='red',pch=2)
    leg = c("True", "Estimated")
    col = c("black","red")
    legend(x=sXX$time[1]*1.1,y=max(sXX$hplus)*.9,legend=leg,col=col,pch=c(1,2))
    
    plot(sXX$time,sXX$hcross,type='l',xlab="Time [s]",ylab="hcross")
    points(sXX$time,hcross[ind0:indn],type='l',col='red',pch=2)
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