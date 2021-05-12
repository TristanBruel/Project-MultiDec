source("data_multiDec.R")
source("specPdgrm.R")
source("multiDec_algebra.R")

library(signal)
library(Matrix)

###########################
##### Inverse Problem #####
###########################

### GAUSSIAN NOISE ###

### Source ###
dec=30
ra=6
t=1302220800

### Signal ###
signal_name="KURODA_TM1_H_resampled.dat"
fs=4096
detectors = c("LHO", "LLO", "VIR", "KAG")

signals = signal_multiDec(dec, ra, t,signal=signal_name, detectors = detectors,
                         actPlot = FALSE);

### Inverse function ###
########################

inverse = function(signals,SNR=50, detectors=c("LHO","LLO","VIR","KAG"),
                          actPlot=TRUE, verbose=TRUE){
  
  m = length(detectors)
  T = fs*signals$duration+1
  F = matrix(0,m*T,2*T)   # contains the antenna responses
  S = rep(0,m*T)   # contains the FT signals

  for (k in 1:m){
    wvf = signals[[k]]
    
    # Time alignment
    wvf$time = wvf$time + time_delay(dec, ra, t, detectors[k]);
    
    ### Noise to get desired SNR ###
    noisy = whiteNoise(wvf,SNR=SNR,main=detectors[k],actPlot = actPlot);

    ### Fourier Transform ###
    S[((k-1)*T+1):(k*T)] = fft(noisy$noisy);
    
    ### Antenna patterns ###
    coeff = grav_response(dec,ra,t,0,detectors[k]);
    A = cbind(diag(rep(coeff$Fplus,T)),diag(rep(coeff$Fcross,T)));
    A = mvfft(A);   # each column replaced by its discrete Fourier transform
    F[((k-1)*T+1):(k*T),1:T] = A[1:T,1:T]
    F[((k-1)*T+1):(k*T),(T+1):(2*T)] = A[1:T,(T+1):(2*T)] 
  }

  ### Resolution : LS estimate ###
  hplus = rep(0,T);
  hcross = rep(0,T);
  
  res = solve(Conj(t(F)) %*% F, Conj(t(F)) %*% S)
  hplus = Re(res[1:T])
  hcross = Re(res[(T+1):(2*T)])
  
  ### Plots ###
  if (actPlot == TRUE){
    gw_filename=paste("Waveforms/",signal_name,sep="")
    sXX = read.table(gw_filename) # V1 time, V2 hplus, V3 hcross
    colnames(sXX) = c ("time","hplus","hcross")
    
    plot(sXX$time,sXX$hplus,type='l',xlab="Time [s]",ylab="hplus")
    points(sXX$time,hplus,type='l',col='red',pch=2)
    leg = c("True", "Estimated")
    col = c("black","red")
    legend(x=sXX$time[1]*1.1,y=max(sXX$hplus)*.9,legend=leg,col=col,pch=c(1,2))
    
    plot(sXX$time,sXX$hcross,type='l',xlab="Time [s]",ylab="hcross")
    points(sXX$time,hcross,type='l',col='red',pch=2)
    leg = c("True", "Estimated")
    col = c("black","red")
    legend(x=sXX$time[1]*1.1,y=max(sXX$hcross)*.9,legend=leg,col=col,pch=c(1,2))
    
    spec=specPdgrm(hplus,sXX$time,l=100,p=90,fs,
                   main=c("Estimated hplus spectrogram : SNR=",SNR));
    spec=specPdgrm(hcross,sXX$time,l=100,p=90,fs,
                   main=c("Estimated hcross spectrogram : SNR=",SNR));
  }
  
  return(list(hplus=hplus,hcross=hcross))
}


########################################################################
whiteNoise = function(waveform, SNR, fs=4096, fcut=0, 
                      actPlot = FALSE, main=NULL, verbose=FALSE){
  ########################################################################
  
  # Compute the constant psd required to get the desired SNR
  n=length(waveform$hoft)
  a=nextpow2(10*n)         #zero padding and rounding to the next power of 2
  n2=2^a
  freq2 = fs*fftfreq(n2)         # two-sided frequency vector
  freq1=freq2[1:int(n2/2)]
  
  vec=rep(0,n2)
  for (i in 1:n){
    vec[n2/4+i]=vec[n2/4+i]+waveform$hoft[i]
  }  
  
  hf=fft(vec)/sqrt(fs);     # normalization wrt the sampling
  
  hf=hf[1:(n2/2)]    # integral performed over positive frequencies
  
  hf=subset(hf,freq1-fcut>0);
  freq1=subset(freq1, freq1-fcut>0);
  
  integrand=abs(hf*Conj(hf))
  p=integrand/fs
  psd=4*trapz(freq1,p)/SNR^2
  
  # Create Noise
  X = rnorm(n, mean=0, sd=1);
  XX = fft(X);
  XXX = XX*sqrt(psd)*sqrt(fs);
  Y = fft(XXX, inverse = TRUE);
  Y = Re(Y)/n;
  
  for (i in 1:n){
    Y[i] = Y[i] + waveform$hoft[i]
  }
  
  if (verbose ==TRUE){
    print(sprintf("Gaussian white noise on detector %s with mean %g and standard deviation %g",
                  main, mean(Y), sd(Y)));
  }
  
  if (actPlot == TRUE){
    plot(waveform$time,Y,type='l',col='black',
         xlab = "Time [s]",ylab="Hoft",main=main)
    points(waveform$time,waveform$hoft,type='l',col='red')
    leg = c("Noisy signal", "Signal")
    col = c("black","red")
    legend ("topleft",legend=leg,cex=.8,col=col,pch=c(1,2))
  }
  
  return(list(noisy=Y,psd=psd))
}