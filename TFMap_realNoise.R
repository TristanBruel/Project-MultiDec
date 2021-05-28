source("data_multiDec.R")
source("timeFreqMaps.R")

library(signal)
library(Matrix)

###########################
##### Inverse Problem #####
###########################

### REAL NOISE ###

### Source ###
dec=30
ra=6
t=1302220800
skyPosition=c(dec,ra,t)

### Signal ###
signal_name="KURODA_TM1_H_resampled.dat"
fs=4096
detectors = c("LHO", "LLO", "VIR", "KAG")
nDet = length(detectors)

signals = signal_multiDec(dec, ra, t,signal=signal_name, detectors = detectors,
                          actPlot = FALSE);
N = fs*signals$duration+1   # number of samples

dist = 1
for (k in 1:nDet){
  snr=compute_SNR(signals[[k]],dist=dist)
  print(sprintf("SNR in detector %s is %s",detectors[k],snr))
}
data = data_multiDec(fs,signals,signals$duration, detectors=detectors,
                     ampl=10/dist,verbose=FALSE,actPlot=FALSE);

l = 100   # interval length to use for each FFT
offset=10   # 90% overlapping
padd = floor(0.02*fs+1);   # samples to ignore at the start and at the end (20ms)

wData = matrix(0,nrow = N, ncol = nDet);
psd = matrix(0,nrow = l/2+1, ncol = nDet);

freq=fs*seq(0,1/2,by=1/l);

for (k in 1:nDet){
  psd[,k]=PSD_fromfiles(freq,1,detectors[k])
  wData[,k] = data[[k]]$y;   # prewhiten data
  wData[,k] = wData[,k]/sqrt((fs/2)*mean(psd[,k]));
}

likelihoods=timeFreqMap(fs,wData,detectors,psd,skyPosition,l,10,padd,
                        freqBand=c(10,Inf),windowType='bartlett',actPlot=TRUE)