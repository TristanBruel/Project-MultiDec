source("data_multiDec.R")
source("inverse.R")

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
nDet = length(detectors)

signals = signal_multiDec(dec,ra,t,signal=signal_name,detectors=detectors,
                          actPlot=FALSE);
N = fs*signals$duration+1;

SNR = 30

psd = matrix(0,nrow = N, ncol = nDet);
wData = matrix(0,nrow = N, ncol = nDet);

for (k in 1:nDet){
  wvf = signals[[k]]
  ### Noise to get desired SNR ###
  noisy = whiteNoise(wvf,SNR=SNR,main=detectors[k],actPlot=FALSE);
  psd[,k] = rep(noisy$psd, N);
  wData[,k] = noisy$noisy/sqrt(psd[,k]);
}

estimates = inverse(fs=fs,wData=wData,detectors=detectors,psd=psd,
                    skyPosition=c(dec,ra,t))
