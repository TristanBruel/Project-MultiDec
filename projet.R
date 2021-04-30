source("specPdgrm.R")
source("data_multiDec.R")
source("functions.R")
source("multiDec_algebra.R")

library(signal)

### Source ###
dec=-65
ra=8
t=1126259462.0

signal_name="KURODA_TM1_H_resampled.dat"
detectors=c("LHO", "LLO", "VIR", "KAG")
signal = signal_multiDec(dec, ra, t, signal=signal_name, detectors=detectors)


fs=4096
filtering_method="prewhiten"
gmode=c("right")

data = data_multiDec(fs,signal,signal$duration, detectors=detectors,
                     ampl=10,verbose=TRUE);


for (k in 1:length(detectors)){
spec=specPdgrm(data[[k]]$y,data[[k]]$t,l=100,p=90,fs,
                main=detectors[k]);
out=findGmodes(spec, l=100, p=90, gmode=gmode, actPlot = TRUE);
}