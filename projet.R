source("specPdgrm.R")
source("data_multiDec.R")
source("functions.R")

library(signal)

signal_name="KURODA_TM1_H_resampled.dat"
signal = signal_multiDec(signal=signal_name)
fs=4096
filtering_method="prewhiten"

wvf_LHO=signal$wvf_LHO
wvf_LLO=signal$wvf_LLO
wvf_VIR=signal$wvf_VIR
duration=signal$duration

data = data_multiDec(fs,duration,wvf_LHO,wvf_LLO,wvf_VIR,ampl=10,verbose=TRUE)

specH=specPdgrm(data$data_H$x,data$data_H$t,l=200,p=90,fs,main="LHO")
specL=specPdgrm(data$data_L$x,data$data_L$t,l=200,p=90,fs,main="LLO")
specV=specPdgrm(data$data_V$x,data$data_V$t,l=200,p=90,fs,main="VIRGO")