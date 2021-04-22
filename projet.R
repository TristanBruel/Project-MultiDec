source("specPdgrm.R")
source("data_multiDec.R")
source("functions.R")

library(signal)

signal = signal_multiDec()
fs=4096

data = data_multiDec(fs,signal$duration,signal$wvf_LHO,signal$wvf_LLO,
                     signal$wvf_VIR,ampl=10,verbose=T)

specH=specPdgrm(data$data_H$x,data$data_H$t,l=200,p=90,fs,main="LHO")
specL=specPdgrm(data$data_L$x,data$data_L$t,l=200,p=90,fs,main="LLO")
specV=specPdgrm(data$data_V$x,data$data_V$t,l=200,p=90,fs,main="VIRGO")