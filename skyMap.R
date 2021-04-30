source("multiDec_algebra.R")
library(Matrix)
library(fields)

detectors = c("LHO", "LLO", "VIR", "KAG")

t = 1302220800  # GPS time 12/04/2021 00:00:00

m = length(detectors)
list_dec = seq(0,360,by=1)
list_ra = seq(0,24,by=1)
res = matrix(0,361,25)
for (k in detectors){
  for (dec in 1:361){
    for (ra in 1:25){
      coeff = grav_response(dec, ra, t, pol=0, k);
      Fplus = coeff$Fplus
      Fcross = coeff$Fcross
      res[dec,ra] = Fplus^2 + Fcross^2
    }
  }
image.plot(list_dec, list_ra, res, main=c("Detector",k),
           xlab="Declination [°]", xaxp=c(0,360,6),
           ylab="Right ascension [hours]", yaxp=c(0,24,6))
}