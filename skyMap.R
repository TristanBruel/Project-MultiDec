source("multiDec_algebra.R")
#source("test_algebra2.R")
library(Matrix)
library(fields)

#detectors = c("LHO", "LLO", "VIR", "KAG", "LIO", "ET1", "ET2", "ET3");
detectors = c("LHO", "LLO", "VIR", "KAG", "LIO");

t = 1302220800  # GPS time 12/04/2021 00:00:00

n = 100

res = matrix(0,2*n,n)
list_dec = seq(-90,by=180/n,length=n);
list_ra = seq(0,by=24/(2*n),length=2*n);

for (k in detectors){
  for (ra in 1:(2*n)){
    for (dec in 1:n){
      F=antenna_patterns(-90+(dec-1)*180/n, (ra-1)*24/(2*n), t, pol=0, k);
      #res[ra,dec] = sqrt(F[1]^2 + F[2]^2);
      res[ra,dec] = F[1];
    }
  }
  filename = sprintf("Plots/skyMaps/Fp_%s.png",k)
  png(filename)
  image.plot(list_ra, list_dec, res, main=paste("Fplus ",k),
             xlab="Right ascension [h]", xaxp=c(0,24,6),
             ylab="Declination [°]", yaxp=c(-90,90,6))
  dev.off()
  
  #filename=sprintf("Plots/skyMaps/Fp_%s.txt",k)
  #write.table(t(res), file=filename, sep=" ", row.names=FALSE, col.names=FALSE)
  }
