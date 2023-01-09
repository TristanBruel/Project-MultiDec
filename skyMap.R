source("multiDec_algebra.R")
library(fields)

detectors = c("LHO", "LLO", "VIR", "KAG", "LAO");
detectors = c("VIR")

t0 = 1325030418  # GPS time 01/01/2022 00:00:00

n = 100

res = matrix(0,2*n,n)
list_dec = seq(-90,by=180/n,length=n);
list_ra = seq(0,by=24/(2*n),length=2*n);

save_dir="./plots/skyMaps/gif/";
dir.create(path=save_dir, showWarnings=FALSE, recursive=TRUE);

for (det in detectors){
  for (dt in 0:24){
    t = t0 + dt*3600
    for (ra in 1:(2*n)){
      for (dec in 1:n){
        F=antenna_patterns(dec=-90+(dec-1)*180/n, ra=(ra-1)*24/(2*n), t=t, pol=0, detectors=det);
        res[ra,dec] = sqrt(F[1]^2 + F[2]^2);
        #res[ra,dec] = F[2];
      }
    }
    filename = sprintf("Feq_%s_h%s.png",det,dt)
    save_path=paste(save_dir, filename, sep='');
    png(save_path)
    imagePlot(list_ra, list_dec, res, main=paste("Feq",det),
              xlab="Right ascension [h]", xaxp=c(0,24,6),
              ylab="Declination [°]", yaxp=c(-90,90,6),
              zlim=c(0.,1.))
    dev.off()
    
    #write.table(t(res), file=filename, sep=" ", row.names=FALSE, col.names=FALSE)
  }
}
