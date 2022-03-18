source("multiDec_algebra.R")

#############################
### Simulation parameters ###
#############################
# Sky position : Galactic center (SIMBAD catalog http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Galactic+Centre)
dec=-29.006
ra=17.761

# Time of arrival at the center of Earth
t0=1325048418

detectors=c("LHO","LLO","VIR","KAG","LAO")
nDet=length(detectors)

# number and size of time steps: event each 30min over 24 hours
time_nb=49L
time_step=0.5

result=array(0,dim=c(nrow=time_nb,ncol=2,nDet))

for (dt in 1:time_nb){
  t = t0+(dt-1)*time_step*3600 # GPS time in seconds
  patterns=antenna_patterns(dec=dec,ra=ra,t=t,detectors=detectors);
  result[dt,,]=t(patterns)
}

filename="galaxyCenter/antenna_patternsHLVKA.txt"
write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)