source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")

#############################
### Simulation parameters ###
#############################
# Sky position : Galactic center (SIMBAD catalog http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Galactic+Centre)
dec=-29.006;
ra=17.761;
skyPosition = c(dec,ra);
dist=8.2; # distance to the galactic center (in kpc)

# Time of arrival at the center of Earth
t0=1325048418;

# Network of detectors
detectors=c("LHO", "LLO", "VIR", "KAG", "LAO")
nDet=length(detectors)

# Waveform to use
signal_name = "s20.0--LS220"
#signal_name = "s15--3D_eqtr"

sig=signal_multiDec(dec=dec,ra=ra,t=t0,fs=4096,signal=signal_name,detectors=detectors,
                    pbOff=TRUE,actPlot=FALSE,verbose=FALSE)

F = antenna_patterns(dec,ra,t0,detectors=detectors)
print(F)

SNRs = rep(0,nDet)
for (d in 1:nDet){
  wvf = list(hoft=sig[[d]], time=sig$time);
  SNR=compute_SNR(wvf, detector=detectors[d], fcut=10, dist=dist, actPlot=FALSE)
  SNRs[d]=SNR
  print(sprintf("SNR for wvf %s in detector %s @ distance %f: %f", 
                signal_name, detectors[d], dist, SNR));
}
SNR_net = sqrt(sum(SNRs**2))
print(sprintf("SNR for wvf %s in network HLVKA @ distance %f: %f", 
              signal_name, dist, SNR_net))