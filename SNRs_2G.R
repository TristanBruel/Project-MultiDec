source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")

library(signal)
library(stringr)


#############################
### Simulation parameters ###
#############################
# Sky position: direction inside the Sagittarius constellation
dec=-16.18;
ra=18.34;
skyPosition=c(dec,ra);
# Time of arrival at the center of Earth
t0=1330480818; #favourable case
#t0=1326628818; #unfavourable case

# List of networks
networks=list(c("LHO","LLO"),c("LHO","LLO","VIR","KAG","LAO"));
network_names=c("HL","HLVKA");

# List of waveforms
signals=c("s11.2--LS220", "s15.0--GShen", "s15.0--LS220", "s15.0--SFHo", 
          "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220",
          "s15--3D_eqtr", "s15--3D_pole");

fs=4096;

# number and size of distance steps
dist_nb=61;
dist_steps=c(0.5,1.0,1.5,1.0,1.0,1.0,2.5,2.5,0.5,0.5);

ind_net=0;
for (detectors in networks){
  ind_net=ind_net+1;
  print(paste("Network of detectors: ",str_c(detectors,collapse=', ')));
  nDet=length(detectors);
  
  ind_sig=0;
  for (signal_name in signals){
    ind_sig = ind_sig+1;
    SNRs <- array(0,c(nrow=dist_nb,ncol=nDet));
    
    wvfs = signal_multiDec(dec=dec,ra=ra,t=t0,fs=fs,signal=signal_name,detectors=detectors,
                           pbOff=TRUE,actPlot=FALSE,verbose=FALSE);
    
    # To use always the same random noise for all waveforms % detectors
    # set.seed(1);
    
    for (d in 1:nDet){
      wvf = list(hoft=wvfs[[d]], time=wvfs$time)
      for (j in 1:dist_nb){
        dist = 0.001*(j==1)+(j-1)*dist_steps[ind_sig]; # avoid dist=0 which creates divergence
        SNR = compute_SNR(wvf = wvf,
                          detector = detectors[d],
                          dist = dist);
        SNRs[j, d] = SNR;
      }
    }
    print(sprintf("SNRs for signal %s in network %s computed.", signal_name, network_names[ind_net]));
    
    save_dir=sprintf("./perf/2G/favourable/%s/", network_names[ind_net]);
    dir.create(path=save_dir, showWarnings=FALSE, recursive=TRUE);
    filename=sprintf("SNRs_%s_%s.txt", signal_name, network_names[ind_net]);
    save_path=paste(save_dir, filename, sep='');
    write.table(SNRs, file=save_path, sep=" ", row.names=FALSE, col.names=detectors);
  }
}