source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")

library(signal)
library(stringr)

#########################
### Statistical model ###
#########################
# Load data to generate model
fits_data = read.table("inputs/1D_simulations/A-A_fits_data_g2.dat", sep = ",");
colnames(fits_data) = c("r", "f");

# Variable variance linear model
mu3 = "~ f + I(f^2) + I(f^3) - 1";
s2 = "~ f + I(f^2) - 1";

Xm  = model.matrix(eval(parse(text=eval(parse(text=mu3)))), fits_data);
Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = FALSE);


#############################
### Simulation parameters ###
#############################
# Sky position: direction inside the Sagittarius constellation
dec=-16.18;
ra=18.34;
skyPosition=c(dec,ra);
# Time of arrival at the center of Earth
t0=1325052478; #favourable case
t0=1325077869; #unfavourable case

# List of networks
networks=list(c("LHO","LLO"), c("LHO","LLO","VIR","KAG","LAO"));
network_names=c("HL", "HLVKA");

# List of waveforms
signals=c("s11.2--LS220", "s15.0--GShen", "s15.0--LS220", "s15.0--SFHo", 
          "s20.0--LS220", "s20.0--SFHo");

fs=4096;
#filtering_method=prewhiten;
filtering_method="spectrum";

# loop over N generation of noisy data
N=100;

# number and size of distance steps
dist_nb=61;
dist_steps=c(1.0,2.0,2.5,2.0,1.5,1.5,2.5,2.5);

ind_net=0;
for (detectors in networks){
  ind_net=ind_net+1;
  print(paste("Network of detectors: ",str_c(detectors,collapse=', ')));
  nDet=length(detectors);
  
  ind_sig=0;
  for (signal_name in signals){
    ind_sig = ind_sig+1;
    result<-array(0,c(nrow=N*dist_nb,ncol=6));
    
    wvfs = signal_multiDec(dec=dec,ra=ra,t=t0,fs=fs,signal=signal_name,detectors=detectors,
                           pbOff=TRUE,actPlot=FALSE,verbose=TRUE);
    true_data = wvfs$true_data;
    startTime = 0.1;  # signal starts 100ms after bounce (t=0)
    L = length(wvfs$time);   # number of samples
    
    l = 400L;   # interval length to use for each FFT
    p = 90L;   # overlapping percentage
    offset = as.integer(round((1-p/100)*l));   # offset between consecutive FTs
    transient = as.integer(0.05*fs+1);   # samples to ignore at the start and end (50ms)
    
    wData = matrix(0,nrow=L,ncol=nDet);
    psd = matrix(0,nrow=l/2+1,ncol=nDet);
    freq = fs*seq(0,1/2,by=1/l);
    
    # To use always the same random noise for all waveforms % detectors
    # set.seed(1);

    for (j in 1:dist_nb){
      dist = 0.001*(j==1)+(j-1)*dist_steps[ind_sig];
      
      # To use always the same random noise for all waveforms and detectors
      # at each distance step
      set.seed(1);
      
      for (i in 1:N){
        d = data_multiDec(fs=fs,wvfs=wvfs,ampl=10/dist,detectors=detectors, 
                          filter=filtering_method, setseed=0,
                          actPlot=FALSE, verbose=FALSE);
        
        for (k in 1:nDet){
          psd[,k] = PSD_fromfiles(freq,1,detectors[k]);
          wData[,k] = d[[k]]$y;   # whitened data
          wData[,k] = wData[,k]/sqrt(mean(psd[,k]));
        }
        
        likelihoods = timeFreqMap(fs=fs,wData=wData,detectors=detectors,psd=psd,
                                  skyPosition=skyPosition,t=t0,
                                  integLength=l,offsetLength=offset,transientLength=transient,
                                  freqLim=c(0,2000),windowType='modifiedHann',startTime=startTime,
                                  logPow=TRUE,actPlot=FALSE,verbose=FALSE);
        r = likelihoods$std;
        
        ### Reduced time window ###
        r2=list();
        r2$t = subset(r$t, r$t<=1.5);
        r2$f = r$f;
        r2$E = r$E[1:length(r2$t),];
        
        out = covpbb_LASSO(r=r2, mod=fit, true_data=true_data, limFreq=c(1000),
                          actPlot=FALSE);
        
        result[i+(j-1)*N,1]=dist;
        result[i+(j-1)*N,2]=out$covpbb[1,1];
        result[i+(j-1)*N,3]=out$covpbb[1,2];
        result[i+(j-1)*N,4]=out$residual[1,1];
        result[i+(j-1)*N,5]=out$residual[1,2];
        result[i+(j-1)*N,6]=out$residual[1,3];
      }
      
      ind1=1+(j-1)*N;
      ind2=N+(j-1)*N;
      print(sprintf("signal %s @ distance: %f kpc. Covpbb mean:%f. Covpbb median: %f",
                    signal_name, dist, mean(result[ind1:ind2,2]), median(result[ind1:ind2,2])));
    
    }
    save_dir=sprintf("./perf/2G/unfavourable/%s/", network_names[ind_net]);
    dir.create(path=save_dir, showWarnings=FALSE, recursive=TRUE);
    filename=sprintf("results_AA_%s_f2_%s.txt", filtering_method, signal_name);
    save_path=paste(save_dir, filename, sep='');
    write.table(result, file=save_path, sep=" ", row.names=FALSE, col.names=FALSE);
  }
}