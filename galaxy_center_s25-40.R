source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")


#########################
### Statistical model ###
#########################
fits_data = read.table("inputs/1D_simulations/A-A_fits_data_g2.dat", sep = ",");# data to generate model
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
# Sky position : Galactic center (SIMBAD catalog http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Galactic+Centre )
dec=-29.006;
ra=17.761;
skyPosition = c(dec,ra);
dist=8.2; # distance to the galactic center (in kpc)

# First arrival time at the center of Earth
t0=1325048418

# List of networks
networks=list(c("LHO","LLO"),c("LHO","LLO","VIR"),c("LHO","LLO","VIR","KAG"),
              c("LHO","LLO","VIR","KAG","LAO"),c("LHO","LLO","VIR","LAO"));
network_names=c("HL","HLV","HLVK","HLVKA","HLVA");

# List of waveforms
signals=c("s25.0--LS220", "s40.0--LS220");

fs=4096;
#filtering_method="prewhiten";
filtering_method="spectrum";

# loop over N generation of noisy data
N=100L;

# number and size of time steps: event each 30min over 24 hours
time_nb=49L;
time_step=0.5;

ind_net=0;
for (detectors in networks){
  ind_net=ind_net+1;
  print(paste("Network of detectors: ",str_c(detectors,collapse=', ')));
  nDet=length(detectors);
  
  for (signal_name in signals){
    result<-array(0,c(nrow=N*time_nb,ncol=6));
    
    for (dt in 1:time_nb){
      
      t = t0+(dt-1)*time_step*3600;  # GPS time in seconds
      
      wvfs = signal_multiDec(dec=dec,ra=ra,t=t,fs=fs,signal=signal_name,detectors=detectors,
                             pbOff=TRUE,actPlot=FALSE,verbose=TRUE);
      
      true_data = wvfs$true_data;
      startTime = 0.1;   # signal starts 100ms after bounce (t=0)
      L = length(wvfs$time);   # number of samples
      
      l = 400L;   # interval length to use for each FFT
      p = 90L;   # overlapping percentage
      offset = as.integer(round((1-p/100)*l));   # offset between consecutive FTs
      transient = as.integer(0.05*fs+1);   # samples to ignore at the start and end (50ms)
      
      wData = matrix(0,nrow = L, ncol = nDet);
      psd = matrix(0,nrow = l/2+1, ncol = nDet);
      freq = fs*seq(0,1/2,by=1/l);
      
      # To use always the same random noise for all waveforms % detectors
      set.seed(1);
      
      for (i in 1:N){
        d = data_multiDec(fs=fs,wvfs=wvfs,ampl=10/dist,detectors=detectors, 
                          filter=filtering_method, setseed=0,
                          actPlot=FALSE, verbose=FALSE);
        
        for (k in 1:nDet){
          psd[,k] = PSD_fromfiles(freq,1,detectors[k],actPlot=FALSE);
          wData[,k] = d[[k]]$y;   # whitened data
          wData[,k] = wData[,k]/sqrt(mean(psd[,k]));
        }
        
        likelihoods = timeFreqMap(fs=fs,wData=wData,detectors=detectors,psd=psd,
                                  skyPosition=skyPosition,t=t,
                                  integLength=l,offsetLength=offset,transientLength=transient,
                                  freqLim=c(0,2000),windowType='modifiedHann',startTime=startTime,
                                  logPow=TRUE,actPlot=FALSE,verbose=FALSE);
        r = likelihoods$std;
        
        ### Reduced time window ###
        r2 = list();
        r2$t = subset(r$t, r$t<=1.5);
        r2$f = r$f;
        r2$E = r$E[1:length(r2$t),];
        
        out = covpbb_LASSO(r=r2, mod=fit, true_data=true_data, limFreq=c(1000),
                           mask_f=c(600,700),
                           actPlot=FALSE);
        
        result[i+(dt-1)*N,1]=time_step*(dt-1);
        result[i+(dt-1)*N,2]=out$covpbb[1,1];
        result[i+(dt-1)*N,3]=out$covpbb[1,2];
        result[i+(dt-1)*N,4]=out$residual[1,1];
        result[i+(dt-1)*N,5]=out$residual[1,2];
        result[i+(dt-1)*N,6]=out$residual[1,3];
      }
      
      ind1=1+(dt-1)*N;
      ind2=N+(dt-1)*N;
      print(sprintf("signal %s @ time: +%s hours. Covpbb mean:%f. Covpbb median: %f",
                    signal_name, time_step*(dt-1), mean(result[ind1:ind2,2]), median(result[ind1:ind2,2])));
    }
    
    save_dir=sprintf("./galaxyCenter/%s/", network_names[ind_net]);
    dir.create(path=save_dir, showWarnings=FALSE, recursive=TRUE);
    filename=sprintf("results_AA_%s_f2_%s.txt", filtering_method, signal_name);
    save_path=paste(save_dir, filename, sep='');
    write.table(result, file=save_path, sep=" ", row.names=FALSE, col.names=FALSE);
  }
}