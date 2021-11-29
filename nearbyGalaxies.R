source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")


#########################
### Statistical model ###
#########################
fits_data = read.table("inputs/A-A_fits_data_g2.dat", sep = ",");# data to generate model
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
# Nearby Galaxies
galaxies_filename = "nearbyGalaxies/nearby_galaxies.csv"
galaxies = read.csv(galaxies_filename, stringsAsFactors=FALSE);

galaxy_names = galaxies$name;
skyPositions = cbind(galaxies$RA,galaxies$DEC)
distances = galaxies$distance

# Waveform
signal_name = c("s20.0--LS220")

# List of networks
networks=list(c("CEH","CEL"),c("ET1","ET2","ET3"),c("CEH","CEL","ET1","ET2","ET3"))
network_names=c("CEH+CEL","ET","CEH+CEL+ET")

fs = 4096
filtering_method = "spectrum"

t0 = 1302220800

# loop over N generation of noisy data
N = 100L

# number and size of time steps: event each hour over 24 hours
time_nb=25L
time_step=1

ind_net=0
for (detectors in networks){
  ind_net=ind_net+1
  print(paste("Network of detectors: ",str_c(detectors,collapse=', ')))
  nDet=length(detectors)
  
  ind_gal = 0
  for (gal in galaxy_names){
    result <- array(0,c(nrow=N*time_nb,ncol=6))
    ind_gal = ind_gal+1
    skyAngles = skyPositions[ind_gal,]
    dist = distances[ind_gal]
    dec = skyAngles[1]
    ra = skyAngles[2]
    skyPosition = c(dec,ra)
    
    for (dt in 1:time_nb){
      
      t = t0+(dt-1)*time_step*3600   # GPS time in seconds
      
      wvfs = signal_multiDec(dec=dec,ra=ra,t=t,fs=fs,signal=signal_name,detectors=detectors,
                             pbOff=TRUE,actPlot=FALSE,verbose=TRUE)
      
      true_data = wvfs$true_data
      startTime = 0.1   # signal starts 100ms after bounce (t=0)
      L = length(wvfs$time)   # number of samples
      
      l = 400L   # interval length to use for each FFT
      p = 90L   # overlapping percentage
      offset = as.integer(round((1-p/100)*l))   # offset between consecutive FTs
      transient = as.integer(0.05*fs+1);   # samples to ignore at the start and end (50ms)
      
      wData = matrix(0,nrow = L, ncol = nDet);
      psd = matrix(0,nrow = l/2+1, ncol = nDet);
      freq = fs*seq(0,1/2,by=1/l);
      
      # To use always the same random noise for all waveforms % detectors
      set.seed(1)
      
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
                                  skyPosition=skyPosition,t=t,
                                  integLength=l,offsetLength=offset,transientLength=transient,
                                  freqLim=c(0,2000),windowType='modifiedHann',startTime=startTime,
                                  logPow=TRUE,actPlot=FALSE,verbose=FALSE);
        r = likelihoods$std
        
        ### Reduced time window ###
        r2 = list();
        r2$t = subset(r$t, r$t<=1.5);
        r2$f = r$f;
        r2$E = r$E[1:length(r2$t),]
        
        out = covpbb_poly(r=r2, mod=fit, true_data=true_data, limFreq=c(1000),
                          actPlot=FALSE);
        
        result[i+(dt-1)*N,1]=time_step*(dt-1)
        result[i+(dt-1)*N,2]=out$covpbb[1,1]
        result[i+(dt-1)*N,3]=out$covpbb[1,2]
        result[i+(dt-1)*N,4]=out$residual[1,1]
        result[i+(dt-1)*N,5]=out$residual[1,2]
        result[i+(dt-1)*N,6]=out$residual[1,3]
      }
      
      ind1=1+(dt-1)*N
      ind2=N+(dt-1)*N
      print(sprintf("galaxy %s @ distance: %f kpc @ time: +%s hours. Covpbb mean:%f. Covpbb median: %f",
                    gal, dist, time_step*(dt-1), mean(result[ind1:ind2,2]), median(result[ind1:ind2,2])))
    }
    
    filename=sprintf("nearbyGalaxies/new/%s/results_AA_%s_f2_%s_%s.txt", 
                     network_names[ind_net], filtering_method, signal_name, gal)
    write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)
  }
}