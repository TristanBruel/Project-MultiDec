source("specPdgrm.R")
source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")

library(signal)

#########################
### Statistical model ###
#########################
fits_data = read.table("inputs/A-A_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("inputs/CoCo_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("inputs/fits_data.dat", sep = ",");# data to generate model
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
# Sky position
dec=60
ra=8
t=1302220800
skyPosition=c(dec,ra)

detectors=c("LHO","LLO","VIR")
detectors=c("ET1","ET2","ET3")
nDet=length(detectors)

signals=c("s11.2--LS220", "s15.0--GShen", "s15.0--LS220", "s15.0--SFHo", 
          "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220")
signals = c("s20.0--LS220")

freq_min=200
distance_step=3

fs=4096
filtering_method="spectrum"

# loop over N generation of noisy data

N=100
#N=1
dist_nb=60
#dist_nb=1
result<-array(0,c(nrow=N*dist_nb,ncol=6))

for (signal_name in signals){
  
  wvfs = signal_multiDec(dec,ra,t,fs,signal=signal_name,detectors=detectors,
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
  
  for (j in 1:dist_nb){
    dist = 1+(j-1)*distance_step
    #dist=3
    for (i in 1:N){
      d = data_multiDec(fs,wvfs,ampl=10/dist,detectors=detectors, 
                      filter=filtering_method, setseed=0,
                      actPlot=FALSE, verbose=FALSE);
      
      for (k in 1:nDet){
        psd[,k] = PSD_fromfiles(freq,1,detectors[k]);
        wData[,k] = d[[k]]$y;   # whitened data
        wData[,k] = wData[,k]/sqrt(mean(psd[,k]));
      }
      
      likelihoods = timeFreqMap(fs,wData,detectors=detectors,psd,skyPosition,t,
                                l,offset,transient,freqBand=c(0,2000),
                                windowType='modifiedHann',startTime=startTime,
                                logPow=TRUE,actPlot=FALSE,verbose=FALSE);
      
      if (nDet==1){
        r = likelihoods$spec
      }else{
        r = likelihoods$std
      }
      
      ### time window : 0 to 1.5s ###
      r2 = list();
      r2$t = subset(r$t, r$t<=1.5);
      r2$f = r$f;
      r2$E = r$E[1:length(r2$t),]

      #r2 = thresh_abs(r2,limit=2.3,actPlot=TRUE);
      
      #max = findGmodes_poly(r2, N=5, setStart=TRUE, m_L=2, initfreq_L=c(freq_min,500), actPlot=TRUE);
      #max2 = findGmodes_LASSO(r2, lambda=1, setStart=FALSE, m_L=2, initfreq_L=c(freq_min,500), actPlot = TRUE);
      
      out = covpbb_poly(r2, mod=fit, setStart=FALSE, m_L=2, initfreq_L=c(freq_min, 500),
                        true_data=true_data, limFreq=c(1000),
                        actPlot=TRUE);
      
      result[i+(j-1)*N,1]=dist
      result[i+(j-1)*N,2]=out$covpbb[1,1]
      result[i+(j-1)*N,3]=out$covpbb[1,2]
      result[i+(j-1)*N,4]=out$residual[1,1]
      result[i+(j-1)*N,5]=out$residual[1,2]
      result[i+(j-1)*N,6]=out$residual[1,3]
    }
    
    ind1=1+(j-1)*N
    ind2=N+(j-1)*N
    print(sprintf("signal %s @ distance: %f kpc. Covpbb mean:%f. Covpbb median: %f",
                  signal_name, dist, mean(result[ind1:ind2,2]), median(result[ind1:ind2,2])))
  
  }
  filename=sprintf("results/multiDec/EinsteinTelescope/results_AA_%s_f2_%s_multiDec.txt", 
                   filtering_method, signal_name)
  #write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)
}