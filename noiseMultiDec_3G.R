source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")

library(signal)

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

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = TRUE);


#############################
### Simulation parameters ###
#############################
# Sky position: direction towards Andromeda Galaxy
dec=41.27;
ra=0.71;
skyPosition=c(dec,ra);
# Time of arrival at the center of Earth
t0=1325113218; #favourable case
t0=1325062818; #unfavourable case

detectors=c("ET1","ET2","ET3","CEH","CEL");

nDet=length(detectors);

signal_name=c("s20.0--LS220");
dist=1e6;

fs=4096;
#filtering_method="prewhiten";
filtering_method="spectrum";

# loop over N generation of noisy data and add signal
N=100;

result<-zeros(N,6);

wvfs = signal_multiDec(dec=dec,ra=ra,t=t0,fs=fs,signal=signal_name,detectors=detectors,
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

# To use always the same random noise realizations
set.seed(1);

#limits = rep(0,N)
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
  r2 = list();
  r2$t = subset(r$t, r$t<=1.5);
  r2$f = r$f;
  r2$E = r$E[1:length(r2$t),];
  
  out = covpbb_LASSO(r=r2, mod=fit, true_data=true_data, limFreq=c(1000),
                     actPlot=FALSE);
  
  result[i,1]=dist;
  result[i,2]=out$covpbb[1,1];
  result[i,3]=out$covpbb[1,2];
  result[i,4]=out$residual[1,1];
  result[i,5]=out$residual[1,2];
  result[i,6]=out$residual[1,3];
}

print(sprintf("signal %s @ distance: %f kpc. Covpbb mean:%f. Covpbb median: %f",
              signal_name, dist, mean(result[1:N,2]), median(result[1:N,2])));

save_dir="./perf/3G/rmsd_favourable/HLVKA/";
dir.create(path=save_dir, showWarnings=FALSE, recursive=TRUE);
filename=sprintf("results_AA_%s_f2_noise.txt", filtering_method);
save_path=paste(save_dir, filename, sep='');
write.table(result, file=save_path, sep=" ", row.names=FALSE, col.names=FALSE)