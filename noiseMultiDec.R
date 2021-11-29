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

# Linear model
mod = lm(r~ 0 + f + I(f^3), data = fits_data);# final model
summary(mod);

# Variable variance linear model
mu1 = "~ f - 1";
mu2 = "~ f + I(f^2) - 1";
mu3 = "~ f + I(f^2) + I(f^3) - 1";
mu4 = "~ f + I(f^3) - 1";
s1 = "~ f - 1"; 
s2 = "~ f + I(f^2) - 1";

Xm  = model.matrix(eval(parse(text=eval(parse(text=mu3)))), fits_data);
Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = TRUE);

############################
### Algorithm parameters ###
############################
# g-mode estimate: starting from the left
gmode = c("left")

#############################
### Simulation parameters ###
#############################
# Sky position
# Nearly best case scenario for Hanford and Livingston
dec=60
ra=8
skyPosition=c(dec,ra)

t0=1302220800

detectors=c("LHO","LLO","VIR", "KAG", "LAO")

nDet=length(detectors)

signal_name=c("s25.0--LS220")
dist=1e6

fs=4096
#filtering_method="prewhiten"
filtering_method="spectrum"

# loop over N generation of noisy data and add signal
N=100

result<-zeros(N,6)

wvfs = signal_multiDec(dec=dec,ra=ra,t=t0,fs=fs,signal=signal_name,detectors=detectors,
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

# To use always the same random noise realizations
set.seed(1)

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
  
  out = covpbb_poly(r=r2, mod=fit, true_data=true_data, limFreq=c(1000),
                    actPlot=FALSE);
  
  result[i,1]=dist
  result[i,2]=out$covpbb[1,1]
  result[i,3]=out$covpbb[1,2]
  result[i,4]=out$residual[1,1]
  result[i,5]=out$residual[1,2]
  result[i,6]=out$residual[1,3]
}

print(sprintf("signal %s @ distance: %f kpc. Covpbb mean:%f. Covpbb median: %f",
                signal_name, dist, mean(result[1:N,2]), median(result[1:N,2])))

filename=sprintf("perf/multiDec/HLVKA/results_AA_%s_f2_noise.txt", 
                 filtering_method)
write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)