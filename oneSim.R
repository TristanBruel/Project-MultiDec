source("data_multiDec.R")
source("functions.R")
source("timeFreqMaps.R")

library(signal)

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
# Sky position
dec=60
ra=8
skyPosition=c(dec,ra)
# Time of arrival at the center of the Earth
t0=1302220800
# Distance of the source
dist=10

detectors=c("LHO","LLO","VIR","KAG","LAO")
nDet=length(detectors)

signal_name="s25.0--LS220"

fs=4096
filtering_method="spectrum"
  
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


set.seed(1)

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
                          logPow=TRUE,actPlot=TRUE,verbose=FALSE);

r = likelihoods$std

### Reduced time window ###
r2=list()
r2$t = subset(r$t, r$t<=1.5);
r2$f = r$f;
r2$E = r$E[1:length(r2$t),]

out = covpbb_poly(r2, mod=fit, true_data=true_data, limFreq=c(1000),
                  actPlot=TRUE);
  
print(sprintf("signal %s @ distance: %f kpc. Covpbb: %f",
                signal_name, dist, out$covpbb[1,1]))