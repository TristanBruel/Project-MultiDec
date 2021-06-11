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
dec=60
ra=8
t=1302220800
skyPosition=c(dec,ra,t)

detectors=c("LHO","LLO","VIR")
nDet=length(detectors)
#signals=c("s11.2--LS220", "s15.0--GShen", "s15.0--LS220", "s15.0--SFHo", 
#          "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220")
signal_name=c("s20.0--LS220")
freq_min=200
distance_step=0.5

fs=4096
filtering_method="prewhiten"

# loop over N generation of noisy data and add signal

N=100
N=1
dist_nb=40
dist_nb=1
result<-array(0,c(nrow=N*dist_nb,ncol=6))

wvfs=signal_multiDec(dec,ra,t,fs,signal_name,pbOff=TRUE,actPlot=FALSE)
true_data=wvfs$true_data
startTime = 0.1   # signal starts 100ms after bounce (t=0)
L = fs*wvfs$duration+1   # number of samples

l = 200   # interval length to use for each FFT
p = 90   # overlapping percentage
offset = (1-p/100)*l   # offset between consecutive FTs
transient = floor(0.05*fs+1);   # samples to ignore at the start and end (50ms)

wData = matrix(0,nrow = L, ncol = nDet);
psd = matrix(0,nrow = l/2+1, ncol = nDet);
freq=fs*seq(0,1/2,by=1/l);

# To use always the same random noise for all waveforms % detectors
set.seed(1)

for (j in 1:dist_nb){
  dist = 1+(j-1)*distance_step
  dist=10
  for (i in 1:N){
    d=data_multiDec(fs,wvfs,ampl=10/dist,detectors=c("LHO","LLO","VIR"), 
                    filter=filtering_method, setseed=0,
                    actPlot=FALSE, verbose=FALSE);
    
    for (k in 1:nDet){
      psd[,k]=PSD_fromfiles(freq,1,detectors[k])
      wData[,k] = d[[k]]$y;   # prewhiten data
      wData[,k] = wData[,k]/sqrt(mean(psd[,k]));
    }
    
    likelihoods=timeFreqMap(fs,wData,detectors,psd,skyPosition,l,offset,transient,
                            freqBand=c(0,2000),windowType='modifiedHann',
                            startTime=startTime,logPow=TRUE,
                            actPlot=TRUE,verbose=FALSE)
    r=likelihoods$std
    
    out = covpbb2(r, mod=fit, l=l, p=p, fs=fs,
                 um_L = 2, dm_L = 0, m_L = 8, initfreq_L = c(200, 500),
                 um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(1000, 1700),
                 gmode = gmode, true_data=true_data, actPlot=FALSE,
                 limFreq = c(1000));
    
    result[i+(j-1)*N,1]=dist
    result[i+(j-1)*N,2]=out$covpbb[1,1]
    result[i+(j-1)*N,3]=out$covpbb[1,2]
    result[i+(j-1)*N,4]=out$residual[1,1]
    result[i+(j-1)*N,5]=out$residual[1,2]
    result[i+(j-1)*N,6]=out$residual[1,3]
  }
  for (k in 1:nDet){
    ind1=1+(j-1)*N
    ind2=N+(j-1)*N
    print(sprintf("signal %s @ distance: %f kpc. Covpbb mean:%f. Covpbb median: %f",
                  signal_name, dist, mean(result[ind1:ind2,2]), median(result[ind1:ind2,2])))
  }
}
for (k in 1:nDet){
  #filename=sprintf("results/multiDec/results_AA_%s_f2_%s_%s.txt", 
  #                 filtering_method, signal_name, "multiDec")
  #write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)
}