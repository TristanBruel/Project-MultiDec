source("data_multiDec.R")
source("timeFreqMaps.R")
source("specPdgrm.R")
source("functions.R")

### Source ###
dec=30
ra=6
t=1302220800
skyPosition=c(dec,ra,t)

### Signal ###
signal_name="KURODA_TM1_H_resampled.dat"
fs=4096

gw_filename=paste("Waveforms/",signal_name,sep="");
sXX = read.table(gw_filename); # V1 time, V2 hplus, V3 hcross
colnames(sXX) = c ("time","hplus","hcross");

detectors = c("LHO", "LLO","VIR","KAG");
nDet = length(detectors);

signals = signal_multiDec(dec, ra, t,signal=signal_name, detectors = detectors,
                          actPlot = FALSE);
N = fs*signals$duration+1   # number of samples

specPlus=specPdgrm(sXX$hplus,l=100,p=90,fs=fs,zoomFreq=c(0,3/4),main='True hplus spectrogram');
gmodeplus=findGmodes(r=specPlus,l=100,p=90,actPlot=T,gmode="left",
                     initfreq_L = c(0,400));
specCross=specPdgrm(sXX$hcross,l=100,p=90,fs=fs,zoomFreq=c(0,3/4),main='True hcross spectrogram');
gmodecross=findGmodes(r=specCross,l=100,p=90,actPlot=T,gmode="left",
                      initfreq_L = c(0,400))
n=length(specPlus$x);
timedata = seq(specPlus$x[1], specPlus$x[n],length=n);
plot(timedata,gmodeplus$maxf_L,col='black',xlab="Time [ms]",ylab="Freq [Hz]",
     main='gmode tracking',panel.first=grid())
points(timedata,gmodecross$maxf_L,col='grey')

distances=c(1,3,5,7,10);
colors=c('red','blue','green','orange','purple');
temp=1

for (dist in distances){
  data = data_multiDec(fs,signals,signals$duration,detectors=detectors,
                       setseed=3,ampl=10/dist,verbose=FALSE,actPlot=FALSE);
  # samples to ignore at the start and end (50ms)
  transient = floor(0.05*fs+1);
  
  wData = matrix(0,nrow=N,ncol=nDet);
  psd = matrix(0,nrow=51,ncol=nDet);
  freq=fs*seq(0,1/2,by=1/100);
  
  for (k in 1:nDet){
    psd[,k]=PSD_fromfiles(freq,1,detectors[k])
    wData[,k] = data[[k]]$y;   # prewhiten data
    wData[,k] = wData[,k]/sqrt(mean(psd[,k]));
  }
  
  likelihoods=timeFreqMap(fs,wData,detectors,psd,skyPosition,100,10,transient,
                          freqBand=c(0,1500),windowType='modifiedHann',actPlot=T);
  r=likelihoods$soft;
  n=length(r$x);
  gmode=findGmodes(r,l=100,p=90,actPlot=T,gmode="left",
                   initfreq_L = c(0,400),um_L=4,dm_L=4)
  timedata = seq(r$x[1],r$x[n],length=n);
  #points(timedata,gmode$maxf_L,col=colors[temp])
  temp=temp+1;
}
#legend ("topleft",legend=c("hplus","hcross",distances),
 #       col=c('black','grey',colors),pch=1)