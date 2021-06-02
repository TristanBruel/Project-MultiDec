library ("stats")
library ("signal")
library ("seewave")
library ("psd")
library ("pracma")
library("plyr")

source("multiDec_algebra.R")


##################################
### Multi detection simulation ###
##################################


########################################################################
signal_multiDec = function(dec=30, ra=6, t=1302220800, 
                           signal="KURODA_TM1_H_resampled.dat",
                           detectors=c("LHO","LLO","VIR"),
                           verbose=FALSE,actPlot=FALSE){
  ######################################################################
  # Inputs:  sky position of the source
  #               declination in °, right ascension in hours
  #          time GPS at which the wave arrives at the center of the Earth
  #          name of the simulated waveform
  #          detectors in which the signal will be measured
  #
  # Outputs: time series for the given detectors
  #          wvf_LHO   time: time vector
  #                    hoft: amplitude measured
  #          wvf_LLO
  #          ...
  ######################################################################
  
  folder="Waveforms/"
  
  gw_filename = paste(folder,signal,sep="");
  sXX = read.table(gw_filename) # V1 time, V2 hplus, V3 hcross
  colnames(sXX) = c ("time","hplus","hcross");
  fs_orig = round(1/(sXX$time[2]-sXX$time[1]));
  n = length(sXX$time);
  
  padding = floor(0.05*fs_orig+1);   # zero padding with 50ms at the start and end

  duration = (n+2*padding-1)/fs_orig;
  
  if (verbose==TRUE){
    print(gw_filename)
    print(sprintf("Duration : %gs", duration))
  }
  
  nDet=length(detectors);
  res=list();
  # antenna responses and time delay for each detector
  F=antenna_patterns(dec,ra,t,0,detectors);
  delays=time_delays(dec,ra,t,detectors);
  # Reset reference position to first detector
  delays=delays-delays[1];
  delayLengths=round(delays*fs);
  
  ### Assign storage ###
  time = rep(0,n+2*padding);
  
  ind1 = 1+padding
  indn = n+padding
  time[ind1:indn] = sXX$time
  t1 = sXX$time[1]   # time at which the GW signal starts
  tn = sXX$time[n]   # time at which the GW signal ends
  time[1:padding] = seq(t1-padding/fs_orig,t1-1/fs_orig,by=1/fs_orig);
  time[(indn+1):(indn+padding)] = seq(tn+1/fs_orig,tn+padding/fs_orig,by=1/fs_orig);

  for (k in 1:nDet){
    signal= rep(0,n+2*padding);
    Fplus = F[k,1];
    Fcross = F[k,2];
    startGW=ind1+delayLengths[k];
    endGW=indn+delayLengths[k];
    signal[startGW:endGW]= Fplus*sXX$hplus + Fcross*sXX$hcross;
    
    wvf=data.frame("time"=time,"hoft"=signal);
    
    # Plot
    if (actPlot == TRUE){
      plot(1000*wvf$time,wvf$hoft,type='l',xlab="Time [ms]",ylab="Hoft",
           main=detectors[k],ylim=c(-3e-22,3e-22))
    }
    
    res$wvf=wvf
    res=rename(res,c("wvf"=sprintf("wvf_%s",detectors[k])))
  }
  res$duration = duration
  
  # Time delays between the arrivals at each detector
  if (verbose==TRUE){
    for (k in 2:nDet){
      print(sprintf("Time shift between %s and %s is %s ms",detectors[k],
                    detectors[1],1000*delays[k]))
    }
  }
  return(res)
}


########################################################################
data_multiDec = function (fs=4096, wvfs, duration, ampl=1, detectors=c("LHO","LLO","VIR"), 
                          filter="prewhiten", setseed=0,
                          actPlot=FALSE, verbose=FALSE){
  ######################################################################
  # Inputs:   fs: sampling frequency
  #           duration: duration (in second) of the output time serie
  #           wvfs: list of dataframes (time=t, hoft=h(t)) that contain 
  #                   the signal waveforms sampled at fs
  #           ampl: multiplication factor to simulate the source distance
  #           detectors: vector of detectors from which data will be extracted
  #           filter: name of the method
  #              "HP" : The fcut parameter is fixed internally (15 Hz)
  #              "spectrum" : the data are whiten in Fourier domain using 
  #                             the noise spectrum estimate
  #              "AR" : AR model
  #              "prewhiten": use the R prewhiten function
  #           setseed: if 0, random seed. Otherwise set to the value
  # 
  # Outputs:  data_H  d$t: time vector
  #                   d$x: noise+signal
  #                   d$y: filtered (noise+signal)
  #           data_L
  #           ...
  ######################################################################
  n=length(wvfs[[1]]$time)
  
  m=length(detectors)
  res=list()
  for (k in 1:m){
    wvf=wvfs[[k]]
    wvf_size=length(wvf$hoft)
    
    if (n<wvf_size){
      print(sprintf("data_generator: signal waveform duration larger than %f", duration))
      return()
    }
    
    # The output vector will be 2 times larger than n
    factor=2
    
    # Add noise
    data=noise_generator(factor,fs, duration, detectors[k], setseed=setseed, 
                         filter=FALSE, actPlot=FALSE, verbose=FALSE)
    
    Y=data$x
    psd=data$psd        # 2 sided PSD
    n_data=length(Y)     # factor x n
    
    if (verbose==TRUE){
      print(sprintf("data_generator:size of the output: %d", n))
      print(sprintf("data_generator:size of the noise : %d", n_data))
      print(sprintf("data_generator:size of the %s signal: %d", detectors[k], wvf_size))
      print(sprintf("data_generator:amplitude of the signal: %f", ampl))
    }
    
    # Signal addition (centered at the middle of the data vector 
    # to avoid filtering leakage at the beginning and end).
    ind1=floor((n_data-wvf_size)/2)
    
    for (i in 1:wvf_size){
      Y[ind1+i]=Y[ind1+i]+ampl*wvf$hoft[i]
    }
    
    # filter the time series if requested
    if (filter != FALSE){
      YY=filtering(Y, fs, filter, psd, verbose)
    }else{
      YY=Y
    }
    
    # generate a time series
    T = seq(1, n_data, by = 1)
    
    # select the original data size
    Tf = seq(1, n, by = 1)
    
    for (i in 1:n){
      Tf[i]=wvf$time[i]
    }
    
    T_wvf=seq(1,wvf_size,by=1)
    
    for (i in 1:wvf_size){
      T_wvf[i]=wvf$time[i]
    }
    
    Yf = seq(1, n, by = 1)
    YYf = seq(1, n, by = 1)
    
    for (i in 1:n){
      Yf[i]=Y[ind1+i]
      YYf[i]=YY[ind1+i]
    }
    
    if (actPlot){
      if (filter == "HP" || filter == "spectrum" || filter == "prewhiten" || filter == "AR"){
        # Plot for LHO
        plot(T, Y, col="black", type="l", pch=1, panel.first = grid(), 
             xlab="Sample",ylab="Hoft",main=detectors[k])
        points(T, YY, col="red", type="l", pch=2);        # (noise + signal) filtered 
        leg = c("noise+signal", "(noise+signal) filtered")
        col = c("black","red")
        legend (x=T[1]*1.1,y=max(Y)*.9,legend=leg,cex=.8,col=col,pch=c(1,2))
        
        plot(Tf, Yf, col="black", type="l", pch=1, panel.first = grid(), 
             xlab="Time [s]",ylab="Hoft", main=detectors[k])
        points(Tf, YYf, col="red", type="l", pch=2);          # (noise + signal) filtered
        points(T_wvf,(wvf$hoft)*ampl,col="green",type="l",pch=3);  # signal only
        
        leg = c("noise", "(noise+signal) filtered", "signal only")
        col = c("black","red","green")
        legend (x=Tf[1]*1.1,y=max(Yf)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))
      
      }else{
        plot (Tf, Yf, type="l", col="black", main=detectors[k])
        legend (x=Tf[1]*1.1, y=max(Yf)*.9, legend="noise+signal")
      }
    }
  res$data=data.frame(t=Tf,x=Yf,y=YYf)
  res=rename(res,c("data"=sprintf("data_%s",detectors[k])))
  }
  
  return(res)
  
}


########################################################################
noise_generator = function (factor,fs, duration, detector, setseed=0,
                            filter=FALSE, actPlot=FALSE, verbose=FALSE){
  ######################################################################
  # Inputs :  fs: sampling frequency
  #           duration: duration (in second) of the output time series
  #           detector: name of the detector whose ASD will be used to generate colored noise
  #           setseed: if 0, random seed. Otherwise set to the value.
  #           filter method:
  #               "HP" : The fcut parameter is fixed internally (15 Hz)
  #               "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #               "AR" : AR model
  #               "prewhiten": use the R prewhiten function
  #
  # Output:   d$t: time vector
  #           d$x: noise
  #           d$y: filtered noise
  #           d$psd: psd
  ######################################################################
  
  if (duration < 10){
    # For 3G detectors we need to use a frequency resolution smaller than 0.1 Hz
    n=factor*(duration*fs+1)
  }
  else{
    n=duration*fs+1
  }
  
  if (verbose==TRUE){
    print(sprintf("noise_generator:size of the noise output vector:%d", n))
  }
  
  # Noise generation
  freq2=fs*fftfreq(n)          # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:floor(n/2)]        # one-sided frequency vector
  
  # Get the 2 sided PSD
  if (detector == "ALIGO"){
    psd=aLIGO_PSD_new(freq2, 2)
  }
  else{
    psd=PSD_fromfiles(freq2, 2, detector, actPlot)
  }
  
  if (setseed >0){
    set.seed(setseed)
  }
  
  X = rnorm(n, mean=0, sd=1);           # Gaussian white noise
  XX = fft(X);                          # FFT computing
  XXX = XX*sqrt(psd)*sqrt(fs);          # Coloring
  Y = fft(XXX, inverse = TRUE);         # FFT inverse
  Y = Re(Y)/n;                          # noise in time domain
  # Note on the normalisation factor: 
  #  - n comes from the FFT and FFT inverse (sqrt(n) each)
  #  - to color properly the noise and keep the amplitude right 
  #    one needs to multiply by sqrt(psd) x sqrt(fs) 
  
  # filter the time series if requested
  if (filter != FALSE){
    YY=filtering(Y, fs, filter, psd, verbose)
  }else{
    YY=Y
  }
  
  if (verbose){
    ss <- std(Y)
    print(sprintf("noise_generator:noise time serie sigma:%g", ss))
    ss <- std(YY)
    print(sprintf("noise_generator:filtered noise time serie sigma:%g", ss))
  }
  
  # generate a time series vector sampled at fs
  Tf = seq(1, n, by = 1)
  
  for (i in 1:n){
    Tf[i]=(i-1)/fs
  }
  
  if (actPlot){
    # Time series
    T = seq(1, n, by = 1)
    ts.plot(Y); # noise only
    points(T, Y, col="black", type="l", pch=1, panel.first = grid())
    points(T, YY, col="red", type="l",pch=2)
    
    legend_str=c("simulated noise", "filtered noise")
    legend (x=0, y=abs(max(Y)), legend=legend_str, col=c("black","red"), pch=c(1,2))   
    
    # spectrum estimated
    psdest <- pspectrum(Y, Y.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, plot=FALSE,verbose=FALSE)
    psdest_filtered <- pspectrum(YY, YY.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, plot=FALSE,verbose=FALSE)
    
    # Fourier transform
    YFT = sqrt(2)*fft(Y)/sqrt(n);
    WFT = sqrt(2)*fft(YY)/sqrt(n);
    ymin=10^(ceiling(log10(min(abs(YFT)[1:floor(n/2)])/sqrt(fs))))
    ymax=10^(ceiling(log10(max(abs(YFT)[1:floor(n/2)])/sqrt(fs))))
    #ymin=1e-24
    #ymax=2e-21
    plot (freq1, abs(YFT)[1:floor(n/2)]/sqrt(fs), log="xy", type="l", xlab="Frequency", ylab="ASD", 
          col="grey", xlim=c(1, fs/2), ylim=c(ymin,ymax), pch=1, panel.first = grid())
    
    lines(fs*psdest$freq, sqrt(psdest$spec)/sqrt(fs), col="blue", pch=2)
    
    lines(freq1, abs(WFT)[1:floor(n/2)]/sqrt(fs), col="black", pch=4)        # factor 2 because FT is 2 sided
    lines(fs*psdest_filtered$freq[1:floor(n/2)],                      # pspectrum is 1 sided
          sqrt(psdest_filtered$spec[1:floor(n/2)])/sqrt(fs), col="green", pch=5)
    
    lines(freq1, sqrt(2*psd[1:floor(n/2)]), col="red", pch=3)         # PSD is 2 sided PSD
    
    legend_str=c("col noise FT", "col noise spectrun", "ASD model", "filtered FT", "filtered spectrum")
    legend ("topright", legend=legend_str, col=c("grey","blue","red","black","green"), pch=c(1,2,3,4,5))   
    
    if (verbose==TRUE){
      s1 <- sqrt(2*trapz(fs*psdest$freq[1:floo(rn/2)], psdest$spec[1:floor(n/2)]/fs))
      print(sprintf("noise_generator:colored noise rms:%g", s1))
      
      s2 <- sqrt(2*trapz(fs*psdest_filtered$freq[1:floor(n/2)], psdest_filtered$spec[1:floor(n/2)]/fs))
      print(sprintf("noise_generator:filtered noise rms:%g", s2))
      
      Sn_min=sqrt(2*min(psd))
      print(sprintf("minimal asd value:%g",Sn_min))
    }
    
  }
  return(list(t=Tf,x=Y,y=YY,psd=psd))
}


########################################################################
PSD_fromfiles=function(f, type, detector, actPlot=FALSE){
  ########################################################################
  # Sensitivity curves for advanced LIGO, advanced Virgo and KAGRA.
  # [Add refs here]
  # f: frequency vector
  # type=1 --> one-sided PSD.
  # type=2 --> two-sided PSD.
  # detector: name of the detector
  
  cutoff=1e-42            # For 2nd generator detectors
  
  if (detector=="aLIGO"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", toupper(detector))
    data=read.table(psd_filename);
    sens=data$V6}   # Design
  
  if ((detector=="KAGRA") || (detector=="KAG")){
    psd_filename=sprintf("PSD/ALIGO_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V6}   # Design
  
  if ((detector=="LHO") || (detector=="LLO") || (detector=="VIR")){
    psd_filename=sprintf("PSD/ALIGO_sensitivity.txt")
    data=read.table(psd_filename);
    sens=data$V6}   # Design
  
  if (exists("sens")==FALSE){
    stop(sprintf("Detector %s is not implemented in this code. You may want to use CE1, CE2, ET_B, ET_C, ET_D, aLIGO, ADV, KAGRA or ALIGO",detector))
  }
  
  n = length(f)
  fmin=f[1]
  if (type==1){
    fmax=f[n]
  } else{
    fmax=abs(f[n/2+1])}
  
  yl=sens[1]
  yr=sens[length(data$V1)]
  
  asd_func = approxfun(x = data$V1, y = sens, method = "linear",
                       yleft=yl, yright=yr, rule = 1, f = 0, ties = "mean");
  
  if (type==1){
    asd = asd_func(f)
    psd = asd*asd
  }else{
    asd = rep(0, n);
    asd_1sided = asd_func(abs(f[1:floor(n/2)]))
    
    if (length(asd_1sided) != floor(n/2)){
      print (sprintf("Warning: ASD vector size %d is different from the frequency vector size %d", 
                     length(asd_1sided), n/2))
    }
    
    
    asd[1]=asd_1sided[1];
    for(i in 2:(floor(n/2))){
      asd[i]=asd_1sided[i];
      
      # Wraparound frequency: f=0 is the first element (i=1), 
      # and all elements are symetric around index n/2+1
      asd[n+2-i]=asd[i]
    }  
    asd[n/2+1]=asd_func(abs(f[floor(n/2)+1]))
    
    # Two sided psd
    asd=asd/sqrt(2);
    psd=asd*asd;
  }
  
  for (i in 1:n){
    if (psd[i]>cutoff){
      psd[i]=cutoff
    }
  }
  
  if (actPlot==TRUE){
    fN=4096
    if (type==1){
      plot(f,psd,log="y",col="blue",xlim=c(1, fN/2),pch=2)
      points(data$V1,sens*sens,col="red",type="l",pch=1)
    }else{
      plot(f,psd,log="y",col="blue",xlim=c(1, fN/2),pch=2)
      points(data$V1,0.5*sens*sens,col="red",type="l",pch=1)
    }
    leg = c(detector,"interpolated")
    col = c("red","blue")
    legend (x=500,y=psd[1]*0.8,legend=leg,cex=.8,col=col,pch=c(1,2))
  }
  
  return(psd)
}


########################################################################
filtering = function(X, fs, method, psd=0, verbose=FALSE){
  ########################################################################
  # data processing of the input vector according to different methods
  # X: input data
  # fs: sampling frequency of X
  # method: filtering method
  #   "HP" : The fcut parameter is fixed internally (15 Hz)
  #   "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #   "AR" : AR model
  #   "prewhiten": use the R prewhiten function
  # psd: PSD required by the AR filering method
  
  # warning: the psd must be the 2 sided PSD. The size of the psd and data vectors must be equal
  if (length(X) != length(psd)){
    print(length(X))
    print(length(psd))
    warning("filtering:filtering::the data and psd vectors must have the same size. Please check")
  }
  
  n=length(X)
  duration=(n-1)/fs
  
  # compute noise sigma 
  freq2=fs*fftfreq(n)          # two-sided frequency vector
  
  s0 <- sqrt(trapz(freq2, psd))
  if (verbose==TRUE){
    print(sprintf("filtering: ASD noise rms: %g", s0))
  }
  
  if (method == "HP"){
    fcut=10
    # filtfilt : zero phase filter (forward& backward)
    myfilter=butter(n=5, W=fcut/(fs/2), type="high")
    Y=filtfilt(filt=myfilter, x=X)} 
  else if (method == "AR"){
    if (psd==0){
      print("Filtering with AR method cannot be performed because noise psd has not been provided")
    }else{
      # generate another noise TS
      X1 = rnorm(n, mean=0, sd=1);            # Gaussian white noise
      XX1 = fft(X1);                          # FFT computing
      XXX1 = XX1*sqrt(psd);                   # Coloring
      Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
      Y1 = Re(Y1)*sqrt(fs)/n;                 # noise in time domain
      
      ar_model <- stats::ar(Y1,order.max=10, aic=FALSE ,method=c("yule-walker"), demean=TRUE);
      b <- stats::filter(x=X, filt=c(1, -ar_model$ar[1], -ar_model$ar[2], -ar_model$ar[3], 
                                     -ar_model$ar[4], -ar_model$ar[5], -ar_model$ar[6], 
                                     -ar_model$ar[7], -ar_model$ar[8], -ar_model$ar[9], 
                                     -ar_model$ar[10]), method="convolution", sides = 1);
      b[1]=b[2]=b[3]=b[4]=b[5]=b[6]=b[7]=b[8]=b[9]=b[10]=b[11]
      Y=b}
  }
  else if (method == "spectrum"){
    if (psd[1]==0){
      print("Filtering with specrum method cannot be performed because noise psd has not been provided")
    }else{
      # generate another noise TS
      X1 = rnorm(n, mean=0, sd=1);          # Gaussian white noise
      XX1 = fft(X1);                        # FFT computing
      XXX1 = XX1*sqrt(psd);                 # Coloring
      Y1 = fft(XXX1, inverse = TRUE);       # FFT inverse
      Y1 = Re(Y1)*sqrt(fs)/n;               # noise in time domain
      
      # compute the PSD
      psdest <- pspectrum(Y1, Y1.frqsamp=fs, ntap.init=6, Nyquist.normalize = TRUE, 
                          plot=FALSE,verbose = FALSE)
      psdwhitening=rep(0, n);
      for(i in 1:(floor(n/2))){
        psdwhitening[i]=psdest$spec[i]
        psdwhitening[n+1-i]=psdest$spec[i]
      }
      a = fft(X)                        # FFT computing and normalization
      b = a/sqrt(psdwhitening)       # whitening
      c = fft(b, inverse = TRUE);       # FFT inverse
      Y = s0*Re(c)/n;                   # Normalization factor of the 2 FFTs
      
      #    myfilter=butter(n=4,W=10/(fs/2),type="high")
      #    YY=filtfilt(filt=myfilter,x=YY)
    }
  }
  else if (method == "prewhiten"){
    # prewhiten
    myts <- ts(X, start=0, end=duration, frequency=fs)
    myts <- prewhiten(myts, AR.max=100, zero.pad="rear", plot=FALSE, verbose=FALSE)
    Y <- myts[['prew_ar']][1:n]}
  else{
    print("No filtering method specify")
    Y=X
  }
  
  return (Y)
}


########################################################################
fftfreq = function(n, d = 1){
  ########################################################################
  # surrogate for the numpy fft.fftfreq function that generates the two sided 
  # frequency vector. Defaults d=1 means sampling frequency is 1. 
  # https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
  #
  # n: samples number
  # d: sample spacing (inverse of the sampling rate). Defaults to 1
  
  if(n%%2 == 0){# n is even
    out = c(seq(0, n/2-1, by = 1), seq(-n/2, -1, by=1)) / (d*n);
  }else{ # n is odd
    out = c(seq(0, (n-1)/2, by = 1), seq(-(n-1)/2, -1, by=1)) / (d*n);
  }
  
  return(out);
}

########################################################################
int = function(n){
  ########################################################################
  # https://stackoverflow.com/questions/31036098/what-is-the-difference-between-int-and-floor-in-python-3
  if(n < 0 ){
    return(-floor(abs(n)))
  }else{
    return(floor(n))
  }
}


########################################################################
compute_SNR = function(wvf, detector="aLIGO", asd=NULL, fcut=0, dist=10, 
                       pbOff=FALSE, actPlot=FALSE){
  ########################################################################
  # Compute the Signal-To-Noise ratio for a given wvf (x$time, x$hoft) 
  # and a given detector
  
  fs=round(1/(wvf$time[2]-wvf$time[1]))
  n=length(wvf$hoft)
  a = nextpow2(2*n)         #zero padding and rounding to the next power of 2
  n2=2^a
  
  # Remove or not 0.100s after the bounce (set hoft values to 0)
  if (pbOff == TRUE){
    ext=which(wvf$time<0.1)
    wvf$hoft[ext]=0
  }
  
  
  freq2 = fs*fftfreq(n2)         # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:floor(n2/2)]       # one-sided frequency vector
  
  # Get the 1 sided PSD
  if (detector == "ALIGO"){
    psd=aLIGO_PSD_new(freq1, 1)
  }else{
    psd=PSD_fromfiles(freq1, 1, detector)
  }
  
  if (!is.null(asd)){
    psd = rep(asd*asd, length(freq1))
  }
  
  vec=rep(0,n2)
  for (i in 1:n){
    vec[n2/4+i]=vec[n2/4+i]+wvf$hoft[i]*10./dist
  }  
  
  hf=fft(vec);          # normalization wrt the sampling
  
  hf=hf[1:(n2/2)]                # The integral is performed over positive freqs
  
  hf=subset(hf,freq1-fcut>0)
  psd=subset(psd,freq1-fcut>0)
  freq1=subset(freq1, freq1-fcut>0)
  
  snr=sqrt(4/fs/n2*sum(abs(hf)^2/psd))
  
  if (actPlot){
    plot (freq1, sqrt(freq1)*abs(hf), log="xy", type="l", xlab="Frequency", ylab="hchar", 
          col="grey", xlim=c(1, fs/2), ylim=c(1e-24,1e-20), pch=1, panel.first = grid())
    points(freq1,sqrt(psd), type="l", col="black",pch=2)
    leg = c("sqrt(fs) x h~(f)", "ASD")
    col = c("grey","black")
    legend (x=1,y=6e-22,legend=leg,cex=.8,col=col,pch=c(1,2))
    title(c("SNR:",snr))
  }
  return(snr)  
}


########################################################################
whiteNoise = function(wvf, SNR, fs=4096, fcut=0, 
                      actPlot = FALSE, main=NULL, verbose=FALSE){
  ########################################################################
  # Create a white gaussian noise that match the desired Signal-To-Noise ratio
  # for a given waveform (x$time, x$hoft)
  # Outputs: noisy: signal+noise
  #          psd: constant psd
  
  # Compute the constant (one-sided) psd required to get the desired SNR
  n=length(wvf$hoft)
  a=nextpow2(2*n)         #zero padding and rounding to the next power of 2
  n2=2^a
  freq2 = fs*fftfreq(n2)         # two-sided frequency vector
  freq1=freq2[1:floor(n2/2)]
  
  vec=rep(0,n2)
  for (i in 1:n){
    vec[n2/4+i]=vec[n2/4+i]+wvf$hoft[i]
  }  
  
  hf=fft(vec);     # normalization wrt the sampling
  
  hf=hf[1:(n2/2)]    # integral performed over positive frequencies
  
  hf=subset(hf,freq1-fcut>0);
  freq1=subset(freq1, freq1-fcut>0);
  
  psd=4/fs/n2*sum(abs(hf)^2)/SNR^2
  
  # Create Noise
  X = rnorm(n, mean=0, sd=1);
  XX = fft(X);
  XXX = XX*sqrt(psd)*sqrt(fs);
  Y = fft(XXX, inverse = TRUE);
  Y = Re(Y)/n;
  
  for (i in 1:n){
    Y[i] = Y[i] + wvf$hoft[i]
  }
  
  if (verbose){
    print(sprintf("Gaussian white noise on detector %s with mean %g and standard deviation %g",
                  main, mean(Y), sd(Y)));
  }
  
  if (actPlot){
    plot(wvf$time,Y,type='l',col='black',
         xlab = "Time",ylab="Hoft",main=main)
    points(wvf$time,wvf$hoft,type='l',col='red')
    leg = c("Noisy signal", "Signal")
    col = c("black","red")
    legend ("topleft",legend=leg,cex=.8,col=col,pch=c(1,2))
  }
  
  return(list(noisy=Y,psd=psd))
}


########################################################################
compute_cor = function(h1,h2,psd=1,fs=4096){
  ########################################################################
  #
  # Compute the average correlation coefficient between two signals
  
  if (length(h1) != length(h2)){
    print("Signals must be of same length")
  }
  T=length(h1);
  freq=fs*fftfreq(2*T)
  freq[1]=0.001
  freq=freq[1:T]   # frequencies used to sample the fourier series
  
  N=100   # number of repetitions to get the average correlation value
  r=rep(0,N);
  for (n in 1:N){
    integrand=(Conj(h1)*h2+h1*Conj(h2))/psd
    inner_p=Re(trapz(freq,integrand))
    
    norm1=Re(trapz(freq,2*abs(h1)^2/psd))
    norm2=Re(trapz(freq,2*abs(h2)^2/psd))
    
    r[n] = inner_p/sqrt(norm1*norm2);
  }
  return(mean(r))
}