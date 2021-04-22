library ("stats")
library ("signal")
library ("seewave")
library ("psd")
library ("pracma")

source("multiDec_algebra.R")

##################################
### Multi detection simulation ###
##################################


########################################################################
signal_multiDec = function(dec=-65, ra=8, t=1293494400, 
                           signal="KURODA_TM1_H_resampled.dat", 
                           verbose=TRUE,actPlot=TRUE){
  ######################################################################
  # Inputs :  sky position of the source
  #               declination in ° and right ascension in hours
  #           time GPS at which the wave arrives at the center of the Earth
  #           name of the (simulated) waveform
  #
  # Outputs : measured time series for the 3 detectors LHO, LLO and VIRGO
  
  folder="Waveforms/"
  fs_orig=4096
  
  gw_filename=paste(folder,signal,sep="")
  sXX = read.table(gw_filename) # V1 time, V2 hplus, V3 hcross
  colnames(sXX) = c ("time","hplus","hcross")
  n = length(sXX$time)
  duration=(n-1)/fs_orig
  
  if (verbose==TRUE){
    print(gw_filename)
    print(sprintf("Number of samples at 24414 Hz: %g", n))
    print(sprintf("Duration : %gs", duration))
  }
  
  # gravitational response for each detector
  coeff_LHO = grav_response(dec,ra,t,0,"LHO")
  coeff_LLO = grav_response(dec,ra,t,0,"LLO")
  coeff_VIR = grav_response(dec,ra,t,0,"VIRGO")
  
  signal_LHO = rep(0,n)
  signal_LLO = rep(0,n)
  signal_VIR = rep(0,n)
  
  for (i in 1:n){
    signal_LHO[i]=coeff_LHO$Fplus*sXX$hplus[i] + coeff_LHO$Fcross*sXX$hcross[i]
    signal_LLO[i]=coeff_LLO$Fplus*sXX$hplus[i] + coeff_LLO$Fcross*sXX$hcross[i]
    signal_VIR[i]=coeff_VIR$Fplus*sXX$hplus[i] + coeff_VIR$Fcross*sXX$hcross[i]
  }
  
  # Time delays between the arrivals at each detector
  if (verbose==TRUE){
    detectors = c("LHO","LLO","VIRGO")
    delays = c(time_delay(dec, ra, t, "LHO"), time_delay(dec, ra, t, "LLO"), 
               time_delay(dec, ra, t, "VIRGO"))
    first = which(delays==max(delays))
    third = which(delays==min(delays))
    second = which(!(c(1,2,3)%in%c(first,third)))
    ## WARNING ## : special cases when 2 detectors measure a signal at the exact same time
    dt_second = delays[first]-delays[second]
    dt_third = delays[first]-delays[third]
  
    print(sprintf("Signal detected first at %s",detectors[first]))
    print(sprintf("Then at %s with a %fs delay",detectors[second],dt_second))
    print(sprintf("And finally at %s with a %fs delay",detectors[third],dt_third))
  }
  
  times_H=sXX$time-time_delay(dec,ra,t,"LHO")
  times_L=sXX$time-time_delay(dec,ra,t,"LLO")
  times_V=sXX$time-time_delay(dec,ra,t,"VIRGO")
  
  wvf_LHO=data.frame("time"=times_H,"hoft"=signal_LHO)
  wvf_LLO=data.frame("time"=times_L,"hoft"=signal_LLO)
  wvf_VIR=data.frame("time"=times_V,"hoft"=signal_VIR)
  
  
  # Plot
  if (actPlot == TRUE){
    plot(wvf_LHO$time,wvf_LHO$hoft,type='l',xlab="Time [s]",ylab="Hoft",main="LHO")
    plot(wvf_LLO$time,wvf_LLO$hoft,type='l',xlab="Time [s]",ylab="Hoft",main="LLO") 
    plot(wvf_VIR$time,wvf_VIR$hoft,type='l',xlab="Time [s]",ylab="Hoft",main="VIRGO") 
  }
  
  return(list(wvf_LHO=wvf_LHO,wvf_LLO=wvf_LLO,wvf_VIR=wvf_VIR,duration=duration))
}


########################################################################
data_multiDec = function (fs=4096, duration, wvf_LHO, wvf_LLO, wvf_VIR, 
                          ampl=0, detector="aLIGO", filter="prewhiten", setseed=0,
                           actPlot=TRUE, verbose=TRUE){
  ########################################################################
  # Inputs:   fs: sampling frequency
  #           duration: duration (in second) of the output time serie
  #           wvf: dataframe (time=time, hoft=h(t)) that contains the signal waveform sampled at fs
  #           ampl: multiplication factor to change the source distance
  #           detector: detector name
  #           filter: name of the method
  #              "HP" : The fcut parameter is fixed internally (15 Hz)
  #              "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #              "AR" : AR model
  #              "prewhiten": use the R prewhiten function
  #           setseed: if 0, random seed. Otherwise set to the value
  #           The wvf is centered in the noise vector
  # 
  # Outputs:  data_H  d$t: time vector
  #                   d$x: noise+signal
  #                   d$y: filtered (noise+signal)
  #           data_L
  #           data_V
  
  n=duration*fs+1
  wvf_size_H=length(wvf_LHO$hoft)
  wvf_size_L=length(wvf_LLO$hoft)
  wvf_size_V=length(wvf_VIR$hoft)
  
  if ((n<wvf_size_H) || (n<wvf_size_L) || (n<wvf_size_V)){
    print(sprintf("data_generator:the signal waveform duration is larger than %f", duration))
    return()
  }
  
  # The output vector will be 2 times larger than n
  factor=2
  
  # Noise LHO
  data=noise_generator(factor,fs, duration, detector, setseed=setseed, filter=FALSE,
                       actPlot=FALSE, verbose=FALSE)
  
  Y_LHO=data$x
  psd_LHO=data$psd           # 2 sided PSD
  n_data=length(Y_LHO)     # factor x n

  # Noise LLO
  data=noise_generator(factor,fs, duration, detector, setseed=setseed, filter=FALSE,
                       actPlot=FALSE, verbose=FALSE)
  Y_LLO=data$x
  psd_LLO=data$psd
  
  # Noise VIRGO
  data=noise_generator(factor,fs, duration, detector, setseed=setseed, filter=FALSE,
                       actPlot=FALSE, verbose=FALSE)
  Y_VIR=data$x
  psd_VIR=data$psd
  
  if (verbose==TRUE){
    print(sprintf("data_generator:size of the output: %d", n))
    print(sprintf("data_generator:size of the noise : %d", n_data))
    print(sprintf("data_generator:size of the LHO signal: %d", wvf_size_H))
    print(sprintf("data_generator:size of the LLO signal: %d", wvf_size_L))
    print(sprintf("data_generator:size of the VIRGO signal: %d", wvf_size_V))
    print(sprintf("data_generator:amplitude of the signal: %f", ampl))
  }
  
  # Signal addition (centered at the middle of the data vector 
  # to avoid filtering leakage at the beginning and end).
  ind1_H=floor((n_data-wvf_size_H)/2)
  ind1_L=floor((n_data-wvf_size_L)/2)
  ind1_V=floor((n_data-wvf_size_V)/2)
  
  for (i in 1:wvf_size_H){
    Y_LHO[ind1_H+i]=Y_LHO[ind1_H+i]+ampl*wvf_LHO$hoft[i]
  }
  for (i in 1:wvf_size_L){
    Y_LLO[ind1_L+i]=Y_LLO[ind1_L+i]+ampl*wvf_LLO$hoft[i]
  }
  for (i in 1:wvf_size_V){
    Y_VIR[ind1_V+i]=Y_VIR[ind1_V+i]+ampl*wvf_VIR$hoft[i]
  }
  
  # filter the time series if requested
  if (filter != FALSE){
    YY_H=filtering(Y_LHO, fs, filter, psd_LHO, verbose)
    YY_L=filtering(Y_LLO, fs, filter, psd_LLO, verbose)
    YY_V=filtering(Y_VIR, fs, filter, psd_VIR, verbose)
  }else{
    YY_H=Y_LHO
    YY_L=Y_LLO
    YY_V=Y_VIR
  }
  
  # generate a time series
  T_H = seq(1, n_data, by = 1)
  T_L = seq(1, n_data, by = 1)
  T_V = seq(1, n_data, by = 1)
  
  # select the original data size
  Tf_H = seq(1, n, by = 1)
  Tf_L = seq(1, n, by = 1)
  Tf_V = seq(1, n, by = 1)
  
  for (i in 1:n){
    Tf_H[i]=wvf_LHO$time[i]
    Tf_L[i]=wvf_LLO$time[i]
    Tf_V[i]=wvf_VIR$time[i]
  }
  
  T_wvf_H=seq(1,wvf_size_H,by=1)
  T_wvf_L=seq(1,wvf_size_L,by=1)
  T_wvf_V=seq(1,wvf_size_V,by=1)
  for (i in 1:wvf_size_H){
    T_wvf_H[i]=wvf_LHO$time[i]
  }
  for (i in 1:wvf_size_L){
    T_wvf_L[i]=wvf_LLO$time[i]
  }
  for (i in 1:wvf_size_V){
    T_wvf_V[i]=wvf_VIR$time[i]
  }
  
  Yf_H = seq(1, n, by = 1)
  Yf_L = seq(1, n, by = 1)
  Yf_V = seq(1, n, by = 1)
  YYf_H = seq(1, n, by = 1)
  YYf_L = seq(1, n, by = 1)
  YYf_V = seq(1, n, by = 1)
  
  for (i in 1:n){
    Yf_H[i]=Y_LHO[ind1_H+i]
    Yf_L[i]=Y_LLO[ind1_L+i]
    Yf_V[i]=Y_VIR[ind1_V+i]
    YYf_H[i]=YY_H[ind1_H+i]
    YYf_L[i]=YY_L[ind1_L+i]
    YYf_V[i]=YY_V[ind1_V+i]
  }
  
  if (actPlot==TRUE){
    if (filter == "HP" || filter == "spectrum" || filter == "prewhiten" || filter == "AR"){
      # Plot for LHO
      plot(T_H, Y_LHO, col="black", type="l", pch=1, panel.first = grid(), 
           xlab="Sample",ylab="Hoft",main="LHO")
      points(T_H, YY_H, col="red", type="l", pch=2);        # (noise + signal) filtered 
      leg = c("noise+signal", "(noise+signal) filtered")
      col = c("black","red")
      legend (x=T_H[1]*1.1,y=max(Y_LHO)*.9,legend=leg,cex=.8,col=col,pch=c(1,2))
      
      plot(Tf_H, Yf_H, col="black", type="l", pch=1, panel.first = grid(), 
           xlab="Time [s]",ylab="Hoft", main="LHO")
      points(Tf_H, YYf_H, col="red", type="l", pch=2);          # (noise + signal) filtered
      points(T_wvf_H,(wvf_LHO$hoft)*ampl,col="green",type="l",pch=3);  # signal only
      
      leg = c("noise", "(noise+signal) filtered", "signal only")
      col = c("black","red","green")
      legend (x=Tf_H[1]*1.1,y=max(Yf_H)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))
      
      # Plot for LLO
      plot(T_L, Y_LLO, col="black", type="l", pch=1, panel.first = grid(), 
           xlab="Sample",ylab="Hoft",main="LLO")
      points(T_L, YY_L, col="red", type="l", pch=2);        # (noise + signal) filtered 
      leg = c("noise+signal", "(noise+signal) filtered")
      col = c("black","red")
      legend (x=T_L[1]*1.1,y=max(Y_LLO)*.9,legend=leg,cex=.8,col=col,pch=c(1,2))
      
      plot(Tf_L, Yf_L, col="black", type="l", pch=1, panel.first = grid(), 
           xlab="Time [s]",ylab="Hoft", main="LLO")
      points(Tf_L, YYf_L, col="red", type="l", pch=2);          # (noise + signal) filtered
      points(T_wvf_L,(wvf_LLO$hoft)*ampl,col="green",type="l",pch=3);  # signal only
      
      leg = c("noise", "(noise+signal) filtered", "signal only")
      col = c("black","red","green")
      legend (x=Tf_L[1]*1.1,y=max(Yf_L)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))
      
      # Plot for VIRGO
      plot(T_V, Y_VIR, col="black", type="l", pch=1, panel.first = grid(), 
           xlab="Sample",ylab="Hoft",main="VIRGO")
      points(T_V, YY_V, col="red", type="l", pch=2);        # (noise + signal) filtered 
      leg = c("noise+signal", "(noise+signal) filtered")
      col = c("black","red")
      legend (x=T_V[1]*1.1,y=max(Y_VIR)*.9,legend=leg,cex=.8,col=col,pch=c(1,2))
      
      plot(Tf_V, Yf_V, col="black", type="l", pch=1, panel.first = grid(), 
           xlab="Time [s]",ylab="Hoft", main="VIRGO")
      points(Tf_V, YYf_V, col="red", type="l", pch=2);          # (noise + signal) filtered
      points(T_wvf_V,(wvf_VIR$hoft)*ampl,col="green",type="l",pch=3);  # signal only
      
      leg = c("noise", "(noise+signal) filtered", "signal only")
      col = c("black","red","green")
      legend (x=Tf_V[1]*1.1,y=max(Yf_V)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))
      
    }else{
      plot (Tf_H, Yf_H, type="l", col="black", main="LHO")
      legend (x=Tf_H[1]*1.1, y=max(Yf_H)*.9, legend="noise+signal")
      
      plot (Tf_L, Yf_L, type="l", col="black", main="LLO")
      legend (x=Tf_L[1]*1.1, y=max(Yf_L)*.9, legend="noise+signal")
      
      plot (Tf_V, Yf_V, type="l", col="black", main="VIRGO")
      legend (x=Tf_V[1]*1.1, y=max(Yf_V)*.9, legend="noise+signal")
    }
  }
  
  data_H = data.frame(t=Tf_H,x=Yf_H,y=YYf_H)
  data_L = data.frame(t=Tf_L,x=Yf_L,y=YYf_L)
  data_V = data.frame(t=Tf_V,x=Yf_V,y=YYf_V)
  return(list(data_H=data_H,data_L=data_L,data_V=data_V))
  
}


########################################################################
noise_generator = function (factor,fs, duration, detector, setseed=0,
                            filter=FALSE, actPlot=TRUE, verbose=TRUE){
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
  ########################################################################  
  
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
  freq1=freq2[1:int(n/2)]        # one-sided frequency vector
  
  # Get the 2 sided PSD
  if (detector == "ALIGO"){
    psd=aLIGO_PSD_new(freq2, 2)
  }
  else{
    psd=PSD_fromfiles(freq2, 2, detector, verbose)
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
  
  if (verbose==TRUE){
    ss <- std(Y)
    print(sprintf("noise_generator:noise time serie sigma:%g", ss))
    ss <- std(YY)
    print(sprintf("noise_generator:filtered noise time serie sigma:%g", ss))
  }
  
  # generate a time series vector sampled at fs
  Tf = seq(1, n, by = 1)
  
  for (i in 1:n){
    Tf[i]=i/fs
  }
  
  if (actPlot==TRUE){
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
    ymin=10^(ceiling(log10(min(abs(YFT)[1:int(n/2)])/sqrt(fs))))
    ymax=10^(ceiling(log10(max(abs(YFT)[1:int(n/2)])/sqrt(fs))))
    #ymin=1e-24
    #ymax=2e-21
    plot (freq1, abs(YFT)[1:int(n/2)]/sqrt(fs), log="xy", type="l", xlab="Frequency", ylab="ASD", 
          col="grey", xlim=c(1, fs/2), ylim=c(ymin,ymax), pch=1, panel.first = grid())
    
    lines(fs*psdest$freq, sqrt(psdest$spec)/sqrt(fs), col="blue", pch=2)
    
    lines(freq1, abs(WFT)[1:int(n/2)]/sqrt(fs), col="black", pch=4)        # factor 2 because FT is 2 sided
    lines(fs*psdest_filtered$freq[1:int(n/2)],                      # pspectrum is 1 sided
          sqrt(psdest_filtered$spec[1:int(n/2)])/sqrt(fs), col="green", pch=5)
    
    lines(freq1, sqrt(2*psd[1:int(n/2)]), col="red", pch=3)         # PSD is 2 sided PSD
    
    legend_str=c("col noise FT", "col noise spectrun", "ASD model", "filtered FT", "filtered spectrum")
    legend (x=100, y=min(abs(tail(YFT,-1)))*50000, legend=legend_str, col=c("grey","blue","red","black","green"), pch=c(1,2,3,4,5))   
    
    if (verbose==TRUE){
      s1 <- sqrt(2*trapz(fs*psdest$freq[1:int(n/2)], psdest$spec[1:int(n/2)]/fs))
      print(sprintf("noise_generator:colored noise rms:%g", s1))
      
      s2 <- sqrt(2*trapz(fs*psdest_filtered$freq[1:int(n/2)], psdest_filtered$spec[1:int(n/2)]/fs))
      print(sprintf("noise_generator:filtered noise rms:%g", s2))
      
      Sn_min=sqrt(2*min(psd))
      print(sprintf("minimal asd value:%g",Sn_min))
    }
    
  }
  return(list(t=Tf,x=Y,y=YY,psd=psd))
}

########################################################################
PSD_fromfiles=function(f, type, detector, verbose=FALSE){
  ########################################################################
  # Sensitivity curves for advanced LIGO, advanced Virgo and KAGRA.
  # [Add refs here]
  # f: frequency vector
  # type=1 --> one-sided PSD.
  # type=2 --> two-sided PSD.
  # detector: name of the detector
  
  cutoff=1e-42            # For 2nd generator detectors
  
  # depending on the detector, the ASD is located in different columns of the data vector
  if (detector=="CE1"){
    psd_filename="PSD/curves_Jan_2020/ce1.txt"
    data=read.table(psd_filename);
    sens=data$V2        
    cutoff=1e-44}   
  
  if (detector=="CE2"){
    psd_filename="PSD/curves_Jan_2020/ce2.txt"
    data=read.table(psd_filename);
    sens=data$V2        
    cutoff=1e-44}   
  
  if (detector=="ET_B"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V2
    cutoff=1e-44}
  
  if (detector=="ET_C"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V2
    cutoff=1e-44}
  
  if (detector=="ET_D"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V4   # HF + LF
    cutoff=1e-44}
  
  if (detector=="ADV"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V7}   # Design
  
  if (detector=="aLIGO"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", toupper(detector))
    data=read.table(psd_filename);
    sens=data$V6}   # Design
  
  if (detector=="aLIGO2"){
    psd_filename="PSD/aLIGODesign.txt"
    data=read.table(psd_filename);
    sens=data$V2}   # Design  
  
  
  if (detector=="KAGRA"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V6}
  
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
    asd_1sided = asd_func(abs(f[1:int(n/2)]))
    
    if (length(asd_1sided) != int(n/2)){
      print (sprintf("Warning: ASD vector size %d is different from the frequency vector size %d", 
                     length(asd_1sided), n/2))
    }
    
    for(i in 1:(int(n/2))){
      asd[i]=asd_1sided[i];
      
      # Wraparound frequency: f=0 is the first element (i=1), 
      # and all elements are symetric around index n/2+1
      if(i>1){
        asd[n+2-i]=asd[i]
      }
    }  
    asd[n/2+1]=asd_func(abs(f[int(n/2)+1]))
    
    # Two sided psd
    asd=asd/sqrt(2);
    psd=asd*asd;
  }
  
  for (i in 1:n){
    if (psd[i]>cutoff){
      psd[i]=cutoff
    }
  }
  
  if (verbose==TRUE){
    fN=4096
    if (type==1){
      plot(f,psd,log="xy",col="blue",xlim=c(1, fN/2),pch=2)
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
  duration=n/fs
  
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
      Y1 = Re(Y1)*sqrt(fs)/n;                   # noise in time domain
      
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
      X1 = rnorm(n, mean=0, sd=1);           # Gaussian white noise
      XX1 = fft(X1);                          # FFT computing
      XXX1 = XX1*sqrt(psd);                   # Coloring
      Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
      Y1 = Re(Y1);                   # noise in time domain
      
      # compute the PSD
      #myts <- ts(Y1, start=0, end=duration, frequency=fs)
      psdest <- pspectrum(Y1, Y1.frqsamp=fs, ntap.init=6, Nyquist.normalize = TRUE, 
                          plot=FALSE,verbose = FALSE)
      psdwhitening=rep(0, n);
      for(i in 1:(int(n/2))){
        psdwhitening[i]=psdest$spec[i]
        psdwhitening[n+1-i]=psdest$spec[i]
      }
      
      a = fft(X)                        # FFT computing and normalization
      b = a/sqrt(psdwhitening)          # whitening
      c = fft(b, inverse = TRUE);       # FFT inverse
      Y = s0*Re(c);                     # Normalisation factor of the 2 FFTs
      
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