source("multiDec_algebra.R")
library(fields)
library(signal)

###########################
### Time Frequency Maps ###
###########################


########################################################################
timeFreqMap = function(fs=4096,wData,detectors=c("LHO","LLo","VIR"),psd,
                       skyPosition,integLength,offsetLength,transientLength,
                       freqBand=c(0,Inf),windowType="hann",
                       verbose=TRUE, actPlot=TRUE,logPow=TRUE){
  ######################################################################
  # Inputs :  fs: sampling frequency
  #           wData: matrix of time-domain whitened data
  #                   one column per detector
  #           detectors: vector of detectors used in the analysis
  #           psd: power spectral densities, one column per detector
  #                 Must be sampled at frequencies fs*[0:1/2, by=1/integLength]
  #                 one column per detector
  #           skyPosition: vector (dec, ra, t)
  #           integLength: size of the data used for each Fourier transform
  #                 Must be a power of 2
  #           offsetLength: number of samples between consecutive FT segments
  #                 Must be a power of 2
  #           transientLength: number of samples to ignore (beginning and end)
  #           freqBand: band of frequencies to include in the maps
  #           windowType: window type to be used in FFTs ('hann' or 'bartlett')
  #
  # Output : time frequency map
  ######################################################################
  
  ### Command line arguments ###
  nDet=length(detectors);
  if (length(wData) < nDet){
    stop("Number of detectors is inconsistent with the data provided")
  }
  nDat=length(wData[,1]);
  if (verbose){
    print(c("Number of samples : ",nDat))
  }
  
  # Vector of one-sided frequencies
  freq1=fs/2*seq(0,1,length=integLength/2+1);
  freqInd=which((freq1>=freqBand[1]) & (freq1<=freqBand[2]));
  inbandFreq=freq1[freqInd];
  nFreqBins=length(inbandFreq);
  psd=psd[freqInd,];
  if (verbose){
    print(c("Number of positive frequency bins : ",nFreqBins))
  }
  
  ### Partition (no overlapping yet) ###
  startInd=seq(transientLength+1,nDat-transientLength-integLength+1,by=integLength);
  stopInd=startInd+integLength-1;
  segmentIndices=cbind(startInd,stopInd);
  
  # Number of segments
  nSeg=length(startInd);
  
  # Fractionnal offset of consecutive segments
  offset=offsetLength/integLength;
  
  # Number of time bins (max covers the case in which there is only one segment)
  nTimeBins=max(nSeg-1,1)/offset;
  if (verbose){
    print(c("Number of time bins",nTimeBins))
    }
  
  ### Detectors information ###
  dec=skyPosition[1];
  ra=skyPosition[2];
  t=skyPosition[3];
  F=antenna_patterns(dec,ra,t,0,detectors);
  if (verbose){
    print("Antenna response matrix : F")
    print(F)
  }
  delays=time_delays(dec,ra,t,detectors);
  
  # Reset reference position to first detector
  delays=delays-delays[1];
  delayLengths=delays*fs;
  intDelayLengths=round(delayLengths);
  residualDelays=(delayLengths-intDelayLengths)/fs;
  
  # Create real time indices
  timeIndices=zeros(nTimeBins,nDet);   # one column per detector
  for (k in 1:nDet){
    timeIndices[,k]=(seq(0,nTimeBins-1)*offsetLength/fs+delays[k])*1000;
  }
  
  ### Storage ###
  TFMapFull=array(dim=c(integLength,nTimeBins,nDet));  # contains all frequencies
  TFMap=array(dim=c(nFreqBins,nTimeBins,nDet));   # contains positive frequencies
  detSpec=array(dim=c(nFreqBins,nTimeBins,nDet));
  
  ### TF maps at zero delay ###
  segStartInd=segmentIndices[1,1];
  segStopInd=segmentIndices[max(nSeg-1,1),2];
  
  # Construct Window
  if (windowType=="hann"){
    window=matrix(hanning(integLength),ncol=1);
  }
  else if (windowType=="bartlett"){
    window=matrix(bartlett(integLength),ncol=1);
  }
  else{
    window=ones(integLength,1);
  }
  window=window/sqrt(mean(window^2));
  windowData=repmat(window,1,max(nSeg-1,1));
  
  # FT the data in loop over detectors
  for (k in 1:nDet){
    if (offset==1){   # means no overlapping
      data=wData[segStartInd:segStopInd,k];
      dataArray=matrix(data,nrow=integLength,ncol=nTimeBins);
      TFMapFull[,,k]=mvfft(windowData*dataArray);
    }
    else{
      for (j in 1:(1/offset)){
        offsetStartInd=segStartInd+integLength*(j-1)*offset;
        offsetStopInd=segStopInd+integLength*(j-1)*offset;
        data=wData[offsetStartInd:offsetStopInd,k];
        dataArray=matrix(data,nrow=integLength);
        timeBins=seq(j,nTimeBins,by=1/offset);
        TFMapFull[,timeBins,k]=mvfft(windowData*dataArray);
      }
    }
    # Extract in-band frequencies
    TFMap[,,k]=TFMapFull[freqInd,,k];
    
    # Save the squared magnitudes
    detSpec[,,k]=Re(TFMap[,,k])^2+Im(TFMap[,,k])^2;
    if(logPow == TRUE){
      detSpec[,,k]=log10(sqrt(detSpec[,,k]));
    }
    if (actPlot){
      image.plot(timeIndices[,k],inbandFreq,t(detSpec[,,k]),
                 xlab="Time [ms]",ylab="Frequency [Hz]",
                 main=c(detectors[k],"Spectrogram"))
    }
  }
  
  ### Compute detector responses ###
  # Noise-weighted responses (one column per detector)
  wFp=repmat(F[,1],nFreqBins,1)/sqrt(psd);
  wFc=repmat(F[,2],nFreqBins,1)/sqrt(psd);
  
  # Convert to Dominant Polarization Frame
  wFpDP=zeros(nFreqBins,nDet);
  wFcDP=zeros(nFreqBins,nDet);
  for (l in 1:nFreqBins){
    DPF=convertToDPF(wFp[l,],wFc[l,]);
    wFpDP[l,]=DPF[,1];
    wFcDP[l,]=DPF[,2];
  }
  
  ### Residual delays ###
  residualDelayPhases=
    exp(1i*2*pi*t(matrix(inbandFreq,nrow=1))%*%matrix(residualDelays,nrow=1));

  ### Construct TF maps for each detector ###
  for (k in 2:nDet){
    ### FFT data ###
    segStartInd=segmentIndices[1,1]+intDelayLengths[k];
    segStopInd=segmentIndices[max(nSeg-1,1),2]+intDelayLengths[k];
    
    if (offset==1){   # means no overlapping
      data=wData[segStartInd:segStopInd,k];
      dataArray=matrix(data,nrow=integLength,ncol=nTimeBins);
      TFMapFull[,,k]=mvfft(windowData*dataArray);
    }
    else{
      for (j in 1:(1/offset)){
        offsetStartInd=segStartInd+integLength*(j-1)*offset;
        offsetStopInd=segStopInd+integLength*(j-1)*offset;
        data=wData[offsetStartInd:offsetStopInd,k];
        dataArray=matrix(data,nrow=integLength);
        timeBins=seq(j,nTimeBins,by=1/offset);
        TFMapFull[,timeBins,k]=mvfft(windowData*dataArray);
      }
    }
    
    # Retain only positive frequencies
    TFMap[,,k]=TFMapFull[freqInd,,k];
    
    # Apply residual time shifts
    TFMap[,,k]=TFMap[,,k]*repmat(matrix(residualDelayPhases[,k]),1,nTimeBins);
  }

  # DPF matrix components
  Mpp=sum(wFpDP*wFpDP)
  Mcc=sum(wFcDP*wFcDP)
  
  # Antenna-response weighted TF maps
  wFpTFMap=zeros(nFreqBins,nTimeBins);
  wFcTFMap=zeros(nFreqBins,nTimeBins);
  
  for (k in 1:nDet){
    wFpTFMap=wFpTFMap+TFMap[,,k]*wFpDP[,k]
    wFcTFMap=wFcTFMap+TFMap[,,k]*wFcDP[,k]
  }
  
  ### Compute likelihood maps ###
  # Plus energy
  plusLikelihood=(1/Mpp)*(Re(wFpTFMap)^2+Im(wFpTFMap)^2);
  
  # Cross energy
  crossLikelihood=(1/Mcc)*(Re(wFcTFMap)^2+Im(wFcTFMap)^2);
  
  # Standard likelihood
  stdLikelihood=plusLikelihood+crossLikelihood;
  
  # log-scale transformation
  if(logPow){
    plusLikelihood=log10(sqrt(plusLikelihood));
    crossLikelihood=log10(sqrt(crossLikelihood));
    stdLikelihood=log10(sqrt(stdLikelihood));
  }
  
  ### Plot Maps ###
  if (actPlot){
    image.plot(timeIndices[,1],inbandFreq,t(plusLikelihood),
               xlab = "Time [ms]", ylab = "Frequency [Hz]", main = "Plus energy")
    image.plot(timeIndices[,1],inbandFreq,t(crossLikelihood),
               xlab = "Time [ms]", ylab = "Frequency [Hz]", main = "Cross energy")
    image.plot(timeIndices[,1],inbandFreq,t(stdLikelihood),
               xlab = "Time [ms]", ylab = "Frequency [Hz]", main = "Standard likelihood")
  }
  return(list(stdLike=stdLikelihood,Ep=plusLikelihood,Ec=crossLikelihood))
}