source("functions.R")
#library(psplinePsd)
library(fields)

###################
### Periodogram ###
###################

########################################################################
pdgrm = function(data, ylog = TRUE, actPlot = TRUE, main = NULL, ...){
######################################################################## 
  
  # This function generates the periodogram
  
  ##############
  ### Inputs ###
  ##############
  
  # data: data set
  # ylog: it produces log-periodogram
  # actPlot: if it is TRUE, it produces the plot
  # main: title for the plot
  
  ### ###
  
  if(class(data) == "data.frame"){
    
    n = dim(data)[1];
    
    X = data$V2;
    
  }else if(class(data) == "numeric"){
    
    n = length(data);
    
    X = data;
    
  }else{
    
    stop("Check format of the data");
    
  }
  
  N    = ifelse( n%%2 == 0, n/2+1, (n + 1)/2);# N=(n+1)/2 (ODD) or N=n/2+1 (EVEN)
  
  freq = seq(0, pi, length = N);
  
  # Frequencies to remove from estimate
  if (n %% 2) {  # Odd length time series
    bFreq = 1;
    
  }else {  # Even length time series
    bFreq = c(1, N);
  }
  
  # Optimal rescaling to prevent numerical issues
  rescale = stats::sd(X);
  auxdata = X / rescale;  # Data now has standard deviation 1
  auxdata = auxdata - mean(auxdata);
  
  N       = length(auxdata)/2 + 1;
  pdgrm   = (abs(stats::fft(auxdata))^2 / (2*pi*length(auxdata)))[1:N];
  pdgrm   = pdgrm * rescale^2;
  
  #pdgrm = (abs(stats::fft(X))^2 / (2 * pi * n))[1:N];
  
  pdgrm = pdgrm[-bFreq];
  
  if(ylog == TRUE){ # log periodogram
    
    pdgrm = log(pdgrm);
    
  }
  
  if(actPlot == TRUE){
    
    if(is.null(main)){
      
      main = c("log10-periodogram");
      
    }else if(is.null(main)){
      
      main = c("Periodogram");
      
    }
    
    graphics::plot.default(freq[-bFreq], pdgrm, type = "l", col = "grey",
                           xlab = "Frequency", ylab = "log PSD", main = paste(main), ...);
    
  }
  
  return(pdgrm);
  
}

###################
### Spectrogram ###
###################

########################################################################
specPdgrm = function(data, time = NULL, l, p, fs = NULL, actPlot = TRUE, 
                     logPow = TRUE, zoomFreq = c(0,1), main  = NULL,
                     method = "ar", y = NULL){
  ########################################################################
  
  # This function produces x, y & z values of the spectrogram
  # If actPlot is TRUE, it also produces the spectrogram
  
  ##############
  ### Inputs ###
  ##############
  
  # data = numeric vector with data (no zeros [specPdgrm])
  # time = numeric vector - time values (optional)
  # l    = interval length
  # p    = overlapping percentage
  # fs   = sampling frequency
  # actPlot  = if it is TRUE, it produces the spectrogram
  # logPow   = if it is TRUE, it works with log-transformation
  # zoomFreq = 
  # main     = optional title for the Spectrogram
  # method   = "ar" uses autoregressive model to estimate PSD (see "spec.ar"
  #            function for more information); "pdgrm" uses periodogram to 
  #            estimate PSD
  # y        = y-axis values
  
  ### ###
  
  if((zoomFreq[1]<0) || (zoomFreq[2]>1) || (zoomFreq[1] > zoomFreq[2])){
    stop("zoomFreq must be a vector c(a,b) with values 0 <= a < b <= 1");
  }
  
  index = ints(n = length(data), l = l, p = p);
  
  R = NULL;
  
  for(i in 1:dim(index)[1]){
    
    auxdata = data[index[i,1]:index[i,2]];
    
    # Optimal rescaling to prevent numerical issues
    rescale = stats::sd(auxdata);
    auxdata = auxdata / rescale;  # Data now has standard deviation 1
    auxdata = auxdata - mean(auxdata);
    
    if(method == "pdgrm"){
      
      aux = pdgrm(auxdata, ylog = FALSE, actPlot = FALSE);
      
    }else if(method == "ar"){
      
      aux = stats::spectrum(auxdata, method = "ar", plot = FALSE)$spec;
      
      #aux = spectrum(auxdata, method = "ar", plot = FALSE);
      #print(aux$method) # Print order of the AR model
      #aux = aux$spec
      
    }else{
      
      stop("Define proper method: 'ar' or 'pdgrm'");
      
    }
    
    R = cbind(R, aux);
    
  }
  
  if(logPow == TRUE){R = log10(sqrt(R));}
  
  index    = round(dim(R)[1] * zoomFreq);
  index    = seq(index[1] + 1, index[2]);
  
  pc = 1-p/100;
  
  i <- 1:dim(R)[2];
  #i <- 1:(length(R)-1); # Rput form psplinePsd
  
  x <- round(l / fs * (i-1) * pc * 1000, 2); # number of intervals
  
  if(!is.null(time)){
    
    if(length(time)==2){
      
      x <- seq(from= time[1], to=time[2], length = length(x));
      
    }else{
      
      x = seq(from= time[1], to=time[length(time)], length = length(x));
      
    }
  }
  
  if(is.null(y)){
    
    y <- seq(1,fs/2,length=dim(R)[1]);
    
  }else{
    
    if(length(y) != dim(R)[1])stop(paste("Length of y must be ", dim(R)[1]));  
    
  }
  z <- t(R);
  #z <- z[, -c(1, dim(z)[2])]; #deleting first and last row in image
  #y <- y[-c(1, dim(z)[2])]; # adjustment
  
  y = y[index]; # zoom
  z = z[, index];
  rownames(z) = NULL;
  
  if(actPlot == TRUE){
    
    if( is.null(main) ){
      
      if( method == "ar"){
        
        main = c("Spectrogram based on Autoregressive model");
        
      }else if(method == "pdgrm"){
        
        main = c("Spectrogram based on periodogram");
        
      }
    }
    
    image.plot(x, y, z, xlab = "Time [s]", ylab = "Frequency [Hz]", main = paste(main)); # yaxt='n'
  }
  
  return(list(x = x, y = y, z = z));
  
}
