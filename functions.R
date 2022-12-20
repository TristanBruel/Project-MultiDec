library ("seewave");
library ("signal")
library(lmvar)
library(glmnet)


########################################################################
findGmodes_LASSO = function(r, lambda = 1,
                            mask_t = NULL, mask_f = c(0.,100.), 
                            actPlot = FALSE, saveToTxt = FALSE){
  ########################################################################
  
  # r : output from specPdgrm
  # lambda : penalization parameter
  # mask_t : time band in which the pixels of maximum intensity 
  #           are not taken into account for the tracking
  # mask_f : frequency band in which the pixels of maximum intensity 
  #           are not taken into account for the tracking
  # actPlot   : logical value to produce plot
  
  x = r$t; y = r$f; E = r$E; # values from spectrogram
  
  # frequency of maximum intensity at each time index
  maxf = y[apply(E, 1, which.max)]; 
  
  maxE = rep(0,length(x))
  for (i in 1:length(x)){
    maxE[i] = E[i, which.max(E[i,])]
  }
    
  # Remove the maxima in the time band mask_t
  if (!is.null(mask_t)){
    t_masked = x<mask_t[1] | x>mask_t[2]
    x = x[t_masked]
    maxf = maxf[t_masked]
    maxE = maxE[t_masked]
  }
  
  # Remove the maxima in the frequency band mask_f
  if (!is.null(mask_f)){
    f_masked = maxf<mask_f[1] | maxf>mask_f[2]
    x = x[f_masked]
    maxf = maxf[f_masked]
    maxE = maxE[f_masked]
  }

  #x = x[maxf>0];
  #maxf = maxf[maxf>0];
  
  if (length(maxf)==0){ # in case all maxima are equal to zero
    return(rep(0,length(r$t)))
  }

  # special case in which all the maxima are the same
  if (length(count(maxf)$freq)==1){
    return(rep(maxf[1],length(r$t)))
  }
  # standard case
  Xmat <- cbind(x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10);
  weights = 10^maxE
  lasso_model = glmnet(Xmat, maxf, alpha=1., lambda=lambda, weights=weights);
  #lasso_model = glmnet(Xmat, maxf, alpha=1., lambda=lambda);
  
  maxf_fit = predict(lasso_model, Xmat);
  
  if(actPlot){
    image.plot(r$t,y,E,xlab="Time [s]",ylab="Frequency [Hz]",
               cex.lab=1.8, cex.axis=1.5)
    points(x, maxf, col='blue',pch=19)
    points(x, maxf_fit, col='black',pch=19)
    leg <- c("Maxima", "First fit")#,"Uncertainty");
    col <- c("blue","black")#,"gray");
    legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
  }
  
  ## Second iteration
  # Compute the rsmd and remove the maxf values further away from the fit
  rmsd = sqrt(sum((maxf_fit-maxf)^2)/length(maxf));
  mask = abs(maxf-maxf_fit)<rmsd;
  #mask = abs(maxf-maxf_fit)<400;
  maxf2 = maxf[mask];
  t2 = x[mask];
  maxE2 = maxE[mask]
  
  # special case in which the fit is too far from the maxima
  if (length(maxf2)==0){
    return(rep(0,length(r$t)))
  }
  
  # special case in which all the maxima are the same
  if (length(count(maxf2)$freq)==1){ 
    return(rep(maxf2[1],length(r$t)))
  }
  
  # Fit again
  Xmat2 <- cbind(t2, t2^2, t2^3, t2^4, t2^5, t2^6, 
                 t2^7, t2^8, t2^9, t2^10);
  weights2 = 10^maxE2
  lasso_model = glmnet(Xmat2, maxf2, alpha=1., lambda=lambda, weights=weights2);
  #lasso_model = glmnet(Xmat2, maxf2, alpha=1., lambda=lambda);
  
  lasso_coef <- coef(lasso_model);
  maxf_fit = 0;
  for (k in 1:length(lasso_coef)){
    maxf_fit = maxf_fit+lasso_coef[k]*r$t^(k-1)
  }
  
  if(actPlot){
    image.plot(r$t,y,E,xlab="Time [s]",ylab="Frequency [Hz]",
               cex.lab=1.8, cex.axis=1.5)
    points(t2, maxf2, col='blue',pch=19)
    points(r$t, maxf_fit, col='black',pch=19)
    leg <- c("Close maxima", "Second fit")#,"Uncertainty");
    col <- c("blue","black")#,"gray");
    legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
    
    image.plot(r$t,y,E,xlab="Time [s]",ylab="Frequency [Hz]",
               cex.lab=1.8, cex.axis=1.5)
    points(x, maxf, col='blue',pch=19)
    points(r$t, maxf_fit, col='black',pch=19)
    leg <- c("Maxima", "Fit")#,"Uncertainty");
    col <- c("blue","black")#,"gray");
    legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
  }
  
  if (saveToTxt){
    dir.create(path='oneSim/stdLike/', showWarnings=FALSE, recursive=TRUE);
    write.table(list(x,maxf),'oneSim/stdLike/maxf.txt', row.names = FALSE, col.names = FALSE);
    write.table(list(r$t,maxf_fit),'oneSim/stdLike/fit.txt', row.names = FALSE, col.names = FALSE);
    print('Text files for maxima tracking saved in ./oneSim/stdLike/')
  }
  
  return(maxf_fit);
}




########################################################################
movf = function(vec, n, f){
  ########################################################################
  
  # This function smoothes the numeric vector "vec". 
  # It applies the "f" function in intervals of length "n".
  #
  # vec : numeric vector.
  # n   : interval length for smoothing.
  # f   : function to smooth vec, for instance, mean or median.
  
  if(n == 1){
    return(vec);
  }
  
  N   = length(vec);
  
  if(N < n){
    
    warning("The length of 'vec' must be greater than 'n'");
    return(vec);
  }
  
  out = rep(NA, N);
  
  if( n %% 2 != 0){ # odd
    aux  = (n+1)/2;
    aux1 = aux - 1; 
    
    for(i in aux:(N-aux+1)){
      out[i] = f(vec[(i - aux1):(i + aux1)]);
    }
    
    index  = seq(1, n - 2, by = 2);
    index1 = seq(N - n + 1 + 2, N, by = 2);
    
    for(i in 1:length(index)){
      
      out[i] = f(vec[1:index[i]]); # calculating first values
      
      out[i + N - aux + 1] = f(vec[index1[i]:N]); # calculating last values
      
    }
    
  }else{
    stop("This version only works with odd n");
  }
  
  return(out);
  
}



########################################################################
ints = function(n, l, p = 00, eq = TRUE){
  ########################################################################
  # Function taken from psplinePsd function
  
  #' This function produces a matrix which each row contains the first and last indexes
  #'  to split a time series of length "n"
  #'
  #' n  = time series length
  #' l  = interval length
  #' p  = overlapping percentage
  #' eq = last interval has length "l"
  
  if(n<=0) stop("n must be a positive integer")
  if( (l<=0) || (l %% 1 != 0)) stop("l must be an even positive integer")
  if( (p<0) || (p>=100)) stop("p must be an integer between 0 and 100")
  if( l >= n ) stop("l must be lower than n")
  
  # This version yields eve interval lengths
  if (l %% 2 != 0) stop("this version of bsplinePsd must have l even")
  
  ovp  = round(l*p/100,0); # number of overlaped points
  a    = l - ovp + 1;
  col1 = c(1, seq(from = a, to = n - l + a -1, by = a-1));
  
  index = cbind(col1, col1 + l-1);
  # checking the lengths
  # index[,2] - index[,1] + 1 == l
  # diff(index[,1]) == a-1
  # diff(index[,2]) == a-1
  
  # fixing last interval
  index[dim(index)[1], 2] = n;
  
  colnames(index) = NULL;
  rownames(index) = NULL;
  
  #cat(paste("The number of data subsets is ",dim(index)[1], sep=""),"\n");
  
  if(eq){
    
    lint = dim(index)[1]
    index[lint, 1] = n - l + 1;
    # index[,2] - index[,1] == l
    
    x = index[lint-1, 2] - index[lint, 1] + 1;
    x = round(100 * x / l, 2);
    
    if(x!=p){
      #cat(paste("Last interval overlaps ", x, "%", sep = ""), "\n");
    }
    
  }else{
    
    aux = index[dim(index)[1], 2] - index[dim(index)[1], 1] + 1;
    
    if(aux %% 2 != 0){
      
      # last interval is been made even
      index[dim(index)[1], 1] = index[dim(index)[1], 1] - 1;
      aux = index[dim(index)[1], 2] - index[dim(index)[1], 1] + 1;
    }
    
    if(aux != l){
      cat(paste("Last interval contains ", aux, " observations", sep = ""), "\n");
    }
  }
  
  return(index);
}



########################################################################
covpbb_LASSO = function(r, mod, movBand = 5, true_data, timeGmode = NULL,
                        limFreq = NULL, mask_t = NULL, mask_f = c(0., 1.), 
                        actPlot = FALSE, saveToTxt = FALSE){
  ########################################################################  
  
  # r   : time-frequency map (output of timeFreqMaps.R)
  # mod : fit model describing the evolution of the ratio with frequency
  # movBand   : define the number of points to smooth the band
  # true_data : simulated ratio time evolution where ratio is M/R^2 (g2 mode)
  # timeGMode : time interval to define g-modes
  # limFreq   : specifies upper threshold (in Hz) for the estimated g-modes
  # actPlot   : logical value to produce plot
  
  true_time = true_data$time; 
  true_ratios = true_data$ratio;
  
  timefreq = r$t;
  
  maxf = findGmodes_LASSO(r, lambda=1,
                          mask_t=mask_t, mask_f=mask_f, 
                          actPlot=actPlot, saveToTxt=saveToTxt);
  
  # g-modes
  if( !is.null(timeGmode)){
    #timeGmode = data0[c(1,length(data0[,1])), 1];
    out = apply(as.matrix(timefreq), 1, 
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
    timefreq = timefreq[out];
    maxf = maxf[out]
  }
  
  if(is.null(limFreq)){
    limFreq = Inf; # case in which there is no threshold for frequencies
  }
  
  out1 = NULL; # to store output - coverage probability
  out2 = NULL; # to store output - stat of residuals
  out3 = NULL; # to store output - chi2

  discFreq = rep(TRUE,length(maxf));   # positions to keep
  for (k in 1:length(maxf)){
    if(maxf[k]>limFreq){
      discFreq[k:length(maxf)]=FALSE;
    }
  }
  
  sfq = sum(discFreq);
  
  if(sfq <= 2){ 
    if(sfq == 0){
      #warning(paste("All frequencies are greater than limFreq", j));
    }else{
      #warning(paste("Only", sfq, "frequency is lower than limFreq", j));
    }
    if(any(class(mod) == "lm")){
      out1 = rbind(out1, c(-1, -1));
      out2 = rbind(out2, c(-1, -1, -1));
      out3 = rbind(out3, c(-1));
    }
    else {
      out1 = rbind(out1, c(0, 1));
      out2 = rbind(out2, c(-2, -2, -2));
      out3 = rbind(out3, c(-2));
    }                
  }
  else { # At least 3 g-modes in maxf are required to generate the covpbb band
    maxf1     = maxf[discFreq]; # discarding frequencies according to limFreq
    timefreq1 = timefreq[discFreq]; # discarding time points
    
    # defining true ratios in the band limit given by limfreq 
    discTime = (true_time >= min(timefreq1)) & (true_time <= max(timefreq1));
    true_time1 = true_time[discTime];  
    true_ratio1 = true_ratios[discTime];
    
    # prediction : pred$fit pred$lwr pred$upr
    new  = data.frame(f = maxf1);
    
    if(any(class(mod) == "lm")){
      pred = predict(mod, new, interval = "prediction"); # predictions
    }
    else if(any(class(mod) == "lmvar")){
      f = maxf1;
      
      ### mu ###
      x = colnames(mod$X_mu)
      if(x[1] == "(Intercept)"){
        
        X_mu = apply(as.matrix(x[-1]), 1, function(y)eval(parse(text = y)));
        colnames(X_mu) = colnames(mod$X_mu)[-1];
        
      }else{
        
        X_mu = apply(as.matrix(x), 1, function(y)eval(parse(text = y)));
        colnames(X_mu) = colnames(mod$X_mu);
      }
      
      ### sigma ###
      x = colnames(mod$X_s);
      if(x[1] == "(Intercept_s)"){
        
        X_s = apply(as.matrix(x[-1]), 1, function(y)eval(parse(text = y)));
        colnames(X_s) = colnames(mod$X_s)[-1];
        
      }else{
        
        X_s = apply(as.matrix(x), 1, function(y)eval(parse(text = y)));
        colnames(X_s) = colnames(mod$X_s);
        
      }
      
      pred = predict(mod, X_mu = X_mu, X_sigma = X_s, 
                     interval = "prediction", sigma = FALSE); # predictions
    }
    
    if (saveToTxt){
      dir.create(path='oneSim/stdLike/', showWarnings=FALSE, recursive=TRUE);
      write.table(list(true_time1,true_ratio1),'oneSim/stdLike/true_ratio.txt', row.names = FALSE, col.names = FALSE);
      write.table(list(timefreq1,pred),'oneSim/stdLike/pred.txt', row.names = FALSE, col.names = FALSE);
      print('Text files for ratio reconstruction saved in ./oneSim/stdLike/')
    }
    
    ### generating band function ### 
    
    # interpolating lower bound for predicted values
    fd = approxfun(x = timefreq1, y = movf(pred[,2],n=movBand,mean), method = "linear",
                   yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
    # interpolating upper bound for predicted values
    fu = approxfun(x = timefreq1, y = movf(pred[,3],n=movBand,mean), method = "linear",
                   yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
    # interpolating point estimates
    fm = approxfun(x = timefreq1, y = movf(pred[,1],n=movBand,mean), method = "linear",
                   yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
    
    # fd & fu use smooth confidence intervals by using "movf"
    aux = cbind(true_ratio1,          # true ratios
                fd(true_time1),  # lower band (using predicted ratios)
                fu(true_time1)); # upper band (using predicted ratios)
    
    if(actPlot){
      yaux = c(true_ratios, pred[,2:3]);
      plot(true_time1, true_ratio1, xlab = "Time [s]", xlim=c(min(true_time1)*.9,max(true_time1)*1.05),
           ylab = "Ratio", ylim = c(min(yaux), 1.3*max(true_ratio1)), type = "n",
           cex.lab=1.8, cex.axis=1.5);
      arrows(timefreq1, pred[,2], timefreq1, pred[,3], code=3, angle=90,
             length=0.05, col="gray",pch=3);
      points(true_time1, true_ratio1, col = "black", pch=1);
      points(timefreq1, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3, pch=2);
      
      leg <- c("Simulation", "Estimation ","Uncertainty");
      col <- c("black","red","gray");
      legend("topleft",legend=leg,cex=1.,col=col,pch=c(1,2,3));
    } # end plot
    
    # discarding the true values which are out of the range of the predicted values
    
    # left side
    disc = which(is.na(aux[,2]));
    
    if(length(disc) != 0){
      aux = aux[-disc,];
    }
    
    # right side
    disc = which(is.na(aux[,3]));
    
    if(length(disc) != 0){
      aux = aux[-disc,];
    }
    
    # testing if the true ratios are inside the bands
    prop = apply(aux, 1, 
                 function(x){
                   if((x[1]>= x[2]) && (x[1] <= x[3])){
                     return(1);
                   }else{
                     return(0);
                   }
                 });
    
    aux[aux[,2]<0,2] = 0; # it replaces negative values in lower limit
    
    l = aux[,3] - aux[,2];
    p = mean(prop);
    out1 = rbind(out1, c(p, median(l))); # covpbb & medBandWidth
    
    # Residual, RMS & precision
    res  = true_ratio1 - fm(true_time1); # true_value - estimate
    res_absres = mean(abs(res))
    res_MSE = mean(res^2)
    res_precision = mean(abs(res)/true_ratio1)
    
    out2 = rbind(out2, c(res_absres, res_MSE, res_precision));

  } # end 'all(discFreq)'
  
  # colnames
  colnames(out1) = c("covpbb", "medBandWidth");
  colnames(out2) = c("absres", "MSE", "precision");
  
  R = list(covpbb = out1, residual = out2);
  
  return(R);
}


########################################################################
compute_SNR = function(wvf, detector="LHO", fcut=0, dist=10, 
                       pbOff=FALSE, actPlot=FALSE){
  ########################################################################
  # Compute the Signal-To-Noise ratio for a given wvf (x$time, x$hoft) 
  # and a given detector
  
  fs=round(1/(wvf$time[2]-wvf$time[1]))
  n=length(wvf$hoft)
  a = nextpow2(2*n)         #zero padding and rounding to the next power of 2
  n2=2^a
  
  # If pbOff remove 0.100s after the bounce (set hoft values to 0)
  if (pbOff){
    ext=which(wvf$time<0.1)
    wvf$hoft[ext]=0
  }
  
  freq2 = fs*fftfreq(n2)       # two-sided frequency vector
  freq2[1]=0.001               # to avoid plotting pb in logscale
  freq1=freq2[1:floor(n2/2)]   # one-sided frequency vector
  
  # Get the 1 sided PSD
  psd=PSD_fromfiles(freq1, 1, detector, actPlot)
  
  vec=rep(0,n2)
  for (i in 1:n){
    vec[n2/4+i]=vec[n2/4+i]+wvf$hoft[i]*10./dist
  }  
  
  hf=fft(vec);
  
  hf=hf[1:(n2/2)]   # The integral is performed over positive freqs
  
  hf=subset(hf,freq1-fcut>0)
  psd=subset(psd,freq1-fcut>0)
  freq1=subset(freq1, freq1-fcut>0)
  
  snr=sqrt(4/fs/n2*sum(abs(hf)^2/psd))
  
  if (actPlot){
    plot (freq1, sqrt(freq1)*abs(hf), log="xy", type="l", 
          xlab="Frequency", ylab="hchar", xlim=c(1, fs/2), ylim=c(1e-24,1e-20), 
          col="grey", pch=1, panel.first = grid(),
          cex.lab=1.8, cex.axis=1.5)
    points(freq1,sqrt(psd), type="l", col="black",pch=2)
    leg = c("sqrt(fs) x h~(f)", "ASD")
    col = c("grey","black")
    legend (x=1,y=6e-22,legend=leg,cex=1.,col=col,pch=c(1,2))
    title(c("SNR:",snr))
  }
  return(snr)  
}