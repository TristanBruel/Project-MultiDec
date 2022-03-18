library ("seewave");
library ("signal")
library(lmvar)
library(glmnet)


########################################################################
findGmodes_LASSO = function(r, lambda = 1, 
                            setStart = FALSE, m_L = 2, initfreq_L = c(0, 300),
                            actPlot = FALSE){
  ########################################################################
  
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes
  
  # r : output from specPdgrm
  # lambda : penalization parameter
  
  # setStart : logical value to constrain the starting point
  
  # actPlot   : logical value to produce plot
  
  x = r$t; y = r$f; z = r$E; # values from spectrogram
  
  maxf = y[apply(z, 1, which.max)]; # maximum frequency at each time index
  
  x = x[maxf>0];
  maxf = maxf[maxf>0];
  
  if (length(maxf)==0){ # in case all maxima are equal to zero
    return(rep(0,length(r$t)))
  }
  
  if (setStart){
    # limiting frequencies according to initfreq in y
    inity = (y >= initfreq_L[1]) & (y <= initfreq_L[2]); # logical vector
    
    initz = z[1:m_L, ]; # m first columns of the spectrogram
    
    if(m_L != 1){
      
      if(initfreq_L[1] != -Inf || initfreq_L[2] != Inf){
        initz[, !inity] = -Inf; # discarding certain freq according to initfreq
      }
      mp = apply(initz, 1, which.max); # maximum in first m columns of spectrogram
      mp = round(median(mp)); # position maximum value
      
    }else if(m_L == 1){
      
      initz = z[1,]; # last column of spectrogram
      
      if(initfreq_L[1] != -Inf || initfreq_L[2] != Inf){
        initz[!inity] = -Inf;
      }
      
      mp = which.max(z[1,]); # position maximum value
      
    }
    
    # we'll force the fit to pass through t0 and f[mp]
    x0 = x[1];
    y0 = y[mp];
    
    Xmat0 <- cbind(x-x0, (x-x0)^2, (x-x0)^3, (x-x0)^4, (x-x0)^5, (x-x0)^6, 
                  (x-x0)^7, (x-x0)^8, (x-x0)^9, (x-x0)^10);
    
    lasso_model = glmnet(Xmat0, maxf-y0, alpha=1, lambda=lambda, intercept=FALSE);
    maxf_fit = predict(lasso_model, Xmat0)+y0;
    
    
    if(actPlot){
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]",
                 cex.lab=1.8, cex.axis=1.5)
      points(x, maxf, col='blue',pch=19)
      points(x, maxf_fit, col='black',pch=19)
      points(x0, y0, col = "green", pch = 19); # blue point in the plot
    }
    
    ## Second iteration
    # Compute the rsmd and remove the maxf values further away from the fit
    rmsd = sqrt(sum((maxf_fit-maxf)^2)/length(maxf))
    indices = abs(maxf-maxf_fit)<rmsd;
    indices = abs(maxf-maxf_fit)<100;
    maxf2 = maxf[indices]-y0;
    t2 = x[indices]-x0;
    
    # special case in which all the maxima are the same
    if (length(count(maxf2)$freq)==1){ 
      return(rep(maxf2[1],length(r$t)))
    }
    
    # Fit again
    Xmat2 <- cbind(t2, t2^2, t2^3, t2^4, t2^5, t2^6, 
                   t2^7, t2^8, t2^9, t2^10);
    
    lasso_model = glmnet(Xmat2, maxf2, alpha=1, lambda=lambda, intercept=FALSE);
    lasso_coef <- coef(lasso_model);
    
    maxf_fit = 0;
    for (k in 2:length(lasso_coef)){
      maxf_fit = maxf_fit+lasso_coef[k]*(r$t-x0)^(k-1)
    }
    maxf_fit = maxf_fit+y0;
    
    if(actPlot){
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]",
                 cex.lab=1.8, cex.axis=1.5)
      points(t2+x0, maxf2+y0, col='blue',pch=19)
      points(r$t, maxf_fit, col='black',pch=19)
      points(x0, y0, col = "green", pch = 19); # green point on the plot
      
      leg <- c("Maxima", "Fit")#,"Uncertainty");
      col <- c("blue","black")#,"gray");
      legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
    }
  }
  
  # case without any constraint on the starting point
  else{
    
    # special case in which all the maxima are the same
    if (length(count(maxf)$freq)==1){
      return(rep(maxf[1],length(r$t)))
    }
    # normal case
    Xmat <- cbind(x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10);
    lasso_model = glmnet(Xmat, maxf, alpha=1, lambda=lambda);
    maxf_fit = predict(lasso_model, Xmat);
    
    if(actPlot){
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]",
                 cex.lab=1.8, cex.axis=1.5)
      points(x, maxf, col='blue',type="p")
      points(x, maxf_fit, col='black',type="p")
      
      leg <- c("Maxima", "Fit")#,"Uncertainty");
      col <- c("blue","black")#,"gray");
      legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
    }
    
    ## Second iteration
    # Compute the rsmd and remove the maxf values further away from the fit
    rmsd = sqrt(sum((maxf_fit-maxf)^2)/length(maxf))
    indices = abs(maxf-maxf_fit)<rmsd;
    maxf2 = maxf[indices];
    t2 = x[indices];
    
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
    
    lasso_model = glmnet(Xmat2, maxf2, alpha=1, lambda=lambda);
    lasso_coef <- coef(lasso_model);
    
    maxf_fit = 0;
    for (k in 1:length(lasso_coef)){
      maxf_fit = maxf_fit+lasso_coef[k]*r$t^(k-1)
    }
    maxf_fit = maxf_fit;
    
    
    if(actPlot){
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]",
                 cex.lab=1.8, cex.axis=1.5)
      #points(t2, maxf2, col='blue',pch=19)
      points(x, maxf, col='blue',pch=19)
      points(r$t, maxf_fit, col='black',pch=19)
      #arrows(r$t, pred[,2], r$t, pred[,3], code=3, angle=90,
      #       length=0.05, col="gray",pch=3);
      
      leg <- c("Maxima", "Fit")#,"Uncertainty");
      col <- c("blue","black")#,"gray");
      legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
    }
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
  
  if(eq == TRUE){
    
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
covpbb_LASSO = function(r, mod, setStart = FALSE, m_L = 5, initfreq_L = c(0, Inf),
                       movBand = 5, timeGmode = NULL, true_data, limFreq = NULL,
                       actPlot = FALSE){
  ########################################################################  
  
  # r   : time-frequency map (output of timeFreqMaps.R)
  # mod : fit model describing the evolution of the ratio with frequency
  
  # setStart   : logical value to constrain starting value
  # m_L        : number of time bins used to calculate starting value
  # initfreq_L : interval to restrict frequency for the starting value
  
  # movBand   : define the number of points to smooth the band
  # timeGMode : time interval to define g-modes
  # true_data : simulated ratio time evolution where ratio is M/R^2 (g2 mode)
  # limFreq   : specifies upper threshold (in Hz) for the estimated g-modes
  # actPlot   : logical value to produce plot
  
  ### ###
  true_time = true_data$time; 
  true_ratios = true_data$ratio;
  
  timefreq = r$t;
  
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
    
    r$t = r$t[out];
    r$E = r$E[out, ]; # row=time, col=freq
  }

  
  maxfs = findGmodes_LASSO(r, lambda=1, setStart=setStart, m_L=m_L, 
                           initfreq_L=initfreq_L, actPlot=actPlot);
  
  if(actPlot){
    image.plot(r$t,r$f,r$E,xlab="Time [s]",ylab="Frequency [Hz]",
               cex.lab=1.8, cex.axis=1.5)
    points(timefreq, maxfs, col='black', pch=19)
    
    leg <- c("g2-mode tracking")
    col <- c("black")
    legend("topleft",legend=leg,cex=1.,col=col,pch=c(19,19));
  }
  
  if(is.null(limFreq)){
    
    limFreq = Inf; # case in which there is no threshold for frequencies
    
  }
  
  out1 = NULL; # to store output - coverage probability
  out2 = NULL; # to store output - stat of residuals
  out3 = NULL; # to store output - chi2
  
  maxf = maxfs;
  
  for(j in limFreq){
    #print(j)
    #discFreq = (maxf < j); # positions to keep 
    discFreq = rep(TRUE,length(maxf));   # positions to keep
    for (k in 1:length(maxf)){
      if(maxf[k]>j){
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
        out1 = rbind(out1, c(0, 1));
        out2 = rbind(out2, c(1, 1, 1));
        out3 = rbind(out3, c(1));
      }
      else {
        out1 = rbind(out1, c(0, 1));
        out2 = rbind(out2, c(1, 1, 1));
        out3 = rbind(out3, c(1));
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
      }else if(any(class(mod) == "lmvar")){
        f    = maxf1;
        
        ### mu ###
        x = colnames(mod$X_mu)
        if(x[1] == "(Intercept)"){
          
          X_mu = apply(as.matrix(x[-1]),1, function(y)eval(parse(text = y)));
          colnames(X_mu) = colnames(mod$X_mu)[-1];
          
        }else{
          
          X_mu = apply(as.matrix(x), 1, function(y)eval(parse(text = y)));
          colnames(X_mu) = colnames(mod$X_mu);
        }
        
        ### sigma ###
        x = colnames(mod$X_s);
        if(x[1] == "(Intercept_s)"){
          
          X_s = apply(as.matrix(x[-1]),1, function(y)eval(parse(text = y)));
          colnames(X_s)  = colnames(mod$X_s)[-1];
          
        }else{
          
          X_s = apply(as.matrix(x), 1, function(y)eval(parse(text = y)));
          colnames(X_s)  = colnames(mod$X_s);
          
        }
        
        pred = predict(mod, X_mu = X_mu, X_sigma = X_s, 
                       interval = "prediction", sigma = FALSE); # predictions
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
      
      if(actPlot == TRUE){
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
    
  } # end loop
  
  # colnames
  colnames(out1) = c("covpbb", "medBandWidth");
  colnames(out2) = c("absres", "MSE", "precision");
  
  R = list(covpbb = out1, residual = out2);
  
  return(R);
}