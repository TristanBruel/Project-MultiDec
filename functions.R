library ("seewave");
library ("signal")
source ("data_multiDec.R")
library(lmvar)

########################################################################
findGmodes = function(r, um_R=0, dm_R = 8, um_L = 8, dm_L = 0, m_R = 8, m_L = 8, 
                      initfreq_R = c(-Inf, Inf), initfreq_L = c(-Inf, Inf),
                      testSlope = FALSE){
########################################################################
  
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes
  
  # r    : output from specPdgrm
  # Paremeters for exploration starting from the right side
  # um_R : up   - to define the neighborhood 
  # dm_R : down - to define the neighborhood
  # m_R  : number of intervals used to calculated starting value
  # initfreq_R: interval to restrict frequency for the starting value
  
  # Paremeters for exploration starting from the left side
  # um_L : up   - to define the neighborhood 
  # dm_L : down - to define the neighborhood
  # m_L  : number of intervals used to calculated starting value
  # initfreq_L: interval to restrict frequency for the starting value
  
  #if(initfreq[1]>=initfreq[2])stop("initfreq[1] must be lower than initfreq[2]");
  
  x = r$x; y = r$y; z = r$z; # values from spectrogram
  
#  for(i in 1:dim(z)[1]){
#    i1=which.max(z[i,]); # position maximum value
#    print(c(i1,y[i1],z[i,i1]))
#    }
  
  ##################
  ### Right side ###
  ##################
  
  # limiting frequencies according to initfreq in y
  inity = (y >= initfreq_R[1]) & (y <= initfreq_R[2]); # logical vector
  
  if(m_R != 1){
    
    zn    = dim(z)[1]; # row number (columns in spectrogram) 
    initz = z[seq(zn, zn - m_R + 1), ]; # m last columns of the spectrogram
    
    if(initfreq_R[1] != -Inf || initfreq_R[2] != Inf){
      initz[, !inity] = -Inf; # discarding certain freq according to initfreq
    }
    
    mp    = apply(initz, 1, which.max); # maximum in last m columns of spectrogram
    mp    = round(median(mp)); # position maximum value
  }else if(m_R == 1){
    
    initz = z[dim(z)[1],]; # last column of spectrogram
    
    if(initfreq_R[1] != -Inf || initfreq_R[2] != Inf){
      initz[!inity] = -Inf;
    }
    
    mp = which.max(z[dim(z)[1],]); # position maximum value
    
  }
  ps = NULL;
  for(i in (dim(z)[1]-1):1){
    
    #print(i);
    ps  = c(ps, mp);
    ind = (mp-dm_R):(mp+um_R);
    ind = ind[ind>0]; # preventing values below 0
    ind = ind[ind<dim(z)[2]]; # preventing from large values
    
    #ind = ind[which(y[ind]<initfreq_R[2])];#preventing indexes to explore region forbidden by initfreq array
    #ind = ind[which(y[ind]>initfreq_R[1])];#preventing indexes to explore region forbidden by initfreq array
    
    mp  = ind[which.max(z[i, ind])];
  }
  
  ps     = c(ps, mp); # adding last mp value
  ps_R   = ps[length(ps):1]; # reordering
  maxf_R = y[ps_R]; # FREQUENCIES FOR LARGEST POWER
  
  #################
  ### Left side ###
  #################

 
  # limiting frequencies according to initfreq in y
  inity = (y >= initfreq_L[1]) & (y <= initfreq_L[2]); # logical vector
  
  initz = z[1:m_L, ]; # m first columns of the spectrogram

  if(m_L != 1){
    
    if(initfreq_L[1] != -Inf || initfreq_L[2] != Inf){
      initz[, !inity] = -Inf; # discarding certain freq according to initfreq
    }
    mp    = apply(initz, 1, which.max); # maximum in first m columns of spectrogram
    mp    = round(median(mp)); # position maximum value
    
  }else if(m_L == 1){
    
    initz = z[1,]; # last column of spectrogram
    
    if(initfreq_L[1] != -Inf || initfreq_L[2] != Inf){
      initz[!inity] = -Inf;
    }
    
    mp = which.max(z[1,]); # position maximum value
    
  }
  
  ps = NULL;
  
  for(i in 2:dim(z)[1]){
    #if (length((mp-dm_L):(mp+um_L)) > 0){
      #print(i);
      ps  = c(ps, mp);
      ind = (mp-dm_L):(mp+um_L);
      ind = ind[ind>0]; # preventing values below 0
      ind = ind[ind<=dim(z)[2]]; # preventing from large values
      mp  = ind[which.max(z[i, ind])];
      #print(c(z[i,ind],mp,y[mp]))
    #}
    #else{
    #  print(c("dm_L+um_L is null",mp,dm_L, um_L))
    #  break
    #}
  } 
  
  ps_L     = c(ps, mp); # adding last mp value
  maxf_L = y[ps_L]; # FREQUENCIES FOR LARGEST POWER
  
  if(testSlope == TRUE){
    
    #####################
    ### Testing Trend ###
    #####################
    
    # H0: beta < 0
    # H1: beta > 0
    
    n_maxf = length(maxf_R);
    x_maxf = 1:n_maxf;
    
    alpha = 0.001;
    
    ### R: starting from the right ###
    
    mod_R = lm(maxf_R~x_maxf);
    fit_R = summary(mod_R);
    
    p_R = pt(coef(fit_R)[2, 3], mod_R$df, lower = FALSE); # pvalue
    
    if((coef(fit_R)[2, 4] > alpha) || (p_R > alpha)){
      # H0:beta=0 is not rejected or H0:beta < 0 is not rejected
      R = FALSE; 
    }else{ # All ok
      # H0:beta=0 is rejected and H0:beta<0 is rejected
      R = TRUE;
    }
    
    #print(p_R); print(coef(fit_R)[2, 4]);
    
    ### L: starting from the left ###
    
    mod_L = lm(maxf_L~x_maxf);
    fit_L = summary(mod_L);
    
    p_L = pt(coef(fit_L)[2, 3], mod_L$df, lower = FALSE); # pvalue
    
    if((coef(fit_L)[2, 4] > alpha) || (p_L > alpha)){
      L = FALSE; 
    }else{
      L = TRUE;
    }
    
    #print(R); print(L);
    
    if( R == FALSE && L == FALSE ){
      cat("None of them work", "\n");
      warning("Both method find estimated g-modes with negative trend");
      maxf = rep(NA, n_maxf);
      
    }else if( R == TRUE && L == FALSE ){ # only Rigth works
      cat("only rigth works", "\n");
      maxf = maxf_R;
      
    }else if( R == FALSE && L == TRUE ){ # only Left works
      
      cat("only left works", "\n");
      maxf = maxf_L;
      
    }else if( R == TRUE && L == TRUE ){ # Both works
      
      cat("Both works", "\n");
      # median of both maxf
      ps   = apply(rbind(ps_L, ps_R), 2, function(x) round(median(x)));
      maxf = y[ps]; 
    }
  }else{
    
    ps   = apply(rbind(ps_L, ps_R), 2, function(x) round(median(x)));
    maxf = y[ps]; 
    
  }
  return(list(maxf_med = maxf, maxf_L = maxf_L, maxf_R = maxf_R));
  
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
  
  if( (n<=0) || (n %% 1 != 0)) stop("n must be a positive integer")
  if( (l<=0) || (l %% 1 != 0)) stop("l must be a positive integer")
  if( (p<0) || (p>=100)) stop("p must be a positive number between 0 and 100")
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
covpbb = function(data, mod, l=200, p=90, fs=16384, movGmode = 11, 
                  um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(-Inf, Inf),
                  um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(-Inf, Inf),
                  gmode = "right",
                  movBand = 5, timeGmode = NULL, 
                  thruth_data, actPlot = FALSE, limFreq = NULL){
########################################################################  
  
  # data     : dataset as matrix (with no zeros [specPdgrm])
  # mod      : model
  # l        : interval length in spectrogram
  # p        : overlaping percentage in spectrogram
  # fs       : sampling frequency
  # movGmode : number of steps to smooth estimated g-modes 
  
  # In 'findGmodes' function
  # Paremeters for exploration starting from the right side
  # um_R : up   - to define the neighborhood 
  # dm_R : down - to define the neighborhood
  # m_R  : number of intervals used to calculated starting value
  # initfreq_R: interval to restrict frequency for the starting value
  
  # Paremeters for exploration starting from the left side
  # um_L : up   - to define the neighborhood 
  # dm_L : down - to define the neighborhood
  # m_L  : number of intervals used to calculated starting value
  # initfreq_L: interval to restrict frequency for the starting value
  
  # gmode = starting side to estimate g-mode: c("right", "left, "median")
  
  # movBand  : define the number of points to smooth the band
  # timeGmode: time interval to define g-modes
  # thruth_data: true ratio and time
  # thruth_data: simulated ratio time evolution where ratio is M/R^2 (g2 mode) or (g3 mode)
  # actPlot  : logical value to produce plot
  # limFreq  : specifies upper threshold (in Hz) for the estimated g-modes
  
  ### ###
  true_time = thruth_data$time; 
  true_ratios = thruth_data$ratio

  # spectrogram
  r = specPdgrm(data$V2, data$V1, l=l, p=p, fs=fs, actPlot=actPlot, logPow=TRUE,
                zoomFreq=c(0,1)); # generating the spectrogram
  
  n = length(data$V1); # number of observations
  
  # starting & ending points of the intervals used in the spectrogram
  index = ints(n=n, l=l, p=p); # from psplinePsd
  
  # centered point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  
  # Adjusting data
  timedata = seq(from=data$V1[1], to=data$V1[n], length = n);
  
  # mean time (centered point) for our g-mode estimates
  timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );
  
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
    
    #maxf     = maxf[out];
    timefreq = timefreq[out];
    
    r$x = r$x[out];
    r$z = r$z[out, ]; # row=time, col=freq
  }
  
  maxfs = findGmodes(r, um_R, dm_R, um_L, dm_L, m_R, m_L, 
                     initfreq_R, initfreq_L, testSlope = FALSE);
  
  maxfs = sapply(maxfs,function(x)movf(x, n=movGmode, mean));# smoothing g-mode estimates
  maxfs = as.data.frame(maxfs);

  if(actPlot == TRUE){
    if(gmode == "left"){
      points(timefreq, maxfs$maxf_L, col='black',type="p")
    }else if (gmode == "right"){
      points(timefreq, maxfs$maxf_R, col='black')
      }else{
        points(timefreq, maxfs$maxf_med, col='black')}
    #arrows(timefreq, maxfs$maxf_L, timefreq, maxfs$maxf_R, code=3, angle=90,
    #       length=0.05, col="gray");
  }
  
  if(is.null(limFreq)){
    
    limFreq = Inf; # case in which there is no threshold for frequencies
    
  }
  
  R = list();
  
  for(gm in gmode){
    
    out1 = NULL; # to store output - coverage probability
    out2 = NULL; # to store output - stat of residuals
    out3 = NULL; # to store output - chi2
    
    if(gm == "left"){
      
      maxf = maxfs$maxf_L;
      
    }else if(gm == "right"){
      
      maxf = maxfs$maxf_R;
      
    }else if(gm == "median"){
      
      maxf = maxfs$maxf_med;
      
    }else{
      
      stop("The gmode must be 'left', 'right' or 'median'");
      
    }  
    
    for(j in limFreq){
      #print(j)
      discFreq = (maxf < j); # positions to keep 
      
      sfq = sum(discFreq);
    
      if(sfq <= 2){ 
        if(sfq == 0){
          warning(paste("All frequencies are greater than limFreq", j));
          print("All frequencies are greater than limFreq")
        }else{
          warning(paste("Only", sfq, "frequency is lower than limFreq", j));
        }
        if(any(class(mod) == "lm")){
          out1 = rbind(out1, c(-1, -1));
          out2 = rbind(out2, c(-1, -1));
          out3 = rbind(out3, c(-1));
        }
        else {
          out1 = rbind(out1, c(-2, -2));
          out2 = rbind(out2, c(-2, -2));
          out3 = rbind(out3, c(-2));
        }                
      }
      else { # At least 3 g-modes (in maxf) are required to generate the covpbb band
        maxf1     = maxf[discFreq];     # discarding frequencies according to limFreq
        timefreq1 = timefreq[discFreq]; # discarding time points

        # defining true ratios in the band limit given by limFreq 
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
          plot(true_time1, true_ratio1, xlab = "Time", xlim=c(min(true_time1)*.9,max(true_time1)*1.05),
               ylab = "Ratio", ylim = c(min(yaux), 1.3*max(true_ratio1)), type = "n");
      #         main = paste("Frequency cutoff", j,"-",gm, "gmode"));
          arrows(timefreq1, pred[,2], timefreq1, pred[,3], code=3, angle=90,
                 length=0.05, col="gray",pch=3);
          points(true_time1, true_ratio1, col = "black", pch=1);
          points(timefreq1, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3, pch=2);
          
          leg <- c("Simulation", "Estimation ","Uncertainty");
          col <- c("black","red","gray");
          legend("topleft",legend=leg,cex=.8,col=col,pch=c(1,2,3));
          
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
                
#        if(actPlot == TRUE){
#          plot(true_time1, res, xlab = "Time",
#               ylab = "Residual", ylim = c(-8e-4, 8e-4), xlim=c(0,max(true_time1)*1.1), type = "n",
#               main = paste("Frequency cutoff", j,"-",gm, "gmode"));
#          points(true_time1, res, col = "black", pch=1);
#        }
      } # end 'all(discFreq)'
      
    } # end loop
    
    # colnames
    colnames(out1) = c("covpbb", "medBandWidth");
    colnames(out2) = c("absres", "MSE", "precision");

    

    
    R[[gm]] = list(covpbb = out1, residual = out2);
    
  }
  
  if(length(gmode) == 1){
    R = R[[1]];
  }
  return(R);
  
}

# THIS FUNCTION MUST BE UPDATED

########################################################################
repcovpbb = function(wvf, duration, ampl, fcut,
                     mod, N, snr = NULL, movGmode = 11, um = 3, dm = 3,
                     movBand = 5, timeGmode = NULL, l=200, p=90, fs=16384, limFreq=NULL,
                     thruth_data){
########################################################################  
  # data : matrix (time, signal)
  # mod  : lm object
  # N    : number of replications
  # ampl : distances (scalar or vector)
  
  if(length(ampl) == 1){
    
    out_cp = NULL;
    out_bl = NULL;
    
    for(i in 1:N){
      noisydata = data_generator(fs, duration, wvf, ampl, filtering = "HP", fcut=fcut, actPlot=FALSE);
      noisydata = data.frame("V1"=noisydata$t,"V2"=noisydata$y);
      if (is.null(limFreq)){
        aux = covpbb(noisydata, mod, l, p, fs, um, dm, 
                     thruth_data=thruth_data, actPlot = FALSE);
      }else{        
        aux = covpbb(noisydata, mod, l, p, fs, um, dm, 
                     thruth_data=thruth_data, actPlot = FALSE);
      }
      out_cp = c(out_cp, aux$covpbb);
      out_bl = c(out_bl, aux$medBandWidth);      
    }
    
    out = list(covpbb = out_cp, medBandWidth = out_bl);
    
  }else if(length(ampl) > 1){
    
    out_cp = list();
    out_bl = list();
    
    for(j in ampl){
      print(j)
      aux_cp = NULL; # to save coverage probabilities
      aux_bl = NULL; # to save band widths
      
      for(i in 1:N){
        
        noisydata = data_generator(fs, duration, wvf, ampl = j, filtering = "HP", fcut=fcut, actPlot=FALSE);
        noisydata = data.frame("V1"=noisydata$t,"V2"=noisydata$y);
        
        if (is.null(limFreq)){
          aux = covpbb(noisydata, mod, l, p, fs, um, dm, 
                       thruth_data=thruth_data, actPlot = FALSE);
        }else{        
          aux = covpbb1(noisydata, mod, l, p, fs, um, dm, 
                        thruth_data=thruth_data, actPlot = FALSE, limFreq);
        }
        
        #print(aux$covpbb)        
        aux_cp = c(aux_cp, aux$covpbb);
        aux_bl = c(aux_bl, aux$medBandWidth);
        
      }
      if(is.null(snr)){
        
        out_cp[[paste(j)]] = aux_cp;
        out_bl[[paste(j)]] = aux_bl;
        
      }else{
        
        out_cp[[paste(snr * j)]] = aux_cp;
        out_bl[[paste(snr * j)]] = aux_bl;  
        
      }
    }
    
    out = list(covpbb = out_cp, medBandWidth = out_bl);
  }
  
  return(out)
  
}



##########################
### BAYESIAN functions ###
##########################

covpbbBayes = function(data, postSamples, l=200, p=90, fs=16384, movGmode = 11, 
                       um = 3, dm = 3, movBand = 5, timeGmode = NULL,
                       thruth_data, actPlot = FALSE, fmean = 0){
  
  # The only difference with 'covpbb' function is how the predicted 
  #  values are calculated.
  
  # data     : dataset as matrix (with no zeros [specPdgrm])
  # mod      : model
  # l        : interval length in spectrogram
  # p        : overlaping percentage in spectrogram
  # movGmode : number of steps to smooth estimated g-modes 
  # um       : define upper neighborhood to find g-modes 
  # dm       : define lower neighborhood to find g-modes
  # movBand  : define the number of points to smooth the band
  # timeGMode: time interval to define g-modes
  # Oberganlinger_data: simulated R and M time evolution
  # fmean: mean of frequencies used to generate the model
  
  # Compute true ratios
  #true_ratios = Obergaunlinguer_data$Mpns / (Obergaunlinguer_data$Rpns^(2));
  true_ratios = thruth_data$x
  
  # spectrogram
  r = specPdgrm(data$V2, data$V1, l = l, p = p, fs = fs, actPlot = FALSE, logPow = T,
                zoomFreq=c(0,1)); # generating the spectrogram
  
  n     = length(data$V1);
  # starting & ending points of the intervals used in the spectrogram
  index = ints(n = n, l = l, p = p); # from psplinePsd 
  
  # centred point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  # Ajusting data / #it can be directly taken from the 'data'
  timedata = seq(from=data$V1[1], to=data$V1[n], length = n);
  
  # mean time (centred point) for our g-mode estimates
  timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );
  
  # g-modes
  
  if( !is.null(timeGmode)){
    
    # g-mode time range
    
    out = apply(as.matrix(timefreq), 1, 
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
    
    timefreq = timefreq[out];
    
    r$x = r$x[out];   # discarding values out of g-mode range
    r$z = r$z[out, ];
    
  }
  
  maxf = findGmodes(r, um=um, dm=dm);
  maxf = movf(maxf, movGmode, median); # smoothing g-mode estimates
  
  cmaxf = maxf - fmean; # centred g-modes (SEE MODEL USED IN JAGS)
  
  xs =  rbind(rep(1, length = length(cmaxf)), cmaxf, cmaxf^2, cmaxf^3); 
  
  # postSamples = (a0, a1, a2, a3, b1); # postSample structure
  postSamples = as.matrix(postSamples); # needed for matrix operations below
  
  # a0+ a1*(fnew[i]-mean(f))+a2*(fnew[i]-mean(f))^2+a3*(fnew[i]-mean(f))^3
  mus = postSamples[,1:4] %*% xs; # means
  
  # standard deviations
  sds = as.matrix(postSamples[,5], ncol=1) %*% t(as.matrix(maxf,ncol=1));
  sds = sqrt(sds);
  
  # simulation from Normal(mus, sds)
  n_maxf = length(maxf);
  X      = apply(cbind(mus, sds), 1, function(x)
    rnorm(n_maxf, mean = x[1:n_maxf], sd = x[-(1:n_maxf)]));
  X      = t(X);
  
  # median, quantiles 0.025 & 0.975
  pred = apply(X, 2, function(x) quantile(x ,probs = c(0.5,0.025,0.975)));
  pred = t(pred);
  
  #pred[pred[,2] < 0, 2] = 0; # discarding negative values in CI
  
  ### generating band function ### 
  
  # interpolating lower bound for predicted values
  fd = approxfun(x = timefreq, y = movf(pred[,2],n=movBand,median), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  # interpolating upper bound for predicted values
  fu = approxfun(x = timefreq, y = movf(pred[,3],n=movBand,median), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  
  # fd & fu use smooth confidence intervals by using "movf"
  
  aux = cbind(true_ratios,                    # true ratios
              fd(Obergaunlinguer_data$time),  # lower band (using predicted ratios)
              fu(Obergaunlinguer_data$time)); # upper band (using predicted ratios)
  
  # discarding the true values which are out of the range of the predicted values
  
  discL = which(is.na(aux[,2])); # left side
  discR = which(is.na(aux[,3])); # right side
  
  disc  = unique(c(discL, discR)); 
  
  if(length(disc) != 0){
    
    aux = aux[-disc, ];
    
    trueRatioTime = Obergaunlinguer_data$time[-disc];
    
  }
  
  out = NULL;
  
  if( !is.null(timeGmode)){
    
    # To test only true values inside g-mode time range
    
    out = apply(as.matrix(trueRatioTime), 1,
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
    
    aux = aux[out, ];
    
  }
  
  if(actPlot == TRUE){
    
    plot(Obergaunlinguer_data$time, true_ratios, 
         ylim=c(min(pred[,2]),max(pred[,3])), xlab = "Time", ylab = "Ratio",
         main = "Bayesian analysis");
    arrows(timefreq, pred[,2], timefreq, pred[,3], 
           code=3, angle=90, length=0.05, col="gray", pch = 3);
    
    points(Obergaunlinguer_data$time, true_ratios, col = "black",pch=1);
    points(timefreq, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3,pch=2);
    
    leg <- c("true ratio", "pred ","pred uncertainty")
    col=c("black","red","gray")
    legend(x=.25,y=0.0037,legend=leg,cex=.8,col=col,pch=c(1,2,3))
    
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
  
  #aux[aux[,2]<0,2] = 0; # it replaces negative values in lower limit
  
  l = aux[,3] - aux[,2];
  p = mean(prop);
  
  return(list(covpbb = p, medBandWidth = median(l)));
  
}


repcovpbbBayes = function(data, postSamples, N, dists, snr = NULL, movGmode = 11, um = 3, dm=3,
                          movBand = 5, timeGmode = NULL, l=200, p=90, fs=16384){
  
  # data : matrix (time, signal)
  # mod  : lm object
  # N    : number of replications
  # dists: distances (scalar or vector)
  
  Fs   = fs;   # 16384; # sampling frequency
  mean = 0.0;
  std  = 1.0;
  n    = dim(data)[1]; # 16384
  
  freq    = Fs*np.fft.fftfreq(n);  # two-sided frequency vector
  psd2    = aLIGO_PSD(freq,2);     # two-sided PSD
  newpsd2 = aLIGO_PSD_new(freq,2); # two-sided PSD              
  
  #ff     = freq[1:int(n/2)];    # one-sided frequency vector 
  #psd    = aLIGO_PSD(ff,1);     # one-sided PSD              
  #newpsd = aLIGO_PSD_new(ff,1); # one-sided PSD 
  
  noisydata = data;
  
  if(length(dists) == 1){
    
    out = NULL;
    
    newdata = data$V2 * dists;
    
    for(i in 1:N){
      
      # Generate noise
      X   = rnorm(n, mean = mean, sd = std); # Gaussian white noise
      XX  = sqrt(Fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
      XXX = XX*sqrt(newpsd2);                # Coloring
      Y   = fft(XXX, inverse = TRUE);        # FFT inverse
      Y   = Re(Y)*sqrt(n);
      
      # noisy data
      noisydata$V2 = newdata + Y; # signal + noise
      hpf = 50; # high pass filter
      noisydata$V2 = ffilter(noisydata$V2, f = fs, from = hpf);
   
      out = c(out, covpbbBayes(noisydata, postSamples, l=l, p=p, fs=fs, movGmode = movGmode, 
                               um = um, dm = dm, movBand = movBand, timeGmode = timeGmode));
      
    }
    
  }else if(length(dists) > 1){
    
    out = list();
    
    for(j in dists){
      
      print(j);
      
      aux = NULL;
      
      newdata = data$V2 * j;
      
      for(i in 1:N){
        
        # Generate noise
        X   = rnorm(n, mean = mean, sd = std); # Gaussian white noise
        XX  = sqrt(Fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
        XXX = XX*sqrt(newpsd2);                # Coloring
        Y   = fft(XXX, inverse = TRUE);        # FFT inverse
        Y   = Re(Y)*sqrt(n);
        
        # noisy data
        noisydata$V2 = newdata + Y; # signal + noise
        hpf = 50; # high pass filter
        noisydata$V2 = ffilter(noisydata$V2, f = fs, from = hpf);
        
        aux = c(aux, covpbbBayes(noisydata, postSamples, l=l, p=p, fs=fs, movGmode = movGmode, 
                                 um = um, dm = dm, movBand = movBand, timeGmode = timeGmode));
        
      }
      
      if(is.null(snr)){
        
        out[[paste(j)]] = aux;
        
      }else{
        
        out[[paste(snr * j)]] = aux;  
        
      }
    }
  }
  
  return(out)
  
}

plotOneSim = function(data, mod = NULL, postSamples=NULL, dist=1, l=200, p=90,
                      fs=16384, movGmode = 11, um = 3, dm = 3, movBand = 5,
                      timeGmode = NULL, filterHz = NULL){
 
  # data: original data. The functions injects coloured noise

  Fs   = fs; # 16384; # sampling frequency
  mean = 0.0;
  std  = 1.0;
  n    = dim(data)[1]; # 16384
 
  freq    = Fs*fftfreq(n);  # two-sided frequency vector
  psd2    = aLIGO_PSD(freq,2);     # two-sided PSD
  newpsd2 = aLIGO_PSD_new(freq,2); # two-sided PSD             
  
  ff     = freq[1:int(n/2)];    # one-sided frequency vector
  psd    = aLIGO_PSD(ff,1);     # one-sided PSD             
  newpsd = aLIGO_PSD_new(ff,1); # one-sided PSD
  
  noisydata = data;
  noisydata$V2 = noisydata$V2 * dist;

  out = NULL;

  # Generate noise
  X   = rnorm(n, mean = mean, sd = std); # Gaussian white noise
  XX  = sqrt(Fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
  XXX = XX*sqrt(newpsd2);                # Coloring
  Y   = fft(XXX, inverse = TRUE);        # FFT inverse
  Y   = Re(Y)*sqrt(n);

  # noisy data
  noisydata$V2 = noisydata$V2 + Y; # signal + noise
 
  if(!is.null(filterHz)){
    
    noisydata$V2 = ffilter(noisydata$V2, f = fs, from = filterHz);
   
  }
 
  r = specPdgrm(noisydata$V2, l = l, p = p, fs = fs, actPlot = FALSE, logPow = T,
                zoomFreq=c(0,1)); # generating the spectrogram
 
  n     = length(data$V1);
  # starting & ending points of the intervals used in the spectrogram
  index = ints(n = n,l = l, p = p); # from psplinePsd
 
  # centred point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  # Ajusting data
  timedata = seq(from=data$V1[1], to=data$V1[n], length = n);
 
  # mean time (centred point) for our g-mode estimates
  timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );

  # g-modes
 
  if( !is.null(timeGmode) ){
   
    #timeGmode = data0[c(1,length(data0[,1])), 1];
   
    out = apply(as.matrix(timefreq), 1,
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
 
    #maxf     = maxf[out];
    timefreq = timefreq[out];
   
  }
 
  r$x = r$x[out];
  r$z = r$z[out, ];
 
  maxf = findGmodes(r, um=um, dm=dm);
  maxf = movf(maxf, movGmode, median); # smoothing g-mode estimates
  points(timefreq,maxf, col = 'green')
 
  if( class(mod) == "lm"){
   
    # prediction
    new   = data.frame(f = maxf);
    pred  = predict(mod, new, interval = "prediction"); # predictions
    title = "Frequentist analysis";
   
  }else if(class(postSamples) == "data.frame"){
   
    # prediction
    aux  = as.matrix(postSamples[,1:3]) %*% t(cbind(maxf, maxf^2, maxf^3));
   
    aux1 = apply(aux,2,median);
   
    sdb  = 1.96 * sqrt(maxf * median(postSamples[,4]));
   
    pred = cbind(apply(aux,2,median), aux1 - sdb, aux1 + sdb);
   
    title = "Bayesian analysis";
  }
 
  ### PLOT ###
 
  yaux = c(true_ratios, pred[,2:3]);
  plot(Obergaunlinguer_data$time, true_ratios, main = paste(title), xlab = "Time",
       ylab = "Ratio", ylim = c(min(yaux), max(yaux)), type = "n");
  arrows(timefreq, pred[,2], timefreq, pred[,3], code=3, angle=90,
         length=0.05, col="gray");
  points(Obergaunlinguer_data$time, true_ratios, col = "black");
  points(timefreq, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3);
 
}
