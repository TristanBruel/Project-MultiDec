library ("seewave");
library ("signal")
library(lmvar)
library(glmnet)

########################################################################
findGmodes = function(r, um_R = 0, dm_R = 8, m_R = 8, um_L = 8, dm_L = 0, m_L = 8, 
                      initfreq_R = c(-Inf, Inf), initfreq_L = c(-Inf, Inf),
                      testSlope = FALSE, actPlot = FALSE, gmode="left"){
  ########################################################################
  
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes
  
  # r    : output from specPdgrm
  # Paremeters for exploration starting from the right side
  # um_R : up   - to define the neighborhood 
  # dm_R : down - to define the neighborhood
  # m_R  : number of intervals used to calculate starting value
  # initfreq_R: interval to restrict frequency for the starting value
  
  # Paremeters for exploration starting from the left side
  # um_L : up   - to define the neighborhood 
  # dm_L : down - to define the neighborhood
  # m_L  : number of intervals used to calculate starting value
  # initfreq_L: interval to restrict frequency for the starting value
  
  x = r$t; y = r$f; z = r$E; # values from spectrogram
  
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
    ps = c(ps, mp);
    ind = (mp-dm_L):(mp+um_L);
    ind = ind[ind>0]; # preventing values below 0
    ind = ind[ind<=dim(z)[2]]; # preventing from large values
    mp = ind[which.max(z[i, ind])];
    
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
  }
  else{
    ps   = apply(rbind(ps_L, ps_R), 2, function(x) round(median(x)));
    maxf = y[ps];
  }
  
  
  if(actPlot){
    image.plot(r$t,r$f,r$E,xlab="Time [s]",ylab="Frequency [Hz]")
    if(gmode == "left"){
      points(r$t, maxf_L, col='black',type="p")
    }else if (gmode == "right"){
      points(r$t, maxf_R, col='black')
    }else{
      points(r$t, maxf, col='black')}
    #arrows(timefreq, maxfs$maxf_L, timefreq, maxfs$maxf_R, code=3, angle=90,
    #       length=0.05, col="gray");
  }
  
  return(list(maxf_med = maxf, maxf_L = maxf_L, maxf_R = maxf_R));
  
}


########################################################################
findGmodes2 = function(r, um_R = 0, dm_R = 8, m_R = 8, um_L = 8, dm_L = 0, m_L = 8,
                       initfreq_R = c(-Inf, Inf), initfreq_L = c(-Inf, Inf),
                       timeGmode = c(0, Inf), gmode = "left",
                       testSlope = FALSE, actPlot = FALSE){
  ########################################################################
  
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes
  
  # r    : output from specPdgrm
  # Paremeters for exploration starting from the right side
  # um_R : up   - to define the neighborhood 
  # dm_R : down - to define the neighborhood
  # m_R  : number of intervals used to calculate starting value
  # initfreq_R: interval to restrict frequency for the starting value
  
  # Paremeters for exploration starting from the left side
  # um_L : up   - to define the neighborhood 
  # dm_L : down - to define the neighborhood
  # m_L  : number of intervals used to calculate starting value
  # initfreq_L: interval to restrict frequency for the starting value
  
  x = r$t; y = r$f; z = r$E; # values from spectrogram

  
  ##################
  ### Right side ###
  ##################
  
  # limiting frequencies according to initfreq in y
  inity = (y >= initfreq_R[1]) & (y <= initfreq_R[2]); # logical vector
  
  if(m_R != 1){
    
    zn = dim(z)[1]; # row number (columns in spectrogram) 
    initz = z[seq(zn, zn - m_R + 1), ]; # m last columns of the spectrogram
    
    if(initfreq_R[1] != -Inf || initfreq_R[2] != Inf){
      initz[, !inity] = -Inf; # discarding certain freq according to initfreq
    }
    
    mp = apply(initz, 1, which.max); # maximum in last m columns of spectrogram
    mp = round(median(mp)); # position maximum value
    
  }else if(m_R == 1){
    
    initz = z[dim(z)[1],]; # last column of spectrogram
    
    if(initfreq_R[1] != -Inf || initfreq_R[2] != Inf){
      initz[!inity] = -Inf;
    }
    
    mp = which.max(z[dim(z)[1],]); # position maximum value
    
  }
  
  ps = NULL;
  for(i in (dim(z)[1]-1):1){
    
    ps = c(ps, mp);
    #ind = (mp-dm_R):(mp+um_R);
    ind = (mp+um_R):(mp-dm_R);
    ind = ind[ind>0]; # preventing values below 0
    ind = ind[ind<dim(z)[2]]; # preventing from large values
    
    mp = ind[which.max(z[i, ind])];
  }
  
  ps = c(ps, mp); # adding last mp value
  ps_R = ps[length(ps):1]; # reordering
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
    mp = apply(initz, 1, which.max); # maximum in first m columns of spectrogram
    mp = round(median(mp)); # position maximum value
    
  }else if(m_L == 1){
    
    initz = z[1,]; # last column of spectrogram
    
    if(initfreq_L[1] != -Inf || initfreq_L[2] != Inf){
      initz[!inity] = -Inf;
    }
    
    mp = which.max(z[1,]); # position maximum value
    
  }
  
  ps = NULL;
  for(i in 2:dim(z)[1]){
    ps = c(ps, mp);
    ind = (mp-dm_L):(mp+um_L);
    ind = ind[ind>0]; # preventing values below 0
    ind = ind[ind<=dim(z)[2]]; # preventing from large values
    mp = ind[which.max(z[i, ind])];
  }
  
  ps_L = c(ps, mp); # adding last mp value
  maxf_L = y[ps_L]; # FREQUENCIES FOR LARGEST POWER
  
  ps   = apply(rbind(ps_L, ps_R), 2, function(x) round(median(x)));
  #ps = ps_L;
  maxf = y[ps];
  
  
  if(actPlot){
    image.plot(r$t,r$f,r$E,xlab="Time [s]",ylab="Frequency [Hz]")
    if(gmode == "left"){
      points(r$t, maxf_L, col='black',type="p")
    }else if (gmode == "right"){
      points(r$t, maxf_R, col='black')
    }else{
      points(r$t, maxf, col='black')}
    #arrows(timefreq, maxfs$maxf_L, timefreq, maxfs$maxf_R, code=3, angle=90,
    #       length=0.05, col="gray");
  }
  
  return(list(maxf_med = maxf, maxf_L = maxf_L, maxf_R = maxf_R));
  
}



########################################################################
findGmodes_poly = function(r, N = 3, 
                           setStart = FALSE, m_L = 2, initfreq_L = c(0, Inf), 
                           actPlot = FALSE){
  ########################################################################
  
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes
  
  # r : output from specPdgrm
  # N : order of the polynomial fit
  
  # Paremeters for exploration starting from the left side
  # m_L  : number of intervals used to calculated starting value
  # initfreq_L: interval to restrict frequency for the starting value
  
  x = r$t; y = r$f; z = r$E; # values from spectrogram
  
  maxf = y[apply(z, 1, which.max)]; # maximum frequency at each time
  
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
    
    model = lm(I(maxf-y0)~ I(x-x0)+I((x-x0)^2)+I((x-x0)^3)-1);
    maxf_fit = predict(model, newdata = list(x)) + y0;
    
    
    if(actPlot){
      image.plot(x,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(x, maxf, col='blue',type="p")
      points(x, maxf_fit, col='black',type="p")
      points(x0, y0, col = "green", pch = 19); # blue point in the plot
    }
    
    # second iteration
    indices = abs(maxf-maxf_fit)<250;
    maxf2 = maxf[indices]-y0;
    t2 = x[indices]-x0;
    
    model = lm(maxf2~ t2+I(t2^2)+I(t2^3)-1);
    coeff = model$coefficients;
    a1 = coeff[1];
    a2 = coeff[2];
    a3 = coeff[3];
    maxf_fit = a1*(x-x0)+a2*(x-x0)^2+a3*(x-x0)^3 + y0;
    
    if(actPlot){
      image.plot(x,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(t2+x0, maxf2+y0, col='blue',type="p")
      points(x, maxf_fit, col='black',type="p")
      points(x0, y0, col = "green", pch = 19); # blue point in the plot
    }
  }
  
  else{
    model = lm(maxf~ stats::poly(r$t,N));
    maxf_fit = model$fitted.values;
    
    if(actPlot){
      image.plot(x,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(x, maxf, col='blue',type="p")
      points(x, maxf_fit, col='black',type="p")
    }
    
    # second iteration
    indices = abs(maxf-maxf_fit)<250;
    maxf2 = maxf[indices];
    t2 = r$t[indices];
    
    model = lm(maxf2~ stats::poly(t2,3));
    new = data.frame(t2 = r$t)
    pred = predict(model, new, interval = "prediction"); # predictions
    maxf_fit = pred[,1];

    if(actPlot){
      image.plot(x,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(t2, maxf2, col='blue',pch=1)
      points(r$t, maxf_fit, col='black',pch=1)
      #arrows(r$t, pred[,2], r$t, pred[,3], code=3, angle=90,
      #       length=0.05, col="gray",pch=3);
      
      leg <- c("Maxima", "Estimation ")#,"Uncertainty");
      col <- c("blue","black")#,"gray");
      legend("topleft",legend=leg,cex=.8,col=col,pch=c(1,1));
    }
  }
  
  return(maxf_fit);
}



########################################################################
findGmodes_LASSO = function(r, lambda =1, 
                            setStart = FALSE, m_L = 2, initfreq_L = c(0, Inf),
                            actPlot = FALSE){
  ########################################################################
  
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes
  
  # r : output from specPdgrm
  # lambda : penalization parameter
  
  x = r$t; y = r$f; z = r$E; # values from spectrogram
  
  maxf = y[apply(z, 1, which.max)]; # maximum frequency at each time
  
  x = x[maxf>0];
  maxf = maxf[maxf>0];
  
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
      image.plot(x,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(x, maxf, col='blue',type="p")
      points(x, maxf_fit, col='black',type="p")
      points(x0, y0, col = "green", pch = 19); # blue point in the plot
    }
    
    # second iteration
    indices = abs(maxf-maxf_fit)<250;
    maxf2 = maxf[indices]-y0;
    t2 = x[indices]-x0;
    
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
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(t2+x0, maxf2+y0, col='blue',type="p")
      points(r$t, maxf_fit, col='black',type="p")
      points(x0, y0, col = "green", pch = 19); # blue point in the plot
    }
  }
  
  else{
    Xmat <- cbind(x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10);
    lasso_model = glmnet(Xmat, maxf, alpha=1, lambda=lambda);
    maxf_fit = predict(lasso_model, Xmat);
    
    if(actPlot){
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(x, maxf, col='blue',type="p")
      points(x, maxf_fit, col='black',type="p")
    }
    
    # second iteration
    indices = abs(maxf-maxf_fit)<250;
    maxf2 = maxf[indices];
    t2 = x[indices];
    
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
      image.plot(r$t,y,z,xlab="Time [s]",ylab="Frequency [Hz]")
      points(t2, maxf2, col='blue',pch=1)
      points(r$t, maxf_fit, col='black',pch=1)
      #arrows(r$t, pred[,2], r$t, pred[,3], code=3, angle=90,
      #       length=0.05, col="gray",pch=3);
      
      leg <- c("Maxima", "Estimation ")#,"Uncertainty");
      col <- c("blue","black")#,"gray");
      legend("topleft",legend=leg,cex=.8,col=col,pch=c(1,1));
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
covpbb = function(r, mod, l=200, p=90, fs=16384, movGmode = 11, 
                  um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(-Inf, Inf),
                  um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(-Inf, Inf),
                  gmode = "right", movBand = 5, timeGmode = NULL, 
                  true_data, actPlot = FALSE, limFreq = NULL){
  ########################################################################  
  
  # r        : output from specPdgrm
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
  # timeGMode: time interval to define g-modes
  # thruth_data: true ratio and time
  # thruth_data: simulated ratio time evolution where ratio is M/R^2 (g2 mode) or (g3 mode)
  # actPlot  : logical value to produce plot
  # limFreq  : specifies upper threshold (in Hz) for the estimated g-modes
  
  ### ###
  true_time = true_data$time; 
  true_ratios = true_data$ratio
  
  nr = length(r$t); # number of time samples in the spectrogram
  duration = r$t[nr]-r$t[1]; # total duration (in s) of the analysed data
  n = duration*fs+1; # number of observations

  # starting & ending points of the intervals used in the spectrogram
  index = ints(n=n, l=l, p=p); # from psplinePsd
  
  # centered point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  
  # Adjusting data
  timedata = seq(from=r$t[1], to=r$t[nr], length = n);
  
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
    
    r$t = r$t[out];
    r$E = r$E[out, ]; # row=time, col=freq
  }
  
  maxfs = findGmodes(r, um_R, dm_R, m_R, um_L, dm_L, m_L, 
                     initfreq_R, initfreq_L, testSlope = FALSE);
  
  maxfs = sapply(maxfs,function(x)movf(x, n=movGmode, mean));# smoothing g-mode estimates
  maxfs = as.data.frame(maxfs);
  
  if(actPlot == TRUE){
    image.plot(r$t,r$f,r$E,xlab="Time [s]",ylab="Frequency [Hz]")
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
      else { # At least 3 g-modes (in maxf) are required to generate the cvvpbb band
        maxf1     = maxf[discFreq];     # discarding frequencies according to limFreq
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


########################################################################
covpbb2 = function(r, mod, l=200, p=90, fs=16384, movGmode = 11, 
                  um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(-Inf, Inf),
                  um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(-Inf, Inf),
                  gmode = "right", movBand = 5, timeGmode = NULL, 
                  true_data, actPlot = FALSE, limFreq = NULL){
  ########################################################################  
  
  # r     : dataset as matrix (with no zeros [specPdgrm])
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
  # timeGMode: time interval to define g-modes
  # thruth_data: true ratio and time
  # thruth_data: simulated ratio time evolution where ratio is M/R^2 (g2 mode) or (g3 mode)
  # actPlot  : logical value to produce plot
  # limFreq  : specifies upper threshold (in Hz) for the estimated g-modes
  
  ### ###
  true_time = true_data$time; 
  true_ratios = true_data$ratio
  
  #nr = length(r$t); # number of time samples in the spectrogram
  #duration = r$t[nr]-r$t[1]; # total duration (in s) of the analysed data
  #n = duration*fs+l; # number of observations
  #print(n)
  
  # starting & ending points of the intervals used in the spectrogram
  #index = ints(n=n, l=l, p=p); # from psplinePsd
  
  # centred point for even "l"
  #mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  
  # Ajusting data
  #timedata = seq(from=r$t[1], to=r$t[nr], length = n);
  
  # mean time (centred point) for our g-mode estimates
  #timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );
  
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
    
    #maxf     = maxf[out];
    timefreq = timefreq[out];
    
    r$t = r$t[out];
    r$E = r$E[out, ]; # row=time, col=freq
  }
  
  maxfs = findGmodes2(r, um_R, dm_R, m_R, um_L, dm_L, m_L,
                     initfreq_R, initfreq_L, testSlope = FALSE);
  

  maxfs = sapply(maxfs,function(x)movf(x, n=movGmode, mean));# smoothing g-mode estimates
  maxfs = as.data.frame(maxfs);
  
  if(actPlot == TRUE){
    image.plot(r$t,r$f,r$E,xlab="Time [s]",ylab="Frequency [Hz]")
    if(gmode == "left"){
      points(timefreq, maxfs$maxf_L, col='black',type="p")
    }else if (gmode == "right"){
      points(timefreq, maxfs$maxf_R, col='black',type="p")
    }else{
      points(timefreq, maxfs$maxf_med, col='black',type="p")}
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
      #discFreq = (maxf < j); # positions to keep 
      discFreq = rep(FALSE,length(maxf));   # positions to keep
      for (k in 1:length(maxf)){
        if(maxf[k]>j){break()}
        discFreq[k]=TRUE;
      }
      
      sfq = sum(discFreq);
      
      if(sfq <= 2){ 
        if(sfq == 0){
          warning(paste("All frequencies are greater than limFreq", j));
          print("All frequencies are greater than limFreq")
        }else{
          warning(paste("Only", sfq, "frequency is lower than limFreq", j));
        }
        if(any(class(mod) == "lm")){
          out1 = rbind(out1, c(0, -1));
          out2 = rbind(out2, c(-1, -1, -1));
          out3 = rbind(out3, c(-1));
        }
        else {
          out1 = rbind(out1, c(0, -2));
          out2 = rbind(out2, c(-2, -2, -2));
          out3 = rbind(out3, c(-2));
        }                
      }
      else { # At least 3 g-modes (in maxf) are required to generate the cvvpbb band
        maxf1     = maxf[discFreq];     # discarding frequencies according to limFreq
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



########################################################################
covpbb_poly = function(r, N=3, mod, setStart = FALSE, m_L = 8, initfreq_L = c(-Inf, Inf),
                   movBand = 5, timeGmode = NULL, true_data, limFreq = NULL,
                   actPlot = FALSE){
  ########################################################################  
  
  # r     : dataset as matrix (with no zeros [specPdgrm])
  # N     : order of the polynomial fit
  # mod   : model

  # Paremeters for exploration starting from the left side
  # m_L  : number of intervals used to calculated starting value
  # initfreq_L: interval to restrict frequency for the starting value
  
  # movBand  : define the number of points to smooth the band
  # timeGMode: time interval to define g-modes
  # thrue_data: simulated ratio time evolution where ratio is M/R^2 (g2 mode) or (g3 mode)
  # actPlot  : logical value to produce plot
  # limFreq  : specifies upper threshold (in Hz) for the estimated g-modes
  
  ### ###
  true_time = true_data$time; 
  true_ratios = true_data$ratio
  
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
  
  
  #maxfs = findGmodes_poly(r, N=N, setStart=setStart, m_L=m_L, initfreq_L=initfreq_L, 
  #                        actPlot=FALSE);
  maxfs = findGmodes_LASSO(r, lambda=1, setStart=setStart, m_L=m_L, initfreq_L=initfreq_L, 
                          actPlot=FALSE);
  
  if(actPlot == TRUE){
    image.plot(r$t,r$f,r$E,xlab="Time [s]",ylab="Frequency [Hz]")
    points(timefreq, maxfs, col='black',type="p")
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
        #print("All frequencies are greater than limFreq")
      }else{
        #warning(paste("Only", sfq, "frequency is lower than limFreq", j));
      }
      if(any(class(mod) == "lm")){
        out1 = rbind(out1, c(0, -1));
        out2 = rbind(out2, c(-1, -1, -1));
        out3 = rbind(out3, c(-1));
      }
      else {
        out1 = rbind(out1, c(0, -2));
        out2 = rbind(out2, c(-2, -2, -2));
        out3 = rbind(out3, c(-2));
      }                
    }
    else { # At least 3 g-modes (in maxf) are required to generate the cvvpbb band
      maxf1     = maxf[discFreq];     # discarding frequencies according to limFreq
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
  
  R = list(covpbb = out1, residual = out2);
  
  return(R);
}


########################################################################
thresh = function(TFMap,threshold=0.05,actPlot=FALSE, main=NULL){
  ######################################################################
  # Inputs :  TFMap: one time-frequency map obtained with timefreqMap
  #           threshold: percentage of pixels to keep in thresholding
  #
  # Output :  TFMap with only the x% more energetic pixel
  ######################################################################
  
  t=TFMap$t;
  f=TFMap$f;
  nbPix=length(t)*length(f);
  selectPix=threshold*nbPix;
  peaks=tail(sort(TFMap$E),selectPix);
  ind=which(TFMap$E %in% peaks);
  
  ### Result ###
  selectMap=zeros(length(t),length(f));
  selectMap[ind]=TFMap$E[ind];
  
  if (actPlot){
    image.plot(t,f,selectMap,xlab="Time [s]",ylab="Frequency [Hz]", main=main)
  }
  return(list(t=t,f=f,E=selectMap))
}


########################################################################
thresh_abs = function(TFMap,limit=1.93,actPlot=FALSE, main=NULL){
  ######################################################################
  # Inputs :  TFMap: one time-frequency map obtained with timefreqMap
  #           threshold: percentage of pixels to keep in thresholding
  #
  # Output :  TFMap with only the x% more energetic pixel
  ######################################################################
  
  t=TFMap$t;
  f=TFMap$f;
  ind=which(TFMap$E >= limit);
  
  ### Result ###
  selectMap=zeros(length(t),length(f));
  selectMap[ind]=TFMap$E[ind];
  
  if (actPlot){
    image.plot(t,f,selectMap,xlab="Time [s]",ylab="Frequency [Hz]", main=main)
  }
  return(list(t=t,f=f,E=selectMap))
}