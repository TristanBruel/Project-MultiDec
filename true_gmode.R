########################################################################
true_gmode = function(TFmap, true_data, movBand=5, actPlot = FALSE){
  ########################################################################  
  
  # r   : time-frequency map (output of timeFreqMaps.R)
  # mod : model describing the evolution of the ratio with frequency
  
  # movBand   : define the number of points to smooth the band
  # true_data : simulated ratio time evolution where ratio is M/R^2 (g2 mode)
  # actPlot   : logical value to produce plot
  
  ###
  true_time = true_data$time; 
  true_ratios = true_data$ratio;
  
  timefreq = TFmap$t;
  
  maxf = TFmap$f[apply(TFmap$E, 1, which.max)]; # maximum frequency at each time index

  
  # defining true ratios in the band limit given by limfreq 
  discTime = (true_time >= min(timefreq)) & (true_time <= max(timefreq));
  true_time1 = true_time[discTime];  
  true_ratio1 = true_ratios[discTime];

  #########################
  ### Statistical model ###
  #########################
  fits_data = read.table("inputs/1D_simulations/A-A_fits_data_g2.dat", sep = ",");# data to generate model
  colnames(fits_data) = c("r", "f");
  
  # Variable variance linear model
  mu3 = "~ r + I(r^2) + I(r^3) - 1";
  s2 = "~ r + I(r^2) - 1";
  
  Xm  = model.matrix(eval(parse(text=eval(parse(text=mu3)))), fits_data);
  Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);
  
  mod = lmvar(fits_data$f, X_mu = Xm, X_sigma = Xs, intercept_mu = FALSE);
  pred = predict(mod, X_mu = Xm, X_sigma = Xs, 
                         interval = "prediction", sigma = FALSE); # predictions

  r = true_ratio1
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

  ### generating band function ### 

  # interpolating lower bound for predicted values
  fd = approxfun(x = true_time1, y = movf(pred[,2],n=movBand,mean), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  # interpolating upper bound for predicted values
  fu = approxfun(x = true_time1, y = movf(pred[,3],n=movBand,mean), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  # interpolating point estimates
  fm = approxfun(x = true_time1, y = movf(pred[,1],n=movBand,mean), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  
  #universal_relation = function(r,b=6.35e5,c=-9.98e7,d=5.7e9){b*r+c*r^2+d*r^3}
  #true_freq = universal_relation(true_ratios)
  
  image.plot(TFmap$t,TFmap$f,TFmap$E,xlab="Time [s]",ylab="Frequency [Hz]",
             cex.lab=1.8, cex.axis=1.5)
  points(timefreq, maxf, col='black', pch=19)
  points(true_time1, pred[,1], col = "blue", cex = pred[,1]/max(pred[,1])+ 0.3, pch=2);
  arrows(true_time1, pred[,2], true_time1, pred[,3], code=3, angle=90,
         length=0.05, col="gray",pch=3);

  leg <- c("Maxima", "True_gmode ","Uncertainty");
  col <- c("black","blue","gray");
  legend("topleft",legend=leg,cex=1.,col=col,pch=c(1,2,3));
  
}
