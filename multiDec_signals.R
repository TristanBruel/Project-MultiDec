### Earth Model WGS-84 params ###

WGS84 = list(a=6378137, b=6356752.314)


### Wave Propagation Frame to Earth Fixed Frame ###

waveProp_to_Earth = function(psi, theta, phi){
  # Input : Sky Angles psi, theta, phi
  # Output : Rotation matrix to fixed earth frame
  Rpsi = matrix(c(cos(psi),-sin(psi),0,sin(psi),cos(psi),0,0,0,1),ncol=3)
  Rtheta = matrix(c(1,0,0,0,-cos(theta),-sin(theta),0,sin(theta),-cos(theta)),ncol=3)
  Rphi = matrix(c(sin(phi),cos(phi),0,-cos(phi),sin(phi),0,0,0,1),ncol=3)
  Rsky = Rpsi%*%(Rtheta%*%Rphi)
  return(Rsky)
}


### Gravitational Wave Tensor in the Earth Fixed Frame ###

grav_tensor = function(psi, theta, phi){
  # Inputs : Sky angles
  # Output : Gravitational wave tensor
  Rsky = waveProp_to_Earth(psi, theta, phi)
  e_wx = Rsky[1,]
  e_wy = Rsky[2,]
  e_wz = Rsky[3,]
  
  e_plus = e_wx%*%t(e_wx) - e_wy%*%t(e_wy)
  e_cross = e_wx%*%t(e_wy) + e_wy%*%t(e_wx)
  
  # h_tensor = h_plus*e_plus + h_cross*e_cross
  
  # Fplus = (e_plus[1,1]-e_plus[2,2])/2
  # Fcross = (e_cross[1,1]-e_cross[2,2])/2
  
  return(list(eplus=e_plus, ecross=e_cross))
}


### Location and Orientation for Detectors in the earth fixed frame ###

location_det = function(l, lambda, h){
  # Inputs : earth model WGS-84 coordinates
  # Output : earth fixed frame coordinates
  a = WGS84$a
  b = WGS84$b
  R = a^2/(sqrt(a^2*cos(l)^2+b^2*sin(l)^2)) #local radius of curvature
  x_e = (R+h)*cos(l)*cos(lambda)
  y_e = (R+h)*cos(l)*sin(lambda)
  z_e = (b^2*R/a^2+h)*sin(l)
  return(c(x_e, y_e, z_e))
}

orientation_arm = function(l, lambda, Psi, omega){
  # Inputs : orientation angles of one arm in the detector frame
  # Outputs : unit orientation vectors in the earth fixed frame
  o1 = -cos(Psi)*cos(omega)*sin(lambda)-sin(Psi)*cos(omega)*cos(lambda)*sin(l)+sin(omega)*cos(l)*cos(lambda)
  o2 = cos(Psi)*cos(omega)*cos(lambda)-sin(Psi)*cos(omega)*sin(lambda)*sin(l)+sin(omega)*cos(l)*sin(lambda)
  o3 = sin(Psi)*cos(omega)*cos(l)+sin(omega)*sin(l)
  return(c(o1,o2,o3))
}

# response_tensor = 0.5*(Omega1%*%t(Omega1)-Omega2%*%t(Omega2))

earthFrame_to_detector = function(l, lambda, h, Psi1, Psi2, omega1, omega2){
  loc = location_det(l, lambda, h)
  arm1 = orientation_arm(l, lambda, Psi1, omega1)
  arm2 = orientation_arm(l, lambda, Psi2, omega2)
  return(list(loc=loc,vectors_arm1=arm1,vectors_arm2=arm2))
}

# Detectors parameters  
det = read.csv("detectors_params.csv", sep=";", stringsAsFactors=FALSE, header=TRUE, dec=",")


### Response for gravitational wave ###
#######################################

grav_response = function(dec, ra, t, pol=0, detector){
  # Inputs : - sky position of the source
  #          - time GMST of emission
  #          - polarization angle (null by default)
  #          - detector name ('LHO', 'LLO' or 'VIRGO')
  #
  # Output : time of arrival at a given detector
  
  if ((detector != "LHO") && (detector != "LLO") && (detector != "VIRGO")){
    stop("detector must be LHO, LLO or VIRGO")
  }
  
  index=which(det$name == detector)
  params = det[index,]
  h = params$height
  l = params$North.latitude
  lambda = params$East.longitude
  Psi1 = params$arm1.polarization
  Psi2 = params$arm2.polarization
  omega1 = params$arm1.angle.to.horizontal
  omega2 = params$arm2.angle.to.horizontal
  detector_frame = earthFrame_to_detector(l, lambda, h, Psi1, Psi2, omega1, omega2)
  
  arm_vec1 = detector_frame$vectors_arm1
  arm_vec2 = detector_frame$vectors_arm2
  
  #Sky angles
  psi = pol
  phi = ra - t
  theta = pi/2 - dec
  
  grav_tsr = grav_tensor(psi,theta,phi)
  eplus = grav_tsr$eplus
  ecross = grav_tsr$ecross
  strain_plus = eplus%*%t(arm_vec1%*%t(arm_vec1)-arm_vec2%*%t(arm_vec2))
  strain_cross = ecross%*%t(arm_vec1%*%t(arm_vec1)-arm_vec2%*%t(arm_vec2))
  
  Fplus = sum(diag(strain_plus))/2
  Fcross = sum(diag(strain_cross))/2
  
  return(list(Fplus=Fplus, Fcross=Fcross))
}


###################
### Time Delays ###
###################

time_delay = function(dec, ra, t, detector){
  # Inputs : - sky position of the source
  #          - time GMST at which the wave arrives at the detector 
  #               (higher order corrections lost here)
  #          - detector name ('LHO', 'LLO' or 'VIRGO')
  #
  # Outputs : - (algebraic) time taken for the wave to reach the center of the earth
  
  if ((detector != "LHO") && (detector != "LLO") && (detector != "VIRGO")){
    stop("detector must be LHO, LLO or VIRGO")
  }
    
  index=which(det$name == detector)
  params = det[index,]
  h = params$height
  l = params$North.latitude
  lambda = params$East.longitude
  Psi1 = params$arm1.polarization
  Psi2 = params$arm2.polarization
  omega1 = params$arm1.angle.to.horizontal
  omega2 = params$arm2.angle.to.horizontal
  detector_frame = earthFrame_to_detector(l, lambda, h, Psi1, Psi2, omega1, omega2)
  
  loc = detector_frame$loc
  
  theta = pi/2 - dec
  phi = ra - t
  # source position unit vector
  e_sky = matrix(c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)))
  
  # projection of the detector location vector onto this unit vector
  dist = -loc%*%e_sky
  return(dist/3e8)
}


### Multi detection simulation ###
##################################


## TO BE UPDATED ##
signal_multiDec = function(dec, ra, t, signal="KURODA_TM1_H_resampled.dat"){
  # Inputs : - sky position of the source
  #          - time GMST at which the wave arrives at the center of the earth
  #          - name of the (simulated) waveform
  #
  # Outputs : detection time series for the 3 detectors LHO, LLO and VIRGO
  
  folder = "Waveforms/"
  
  gw_filename=paste(folder,signal,sep="")
  sXX = read.table(gw_filename); # V1 time, V2 hplus, V3 hcross
  colnames(sXX) = c ("time","hplus","hcross");
  n = length(sXX$time)
  
  coeff_LHO = grav_response(dec,ra,t,0,"LHO")
  coeff_LLO = grav_response(dec,ra,t,0,"LLO")
  coeff_VIRGO = grav_response(dec,ra,t,0,"VIRGO")
  
  signal_LHO = rep(0,n)
  signal_LLO = rep(0,n)
  signal_VIRGO = rep(0,n)
  
  for (i in 1:n){
    signal_LHO[i]=coeff_LHO$Fplus*sXX$hplus[i] + coeff_LHO$Fcross*sXX$hcross[i]
    signal_LLO[i]=coeff_LLO$Fplus*sXX$hplus[i] + coeff_LLO$Fcross*sXX$hcross[i]
    signal_VIRGO[i]=coeff_VIRGO$Fplus*sXX$hplus[i] + coeff_VIRGO$Fcross*sXX$hcross[i]
  }
  
  plot(sXX$time, signal_LHO, type='l', xlab = "Time [s]", ylab = "Hoft_LHO") 
  plot(sXX$time, signal_LLO, type='l', xlab = "Time [s]", ylab = "Hoft_LLO") 
  plot(sXX$time, signal_VIRGO, type='l', xlab = "Time [s]", ylab = "Hoft_VIRGO") 
  
  # at this point the 3 signals are artificially synchronized (the gw arrives at the 3 detectors simultaneously)
  # no noise added yet
  
  # return(list(s_LHO=signal_LHO,s_LLO=signal_LLO,s_VIRGO=signal_VIRGO))
}