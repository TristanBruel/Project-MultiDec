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

orientation_arm = function(l, lambda, az, polar){
  # Inputs : orientation angles of one arm in the detector frame
  # Outputs : unit orientation vectors in the earth fixed frame
  o1 = -cos(az)*cos(polar)*sin(lambda)-sin(az)*cos(polar)*cos(lambda)*sin(l)+sin(polar)*cos(l)*cos(lambda)
  o2 = cos(az)*cos(polar)*cos(lambda)-sin(az)*cos(polar)*sin(lambda)*sin(l)+sin(polar)*cos(l)*sin(lambda)
  o3 = sin(az)*cos(polar)*cos(l)+sin(polar)*sin(l)
  return(c(o1,o2,o3))
}

# response_tensor = 0.5*(Omega1%*%t(Omega1)-Omega2%*%t(Omega2))

earthFrame_to_detector = function(l, lambda, h, az1, az2, polar1, polar2){
  loc = location_det(l, lambda, h)
  arm1 = orientation_arm(l, lambda, az1, polar1)
  arm2 = orientation_arm(l, lambda, az2, polar2)
  return(list(loc=loc,vectors_arm1=arm1,vectors_arm2=arm2))
}

# Detectors parameters  
dets = read.csv("detectors_params.csv", sep=";", stringsAsFactors=FALSE, header=TRUE)


### Response for gravitational wave ###
#######################################

grav_response = function(dec, ra, t, pol=0, detector){
  # Inputs : - sky position of the source
  #               declination in ° and right ascension in hours
  #          - time GPS of arrival at the detector
  #          - polarization angle (null by default)
  #          - detector name ('LHO', 'LLO' or 'VIRGO')
  #
  # Output : time of arrival at a given detector
  
  if (!(detector %in% c("LHO","LLO","VIRGO"))){
    stop("detector must be LHO, LLO or VIRGO")
  }
  
  index=which(dets$name == detector)
  params = dets[index,]
  h = params$height
  l = params$North.latitude
  lambda = params$East.longitude
  az1 = params$arm1.azimuth
  az2 = params$arm2.azimuth
  polar1 = params$arm1.polar
  polar2 = params$arm2.polar
  detector_frame = earthFrame_to_detector(l, lambda, h, az1, az2, polar1, polar2)
  
  arm_vec1 = detector_frame$vectors_arm1
  arm_vec2 = detector_frame$vectors_arm2
  
  #Sky angles
  psi = pol
  phi = ra*(2*pi/24) - GPS_to_GMST(t)
  theta = pi/2 - dec*(2*pi/360)
  
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
  #               declination in ° and right ascension in hours
  #          - time GPS at which the wave arrives at the detector 
  #               (higher order corrections lost here)
  #          - detector name ('LHO', 'LLO' or 'VIRGO')
  #
  # Outputs : - (algebraic) time taken for the wave to reach the center of the earth
  
  if ((detector != "LHO") && (detector != "LLO") && (detector != "VIRGO")){
    stop("detector must be LHO, LLO or VIRGO")
  }
    
  index=which(dets$name == detector)
  params = dets[index,]
  h = params$height
  l = params$North.latitude
  lambda = params$East.longitude
  az1 = params$arm1.azimuth
  az2 = params$arm2.azimuth
  polar1 = params$arm1.polar
  polar2 = params$arm2.polar
  detector_frame = earthFrame_to_detector(l, lambda, h, az1, az2, polar1, polar2)
  
  loc = detector_frame$loc
  
  theta = pi/2 - dec*(2*pi/360)
  phi = ra*(2*pi/24) - GPS_to_GMST(t)
  # source position unit vector
  e_sky = matrix(c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)))
  
  # projection of the detector location vector onto this unit vector
  dist = -loc%*%e_sky
  return(dist[[1]]/299792458)
}


# Function taken from Adrian's work

GPS_to_GMST = function(time){
  # Input : time GPS (in s)
  # Output : time GMST (in rad)
  
  w_E = 2*pi*(1/365.2425+1)/86400 # Earth sidereal rotation speed (rad/s)
  gps2000 = 630720013.0 # GPS time for 01/01/2000 00:00:00 UTC
  s0  = (6.0 + 39.0/60.0 + 51.251406103947375/3600)*pi/12 # Sidereal time at gps2000 (rad)
  
  GMST = (w_E*(time-gps2000)+s0)%%(2*pi)
  
  return(GMST)
}

