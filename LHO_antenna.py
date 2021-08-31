# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
from numpy import cos,sin,pi

def eplus(psi,theta,phi):
    e11 = -2*cos(psi)**2*cos(theta)**2*cos(phi)**2-4*cos(psi)*sin(phi)*sin(psi)*cos(theta)*cos(phi)-2*cos(psi)**2*cos(phi)**2+cos(theta)**2*cos(phi)**2+2*cos(psi)**2+cos(phi)**2-1 
    e13 = sin(theta)*(2*cos(psi)**2*cos(theta)*cos(phi)+2*cos(psi)*sin(phi)*sin(psi)-cos(theta)*cos(phi))
    
    e22 = 2*cos(psi)**2*cos(theta)**2*cos(phi)**2+4*cos(psi)*sin(phi)*sin(psi)*cos(theta)*cos(phi)-2*cos(psi)**2*cos(theta)**2+2*cos(psi)**2*cos(phi)**2-cos(theta)**2*cos(phi)**2+cos(theta)**2-cos(phi)**2
    
    e31 = sin(theta)*(2*cos(psi)**2*cos(theta)*cos(phi)+2*cos(psi)*sin(phi)*sin(psi)-cos(theta)*cos(phi))
    return(e11,e13,e22,e31)
    
def ecross(psi,theta,phi):
    e13 = sin(theta)*(-2*cos(psi)*sin(psi)*cos(theta)*cos(phi)+2*cos(psi)**2*sin(phi)-sin(phi))
    
    e33 = 2*sin(psi)*sin(theta)**2*cos(psi)
    
    return(e13,e33)
    
def Fplus(psi,theta,phi):
    res = -cos(theta)**2*cos(2*psi)*cos(2*phi)/2-cos(2*psi)*cos(2*phi)/2-cos(theta)*sin(2*psi)*sin(2*phi)
    return(res)
    
def Fcross(psi,theta,phi):
    res = cos(theta)**2*sin(2*psi)*cos(2*phi)/2+sin(2*psi)*cos(2*phi)/2-cos(theta)*sin(2*phi)*cos(2*psi)
    return(res)
    
def ResponseLHO(psi,dec,ra,GMST):
    alpha = -ra+GMST
    Fp = 0.3104536563*cos(psi)*sin(psi)*sin(dec)\
    -0.1552268269*sin(alpha)*cos(alpha)*cos(dec)**2\
    -0.6209073106*cos(psi)**2*sin(alpha)*cos(alpha)\
    +0.4947780908*cos(dec)*sin(dec)*cos(alpha)\
    -0.9895561806*cos(dec)*sin(dec)*cos(alpha)*cos(psi)**2\
    +0.4559956647*cos(dec)*sin(alpha)*sin(dec)\
    -1.424276333*cos(psi)**2*cos(alpha)**2*cos(dec)**2\
    -1.424276335*cos(psi)**2+0.4928680670*cos(psi)**2*cos(dec)**2\
    +0.9895561807*sin(psi)*cos(dec)*cos(psi)*sin(alpha)\
    -0.6209073102*sin(psi)*sin(dec)*cos(alpha)**2*cos(psi)\
    +0.3104536568*sin(alpha)*cos(alpha)\
    +2.848552664*cos(psi)**2*cos(alpha)**2\
    -2.848552664*cos(psi)*sin(alpha)*sin(psi)*sin(dec)*cos(alpha)\
    +0.7121381649*cos(alpha)**2*cos(dec)**2\
    +0.3104536565*cos(psi)**2*sin(alpha)*cos(alpha)*cos(dec)**2+0.7121381662\
    -0.9119913309*sin(psi)*cos(dec)*cos(psi)*cos(alpha)\
    -1.424276333*cos(alpha)**2-0.2464340348*cos(dec)**2\
    -0.9119913309*cos(dec)*sin(alpha)*sin(dec)*cos(psi)**2
    
    Fc = -0.4947780908*cos(dec)*sin(alpha)\
    -0.3104536569*cos(psi)*sin(alpha)*sin(psi)*cos(alpha)*cos(dec)**2\
    +1.424276334*cos(psi)*sin(psi)\
    +0.4559956650*cos(dec)*cos(alpha)\
    -0.4928680678*cos(psi)*sin(psi)*cos(dec)**2-0.1552268269*sin(dec)\
    +0.3104536564*cos(psi)**2*sin(dec)\
    -2.848552664*cos(psi)**2*sin(alpha)*sin(dec)*cos(alpha)\
    +0.6209073092*cos(psi)*sin(alpha)*sin(psi)*cos(alpha)\
    +1.424276334*cos(psi)*sin(psi)*cos(alpha)**2*cos(dec)**2\
    +1.424276334*sin(dec)*cos(alpha)*sin(alpha)\
    +0.9895561806*cos(psi)*cos(dec)*sin(psi)*sin(dec)*cos(alpha)\
    -0.6209073097*cos(psi)**2*sin(dec)*cos(alpha)**2\
    +0.9119913300*cos(psi)*cos(dec)*sin(psi)*sin(alpha)*sin(dec)\
    +0.3104536564*sin(dec)*cos(alpha)**2\
    -2.848552664*cos(psi)*sin(psi)*cos(alpha)**2\
    +0.9895561807*cos(psi)**2*cos(dec)*sin(alpha)\
    -0.9119913310*cos(psi)**2*cos(dec)*cos(alpha)
    
    return(Fp,Fc)
