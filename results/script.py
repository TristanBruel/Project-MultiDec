import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os

def minmax(a, ind_mean, ind_var):
    N,Y=a.shape
    
    out=np.zeros((N,2))    
    for i in range(N):
        out[i,0]=max(0, a[i,ind_mean]-a[i,ind_var])
        out[i,1]=min(1, a[i,ind_mean]+a[i,ind_var])
    return out

# 0 dist
# 1 covpbb                
# 2 medbandwidth
# 3 absolute residual (mean)
# 4 MSE (mean)
# 5 precision (mean)

#name="s11.2--LS220"
#name="s15.0--GShen"
#name="s15.0--LS220"
#name="s15.0--SFHo"
name="s20.0--LS220"
#name="s25.0--LS220"
#name="s40.0--LS220"
#detectors=["LHO","LLO","VIR"]
detectors=['multiDec']
#method="singleDec"
method="multiDec"

for det in detectors:
    title=name + ' ' + det
    filename=method+'/results_AA_prewhiten_f2_' + name + '_' + det + '.txt'
    
    a= np.loadtxt(filename, dtype='f', delimiter=' ')
    dist=np.unique(a[:,0]) 
    dist_nb=len(dist)
    print('Number of distances: ', dist_nb)
    
    try:
        os.stat(name)
    except:
        os.mkdir(name)       
    
    covpbb=np.zeros((dist_nb,6),dtype=float)
    medbw=np.zeros((dist_nb,6),dtype=float)
    MSE=np.zeros((dist_nb,6),dtype=float)
    prec=np.zeros((dist_nb,6),dtype=float)
    
    col=['b', 'r']
    colmap=['viridis', 'jet']
    mar=['o', 'x']
    
    
    for i in range(dist_nb):
        ext=np.where((a[:,0] == dist[i]))
    
        covpbb[i,0]=np.median(a[ext,1])
        covpbb[i,1]=np.mean(a[ext,1])
        covpbb[i,2]=np.median(abs(a[ext,1]-covpbb[i,0]))
        covpbb[i,3]=np.std(a[ext,1])
        covpbb[i,4]=np.percentile(a[ext,1],95)
        covpbb[i,5]=1.96*np.std(a[ext,1])/np.sqrt(len(ext))
    
        medbw[i,0]=np.median(a[ext,2])
        medbw[i,1]=np.mean(a[ext,2])
        medbw[i,2]=np.median(abs(a[ext,2]-medbw[i,0]))
        medbw[i,3]=np.std(a[ext,2])
        medbw[i,4]=np.percentile(a[ext,2],95)
        medbw[i,5]=1.96*np.std(a[ext,2])/np.sqrt(len(ext))
    
        MSE[i,0]=np.median(a[ext,4])
        MSE[i,1]=np.mean(a[ext,4])
        MSE[i,2]=np.median(abs(a[ext,4]-MSE[i,0]))
        MSE[i,3]=np.std(a[ext,4])
        MSE[i,4]=np.percentile(a[ext,4],95)
        MSE[i,5]=1.96*np.std(a[ext,4])/np.sqrt(len(ext))
    
        prec[i,0]=np.median(a[ext,5])
        prec[i,1]=np.mean(a[ext,5])
        prec[i,2]=np.median(abs(a[ext,5]-prec[i,0]))
        prec[i,3]=np.std(a[ext,5])
        prec[i,4]=np.percentile(a[ext,5],95)
        prec[i,5]=1.96*np.std(a[ext,5])/np.sqrt(len(ext))
    
    
        
    #plt.figure(figsize=(4, 3))
    plt.figure()
    
    plt.plot(dist, covpbb[:,0], 'b', label='Median +/- MAD')
    a=minmax(covpbb,0,2)
    plt.plot(dist,a,'--b')
    plt.fill_between(dist, a[:,0], a[:,1], alpha=0.05, facecolor='b')
    
    #plt.plot(dist, covpbb[:,1], 'r', label='Mean +/- std')
    #a=minmax(covpbb,1,3)
    #plt.plot(dist,a,'--r')
    #plt.fill_between(dist, a[:,0], a[:,1], alpha=0.05, facecolor='r')
    
    #plt.plot(dist, covpbb[:,1], 'g', label='Mean +/- 95%CI')
    #a=minmax(covpbb,1,5)
    #plt.plot(dist,a,'--g')
    #plt.fill_between(dist, a[:,0], a[:,1], alpha=0.05, facecolor='g')
    
    #plt.plot(dist, prec[:,0], 'b')
    #a=minmax(prec,0,2)
    #plt.plot(dist,a,'--b')
    #plt.fill_between(dist, a[:,0], a[:,1], alpha=0.05, facecolor='b')
    
    #plt.plot(dist, prec[:,1], 'r')
    #a=minmax(prec,1,3)
    #plt.plot(dist,a,'--r')
    #plt.fill_between(dist, a[:,0], a[:,1], alpha=0.05, facecolor='r')
    
    #plt.plot(dist, prec[:,1], 'g')
    #a=minmax(prec,1,5)
    #plt.plot(dist,a,'--g')
    #plt.fill_between(dist, a[:,0], a[:,1], alpha=0.05, facecolor='g')
    
    
    plt.xlabel('Distance [kpc]')
    plt.ylabel('Coverage probability')
    plt.ylim([0, 1.05])
    plt.title(title)
    plt.grid(True)
    plt.legend(loc='best')
    #fig_name=sig + '/' + sig + '_covpbb_prec' + det + '_b200.png';
    #plt.savefig(fig_name)
    
    plt.show()
