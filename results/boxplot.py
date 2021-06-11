import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os
import pylab

def minmax(a, ind_mean, ind_var):
    N,Y=a.shape
    
    out=np.zeros((N,2))    
    for i in range(N):
        out[i,0]=max(0, a[i,ind_mean]-a[i,ind_var])
        out[i,1]=min(1, a[i,ind_mean]+a[i,ind_var])
    return out

def buildbox(a,index):
    dist=np.unique(a[:,0]) 
    dist_nb=len(dist)
    #print('Number of distances: ', dist_nb)
    data = []
    box_name = []
    q1 = []
    q2 = []
    q3 = []
    for i in range(dist_nb):
        ext=np.where((a[:,0] == dist[i]))
        data_boxi = a[ext[0][0]:ext[0][-1],index]
        data.append(data_boxi)
        if (i%5 == 0):
            box_name.append(str(dist[i]))
        else:
            box_name.append('')
        q1.append(np.percentile(data_boxi, 25))
        q2.append(np.percentile(data_boxi, 50))
        q3.append(np.percentile(data_boxi, 75))
    
    return data,box_name,q1,q2,q3    


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
detectors=["multiDec"]
#method="singleDec"
method="multiDec"

try:
    os.stat(name)
except:
    os.mkdir(name)       

for det in detectors:
    title= name + ' ' + det
    filename=method+'/results_AA_prewhiten_f2_' + name + '_' + det + '.txt'

    a= np.loadtxt(filename, dtype='f', delimiter=' ')
    dist=np.unique(a[:,0]) 
    dist_nb=len(dist)
    print('Number of distances: ', dist_nb)
    
    # coverage
    data,box_name,q1,q2,q3=buildbox(a,1)

    #plt.figure(figsize=(4, 3))
    plt.figure()
    blue_cross = dict(markerfacecolor='b', marker='+')
    bx=plt.boxplot(data,flierprops=blue_cross)

    pylab.xticks(range(1,dist_nb+1), box_name)
    
    plt.xlabel('Distance [kpc]')
    plt.ylabel('Coverage probability')
    plt.title(title)
    plt.grid(False)
    plt.xlim((1,30))
    plt.ylim((0,1.05))
    fig_name=name + '/' + name + '_covpbb_boxplot_' + det + '.png';
    plt.savefig(fig_name)
    fig_name=name + '/' + name + '_covpbb_boxplot_' + det + '.eps';
    plt.savefig(fig_name)
    fig_name=name + '/' + name + '_covpbb_boxplot_' + det + '.pdf';
    plt.savefig(fig_name)
    
    
    #plt.figure(figsize=(4, 3))
    #plt.figure()
    #plt.plot(dist, q2, 'b', label='Median')
    #plt.fill_between(dist, q1, q3, alpha=0.05, facecolor='b', label="IQR")
    
    #plt.xlabel('Distance [kpc]')
    #plt.ylabel('coverage')
    #plt.title(title)
    #plt.grid(True)
    #plt.legend(loc='best')

    # delta
    data,box_name,q1,q2,q3=buildbox(a,5)
    
    #plt.figure(figsize=(4, 3))
    plt.figure()
    blue_cross = dict(markerfacecolor='b', marker='+')
    bx=plt.boxplot(data, flierprops=blue_cross)
    pylab.xticks(range(1,dist_nb+1), box_name)
    
    plt.xlabel('Distance [kpc]')
    plt.ylabel('$\Delta$')
    plt.title(title)
    plt.grid(False)
    plt.ylim((0,3.5))
    plt.xlim((1,30))
    fig_name=name + '/' + name + '_error_boxplot_' + det + '.png';
    plt.savefig(fig_name)
    fig_name=name + '/' + name + '_error_boxplot_' + det + '.eps';
    plt.savefig(fig_name)
    fig_name=name + '/' + name + '_error_boxplot_' + det + '.pdf';
    plt.savefig(fig_name)
    
    
    #plt.figure(figsize=(4, 3))
    #plt.figure()
    #plt.plot(dist, q2, 'b', label='Median')
    #plt.fill_between(dist, q1, q3, alpha=0.05, facecolor='b', label="IQR")
    
    #plt.xlabel('Distance [kpc]')
    #plt.ylabel('$\Delta$')
    #plt.title(title)
    #plt.grid(True)
    #plt.legend(loc='best')

    plt.show()

