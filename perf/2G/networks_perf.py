import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage.filters import gaussian_filter1d

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
    return dist,q1,q2,q3

def find_chardist(dist,q2, threshold, type):
    x=-1
    hs=(dist[1]-dist[0])/2
    dist_nb=len(q2)
    if type == 1:
        test=0
        for i in range(dist_nb):
            if q2[i]>threshold:
                test=1
            if (q2[i]<threshold and test==1) :
                x=dist[i]-hs
                test=0
    else:
        test=0
        for i in range(dist_nb):
            if q2[i]<threshold:
                test=1
            if (q2[i]>threshold and test==1) :
                x=dist[i]-hs
                test=0
        
    return x

def smooth(y, sigma):
    n_rows=len(y)
    n_cols=len(y[0])
    y_smooth=np.zeros((n_rows,n_cols))    
    for i in range(n_cols):
        y_smooth[:,i] = gaussian_filter1d(y[:,i], sigma=sigma)
    return y_smooth

# 0 dist
# 1 covpbb                
# 2 medbandwidth
# 3 absolute residual (mean)
# 4 MSE (mean)
# 5 precision (mean)

folder = 'favourable/'

signals = ["s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", 
            "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220"]
signal_names = ["s11", "s15", "s15S", "s15G", "s20", "s20S",
                "s25", "s40", "no signal"]

# signals = ["s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", 
#             "s20.0--LS220", "s20.0--SFHo"]
# signal_names = ["s11", "s15", "s15S", "s15G", "s20", "s20S", "no signal"]

filt = "spectrum"

dist_nb = 60
dist_max = np.zeros(dist_nb)

networks = ["HL", "HLVKA"]
markers = ['--', '-']
ind_net = 0

plt.figure()
for net in networks:
    
    inputfolder = folder + net
    
    sig_nb=np.size(signal_names)
    qq1=np.zeros((dist_nb,sig_nb))
    qq2=np.zeros((dist_nb,sig_nb))
    qq3=np.zeros((dist_nb,sig_nb))
    dd=np.zeros((dist_nb,sig_nb))
    
    quantity="coverage"
    #quantity="delta"
    
    if quantity == "coverage":
        index=1
    else:
        index=5
    
# =============================================================================
#     if quantity == "coverage":
#         threshold=q95
#     else:
#         threshold=q5
# =============================================================================
    
    ind=0
    for sig in signals:
        filename= inputfolder + '/results_AA_' + filt + '_f2_' + sig + '.txt'
        
        a= np.loadtxt(filename, dtype='f', delimiter=' ')
        dist=np.unique(a[:,0]) 
        dist_nb=len(dist)
        
        dist,q1,q2,q3=buildbox(a,index)
        qq1[:,ind]=q1
        qq2[:,ind]=q2
        qq3[:,ind]=q3
        dd[:,ind]=dist
        
        # dist_max is the list that contains the largest distance
        if max(dist)>max(dist_max):
            dist_max = dist
            
        ind=ind+1
    
        
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.turbo(np.linspace(0, 1, sig_nb))))
    lineObjects = plt.plot(dd, smooth(qq2,3), markers[ind_net])
    #plt.plot(dd[:,0:sig_nb-1], qq2[:,0:sig_nb-1], marker="+",linestyle="")
        
    ind_net += 1


# Noise results    
filename= folder + 'HLVKA/results_AA_' + filt + '_f2_noise.txt'
a= np.loadtxt(filename, dtype='f', delimiter=' ')

# percentile
q50=np.percentile(a[:,index], 50)
q95=np.percentile(a[:,index], 95)
q5=np.percentile(a[:,index], 5)

qq1[:,sig_nb-1]=q5
qq2[:,sig_nb-1]=q50
qq3[:,sig_nb-1]=q95
dd[:,sig_nb-1]=0

# blue-ish area corresponding to the noise reponse for the latter network
plt.fill_between(dist_max, qq1[:,sig_nb-1], qq3[:,sig_nb-1], alpha=0.05, facecolor='b')


if quantity == "coverage":
    plt.ylim((0,1.02))
    plt.xlim((1,max(dist_max)))
    #plt.xlim((1,60))
    plt.ylabel('Coverage')
    plt.legend(iter(lineObjects), signal_names, loc='lower right')
else:
    plt.ylim((0.05,1.05*qq2[0,sig_nb-1]))
    plt.xlim((1,max(dist_max)))
    #plt.xlim((1,50))
    plt.ylabel('$\Delta$')
    plt.legend(iter(lineObjects), signal_names, loc='upper right')

plt.xlabel('Distance [kpc]')

plt.grid(True)

fig_name = 'HLvsHLVKA_favourable.png';
plt.savefig(fig_name)
fig_name = 'HLvsHLVKA_favourable.pdf';
plt.savefig(fig_name)

plt.show()    

