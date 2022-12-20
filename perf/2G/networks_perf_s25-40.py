import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter1d

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

folder = 'unfavourable'

signals = ["s25.0--LS220", "s40.0--LS220"]
sig_nb=np.size(signals)
signal_names = ["s25", "s40", "no signal"]

filt = "spectrum"

dist_nb = 61
dist_max = np.zeros(dist_nb)

networks = ["HL", "HLVKA"]
ls = ["--", "-"]
markers = ['+', 'x']

### Plotting parameters ###
# CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
#                   '#f781bf', '#a65628', '#984ea3',
#                   '#999999', '#e41a1c', '#dede00']
CB_color_cycle = ['#377eb8', '#ff7f00']

fs = 20
lw=3

plt.rc('font', family='serif')
plt.rc('lines', linewidth=lw)
custom_cycler = plt.cycler(color=CB_color_cycle)
plt.rc('axes', prop_cycle=custom_cycler)

fig, ax = plt.subplots(1, 1, figsize=(10, 7))
ind_net = 0
for net in networks:
    plt.gca().set_prop_cycle(custom_cycler)
    
    inputfolder = folder + '/' + net
    
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
            
        ind+=1
    
    lineObjects = ax.plot(dd, smooth(qq2,1), linestyle=ls[ind_net])
    plt.gca().set_prop_cycle(custom_cycler)
    ax.plot(dd[:,0:sig_nb], qq2[:,0:sig_nb], marker=markers[ind_net], linestyle="", ms=7)
    
    ind_net += 1
        

# Noise results
filename= inputfolder + '/results_AA_' + filt + '_f2_noise.txt'
a= np.loadtxt(filename, dtype='f', delimiter=' ')

# Percentiles
q50=np.percentile(a[:,index], 50)
q95=np.percentile(a[:,index], 95)
q5=np.percentile(a[:,index], 5)

lineObjects += ax.plot(dist_max, q50*np.ones(len(dist_max)), 'k')

# blue-ish area corresponding to the noise reponse for the last network
ax.fill_between(dist_max, q5, q95, alpha=0.05, facecolor='b')

# Add distances from know objects
ax.axvline(8.2, color='k', linestyle='--', lw=lw-1)
ax.text(6, 1.04, 'SgrA*', fontsize=0.8*fs)
ax.axvline(50, color='k', linestyle='--', lw=lw-1)
ax.text(47, 1.04, 'LMC', fontsize=0.8*fs)

# Set labels and legend
ax.set_xlim(1,100)
ax.set_xlabel('Distance [kpc]', fontsize=fs)

if quantity == "coverage":
    ax.set_ylim((0,1.02))
    ax.set_ylabel('Coverage', fontsize=fs)
    ax.legend(iter(lineObjects), signal_names, loc='upper right', fontsize=0.8*fs)
else:
    ax.set_ylim((0.05,1.05*qq2[0,sig_nb-1]))
    ax.set_ylabel('$\Delta$')
    ax.legend(iter(lineObjects), signal_names, loc='upper right', fontsize=0.8*fs)

ax.tick_params(labelsize=0.8*fs)

ax.grid(True)

fig_name = 'perfHLVKA_'+folder+'_s25-40.png'
plt.savefig(fig_name, bbox_inches='tight')
fig_name = 'perfHLVKA_'+folder+'_s25-40.pdf'
plt.savefig(fig_name, bbox_inches='tight')

plt.show()    

