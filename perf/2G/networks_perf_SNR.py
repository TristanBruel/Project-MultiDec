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

folder = 'favourable'

# signals = ["s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", 
#             "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220"]
# sig_nb=np.size(signals)
# signal_names = ["s11", "s15", "s15S", "s15G", "s20", "s20S",
#                 "s25", "s40", "no signal"]

signals = ["s15--3D_eqtr", "s15--3D_pole", "s11.2--LS220", "s15.0--LS220", 
           "s15.0--SFHo", "s15.0--GShen", "s20.0--LS220", "s20.0--SFHo",
           "s25.0--LS220", "s40.0--LS220"]
sig_nb=np.size(signals)
signal_names = ["s15--3De", "s15--3Dp", "s11", "s15", 
                "s15S", "s15G", "s20", "s20S", 
                "s25", "s40"]

filt = "spectrum"

dist_nb = 61
dist_max = np.zeros(dist_nb)

networks = ["HL", "HLVKA"]
ls = ["--", "-"]
markers = ['+', 'x']

networks = ['HLVKA']
ls = ['-']
markers = ['x']

### Plotting parameters ###
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00', '#000000']

fs = 20
lw=3

plt.rc('font', family='serif')
plt.rc('lines', linewidth=lw)
custom_cycler = plt.cycler(color=CB_color_cycle)
plt.rc('axes', prop_cycle=custom_cycler)

fig, ax = plt.subplots(1, 1, figsize=(10, 7))
ind_net = 0
for net in networks:
    print('Network', net)
    plt.gca().set_prop_cycle(custom_cycler)
    
    inputfolder = folder + '/' + net
    
    qq1=np.zeros((dist_nb,sig_nb))
    qq2=np.zeros((dist_nb,sig_nb))
    qq3=np.zeros((dist_nb,sig_nb))
    dd=np.zeros((dist_nb,sig_nb))
    SNRs=np.zeros((dist_nb,sig_nb))
    
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
        
        SNR_file= inputfolder + '/SNRs_' + sig + '_' + net + '.txt'
        SNR_sig = np.loadtxt(SNR_file, skiprows=1)
        #SNRs[:,ind] = SNR_sig[:,0] # SNR in LIGO Hanford
        #SNRs[:,ind] = np.max(SNR_sig, axis=1) # max SNR in the network
        SNRs[:,ind] = np.sqrt(np.sum(SNR_sig**2, axis=1)) # quadratic SNR
        
        # dist_max is the list that contains the largest distance
        if max(dist)>max(dist_max):
            dist_max = dist
            
        ind+=1
    
    lineObjects = ax.plot(SNRs, smooth(qq2,1), linestyle=ls[ind_net])
    plt.gca().set_prop_cycle(custom_cycler)
    ax.plot(SNRs, qq2, marker=markers[ind_net], linestyle="", ms=7)
    
    ind_net += 1
    

# Set labels and legend
ax.set_xlim(0,100)
#ax.set_xlabel('SNR in LIGO Hanford', fontsize=fs)
#ax.set_xlabel('Maximum SNR in network HLVKA', fontsize=fs)
ax.set_xlabel('Quadratic SNR in network HLVKA', fontsize=fs)

if quantity == "coverage":
    ax.set_ylim((0,1.02))
    ax.set_ylabel('Coverage', fontsize=fs)
    ax.legend(iter(lineObjects), signal_names, loc='lower right', fontsize=0.8*fs)
else:
    ax.set_ylim((0.05,1.05*qq2[0,sig_nb-1]))
    ax.set_ylabel('$\Delta$')
    ax.legend(iter(lineObjects), signal_names, loc='upper right', fontsize=0.8*fs)

ax.tick_params(labelsize=0.8*fs)

ax.grid(True)

fig_name = 'perfHLVKA_SNR_'+folder
plt.savefig(fig_name+'.png', bbox_inches='tight')
plt.savefig(fig_name+'.pdf', bbox_inches='tight')

plt.show()    

