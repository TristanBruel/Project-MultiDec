import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter1d


def buildbox(a,index):
    time=np.unique(a[:,0]) 
    time_nb=len(time)
    #print('Number of time indices: ', time_nb)
    
    data = []
    box_name = []
    q1 = []
    q2 = []
    q3 = []
    for i in range(time_nb):
        ext=np.where((a[:,0] == time[i]))
        data_boxi = a[ext[0][0]:ext[0][-1],index]
        data.append(data_boxi)
        if (i%5 == 0):
            box_name.append(str(time[i]))
        else:
            box_name.append('')
        q1.append(np.percentile(data_boxi, 25))
        q2.append(np.percentile(data_boxi, 50))
        q3.append(np.percentile(data_boxi, 75))
    return time,q1,q2,q3



# 0 time
# 1 covpbb                
# 2 medbandwidth
# 3 absolute residual (mean)
# 4 MSE (mean)
# 5 precision (mean)

networks = ["HL","HLV","HLVK","HLVA","HLVKA"]
net_nb = len(networks)

signals = ["s15--3D_eqtr", "s15--3D_pole", "s11.2--LS220", "s15.0--LS220",
           "s15.0--SFHo", "s15.0--GShen", "s20.0--LS220", "s20.0--SFHo",
           "s25.0--LS220", "s40.0--LS220"]

signal_nb = len(signals)

#filtering = "prewhiten"
filtering = "spectrum"

inputfolder = "./"

### Plotting parameters ###
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
CB_marker_cycle = ['solid', 'dashed', 'dotted',
                   'solid', 'dashed', 'dotted',
                   'solid', 'dashed', 'dotted']

fs = 20
lw=3

plt.rc('font', family='serif')
plt.rc('lines', linewidth=lw)
#custom_cycler = (plt.cycler(color=CB_color_cycle)+plt.cycler(linestyle=CB_marker_cycle))
custom_cycler = plt.cycler(color=CB_color_cycle)
plt.rc('axes', prop_cycle=custom_cycler)

for sig in signals:
    print("Waveform: ",sig)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    time_nb=49

    qq1=np.zeros((time_nb,net_nb))
    qq2=np.zeros((time_nb,net_nb))
    qq3=np.zeros((time_nb,net_nb))
    tt=np.zeros((time_nb,net_nb))
        
    quantity="coverage"
    #quantity="delta"
    
    if quantity == "coverage":
        index=1
    else:
        index=5


    ind=0
    for net in networks:
        filename=inputfolder + net + '/results_AA_' + filtering + '_f2_' + sig + '.txt'
        #filename=inputfolder + net + '/results_AA_' + filtering + '_f2_' + sig + '_8.2kpc.txt'
        
        a= np.loadtxt(filename, dtype='f', delimiter=' ')
        
        time,q1,q2,q3=buildbox(a,index)
        qq1[:,ind]=q1
        qq2[:,ind]=q2
        qq3[:,ind]=q3
        tt[:,ind]=time
        
        color = next(ax._get_lines.prop_cycler)['color']
        smooth = gaussian_filter1d(q2, sigma=1)
        ax.plot(time, smooth, label=net, color=color)
        ax.plot(time, q2, marker="x", linestyle="", color=color, ms=6)
        
        print('Network', net)
        print('Mean over 24h: ')
        print(np.mean(q2[:-1]))
        
        print('Fraction over 0.8: ')
        #print(net, sum(x > 0.8 for x in smooth)/len(smooth))
        print(sum(x >= 0.8 for x in q2[:-1])/len(q2[:-1]))
        
        print('\n')
        
        ind +=1
        
            
    time_max=time[time_nb-1]
    
    # Set labels and legend
    if quantity == "coverage":
        ax.set_ylim((0,1.02))
        ax.set_ylabel('Coverage', fontsize=fs)
        ax.legend(loc='lower right', fontsize=0.8*fs)
    else:
        ax.set_ylim((0.05,1.05*qq2[0,net_nb-1]))
        ax.set_ylabel('$\Delta$', fontsize=fs)
        ax.legend(loc='lower right', fontsize=0.8*fs)
        
    ax.set_xlabel('Time [hours]', fontsize=fs)
    ax.set_xlim((1,time_max))
    ax.set_xticks(np.arange(0,30,6))
    ax.tick_params(labelsize=0.8*fs)
    
    ax.grid(True)
    
    if sig == 's20.0--LS220' or sig == 's15--3D_eqtr':
	    fig_name = sig+'_galacenter';
	    plt.savefig(fig_name+'.png', bbox_inches='tight')
	    plt.savefig(fig_name+'.pdf', bbox_inches='tight')
    
    plt.show()
