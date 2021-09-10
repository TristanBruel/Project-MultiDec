import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.ndimage.filters import gaussian_filter1d


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

networks = ["HL","HLV","HLVI"]
#networks = ["HLV","HLVK"]

net_nb = len(networks)

signal = ["s20.0--LS220"]

signal_nb = len(signal)

#filtering = "prewhiten"
filtering = "spectrum"


for sig in signal:
    plt.figure()
    #colors = ["blue", "red", "orange", "green"]
    colors = sns.color_palette(None, net_nb)
    
    time_nb=96

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
        filename='results_AA_' + filtering + '_f2_' + sig + '_multiDec' + net + '.txt'
        
        a= np.loadtxt(filename, dtype='f', delimiter=' ')
        
        time,q1,q2,q3=buildbox(a,index)
        qq1[:,ind]=q1
        qq2[:,ind]=q2
        qq3[:,ind]=q3
        tt[:,ind]=time
        
        plt.plot(time, gaussian_filter1d(q2, sigma=1), color=colors[ind],label=net)
        #plt.plot(time, q2, marker="+",linestyle="", color=colors[ind])
        
        ind += 1
        print(net, np.mean(q2))
        
    # det_nb=len(net)
    # Feq=np.zeros((time_nb,det_nb))
    
    # for det in range(det_nb):
    #     F=np.unique(a[:,(6+det)])
    #     indexes = np.unique(a[:,(6+det)], return_index=True)[1]
    #     F=[a[index,(6+det)] for index in sorted(indexes)]
    #     Feq[:,det]=F
    
    #     plt.plot(tt[:,ind-1], abs(Feq[:,det]), marker="+",linestyle="",color=colors[det],label=net[det])
    #     plt.legend()
        
        
            
    time_max=time[time_nb-1]
    
    if quantity == "coverage":
        plt.ylim((0,1.02))
        plt.xlim((1,time_max))
        plt.ylabel('Coverage')
        plt.legend()
    else:
        plt.ylim((0.05,1.05*qq2[0,net_nb-1]))
        plt.xlim((1,time_max))
        plt.ylabel('$\Delta$')
        
    
    plt.xlabel('Time [hours]')
    
    plt.grid(True)
    
    plt.title('ccSN in galactic center')
    
    plt.show()
