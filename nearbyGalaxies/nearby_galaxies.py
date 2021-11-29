import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import seaborn as sns
import cartopy.crs as ccrs


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


### Import galaxy sky locations and names ###
with open('nearby_galaxies.csv') as csv_gal:
    galaxies = csv.reader(csv_gal)
    names = []
    RA = np.array([])
    DEC = np.array([])
    distances = np.array([])
    print("Columns in csv file: " + csv_gal.readline())
    for row in galaxies:
        #print(row)
        names.append(row[1])
        RA = np.concatenate((RA,[float(row[2])]))
        DEC = np.concatenate((DEC,[float(row[3])]))
        distances = np.concatenate((distances,[float(row[4])]))


galaxy_nb = len(names)

signals = ["s20.0--LS220"]

filt = "spectrum"

networks = ["CEH+CEL", "ET", "CEH+CEL+ET"]

for net in networks:
    print("Network: "+net)
    for sig in signals:
        # plt.figure()
        colors = sns.color_palette(None, galaxy_nb)
        
        time_nb=25
    
        qq1=np.zeros((time_nb,galaxy_nb))
        qq2=np.zeros((time_nb,galaxy_nb))
        qq3=np.zeros((time_nb,galaxy_nb))
        tt=np.zeros((time_nb,galaxy_nb))
        
        quantity="coverage"
        #quantity="delta"
        
        if quantity == "coverage":
            index=1
        else:
            index=5
    
        ind=0
        
        for gal in names:    
    
            filename=net+'/results_AA_'+filt+'_f2_'+sig+'_'+gal+'.txt'
            
            a= np.loadtxt(filename, dtype='f', delimiter=' ')
            
            time,q1,q2,q3=buildbox(a,index)
            qq1[:,ind]=q1
            qq2[:,ind]=q2
            qq3[:,ind]=q3
            tt[:,ind]=time
            
            ind += 1
    
    ### Projection on the sky ###
    average_cov = np.mean(qq2[:-1,],0) # remove the last item which corresponds to the 24th hour
    print("Galaxy name, Average coverage")
    for i in range(len(names)):
        print(names[i],average_cov[i])
    
    deg = 180/12
    ra = RA*deg
    dec = DEC
    
    
    fig = plt.figure(figsize=(10,5))
    
    ax = fig.add_subplot(1,1,1, projection=ccrs.Mollweide())
    ax.set_global()
    
    cmap = plt.get_cmap("plasma")
    
    plt.scatter(-ra, dec, marker = 'o', color=cmap(average_cov), cmap=cmap,
                transform=ccrs.PlateCarree())
    
    sm = cm.ScalarMappable(cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, orientation='horizontal',fraction=0.05,pad=0.05)
    
    ax.gridlines()
    
    fig_name = net+'nearbygala.png';
    plt.savefig(fig_name)
    fig_name = net+'nearbygala.pdf';
    plt.savefig(fig_name)
    
    plt.show()
