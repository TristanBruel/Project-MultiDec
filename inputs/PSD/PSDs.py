# Standard python numerical analysis imports:
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns


# Get ET spectrum
ET_D_psd=np.loadtxt('ET_D_sensitivity.txt', skiprows=0)

# Get CE spectrum
CE1_psd=np.loadtxt('curves_Jan_2020/ce1.txt', skiprows=0)

# Get CE spectrum
CE2_psd=np.loadtxt('curves_Jan_2020/ce2.txt', skiprows=0)

# Get aLIGO spectrum
AL_psd=np.loadtxt('ALIGO_sensitivity.txt', skiprows=7)

# Get ADV spectrum
ADV_psd=np.loadtxt('AVIRGO_sensitivity.txt', skiprows=7)

# Get KAG spectrum
KAG_psd=np.loadtxt('KAGRA_sensitivity.txt', skiprows=7)


#########################################################################

# plot the spectra
fmin=10
fmax=3000

fig1 = plt.figure()
colors = sns.color_palette(None, 6)
plt.loglog(KAG_psd[:,0],KAG_psd[:,5],color=colors[0],label='KAGRA design')
plt.loglog(ADV_psd[:,0],ADV_psd[:,6],color=colors[1],label='Virgo design')
plt.loglog(AL_psd[:,0],AL_psd[:,6],color=colors[2],label='aLIGO design')
plt.loglog(ET_D_psd[:,0],ET_D_psd[:,3],color=colors[3],label='ET')
plt.loglog(CE1_psd[:,0],CE1_psd[:,1],color=colors[4],label='CE stage 1')
#plt.loglog(CE1_psd[:,0],4*CE1_psd[:,1],color=colors[5],label='CE 20km')
plt.loglog(CE2_psd[:,0],CE2_psd[:,1],color=colors[5],label='CE stage 2')


#plt.axis([1, fmax, 1e-25, 4e-20])

plt.grid(True)
plt.ylabel('ASD [strain/sqrt(Hz)]')
plt.xlabel('Frequency [Hz]')
plt.legend(loc='upper right')

fig_name = 'PSDs_2.png';
plt.savefig(fig_name)
fig_name = 'PSDs_2.pdf';
plt.savefig(fig_name)
    
plt.show()
