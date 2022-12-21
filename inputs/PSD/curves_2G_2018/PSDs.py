#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Standard python numerical analysis imports:
import numpy as np
import matplotlib.pyplot as plt


# Get ET spectrum
ET_D_psd=np.loadtxt('ET_D_sensitivity.txt', skiprows=0)

# Get CE spectrum
#CE1_psd=np.loadtxt('curves_Jan_2020/ce1.txt', skiprows=0)
CE1_psd=np.loadtxt('ce_strain/cosmic_explorer_20km.txt', skiprows=0)

# Get CE spectrum
#CE2_psd=np.loadtxt('curves_Jan_2020/ce2.txt', skiprows=0)
CE2_psd=np.loadtxt('ce_strain/cosmic_explorer.txt', skiprows=0)

# Get aLIGO spectrum
AL_psd=np.loadtxt('ALIGO_sensitivity.txt', skiprows=7)

# Get ADV spectrum
ADV_psd=np.loadtxt('AVIRGO_sensitivity.txt', skiprows=7)

# Get KAG spectrum
KAG_psd=np.loadtxt('KAGRA_sensitivity.txt', skiprows=7)


#########################################################################
fmin = 5
fmax = 5000

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
custom_cycler = plt.cycler(color=CB_color_cycle)
plt.rc('axes', prop_cycle=custom_cycler)

fig, ax = plt.subplots(1, 1, figsize=(10, 7))

plt.loglog(AL_psd[:,0],AL_psd[:,5], label='aLIGO design')
plt.loglog(ADV_psd[:,0],ADV_psd[:,5], label='Virgo late high')
plt.loglog(KAG_psd[:,0],KAG_psd[:,5], label='KAGRA design')
plt.loglog(ET_D_psd[:,0],ET_D_psd[:,3], label='ET')
plt.loglog(CE1_psd[:,0],CE1_psd[:,1], label='CE 20km')
plt.loglog(CE2_psd[:,0],CE2_psd[:,1], label='CE 40km') 

# Set labels and legend
ax.grid(True)
ax.set_xlabel('Frequency [Hz]', fontsize=fs)
ax.set_xlim(xmin=fmin, xmax=fmax)
ax.set_ylabel('ASD [strain/sqrt(Hz)]', fontsize=fs)
ax.set_ylim(ymin=1e-25, ymax=1e-21)
ax.legend(loc='upper right', fontsize=0.8*fs)
ax.tick_params(labelsize=0.8*fs)

fig_name = 'PSDs.png';
plt.savefig(fig_name, bbox_inches='tight')
fig_name = 'PSDs.pdf';
plt.savefig(fig_name, bbox_inches='tight')
    
plt.show()
