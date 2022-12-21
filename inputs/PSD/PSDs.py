#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Standard python numerical analysis imports:
import numpy as np
import matplotlib.pyplot as plt


# Get ET spectrum
ET_D_asd=np.loadtxt('ET_D_sensitivity.txt', skiprows=0)

# Get CE spectrum
#CE1_asd=np.loadtxt('curves_Jan_2020/ce1.txt', skiprows=0)
CE1_asd=np.loadtxt('ce_strain/cosmic_explorer_20km.txt', skiprows=0)

# Get CE spectrum
#CE2_asd=np.loadtxt('curves_Jan_2020/ce2.txt', skiprows=0)
CE2_asd=np.loadtxt('ce_strain/cosmic_explorer.txt', skiprows=0)

# Get aLIGO spectrum
AL_asd=np.loadtxt('AplusDesign.txt', skiprows=0)

# Get ADV spectrum
ADV_asd=np.loadtxt('avirgo_O5high_NEW.txt', skiprows=0)

# Get KAG spectrum
KAG_asd=np.loadtxt('kagra_128Mpc.txt', skiprows=0)


#########################################################################
fmin = 10
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

plt.loglog(AL_asd[:,0],AL_asd[:,1], label='A+ design')
plt.loglog(ADV_asd[:,0],ADV_asd[:,1], label='AdV+ phase 2')
plt.loglog(KAG_asd[:,0],KAG_asd[:,1], label='KAGRA design')
plt.loglog(ET_D_asd[:,0],ET_D_asd[:,3], label='ET')
plt.loglog(CE1_asd[:,0],CE1_asd[:,1], label='CE 20km')
plt.loglog(CE2_asd[:,0],CE2_asd[:,1], label='CE 40km') 

# Set labels and legend
ax.grid(True)
ax.set_xlabel('Frequency [Hz]', fontsize=fs)
ax.set_xlim(xmin=fmin, xmax=fmax)
ax.set_ylabel(r'Strain noise $[1/ \sqrt{Hz}]$', fontsize=fs)
ax.set_ylim(ymin=1e-25, ymax=1e-21)
ax.legend(loc='upper right', fontsize=0.8*fs)
ax.tick_params(labelsize=0.8*fs)

fig_name = 'asds';
plt.savefig(fig_name+'.png', bbox_inches='tight')
plt.savefig(fig_name+'.pdf', bbox_inches='tight')
    
plt.show()
