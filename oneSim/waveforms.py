#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

wvfs_folder = '../inputs/2D_simulations/waveforms/'

wvf = 'gw_s20.0--LS220--GravA'

t_bounce = 0.185
t_bounce = 0.2152

gw_filename = wvfs_folder+wvf+'.dat'
waveform = np.loadtxt(gw_filename)

time = waveform[:,0]-t_bounce
hoft = waveform[:,1]

sample_freq = int(1/(time[1]-time[0]))
resamp_factor = int(sample_freq / 4096)

resampled_hoft, resampled_time = signal.resample(hoft, num=int(len(hoft)/resamp_factor),
                                                 t=time)

### Plotting parameters ###
fs = 24
lw = 3

plt.rc('font', family='serif')
plt.rc('lines', linewidth=lw)

### Plot the time-frequency map of standard likelihood ###
fig, ax = plt.subplots(1, 1, figsize=(20, 7))

#ax.plot(time, hoft, color='k')
ax.plot(resampled_time, resampled_hoft, color='k')

# Set labels and legend
ax.set_xlabel('Time after bounce [s]', fontsize=fs)
ax.set_xticks([0,0.5,1,1.5])
ax.set_ylabel('Strain amplitude', fontsize=fs)
ax.set_ylim(-5.5e-22,5.5e-22)

ax.tick_params(labelsize=0.8*fs)
ax.yaxis.offsetText.set_fontsize(0.8*fs)
ax.grid()
#ax.set_title(wvf, fontsize=fs)

fig_name = 'waveform.png';
# plt.savefig(fig_name, bbox_inches='tight')
    
plt.show()