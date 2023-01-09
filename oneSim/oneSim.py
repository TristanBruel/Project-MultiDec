#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

folder = './stdLike/'

times = np.loadtxt(folder+'times.txt')
freqs = np.loadtxt(folder+'freqs.txt')
energies = np.loadtxt(folder+'energies.txt')

# t_bounce = 0.185
# t_bounce = 0.2152

### Plotting parameters ###
fs = 20
lw = 3

plt.rc('font', family='serif')
plt.rc('lines', linewidth=lw)
cmap = 'inferno'

### Plot the time-frequency map of standard likelihood ###
fig, ax = plt.subplots(1, 1, figsize=(10, 7))
X, Y = np.meshgrid(times, freqs)
im = ax.pcolormesh(X, Y, energies, 
                   norm=mpl.colors.Normalize(),
                   cmap=cmap)

Emin = int(energies.min())
Emax = int(energies.max())+1
nticks = 2*(Emax-Emin)+1
colorbar_ticks = np.linspace(Emin, Emax, nticks)
cb = fig.colorbar(im, ticks=colorbar_ticks)
cb.set_label('log(standard likelihood)', fontsize=fs)
cb.ax.set_yticklabels(['{:.1f}'.format(i) for i in cb.ax.get_yticks()])
cb.ax.tick_params(labelsize=0.8*fs)

# Set labels and legend
ax.set_xlabel('Time after bounce [s]', fontsize=fs)
# ax.set_xticks([0,0.5,1,1.5])
ax.set_ylabel('Frequency [Hz]', fontsize=fs)

ax.tick_params(labelsize=0.8*fs)
ax.yaxis.offsetText.set_fontsize(0.8*fs)

fig_name = folder+'stdLike'
plt.savefig(fig_name+'.png', bbox_inches='tight')
plt.savefig(fig_name+'.pdf', bbox_inches='tight', format='pdf')

#%%
### Plot the tracking on the TF map ###
maxf = np.loadtxt(folder+'maxf.txt')
fit = np.loadtxt(folder+'fit.txt')

fig, ax = plt.subplots(1, 1, figsize=(10, 7))
X, Y = np.meshgrid(times, freqs)
im = ax.pcolormesh(X, Y, energies, 
                   norm=mpl.colors.Normalize(),
                   cmap=cmap)

ax.scatter(maxf[:,0], maxf[:,1], c='b', label='Maxima', s=20)
ax.plot(fit[:,0], fit[:,1], c='k', lw=4, label='Fit')

cb = fig.colorbar(im, ticks=colorbar_ticks)
cb.set_label('log(standard likelihood)', fontsize=fs)
cb.ax.set_yticklabels(['{:.1f}'.format(i) for i in cb.ax.get_yticks()])
cb.ax.tick_params(labelsize=0.8*fs)

# Set labels and legend
ax.set_xlabel('Time after bounce [s]', fontsize=fs)
ax.set_ylabel('Frequency [Hz]', fontsize=fs)
ax.legend(loc='upper left', fontsize=0.8*fs)

ax.tick_params(labelsize=0.8*fs)
ax.yaxis.offsetText.set_fontsize(0.8*fs)

fig_name = folder+'tracking'
plt.savefig(fig_name+'.png', bbox_inches='tight')
plt.savefig(fig_name+'.pdf', bbox_inches='tight', format='pdf')
    
plt.show()

#%%
### Plot the ratio estimated and the true one ###
true_ratio = np.loadtxt(folder+'true_ratio.txt')
pred = np.loadtxt(folder+'pred.txt')

fig, ax = plt.subplots(1, 1, figsize=(10, 7))

ax.errorbar(pred[:,0], pred[:,1], 
            yerr=np.array([pred[:,1]-pred[:,2], pred[:,3]-pred[:,1]]), 
            c='k', ecolor='gray', capsize=10, fmt='o',
            label='Inferred ratio')
ax.scatter(true_ratio[:,0], true_ratio[:,1], c='r', marker='D', zorder=3,
           label='True ratio from CCSN simulation')

# Set labels and legend
ax.set_xlabel('Time after bounce [s]', fontsize=fs)
ax.set_ylabel('Frequency [Hz]', fontsize=fs)
ax.legend(loc='upper left', fontsize=0.8*fs)
ax.grid()

ax.tick_params(labelsize=0.8*fs)
ax.yaxis.offsetText.set_fontsize(0.8*fs)

fig_name = folder+'ratio'
plt.savefig(fig_name+'.png', bbox_inches='tight')
plt.savefig(fig_name+'.pdf', bbox_inches='tight', format='pdf')
    
plt.show()