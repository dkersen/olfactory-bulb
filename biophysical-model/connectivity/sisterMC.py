# generates distribution of shared GC connectivity for sister MCs
import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('sisterMC.mat')
total = mat['total']
fig = plt.figure(figsize=(5, 3.5))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8]) # main axes
ax.hist(total, bins = 35, normed=False, alpha=1, color='blue')
ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax.set_xticklabels(['0','0.25', '0.5', '0.75','1'])
ax.set_yticks([0, 5538.1,11076])
ax.set_yticklabels(['0','0.08','0.16'])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()
