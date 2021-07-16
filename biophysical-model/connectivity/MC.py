# Generates distribution of MC-to-GC connectivity for all MCs

import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('MCcount.mat')
MCcount = mat['MCcount']
loc, scale=stats.expon.fit(MCcount)

fig = plt.figure(figsize=(10,3.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
(n,binners,patches)=ax.hist(MCcount, bins=30, normed=True, alpha=1, color='red')
xmin, xmax = plt.xlim()
x = np.linspace(0, xmax, 100)
p = stats.expon.pdf(x, scale=1225.83) 
ax.plot(x, p, 'k', linewidth=5)
ax.set_xticks([0, 3000, 6000])
ax.set_xticklabels(['0','3000','6000'])
ax.set_yticks([0, 0.0009369144284821987])
ax.set_yticklabels(['0', '0.2'])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()

