# generates distribution of MC-to-GC connectivity for type I MCs

import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('MCcount1.mat')
MCcount1 = mat['MCcount1']
loc,scale=stats.expon.fit(MCcount1)

fig = plt.figure(figsize=(5, 3.5))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8]) # main axes
(n, binners, patches) = ax.hist(MCcount1, bins=15, normed=True, alpha=1, color='darkred')
xmin, xmax = plt.xlim()
x = np.linspace(0, xmax, 100)
p = stats.expon.pdf(x,scale=1426) 
ax.plot(x, p, 'k', linewidth=5)
ax.set_xticks([0, 3000, 6000])
ax.set_xticklabels(['0','3000','6000'])
ax.set_yticks([0,0.0008209538702111023])
ax.set_yticklabels(['0', '0.35'])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

