# generates distribution of MC-to-GC connectivity for type II MCs

import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('MCcount2.mat')
MCcount2 = mat['MCcount2']
loc,scale=stats.expon.fit(MCcount2)

fig = plt.figure(figsize=(5, 3.5))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8]) # main axes
(n, binners, patches) =ax.hist(MCcount2, bins=25, normed=True, alpha=1, color='darkred')
xmin, xmax = plt.xlim()
x = np.linspace(0, xmax, 100)
p = stats.expon.pdf(x,scale=scale) 
ax.plot(x, p, 'k', linewidth=5)
ax.set_xticks([0, 2500,5000])
ax.set_xticklabels(['0','2500','5000'])
ax.set_yticks([0, 0.0016393442622950822])
ax.set_yticklabels(['0','0.3'])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

