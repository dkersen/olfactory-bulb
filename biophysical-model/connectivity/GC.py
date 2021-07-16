# Generates distribution of GC-to-MC connectivity

import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('GCcount.mat')
GCcount = mat['GCcount']
ae, loce, scalee =stats.skewnorm.fit(GCcount)

fig = plt.figure(figsize=(10,3.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
(n,binners,patches)=ax.hist(GCcount, bins=37, normed=True, alpha=1, color='g')
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = stats.skewnorm.pdf(x,ae, loce, scalee) 
ax.plot(x, p, 'k', linewidth=5)
ax.set_xticks([0,100,200,300])
ax.set_xticklabels(['0','100','200','300'])
ax.set_yticks([0, 0.01027777777777778])
ax.set_yticklabels(['0','0.08'])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
