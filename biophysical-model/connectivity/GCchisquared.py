# performs chi-square test for fitting GC data

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
(n, bins2, patches) = ax.hist(GCcount, bins=37, normed=True, alpha=1, color='g')
freq = np.zeros((len(n),))

for i in range(len(bins2)-1):
    freq[i] = stats.skewnorm.cdf(x=bins2[i+1],a=ae, loc=loce, scale=scalee) - stats.skewnorm.cdf(x=bins2[i],a=ae, loc=loce, scale=scalee)
n = n*1/sum(n)*len(GCcount)
freq = freq * len(GCcount)

(h,p) = stats.chisquare(n,freq)
