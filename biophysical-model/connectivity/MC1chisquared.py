# performs chi-square test for fitting type 1 MC data

import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('MCcount1.mat')
MCcount1 = mat['MCcount1']
loc, scale=stats.expon.fit(MCcount1)

fig = plt.figure(figsize=(10,3.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
(n, bins2, patches) = ax.hist(MCcount1, bins=15, normed=False, alpha=1, color='red')
freq = np.zeros((len(n),))

for i in range(len(bins2)-1):
    freq[i] = stats.expon.cdf(x=bins2[i+1],scale=1426) - stats.expon.cdf(x=bins2[i],scale=1426)
n = n*1/sum(n)*len(MCcount1)
freq = freq * len(MCcount1)

(h,p) = stats.chisquare(n,freq)
