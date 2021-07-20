# performs chi-square test for fitting type 2 MC data

import scipy.io
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os

mat = scipy.io.loadmat('MCcount2.mat')
MCcount2 = mat['MCcount2']
loc, scale=stats.expon.fit(MCcount2)

fig = plt.figure(figsize=(10,3.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
(n, bins2, patches) = ax.hist(MCcount2, bins=25, normed=False, alpha=1, color='red')
freq = np.zeros((len(n),))

for i in range(len(bins2)-1):
    freq[i] = stats.expon.cdf(x=bins2[i+1],scale=818.7) - stats.expon.cdf(x=bins2[i],scale=818.7)
n = n*1/sum(n)*len(MCcount2)
freq = freq * len(MCcount2)

(h,p) = stats.chisquare(n,freq)
