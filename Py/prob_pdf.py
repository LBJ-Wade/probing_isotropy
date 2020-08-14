from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib
import math
import numpy as np
from scipy import stats




p=np.genfromtxt("/Users/behnam/work/Union2.1/25_prob_PDF/90/122.txt",unpack=True)
weights = np.ones_like(p)/len(p)

density = stats.kde.gaussian_kde(p)
#density.covariance_factor = lambda : .5
#density._compute_covariance()
x = np.arange(0., 1.1, .1)

plt.figure(figsize=(10,8))
plt.plot(x, density(x),'r',linewidth=3)
plt.hist(p,rwidth=0.9,weights=weights,color='#A81450',bins=10,alpha=.4)

plt.tick_params(labelsize=20,pad=10)



plt.xlabel('P',fontsize=20)
plt.ylabel('PDF',fontsize=20)
plt.grid(True)
plt.show()

#plt.plot(density)
#plt.show()
