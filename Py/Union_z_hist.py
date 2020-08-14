import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("/Users/behnam/work/Union2.1/25_prob_PDF/30/99.txt")



plt.hist(data,rwidth=0.9,color='orange',bins=50)


plt.tick_params(labelsize=20,pad=10)

plt.xlabel("P",fontsize=20)
plt.ylabel('N',fontsize=20)

plt.tight_layout(pad=0.9)

#plt.xlimit(0,0.2)

plt.show()
