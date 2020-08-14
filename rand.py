import numpy as np
import matplotlib.pyplot as plt
a=np.random.normal(0.3,0.1,10000)
b=np.random.normal(0,0.11,10000)

c=np.random.random(10000)
plt.hist(c,bins=200)
plt.show()