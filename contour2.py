import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from pylab import figure, show, rand

j=2
omeg_n_1, delta_n_1=np.genfromtxt("/Users/behnam/work/python_SN/in_1_sigma_north.txt",unpack=True)
omeg_n_2, delta_n_2=np.genfromtxt("/Users/behnam/work/python_SN/in_2_sigma_north.txt",unpack=True)
omeg_s_1, delta_s_1=np.genfromtxt("/Users/behnam/work/python_SN/in_1_sigma_south.txt",unpack=True)
omeg_s_2, delta_s_2=np.genfromtxt("/Users/behnam/work/python_SN/in_2_sigma_south.txt",unpack=True)

#***************************************** 1sigma north
cov_n_1 = np.cov(omeg_n_1, delta_n_1)
lambda_n_1, v_n_1 = np.linalg.eig(cov_n_1)
lambda_n_1 = np.sqrt(lambda_n_1)


ell = Ellipse(xy=(np.mean(omeg_n_1), np.mean(delta_n_1)),width=lambda_n_1[0]*j*2, height=lambda_n_1[1]*j*2,angle=180-np.rad2deg(np.arccos(v_n_1[0, 0]))
,color='blue',linestyle='solid',linewidth=3)


fig = plt.figure(figsize=(11,11))
ax = plt.subplot(111, aspect='equal')


ell.set_facecolor('none')
ax.add_artist(ell)
plt.scatter(omeg_n_1, delta_n_1,color='white')

#****************************************** end of 1sigma north

#****************************************** 2sigma north
cov_n_2 = np.cov(omeg_n_2, delta_n_2)
lambda_n_2, v_n_2 = np.linalg.eig(cov_n_2)
lambda_n_2 = np.sqrt(lambda_n_2)


ell_n_2 = Ellipse(xy=(np.mean(omeg_n_2), np.mean(delta_n_2)),width=lambda_n_2[0]*j*2, height=lambda_n_2[1]*j*2,angle=180-np.rad2deg(np.arccos(v_n_2[0, 0]))
,color='blue',linestyle='dashed',linewidth=3)

ell_n_2.set_facecolor('none')
ax.add_artist(ell_n_2)
plt.scatter(omeg_n_2, delta_n_2,color='white')

#****************************************** end of 2sigma north



j=1
#****************************************** 1sigma south
cov_s_1 = np.cov(omeg_s_1, delta_s_1)
lambda_s_1, v_s_1 = np.linalg.eig(cov_s_1)
lambda_s_1 = np.sqrt(lambda_s_1)


ell_s_1 = Ellipse(xy=(np.mean(omeg_s_1), np.mean(delta_s_1)),width=lambda_n_2[0]*j*2, height=lambda_n_2[1]*j*2,angle=180-np.rad2deg(np.arccos(v_n_2[0, 0]))
,color='red',linestyle='solid',linewidth=3)

ell_s_1.set_facecolor('none')
ax.add_artist(ell_s_1)
plt.scatter(omeg_s_1, delta_s_1)#,color='white')

#****************************************** end of 1sigma south



plt.show()


