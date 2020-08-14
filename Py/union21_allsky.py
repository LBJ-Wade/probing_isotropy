#plot union2.1 all sky

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib
import math


def galat(alpha,delta):
    return np.arcsin(np.cos(delta)*np.cos(0.4735)*np.cos(alpha-3.3660)+np.sin(delta)*np.sin(0.4735))

def galong(alpha,delta,gal_lat):
    return np.arctan2((np.sin(delta)-np.sin(gal_lat)*np.sin(0.4735)),(np.cos(delta)*np.sin(alpha-3.3660)*np.cos(0.4735)))+0.5747


galat=np.vectorize(galat) 
galong=np.vectorize(galong) 

alpha = np.linspace(0,2*np.pi,1000)

equatorlat=galat(alpha,0)
equatorlong=galong(alpha,0,equatorlat)

for i in range(len(equatorlat)):
    if equatorlong[i]<np.pi:
        equatorlong[i]=-equatorlong[i]
    else:
        equatorlong[i]=2*np.pi-equatorlong[i]
        
        
num=580
l=np.zeros(num)
b=np.zeros(num)

z,dm,dm_err,longt,lat=np.genfromtxt("z_dm_dm-err_long_lat.txt",unpack=True)

for i in range(num):
    b[i]=lat[i]*math.pi/180.0
    if longt[i]<180:
        l[i]=-(longt[i]*math.pi/180.0)
    else:
        l[i]=2*math.pi-(longt[i]*math.pi/180.0)
    
    
    
fig = plt.figure(figsize=(13,8))
#fig = plt.figure()
ax = fig.add_subplot(111, projection="mollweide")
plt.setp(ax.get_yticklabels(), visible=False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True)
ax.grid(color='black', linestyle='-', linewidth=1)
im=ax.scatter(l,b,c=z,s=80,edgecolors='none')
cbar=fig.colorbar(im,cax=None,ticks=np.arange(0.0,1.5,0.2),shrink=0.7,orientation="horizontal",pad=0.02)
cbar.set_label("Redshift", size=25)
cbar.ax.tick_params(labelsize=25) 

imequator=ax.scatter(equatorlong,equatorlat,s=10,marker='.',facecolor='None',lw=0.8,color="black")

ax.text(0,-1.5,"(0,-90)",fontsize=22)
ax.text(0,1.4,"(0,90)",fontsize=22)
ax.text(-math.pi/2.0,0,"(90,0)",fontsize=22)
ax.text(math.pi/2.0,0,"(270,0)",fontsize=22)
ax.text(0,0,"(0,0)",fontsize=22)
ax.text(-0.65,0.4,"Celestial Equator",fontsize=16,rotation='65')

plt.tight_layout(pad=0.9)

plt.show()
