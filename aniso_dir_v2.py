import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math

#dipole,-dipole, theta90,theta60,quadrapole,-quadrapole,octapole,-octapole,,theta30
#color=['black','black','blue','red','magenta','magenta','brown','brown','green']


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



longt=[263.99,83.99,33.7,112.5,238.5,58.5,239.0,59,67.5,0,56.2,101]
lat=[48.26,-48.26,-19.5,-9.6,76.6,-76.6,64.3,-64.3,-66.4,-30,-41.8,-41.8]

l=np.zeros(len(longt))
b=np.zeros(len(lat))

for i in range(len(lat)):
    b[i]=lat[i]*math.pi/180.0
    if longt[i]<180:
        l[i]=-(longt[i]*math.pi/180.0)
    else:
        l[i]=2*math.pi-(longt[i]*math.pi/180.0)




fig = plt.figure(figsize=(12,6.5))
ax = fig.add_subplot(111, projection="mollweide")
plt.setp(ax.get_yticklabels(), visible=False)
plt.setp(ax.get_xticklabels(), visible=False)

plt.grid(True)
ax.grid(color='black', linewidth=1)


#im=ax.scatter(l,b,c=color,s=500,edgecolors='none')
for i in range(0,2):
    im=ax.scatter(l[i],b[i],c='#696969',s=400,marker='+', linewidth='6')
    
for i in range(4,6):
    im=ax.scatter(l[i],b[i],c='magenta',s=500,edgecolors='none',marker='^')    

for i in range(6,8):
    im=ax.scatter(l[i],b[i],c='brown',s=500,edgecolors='none',marker='d') 

im=ax.scatter(l[2],b[2],c='b',s=600,edgecolors='none',marker='*')
im=ax.scatter(l[3],b[3],c='r',s=600,edgecolors='none',marker='*')
im=ax.scatter(l[8],b[8],c='g',s=600,edgecolors='none',marker='*')

im=ax.scatter(l[9],b[9],c='b',s=300,edgecolors='none')
im=ax.scatter(l[10],b[10],c='r',s=300,edgecolors='none')
im=ax.scatter(l[11],b[11],c='g',s=300,edgecolors='none')


imequator=ax.scatter(equatorlong,equatorlat,s=10,marker='.',facecolor='None',lw=0.8,color="black",zorder=0)

ax.text(-0.5,0.8,"Celestial Equator",fontsize=16,rotation='60')


ax.annotate('CDP', (l[1]-1.05,b[1]-0.2),fontsize=30,color="#696969")
ax.annotate('CDP', (l[0]+0.1,b[0]),fontsize=30,color="#696969")
ax.annotate('COP', (l[6]+0.1,b[6]-0.05),fontsize=30,color='brown')
ax.annotate('COP', (l[7]+0.1,b[7]),fontsize=30,color='brown')
ax.annotate('CQP', (l[4]-1.3,b[4]-0.3),fontsize=30,color='magenta')
ax.annotate('CQP', (l[5]+0.16,b[5]),fontsize=30,color='magenta')

#ax.annotate('$90^{\circ}$', (l[2]+0.1,b[2]),fontsize=25,color='blue')
#ax.annotate('$60^{\circ}$', (l[3]+0.1,b[3]),fontsize=25,color='red')
#ax.annotate('$30^{\circ}$', (l[8]-0.7,b[8]-0.05),fontsize=25,color='green')

ax.annotate(r'$\frac{\pi}{2}$', (l[2]+0.1,b[2]),fontsize=50,color='blue')
ax.annotate(r'$\frac{\pi}{3}$', (l[3]+0.1,b[3]),fontsize=50,color='red')
ax.annotate(r'$\frac{\pi}{6}$', (l[8]-0.6,b[8]-0.05),fontsize=50,color='green')

ax.annotate(r'$\frac{\pi}{2}$', (l[9]+0.1,b[9]),fontsize=50,color='blue')
ax.annotate(r'$\frac{\pi}{3}$', (l[10]+0.1,b[10]),fontsize=50,color='red')
ax.annotate(r'$\frac{\pi}{6}$', (l[11]+0.15,b[11]+0.09),fontsize=50,color='green')


ax.text(0.4,-1.5,"(0,-90)",fontsize=30)
ax.text(0,1.4,"(0,90)",fontsize=30)
ax.text(-math.pi/2.0,0,"(90,0)",fontsize=30)
ax.text(math.pi/2.0,0,"(270,0)",fontsize=30)
ax.text(0,0,"(0,0)",fontsize=30)

plt.tight_layout(pad=0.9)

plt.show()
