import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_agg
from math import pi

#array between 0 and 360 deg
RA = np.random.random(10000)*360
#array between -45 and 90 degrees. By construction!
DEC= np.random.random(10000)*135-45

fig = plt.Figure((10, 4.5))
ax = fig.add_subplot(111,projection='mollweide')
ax.grid(True)
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

ax.set_xticklabels(np.arange(30,331,30))
hist,xedges,yedges = np.histogram2d(DEC,RA,bins=[90,180],range=[[-90,90],[0,360]])
#TO RECOVER THE EXPECTED BEHAVIOUR, I HAVE TO CHANGE -90 FOR -80 IN THE PREVIOUS LINE:
#hist,xedges,yedges = np.histogram2d(DEC,RA,bins=[90,180],range=[[-80,90],[0,360]])
#I DO NOT WHY!

extent = (-pi,pi,-pi/2.,pi/2.)
image = ax.imshow(hist,extent=extent,clip_on=False,aspect=0.5,origin='lower')

cb = fig.colorbar(image, orientation='horizontal')
canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)

fig.canvas.print_figure("image1.png")