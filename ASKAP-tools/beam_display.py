import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')
import os
from os import path
from matplotlib import colors

pitchangle=0.9
d=pitchangle


twilightcmap=plt.cm.twilight

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rc('font', size=12)          # controls default text sizes
fig_directory='/Users/emma/OneDrive/PhD/thesis/Figures/'


def main():
	beamnumber=0
	beamcenters=np.empty(shape=(6,6,3))

	#these are the centre coordinates of beam zero for closepack
	x0=3*d - d/(2*np.sqrt(2))
	y0=-1*d/(2*np.sqrt(2)) - np.sqrt(3)*d

	fig,axs=plt.subplots(1,2,sharey=True)
	axs[0].set_xlim(-3.5,3.5)
	axs[0].set_ylim(-3.5,3.5)
	axs[0].set_aspect("equal")
	for j in range(0,6):
		for i in range(0,6):
			y=y0 + j*d*(np.sqrt(3)/2)
			if j==0 or j==2 or j==4:
				x=x0 - i*d 
			else: 
				x=x0 - (i+0.5)*d 
			
			ellipse=Ellipse(xy=(x,y),width=d,height=d,angle=0,edgecolor=twilightcmap(0.5),fc='none',lw=1)
			axs[0].add_patch(ellipse)
			axs[0].annotate(text=str(beamnumber),xy=(x,y),color=twilightcmap(0.3),ha='center',va='center')
			beamnumber+=1

	#px=x0 - (5.5)*d 
	#py=y0+ 5.8*d*(np.sqrt(3)/2)
	#rectangle=Rectangle(xy=(px,py),width=d,height=0.05,color=twilightcmap(0.7))
	#axs[0].add_patch(rectangle)
	#axs[0].annotate(text='Pitch',xy=(px,py+0.12),color=twilightcmap(0.5))
	axs[0].set_xlabel('x offset (degrees)')
	axs[0].set_ylabel('y offset (degrees)')

	x0=+2.5*d
	y0=-2.5*d
	axs[1].set_xlim(-3.5,3.5)
	axs[1].set_ylim(-3.5,3.5)

	square_beams=[31,30,29,28,27,26,32,13,12,11,10,25,33,14,2,3,9,24,34,15,0,1,8,23,35,4,5,6,7,22,16,17,18,19,20,21]
	k=0
	for j in range(0,6):
		for i in range(0,6):
			x=x0-(i*d)
			y=y0+(j*d)
			ellipse=Ellipse(xy=(x,y),width=d,height=d,angle=0,edgecolor=twilightcmap(0.5),fc='none',lw=1)
			axs[1].add_patch(ellipse)
			axs[1].annotate(text=str(square_beams[k]),xy=(x,y),color=twilightcmap(0.3),ha='center',va='center')
			k+=1

	axs[1].set_aspect("equal")
	axs[1].set_xlabel('x offset (degrees)')
	axs[1].tick_params(axis='y',which='both',left=False)

	#px=-2.5*d
	#py=d*3.2
	#rectangle=Rectangle(xy=(px,py),width=d,height=0.05,color=twilightcmap(0.7))
	#axs[1].add_patch(rectangle)
	#axs[1].annotate(text='Pitch',xy=(px,py+0.12),color=twilightcmap(0.5))


	axs[0].set_title('closepack36')
	axs[1].set_title('square_6x6')

	plt.subplots_adjust(left=0.1, bottom=0.01, right=0.99, top=0.99, wspace=0, hspace=0.01)

	plt.savefig(fig_directory+'askapbeams.png',transparent='True',dpi=600)
	#plt.show()


if __name__ == "__main__":
	main()
