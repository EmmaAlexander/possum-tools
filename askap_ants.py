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

cmap=plt.cm.twilight


plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rc('font', size=10)          # controls default text sizes
fig_directory='/Users/emma/OneDrive/PhD/thesis/Figures/'
plt.rc('legend', fontsize=8)    # legend fontsize


def main():

	ants=np.loadtxt('askap_ants.txt',dtype={'names': ('ak', 'x', 'y', 'el'),
                     'formats': ('S1', 'f4', 'f4', 'f4')})
	x=ants['x']/1000 - ants[29]['x']/1000
	y=ants['y']/1000 - ants[29]['y']/1000


	ants_16=np.loadtxt('askap_ants_16.txt',dtype={'names': ('ak', 'x', 'y', 'el'),
                     'formats': ('S1', 'f4', 'f4', 'f4')})
	x_16=ants_16['x']/1000 - ants[29]['x']/1000
	y_16=ants_16['y']/1000 - ants[29]['y']/1000

	fig=plt.figure(figsize=(6,6))
	ax1=plt.subplot(221)
	ax2=plt.subplot(222)
	ax3=plt.subplot(223)
	ax4=plt.subplot(224)
	#add in the ants
	ax1.scatter(x,y,s=10,color=cmap(0.4),zorder=80,label='ASKAP-36')
	ax1.scatter(x_16,y_16,marker='x',s=10,color=cmap(0.7),zorder=90,label='ASKAP-16')

	#lines from ants to origin
	for i in range(0,36):
		x_values=[0,x[i]]
		y_values=[0,y[i]]
		ax1.plot(x_values, y_values,color=cmap(0.1),zorder=i)

	#highlight ak02 and ak05
	
	for i in range(1,4):			
		ellipse=Ellipse(xy=(0,0),width=i*2,height=i*2,angle=0,edgecolor=cmap(0.9),fc='none',lw=1,zorder=50)
		ax1.add_patch(ellipse)

	ax1.set_aspect("equal")
	ax1.set_xlim(-3.5,3.5)
	ax1.set_ylim(-3.5,3.5)
	ax1.set_xlabel('Eastern distance from ak01 [km]')
	ax1.set_ylabel('Northern distance from ak01 [km]')
	ax1.set_title('a) ASKAP antenna positions')
	#axs[0].legend(loc='upper right',framealpha=1,facecolor='w')

	#okay now zoom in 
	#add in the ants
	ax2.scatter(x,y,s=10,color=cmap(0.4),zorder=80)
	ax2.scatter(x_16,y_16,s=10,marker='x',color=cmap(0.7),zorder=90)

	#lines from ants to origin
	for i in range(0,36):
		x_values=[0,x[i]]
		y_values=[0,y[i]]
		ax2.plot(x_values, y_values,color=cmap(0.1),zorder=i)

	#highlight ak02 and ak05
	#axs[1].scatter(x[28],y[28],color=cmap(0.75),zorder=100)
	#axs[1].scatter(x[14],y[14],color=cmap(0.75),zorder=101)
	
	for i in range(1,4):			
		ellipse=Ellipse(xy=(0,0),width=i*2,height=i*2,angle=0,edgecolor=cmap(0.9),fc='none',lw=1,zorder=1)
		ax2.add_patch(ellipse)

	ax2.set_aspect("equal")
	ax2.set_xlim(-1.01,1.01)
	ax2.set_ylim(-1.01,1.01)
	ax2.set_xlabel('Eastern distance from ak01 [km]')
	ax2.set_ylabel('Northern distance from ak01 [km]')
	ax2.set_title('b) Antenna positions (zoom)')

	#plt.savefig(fig_directory+'askapants.png',transparent='True',dpi=600)
	#plt.show()


	#snapshot uv coverage
	#fig2, axs = plt.subplots(1,2)

	bl=[]

	for i in range(0,36):
		for j in range(0,36):
			if i!=j:
				u=x[i]-x[j]
				v=y[i]-y[j]
				r=np.sqrt((u*u)+(v*v))
				bl.append(r)

				ax3.scatter(u,v,marker='o',s=10,color=cmap(0.4))

	ax3.set_xlabel('u [km]')
	ax3.set_ylabel('v [km]')
	ax3.set_title('c) Snapshot uv coverage')

	bl=np.array(bl)

	bl_16=[]
	for i in range(0,16):
		for j in range(0,16):
			if i!=j:
				u=x_16[i]-x_16[j]
				v=y_16[i]-y_16[j]
				r=np.sqrt((u*u)+(v*v))
				bl_16.append(r)
				ax3.scatter(u,v,marker='o',s=10,color=cmap(0.7))
	ax3.set_aspect("equal")

	bl_16=np.array(bl_16)

	print(np.nanmin(bl))
	print(np.nanmin(bl_16))

	print(np.nanmax(bl))
	print(np.nanmax(bl_16))

	bins_array=np.arange(start=0,stop=6.6,step=0.2)
	kwargs = dict(histtype='stepfilled',bins=bins_array,alpha=1,linewidth=1)

	ax4.hist(np.ndarray.flatten(bl),**kwargs,color=cmap(0.4),label='ASKAP-36')
	ax4.hist(np.ndarray.flatten(bl_16),**kwargs,color=cmap(0.7),label='ASKAP-16')
	ax4.set_xlim(0,6.4)
	ax4.set_aspect(1./ax4.get_data_ratio(),adjustable='box')
	ax4.set_xlabel('uv distance [km]')
	ax4.set_ylabel('No. baselines')
	ax4.legend(loc='upper right',framealpha=1,facecolor='w')
	ax4.set_title('d) Distribution of baseline lengths')

	plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.4, hspace=0.4)
	#plt.savefig(fig_directory+'askapbaselines.png',transparent='True',dpi=600)
	plt.show()


if __name__ == "__main__":
	main()