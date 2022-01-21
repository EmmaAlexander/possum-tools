from mpl_toolkits.mplot3d import Axes3D 
import math as math
import matplotlib.pyplot as plt
import os
import argparse
import warnings


c=299792458. #m/s

def main(args):

	d=args.distance
	pre_d=d/4. #m, before rotation
	post_d=d/4. #m, after rotation
	n=args.nelec
	lamda=args.wavelength
	amp=args.amplitude

	v=args.verbose
	g=args.grid

	if args.bfield == None and args.rotation != None:
		#specified angle
		angle =args.rotation
		# calculate magnetic field
		angle_r=(angle/360.)*(2*math.pi)
		B=angle_r/(n*d*lamda*lamda)
	elif args.bfield != None and args.rotation == None:
		#specified magnetic field
		B=args.bfield
		# calculate rotation angle 
		angle_r=B*n*d*lamda*lamda
		angle=360.*(angle_r/(2*math.pi))
	elif args.bfield == None and args.rotation == None:
		# neither magnetic field or rotation angle given
		print("Warning: neither B or angle specified.")
		print("Using angle of 360 degrees.")
		angle=360
		angle_r=2*math.pi
		B=angle_r/(n*d*lamda*lamda)	
	elif args.bfield != None and args.rotation != None:
		# both have been given 
		angle=args.rotation
		print("Warning: both B or angle specified.")
		print("Using angle of {} degrees (B input ignored).".format(angle))
		angle_r=(angle/360.)*(2*math.pi)
		B=angle_r/(n*d*lamda*lamda)		
	else:
		print("Something has gone wrong...")

	#RM=2.62*n*d*B*math.pow(10,-13)
	alpha=(2.*math.pi*c)/lamda
	beta=B*n*c*lamda*lamda


	if v:
		print("B = {:.3} T".format(B))
		print("Electron density = {:.3} m^-3".format(n))
		print("Amplitude = {:.3} m".format(amp))
		print("Wavelength = {:.3} m".format(lamda))
		print("Rotating distance = {:.3} m".format(d))
		print("Rotation angle = {:.1f} degrees".format(angle))


	peaks_z=np.arange(0,d,0.01)
	peaks_t=np.divide(peaks_z,c)
	peaks_x=amp*np.cos(np.multiply(beta,peaks_t))
	peaks_y=amp*np.sin(np.multiply(beta,peaks_t))

	plt.rcParams['legend.fontsize'] = 10
	fig = plt.figure()
	ax = fig.gca(projection='3d',aspect='auto')

	plt.plot(peaks_z,peaks_y,peaks_x,color='grey')
	plt.plot(peaks_z,-1*peaks_y,-1*peaks_x,color='k')

	# Prepare arrays x, y, z
	z=np.arange(0,d,0.0001)
	t=np.divide(z,c)
	x=amp*np.multiply(np.sin(np.multiply(alpha,t)),np.cos(np.multiply(beta,t)))
	y=amp*np.multiply(np.sin(np.multiply(alpha,t)),np.sin(np.multiply(beta,t)))
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore")
		# ignoring warnings that come up with arctan
		phi=np.arctan(np.divide(y,x)) #range is -pi/2 to pi/2

	#make the arrays for before and after
	pre_z=np.arange(-1*pre_d,0,0.0001)
	pre_t=np.divide(pre_z,c)
	pre_x=amp*np.multiply(np.sin(np.multiply(alpha,pre_t)),1) #because B=0,beta=0, cos(...)=1
	pre_y=amp*np.multiply(np.sin(np.multiply(alpha,pre_t)),0) #because B=0,beta=0, sin(...)=0
	pre_phi=np.arctan(np.divide(pre_y,pre_x)) #range is -pi/2 to pi/2
	plt.plot([-1*pre_d,0],[0,0],[amp,amp],color='grey')
	plt.plot([-1*pre_d,0],[0,0],[-1*amp,-1*amp],color='k')

	final_phi=phi[-1]
	post_z=np.arange(d,d+post_d,0.0001)
	post_t=np.divide(post_z,c)
	post_x=amp*np.multiply(np.sin(np.multiply(alpha,post_t)),np.cos(angle_r)) #because B=0,beta=0, cos(...)=1
	post_y=amp*np.multiply(np.sin(np.multiply(alpha,post_t)),np.sin(angle_r)) #because B=0,beta=0, sin(...)=0
	post_phi=np.arctan(np.divide(post_y,post_x))
	plt.plot([d,d+post_d],[amp*np.sin(angle_r),amp*np.sin(angle_r)],[amp*np.cos(angle_r),amp*np.cos(angle_r)],color='grey')
	plt.plot([d,d+post_d],[-1*amp*np.sin(angle_r),-1*amp*np.sin(angle_r)],[-1*amp*np.cos(angle_r),-1*amp*np.cos(angle_r)],color='k')

	#mlab.plot3d(x, y, z, t, tube_radius=0.025, colormap='Spectral')
	#ax.plot(z, y, x, label='parametric curve')
	#make axes
	#plt.plot([0,d],[0,0],[0,0],color='k')
	#plt.plot([0,0],[-amp,amp],[0,0],color='k')
	#plt.plot([0,0],[0,0],[-amp,amp],color='k')
	#plt.plot([d,d],[-amp,amp],[0,0],color='k')
	#plt.plot([d,d],[0,0],[-amp,amp],color='k')

	#plt.plot([d,d],[0,amp*np.sin(angle_r)],[0,amp*np.cos(angle_r)],color='r')
	#plt.plot([0,0],[0,0],[0,amp],color='r')

	#ax.grid(False)
	if g==False:
		#remove grid
		plt.axis('off')

	ax.scatter3D(z,y,x, c=phi, cmap='twilight',marker='.',s=1,vmin=-0.5*math.pi,vmax=0.5*math.pi)
	ax.scatter3D(pre_z,pre_y,pre_x, c=pre_phi, cmap='twilight',marker='.',s=1,vmin=-0.5*math.pi,vmax=0.5*math.pi)
	ax.scatter3D(post_z,post_y,post_x, c=post_phi, cmap='twilight',marker='.',s=1,vmin=-0.5*math.pi,vmax=0.5*math.pi)

	#axisEqual3D(ax)

	lim_max=1.2*np.nanmax([x,y,z])
	lim_min=1.2*np.nanmin([x,y,z])
	#ax.legend()
	ax.set_xlim(np.array([lim_min,lim_max]))
	ax.set_ylim(np.array([lim_min,lim_max]))
	ax.set_zlim(np.array([lim_min,lim_max]))

	ax.view_init(elev=args.theta, azim=args.phi)

	plt.tight_layout()
	if args.savefig!=None:
		plt.savefig(args.savefig,dpi=600,transparent='true')
	plt.show()

if __name__ == "__main__":

	ap = argparse.ArgumentParser(description="Make an EM wave with Faraday Rotation.")
	
	ap.add_argument("-s","--savefig",help='Name of image file to save',
		default=None,type=str)

	ap.add_argument("-a","--amplitude",help='Amplitude of wave (default = 1m)',
		default=1.,type=float)
	ap.add_argument("-w","--wavelength",help='Wavelength (default = 0.2m)',
		default=0.2,type=float)	
	ap.add_argument("-d","--distance",help='Faraday rotating path length (default = 5m)',
		default=5.,type=float)	

	ap.add_argument("-b","--bfield",help='Magnetic field strength (no default)',
		default=None,type=float)
	ap.add_argument("-r","--rotation",help='Desired angle of rotation (no default)',
		default=None,type=float)
	ap.add_argument("-n","--nelec",help='n_e (default = 2 m^-3)',
		default=2.,type=float)
	ap.add_argument("-t","--theta",help='viewing angle theta (default=45 degrees)',
		default=45,type=float)
	ap.add_argument("-p","--phi",help='viewing angle phi (default=-45 degrees)',
		default=-45,type=float)
	ap.add_argument("-g","--grid",help='Show axis grid (default=False)', default=False, action='store_true')
	ap.add_argument("-v","--verbose",help='See more detail', default=False, action='store_true')	

	args = ap.parse_args()

	main(args)

