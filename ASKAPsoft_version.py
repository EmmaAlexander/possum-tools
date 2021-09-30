import glob
from astropy.io import fits
from functions import fitsopen
import numpy as np


directory='/Volumes/TARDIS/Work/askap/fullfields/'
files=glob.glob(directory+'*.fits')

SBs= []
versions=[]
filenames=[]

for file in files:
	#print(file)
	data,header=fitsopen(file)
	filename=file.split('/')[-1]
	filenames.append(filename)
	try:
		SBs.append(header['SBID'])
	except:
		print("No SB found for file {}".format(file))
		SBs.append('-')
	try:
		version = header['HISTORY'][0].split()[-1]
		versions.append(version)
	except:
		print('No file history found')
		versions.append('-')

SBs=np.asarray(SBs)
versions=np.asarray(versions)
filenames=np.asarray(filenames)

together=np.vstack((SBs,versions,filenames)).T

np.savetxt('ASKAPsoftversions.txt',together,fmt='%s')                                                                                                                                                                  