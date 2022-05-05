#converts a list of RA and Dec strings into decimal degrees
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import argparse


def main(args):

	sources=np.loadtxt(args.infile,dtype='str')
	ra=sources[:,0]
	dec=sources[:,1]
	coords=SkyCoord(ra,dec)
	out=np.stack([coords.ra.deg,coords.dec.deg],axis=1)
	np.savetxt(args.outfile,out,fmt='%3.5f')

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument('infile',type=str)
	ap.add_argument('outfile',type=str)
	args = ap.parse_args()
	main(args)
