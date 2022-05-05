from functions import *
import argparse
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import warnings
warnings.filterwarnings('ignore')
import os
from os import path

dpi=300

#possum_freq=1.367490741
#emu_freq=0.94349
beamlogs='/Users/emma/OneDrive/postdoc/beamlogs/'
beamfile='beamlogConvolve-natural.image.restored.i.SB10083.contcube.EMU_2154-55.beam00.txt'

def main():
	beams=np.loadtxt(beamlogs+beamfile)
	print(beams.shape)
if __name__ == "__main__":
	main()