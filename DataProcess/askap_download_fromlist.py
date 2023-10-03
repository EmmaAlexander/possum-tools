import numpy as np
import os
import wget

dataloc='/share/nas2/ela/ASKAP/data/'
# the file "to download" comes from the CASDA data access portal
filelist=np.genfromtxt('todownload.txt',dtype='str')

for url in filelist:
    wget.download(url,out=dataloc)
    if len(url)<150:
        #checksum
        fname=url.split('/')[-1]
    elif "islands" in url:
        fname="selavy-image"+url.split("filename")[1].split("image")[1].split("xml")[0]+"xml"
        print(fname) 
    else:
        fname="image"+url.split("filename")[1].split("image")[1].split("fits")[0]+"fits"
    fname=dataloc+fname
    print(fname)
    os.system(("wget -O {} '{}'".format(fname,url)))