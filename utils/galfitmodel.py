import os
import subprocess
from secrets import token_hex

import numpy as np

from astropy.io import fits

def sersic2d(x: float, y: float, xaxis: int, yaxis: int, ar: float, pa: float, re: float, n: float, mag: float, mag_zpt: float=28.086519392283982) -> np.ndarray:
    """
    2D Sersic profile

    Parameters
    ----------
    x : float
        x-centroid
    y : float
        y-centroid
    xaxis : int
        full image x-axis
    yaxis : int
        full image y-axis
    ar : float
        axis ratio
    pa : float
        position angle
    re : float
        effective radius
    n : float
        Sersic index
    mag : float
        total AB magnitude
    """
    token = token_hex(16)
    out_path = f'temp_{token}.fits'
    zero_point = mag_zpt

    fdme = _galfitconfig(xaxis,yaxis,zero_point,out_path) + _modelconfig(x,y,ar,pa,re,n,mag)
    fdme_path = out_path.replace('.fits','.galfit')
    with open(fdme_path, 'w') as f:
        f.write(fdme)

    #os.chdir(os.getcwd())
    cmd = f"galfit {fdme_path}"
    child = subprocess.call(['/bin/zsh',"-i", "-c", cmd], stdout=subprocess.DEVNULL)

    with fits.open(out_path) as hdul:
        img = np.float64(hdul[0].data)
    os.remove(fdme_path)
    os.remove(out_path)

    return img

def _modelconfig(x,y,ar,pa,re,n,mag):
    return f"""0) sersic             # Object type
 1) {x} {y}  1 1    # position x, y        [pixel]
 3) {mag}      1       # total magnitude    
 4) {re}       1       #     R_e              [Pixels]
 5) {n}       1       # Sersic exponent (deVauc=4, expdisk=1)  
 9) {ar}  1       # axis ratio (b/a)   
10) {pa}      1       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
"""

def _galfitconfig(xaxis,yaxis,zero_point,out_path) -> str:
    return f"""A)                     # Input data image (FITS file)
B) {out_path}       # Output data image block
C) none                # Sigma image name (made from data if blank or "none") 
D)                     # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none                # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1 {xaxis} 1 {yaxis} # Image region to fit (xmin xmax ymin ymax)
I) 100    100          # Size of the convolution box (x y)
J) {zero_point}   # Magnitude photometric zeropoint 
K) 0.030  0.030        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 1                   # Options: 0=normal run; 1,2=make model/imgblock & quit

"""
#test
if __name__ == '__main__':
    print(_galfitconfig(10,10,'test.fits') + _modelconfig(10,10,0.5,0,10,1,20))