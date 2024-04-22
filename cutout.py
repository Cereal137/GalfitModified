import os
import shutil
import argparse
from typing import Union
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, stretch
from astropy.table import Table
from photutils.segmentation import SegmentationImage

from utils.segmentcut import SegmCut

prep_base = './io/prep/'
sample_base_base = './io/sample/'

def cutout_worker(id: str, x_c: float, y_c: float, label: Union[int, list], size: int, 
                  tab_img: Table, segm: SegmentationImage, sample_base: str, 
                  overwrite=False, max_size: int = 270) -> None:

    if isinstance(label, list):
        pass
    else: # convert to list
        label = [label]
    if 0 in label: # avoid label 0
        raise ValueError('Label 0 is not allowed for a target')
    
    sample_dir = sample_base + str(id) + '/'
    if os.path.exists(sample_dir):
        if overwrite:
            shutil.rmtree(sample_dir)
            os.mkdir(sample_dir)
        else:
            raise ValueError('Directory already exists: ' + sample_dir)
    else:
        os.mkdir(sample_dir)

    segm_cut = Cutout2D(segm.data, (x_c, y_c), (size, size)).data # cutout(segm.data, x_c, y_c, size)   
    areas = segm.areas
    sc = SegmCut(segm_cut, areas, label, verbose=False)
    valid = True

    fig, ax = plt.subplots(len(tab_img),3,figsize=(8,20))

    raw_path = sample_dir + f'raw_{id}.fits'
    hdr_raw = fits.Header()
    hdr_raw['EXTNAME'] = 'SEGMAP'
    fits.writeto(raw_path, segm_cut, header=hdr_raw)

    for i,bandrow in enumerate(tab_img):
        band = bandrow['band']
        img_band_path = bandrow['img_path']
        with fits.open(img_band_path) as hdul:
            w = WCS(hdul['SCI'].header)
            cut = Cutout2D(hdul['SCI'].data, (x_c, y_c), (size, size), wcs=w)
            sci_cut, w_cut = cut.data, cut.wcs
            err_cut = Cutout2D(hdul['ERR'].data, (x_c, y_c), (size, size)).data # cutout(hdul['ERR'].data, x_c, y_c, size)
        try:
            scimap, errmap, bpmask, bkg1d = sc.gen_cutout(sci_cut, err_cut, hist=True)
        except ValueError:
            print('Error: ' + sample_dir + ' is not a valid cutout')
            valid = False

        if bpmask[size//2,size//2]==1:
            print('Error: ' + sample_dir + ' have bad pixel in the center, which may suggest a bad cutout in ' + band)
            valid = False

        if size > max_size:
            print('Error: ' + sample_dir + f' exceeds the maximum size: {max_size}, which may suggest a bad segmentation')
            valid = False
        
        if not valid:
            shutil.rmtree(sample_dir)
            plt.close(fig)
            return 0
        
        sci_path = sample_dir + 'sci_' + band + '.fits'
        hdr_sci = w_cut.to_header()
        hdr_sci['EXTNAME'] = 'SCIENCE'
        hdr_sci['xc'] = x_c
        hdr_sci['yc'] = y_c
        fits.writeto(sci_path, scimap, header=hdr_sci)

        err_path = sample_dir + 'err_' + band + '.fits'
        hdr_sig = fits.Header()
        hdr_sig['EXTNAME'] = 'SIGMA'
        hdr_sig['xc'] = x_c
        hdr_sig['yc'] = y_c
        fits.writeto(err_path, errmap, header=hdr_sig)

        hdr_err = fits.Header()
        hdr_err['EXTNAME'] = 'ERR_'+band
        hdr_err['xc'] = x_c
        hdr_err['yc'] = y_c
        hdr_sci['EXTNAME'] = 'SCI_'+band
        fits.append(raw_path, sci_cut, header=hdr_sci)
        fits.append(raw_path, err_cut, header=hdr_err)
            
        bpmask_path = sample_dir + 'bpmask_' + band + '.fits'
        hdr_bpmask = fits.Header()
        hdr_bpmask['EXTNAME'] = 'BPMASK'
        hdr_bpmask['xc'] = x_c
        hdr_bpmask['yc'] = y_c
        fits.writeto(bpmask_path, bpmask.astype(int), header=hdr_bpmask)

        ax[i,0].imshow(scimap, origin='lower', cmap='viridis', norm=ImageNormalize(stretch=stretch.LogStretch()))
        ax[i,0].set_title('sci')
        ax[i,1].imshow(bpmask, origin='lower', cmap='gray', vmin=0, vmax=1)
        ax[i,1].set_title('bpmask')
        ax[i,2].hist(bkg1d, bins=50, histtype='step', color='k')
        ax[i,2].set_title('background')
        
        
    plt.tight_layout()
    plt.savefig(sample_dir + f'mask_{id}.png')
    plt.close(fig)

    return 1
        

def main():
    parser = argparse.ArgumentParser(description='Cutout images')
    parser.add_argument('-i','--imgname',  help='image name to be cut', default='none')
    parser.add_argument('-o','--overwrite',  help='overwrite existing files', action='store_true')
    parser.add_argument('-n','--nproc',  help='number of processes', type=int, default=1)

    args = parser.parse_args()

    img_name = input("Enter image name: ") if args.imgname == 'none' else args.imgname

    # load info from pregalfit.fits
    prep_dir = prep_base + img_name + '/'
    pregalfit_path = prep_dir + 'pregalfit.fits'
    with fits.open(pregalfit_path) as hdul:
        tab_img = Table(hdul['IMGLIST'].data)
        tab_ini = Table(hdul['INIPARAM'].data)
        segm = SegmentationImage(hdul['SEGMAP'].data)
    areas = segm.areas

    sample_base = sample_base_base + img_name + '/'
    if args.overwrite:
        print('Overwriting existing directory: ' + sample_base)
        shutil.rmtree(sample_base)
        os.mkdir(sample_base)

    nproc = args.nproc if args.nproc <= os.cpu_count() else os.cpu_count()
    if nproc==1:
        for row in tab_ini:
            id = str(row['ID'])
            x_c = row['xcentroid']
            y_c = row['ycentroid']
            label = int(row['label'])
            size = 2*int(7*row['kron_radius'])
            cutout_worker(id, x_c, y_c, label, size,
                          tab_img, segm, sample_base, args.overwrite)
    else:
        from multiprocessing import Pool
        with Pool(nproc) as p:
            p.starmap(cutout_worker, [(str(id), x_c, y_c, label, 2*int(7*rad),
                                        tab_img, segm, sample_base) 
                                        for (id, x_c, y_c, label, rad) in tab_ini[['ID','xcentroid','ycentroid','label','kron_radius']]])
    print('Cutout Successfully Done!')

if __name__ == '__main__':
    main()