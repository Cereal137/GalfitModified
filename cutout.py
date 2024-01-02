import os
import shutil
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from photutils.segmentation import SegmentationImage

prep_base = './io/prep/'
sample_base_base = './io/sample/'

def cutout(img, x, y, size):
    """
    Cutout image
    """
    x, y = int(x), int(y)
    return img[y-size:y+size, x-size:x+size]

def dilate_mask(mask, size=5):
    """
    Dilate mask
    """
    from scipy.ndimage import binary_dilation
    return binary_dilation(mask, structure=np.ones((size, size)))

def estimate_local_background(data1d):
    """
    Estimate local background
    """
    return sigma_clipped_stats(data1d, sigma=3.0)


def main():
    parser = argparse.ArgumentParser(description='Cutout images')
    parser.add_argument('-i','--imgname',  help='image name to be cut', default='none')

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
    for row in tab_ini:
        seq = row['EGS-ID']
        sample_dir = sample_base + 'EGS-' + str(seq) + '/'
        if os.path.exists(sample_dir):
            print('Warning: ' + sample_dir + ' exists')
            print('Remove existing directory: ' + sample_dir)
            shutil.rmtree(sample_dir)
        os.mkdir(sample_dir)

        x_c = row['xcentroid']
        y_c = row['ycentroid']
        label = row['label']
        srcarea = areas[label-1]
        size = int(np.sqrt(srcarea/3.14)*3.7)
        #print('Size: ' + str(size))

        segm_cut = cutout(segm.data, x_c, y_c, size)
        # growing size mask
        mask_all = np.zeros_like(segm_cut, dtype=bool) # mask all segmentations
        for label_mask in np.unique(segm_cut):
            if label_mask == 0:
                continue
                
            mask_label = (segm_cut == label_mask)
            dilation_size = int(min(8,areas[label_mask-1]/3.14))
            #print('Masking label: ' + str(label_mask)+ ' with dilation size: ' + str(dilation_size))
            mask_label_dilated = dilate_mask(mask_label, dilation_size)
            mask_all = mask_all | mask_label_dilated # expand mask_all
            if label_mask == label:
                mask_source = mask_label_dilated

        fig, ax = plt.subplots(len(tab_img),2,figsize=(8,20))

        for i,bandrow in enumerate(tab_img):
            band = bandrow['band']
            img_band_path = bandrow['img_path']
            #print(f'{band}: Processing with ' +img_band_path)
            with fits.open(img_band_path) as hdul:
                sci = hdul['SCI'].data
                err = hdul['ERR'].data

            sci_cut = cutout(sci, x_c, y_c, size)
            err_cut = cutout(err, x_c, y_c, size)
            err_cut[err_cut == np.inf] = np.nan
            
            mask_nan = np.isnan(err_cut)
            sci_cut[mask_nan] = np.nan
            
            bkg_mean, _, bkg_std = estimate_local_background(sci_cut[~mask_all].flatten())
            ax[i,0].hist(sci_cut[~mask_all].flatten(), bins=np.linspace(-.05,.05,100), label='Background', alpha=.5)
            ax[i,0].set_title(f'{band}: std={bkg_std:.3f}')

            #print('Subtracting with mean background: ' + str(bkg_mean))
            sci_bksub = sci_cut - bkg_mean
            sci_path = sample_dir + 'sci_' + band + '.fits'
            hdr_sci = fits.Header()
            hdr_sci['EXTNAME'] = 'SCI_BKSUB'
            hdr_sci['xc'] = x_c
            hdr_sci['yc'] = y_c
            fits.writeto(sci_path, sci_bksub, header=hdr_sci)

            err_mean = np.nanmean(err_cut[~mask_all])
            #print('Scaling error to the level of background std: ' + str(bkg_std))
            scale_factor = bkg_std/err_mean
            err_scaled = np.choose(mask_all, (err_cut*scale_factor,err_cut))
            err_path = sample_dir + 'err_' + band + '.fits'
            hdr_err = fits.Header()
            hdr_err['EXTNAME'] = 'ERR'
            hdr_err['xc'] = x_c
            hdr_err['yc'] = y_c
            fits.writeto(err_path, err_scaled, header=hdr_err)

            bpmask = ((mask_all | (sci_cut> 7.*bkg_std)) & ~mask_source) | mask_nan # (sci_cut> 7.*bkg_std) is for masking contaminations that are not recognized by segmentation
            bpmask_path = sample_dir + 'bpmask_' + band + '.fits'
            hdr_bpmask = fits.Header()
            hdr_bpmask['EXTNAME'] = 'BPMASK'
            hdr_bpmask['xc'] = x_c
            hdr_bpmask['yc'] = y_c
            fits.writeto(bpmask_path, bpmask.astype(int), header=hdr_bpmask)

            ax[i,1].imshow(bpmask, origin='lower', cmap='gray', vmin=0, vmax=1)
            ax[i,1].set_title('bpmask')
        plt.tight_layout()
        plt.savefig(sample_dir + f'mask_EGS-{seq}.png')
    pass

if __name__ == '__main__':
    main()