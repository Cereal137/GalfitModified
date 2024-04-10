import os
import shutil
import argparse
import warnings
warnings.filterwarnings('ignore')

from tqdm import tqdm
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from photutils.segmentation import SegmentationImage

from utils.segmentcut import SegmCut, cutout

prep_base = './io/prep/'
sample_base_base = './io/sample/'

def main():
    parser = argparse.ArgumentParser(description='Cutout images')
    parser.add_argument('-i','--imgname',  help='image name to be cut', default='none')
    parser.add_argument('-o','--overwrite',  help='overwrite existing files', action='store_true')

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

    for row in tab_ini:
        id = row['ID']
        sample_dir = sample_base + str(id) + '/'
        if os.path.exists(sample_dir):
            print('Warning: ' + sample_dir + ' exists')
            if not args.overwrite:
                print('Skip by default [Overwrite==False]')
                continue
        else:
            os.mkdir(sample_dir)

        print('Cutting out ' + sample_dir + ' ...', end=' \r', flush=True)
        x_c = row['xcentroid']
        y_c = row['ycentroid']
        label = int(row['label'])
        """
        srcarea = areas[label-1]
        size = int(np.sqrt(srcarea/3.14)*3.7)
        """
        size = int(7*row['kron_radius'])
        #print('Size: ' + str(size))

        segm_cut = Cutout2D(segm.data, (x_c, y_c), (size*2, size*2)).data # cutout(segm.data, x_c, y_c, size)
        try:
            sc = SegmCut(segm_cut, size, areas, [label], verbose=False)
        except ValueError:
            print('Failed to cut out ' + sample_dir + ' ...\n')
            continue


        fig, ax = plt.subplots(len(tab_img),2,figsize=(8,20))

        for i,bandrow in enumerate(tab_img):
            band = bandrow['band']
            img_band_path = bandrow['img_path']
            with fits.open(img_band_path) as hdul:
                w = WCS(hdul['SCI_BKSUB'].header)
                cut = Cutout2D(hdul['SCI_BKSUB'].data, (x_c, y_c), (size*2, size*2), wcs=w)
                sci_cut, w_cut = cut.data, cut.wcs
                #sci_cut = cutout(hdul['SCI_BKSUB'].data, x_c, y_c, size)
                err_cut = Cutout2D(hdul['ERR'].data, (x_c, y_c), (size*2, size*2)).data # cutout(hdul['ERR'].data, x_c, y_c, size)
            scimap, errmap, bpmask = sc.gen_cutout(sci_cut, err_cut)
            
            sci_path = sample_dir + 'sci_' + band + '.fits'
            hdr_sci = w_cut.to_header()
            hdr_sci['EXTNAME'] = 'SCI_BKSUB'
            hdr_sci['xc'] = x_c
            hdr_sci['yc'] = y_c
            fits.writeto(sci_path, scimap, header=hdr_sci)

            err_path = sample_dir + 'err_' + band + '.fits'
            hdr_err = fits.Header()
            hdr_err['EXTNAME'] = 'ERR'
            hdr_err['xc'] = x_c
            hdr_err['yc'] = y_c
            fits.writeto(err_path, errmap, header=hdr_err)

            err0_path = sample_dir + 'err0_' + band + '.fits'
            hdr_err0 = fits.Header()
            hdr_err0['EXTNAME'] = 'ERR0'
            hdr_err0['xc'] = x_c
            hdr_err0['yc'] = y_c
            fits.writeto(err0_path, err_cut, header=hdr_err0)
            
            bpmask_path = sample_dir + 'bpmask_' + band + '.fits'
            hdr_bpmask = fits.Header()
            hdr_bpmask['EXTNAME'] = 'BPMASK'
            hdr_bpmask['xc'] = x_c
            hdr_bpmask['yc'] = y_c
            fits.writeto(bpmask_path, bpmask.astype(int), header=hdr_bpmask)

            ax[i,0].imshow(scimap, origin='lower', cmap='gray')
            ax[i,0].set_title('sci_' + band)

            ax[i,1].imshow(bpmask, origin='lower', cmap='gray', vmin=0, vmax=1)
            ax[i,1].set_title('bpmask')
        plt.tight_layout()
        plt.savefig(sample_dir + f'mask_{id}.png')
        
    print('Cutout Successfully Done!')

if __name__ == '__main__':
    main()