import glob
import argparse
import re
import os

import numpy as np

from astropy.io import fits
from astropy.table import Table

from galfitclass import GalfitClass, SersicComponent

prep_base = './io/prep/'

band_labels = ['F115W','F150W','F200W','F277W','F356W','F410M','F444W']
band_wavelengths = [115.40, 150.00, 198.80, 276.10, 356.80, 408.20, 440.40]
psf_list = [f'../GalfitPy3/io/psf/{band.lower()}_psf.fits' for band in band_labels]
constraint_path = './constraint.txt'

# custom
def calc_zpt(x: int, y: int)-> list:
    """
    Calculate zeropoint
    args:
        x: int
            x position

        y: int
            y position
    """
    CR = {'F115W':[0.914,0.908,0.862,0.813,0.959,0.897,0.995,0.831],
            'F150W':[1.005,0.966,0.932,0.855,1.000,0.968,0.988,0.865],
            'F200W':[0.865,0.889,0.800,0.791,0.873,0.829,0.901,0.812],
            'F277W':[1.107,1.107,1.107,1.107,1.000,1.000,1.000,1.000],
            'F356W':[1.088,1.088,1.088,1.088,1.000,1.000,1.000,1.000],
            'F410M':[1.031,1.031,1.031,1.031,1.042,1.042,1.042,1.042],
            'F444W':[1.036,1.036,1.036,1.036,1.096,1.096,1.096,1.096]
            }
    CR0 = {'F115W':[1,1,1,1,1,1,1,1],
            'F150W':[1,1,1,1,1,1,1,1],
            'F200W':[1,1,1,1,1,1,1,1],
            'F277W':[1,1,1,1,1,1,1,1],
            'F356W':[1,1,1,1,1,1,1,1],
            'F410M':[1,1,1,1,1,1,1,1],
            'F444W':[1,1,1,1,1,1,1,1]
            }
    
    let = x<5400 #letter B for left half, A for right half
    index = let << 2
    if let:
        index = index | ((x>2350) << 1) | (y>2300)
    else:
        index = index | ((x<8140) << 1) | (y<2400)
            
    #for updated version of PHOTMJSRs
    #print('zeropoints are adopted from jwst_0995.pmap\n')
    #return [28.086519392283982-2.5*np.log10(CR0[band][index]) for band in self.bands]

    #for ealier version of PHOTMJSRs
    print('Using magnitude zeropoints for F150W from Boyer et al. (2022), and for the other six filters, we use the zero points from Brammer (2022).')
    return [28.086519392283982-2.5*np.log10(CR[str(band)][index]) for band in CR.keys()]


def main():
    parser = argparse.ArgumentParser(description='Cutout images')
    parser.add_argument('-i','--imgname',  help='image name to be cut', default='none')
    parser.add_argument('-o','--overwrite',  help='overwrite existing files', action='store_true')

    args = parser.parse_args()

    img_name = input("Enter image name: ") if args.imgname == 'none' else args.imgname

    pregalfit_path = prep_base + img_name + '/pregalfit.fits'
    with fits.open(pregalfit_path) as hdul:
        tab_ini = Table(hdul['INIPARAM'].data)

    sample_base = './io/sample/' + img_name + '/'
    sample_dir_list = glob.glob(sample_base + 'EGS-*/')
    sample_dir_list.sort()

    for sample_dir in sample_dir_list:
        seq = re.match(r'.*EGS-(\d+)/', sample_dir).group(1)
        seq = int(seq)
        print('\nProcessing: EGS-' + str(seq))

        img_list = [sample_dir + 'sci_' + band + '.fits' for band in band_labels]
        err_list = [sample_dir + 'err_' + band + '.fits' for band in band_labels]
        bpmask_list = [sample_dir + 'bpmask_' + band + '.fits' for band in band_labels]

        x_c = fits.getheader(img_list[0])['xc']
        y_c = fits.getheader(img_list[0])['yc']
        zpts = calc_zpt(x_c, y_c)
        dict_band = dict(zip(band_labels, zip(band_wavelengths, zpts)))

        # create a galfit object
        gal_save_path = sample_dir + 'output.fits'
        if os.path.exists(gal_save_path):
            print(f'Warning: output {gal_save_path} already exists')
            if args.overwrite:
                print(f'Overwriting output file: {gal_save_path}')
            else:
                print(f'Skipping {seq} as default [overwrite==False]')
                continue

        gal_obj = GalfitClass(dict_band, gal_save_path, plate_scale=(0.03,0.03))

        gal_obj.load_data(img_list, err_list, bpmask_list)
        gal_obj.load_psf(psf_list)
        gal_obj.load_constraint(constraint_path)

        # add a sersic component
        sersic_comp = SersicComponent(gal_obj.nbands, skip=False)
        sersic_comp.config_x(1, [gal_obj.image_size[0]/2+1])
        sersic_comp.config_y(1, [gal_obj.image_size[1]/2+1])

        row_ini = tab_ini[tab_ini['EGS-ID']==seq]
        mag_ini_list = [row_ini[f'KronPhot_{band}_mag'][0] for band in band_labels]
        sersic_comp.config_mag(7, mag_ini_list)
        sersic_comp.config_n(2, [2.0]*2)
        sersic_comp.config_re(2, [5.0]*2)
        sersic_comp.config_q(1, [row_ini['eccentricity'][0]])
        sersic_comp.config_pa(1, [row_ini['orientation'][0]])

        gal_obj.add_component(sersic_comp)

        # generate feedme file
        feedme = gal_obj.genstr_feedme()
        feedme_path = sample_dir + f'{seq}.galfit'
        with open(feedme_path, 'w') as f:
            f.write(feedme)

        gal_obj.run(feedme_path)

if __name__ == '__main__':
    main()