import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

from photutils.segmentation import SourceCatalog, SegmentationImage

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--imgname', help='image name to be checked', default='none')
    parser.add_argument('-r','--reference', help='reference catalog path', default='none')

    args = parser.parse_args()
    img_name = input("Enter image name: ") if args.imgname == 'none' else args.imgname
    ref_path = input("Enter reference catalog path: ") if args.reference == 'none' else args.reference

    with fits.open(f'io/image/{img_name}/merged.fits') as hdul:
        w = WCS(hdul[0].header)
        image_merged_conv = hdul['SCI_BKSUB_CONV'].data
        segm_deblend = SegmentationImage(hdul['SEGMAP'].data)

    cat = SourceCatalog(image_merged_conv, segm_deblend, wcs=w, convolved_data=image_merged_conv)
    skycoord_list = cat.sky_centroid_icrs
    label_list = cat.labels

    tab_ref = Table.read(ref_path)
    print(tab_ref.colnames)
    key_ra = input("Enter RA key: ")
    key_dec = input("Enter Dec key: ")
    key_id = input("Enter ID key: ")
    ref_list = SkyCoord(ra=tab_ref[key_ra]*u.degree, dec=tab_ref[key_dec]*u.degree, frame='fk5').transform_to('icrs')

    max_sep = 1. # arcsec
    idx1, idx2, _, _ = skycoord_list.search_around_sky(ref_list, max_sep*u.arcsec)
    idx1_unique, pos_unique = np.unique(idx1, return_index=True)
    tab_ref_matched = tab_ref[idx1_unique]
    ra_matched = tab_ref_matched[key_ra]
    dec_matched = tab_ref_matched[key_dec]
    skypos_matched = SkyCoord(ra=ra_matched*u.degree, dec=dec_matched*u.degree, frame='fk5').transform_to('icrs')
    pixpos_matched = skypos_matched.to_pixel(w)
    label_matched = segm_deblend.data[pixpos_matched[1].astype(int), pixpos_matched[0].astype(int)]

    match_dir = f'./io/match/{img_name}/'
    # crossmatch check in csv format, convienient for human check
    tab_ref_matched.keep_columns([key_ra, key_dec, key_id])
    tab_ref_matched.rename_column(key_ra, 'RA')
    tab_ref_matched.rename_column(key_dec, 'Dec')
    tab_ref_matched.rename_column(key_id, 'ID')
    
    tab_ref_matched.add_column(label_matched, name='label')
    tab_ref_matched.write(match_dir + f'tab_ref_matched_{img_name}.csv', format='ascii.csv', overwrite=True)
    print(f'Crossmatch result has been successfully saved in {match_dir}tab_ref_matched_{img_name}.csv')


if __name__ == '__main__':
    main()