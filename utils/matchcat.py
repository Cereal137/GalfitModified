import os
import glob
import matplotlib.pyplot as plt
import numpy as np

import warnings
warnings.filterwarnings('ignore')

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from photutils.segmentation import SourceCatalog, SegmentationImage

match_base = './io/match/'
image_base = './io/image/'

class MatchCat:
    def __init__(self, img_name :str, match_dir :str = None, verbose :bool = False):
        """
        Class to store matched catalog and images for a given image
        
        Parameters
        ----------
        img_name : str
            Name of the image
        """
        self.img_name = img_name
        self.match_dir = match_base + img_name + '/' if match_dir is None else match_dir
        if verbose:
            print(f'Image name: {self.img_name}')
            print(f'Match directory: {self.match_dir}')
        img_dir = image_base + img_name +'/'
        if not os.path.exists(img_dir):
            raise FileNotFoundError(f'Image directory {img_dir} does not exist')

        self.img_path = img_dir + 'merged.fits'
        self.tab_path = self.match_dir + f'tab_ref_matched_{img_name}.csv'
        if not os.path.exists(self.tab_path):
            raise FileNotFoundError(f'Matched catalog {self.tab_path} does not exist')

    def _get_hst(self, hst_path :str):
        """
        Get HST image and WCS

        Parameters
        ----------
        hst_path : str
            Path to the HST image
        """
        with fits.open(hst_path) as hdul:
            #hdul.info()
            hst_data = hdul[1].data
            hst_wcs = WCS(hdul[1].header)
        return hst_data, hst_wcs
    
    def _get_jwst(self):
        """
        Get JWST convolved image, WCS and segmentation map
        """
        with fits.open(self.img_path) as hdul:
            image_merged_conv = hdul['SCI_BKSUB_CONV'].data
            segmap = hdul['SEGMAP'].data
            img_wcs = WCS(hdul[0].header)
        return image_merged_conv, img_wcs, segmap
    
    def interactive_match(self, hst_path :str):
        """
        Interactive matching of objects in HST and jwst images
        
        Parameters
        ----------
        img_name : str
            Name of the image
        hst_path : str
            Path to the HST image
        """

        tab_ref_matched = Table.read(self.tab_path, format='ascii.csv')

        hst_data, hst_wcs = self._get_hst(hst_path)
        image_merged_conv, img_wcs, segmap = self._get_jwst()

        segm = SegmentationImage(segmap)
        srccat = SourceCatalog(image_merged_conv, segm)
        xcen, ycen = srccat.xcentroid, srccat.ycentroid

        i = 0
        while i <len(tab_ref_matched):
            print(f'Object {i+1} of {len(tab_ref_matched)}')

            hst_x, hst_y = hst_wcs.all_world2pix(tab_ref_matched['RA'][i], tab_ref_matched['Dec'][i], 0)
            hst_x, hst_y = int(hst_x), int(hst_y)

            img_x, img_y = img_wcs.all_world2pix(tab_ref_matched['RA'][i], tab_ref_matched['Dec'][i], 0)
            img_x, img_y = int(img_x), int(img_y)

            hst_size = 50
            img_size = 50
            segm_cutout = segmap[img_y-img_size:img_y+img_size, img_x-img_size:img_x+img_size]
            labels = np.unique(segm_cutout)

            fig, ax = plt.subplots(1, 2, figsize=(10, 5))
            im_hst = ax[0].imshow(hst_data[hst_y-hst_size:hst_y+hst_size, hst_x-hst_size:hst_x+hst_size], cmap='gray')
            ax[0].set_title('HST')
            plt.colorbar(im_hst, ax=ax[0])
            im_jwst = ax[1].imshow(image_merged_conv[img_y-img_size:img_y+img_size, img_x-img_size:img_x+img_size], cmap='gray')
            ax[1].set_title('JWST')
            plt.colorbar(im_jwst, ax=ax[1])
            for label in labels[1:]:
                ylabel, xlabel = xcen[label-1], ycen[label-1]
                ax[1].text(ylabel-img_x+img_size, xlabel-img_y+img_size, label, color='red')
            fig.suptitle(f'ID: {tab_ref_matched["ID"][i]}')
            plt.show(block=False)

            goback = input('Go back? (n): ')
            if goback == 'y':
                i -= 2 if i > 0 else 1
                
            i += 1
            plt.close()

    def get_image(self, hst_path: str, i: int, img_size: int = 50):
        """
        Get image of object i in HST and JWST images

        Parameters
        ----------
        hst_path : str
            Path to the HST image
        i : int
            Index of the object
        """

        hst_data, hst_wcs = self._get_hst(hst_path)
        image_merged_conv, img_wcs, segmap = self._get_jwst()

        tab_ref_matched = Table.read(self.tab_path, format='ascii.csv')

        print(f'Object {i+1} of {len(tab_ref_matched)}: {tab_ref_matched["ID"][i]}')

        hst_x, hst_y = hst_wcs.all_world2pix(tab_ref_matched['RA'][i], tab_ref_matched['Dec'][i], 0)
        hst_x, hst_y = int(hst_x), int(hst_y)

        img_x, img_y = img_wcs.all_world2pix(tab_ref_matched['RA'][i], tab_ref_matched['Dec'][i], 0)
        img_x, img_y = int(img_x), int(img_y)

        hst_size = img_size # same size for HST and JWST
        
        hst_cutout = hst_data[hst_y-hst_size:hst_y+hst_size, hst_x-hst_size:hst_x+hst_size]
        jwst_cutout = image_merged_conv[img_y-img_size:img_y+img_size, img_x-img_size:img_x+img_size]
        segm_cutout = segmap[img_y-img_size:img_y+img_size, img_x-img_size:img_x+img_size]

        return hst_cutout, jwst_cutout, segm_cutout#, hst_x, hst_y, img_x, img_y

    def nmatch(self) -> int:
        """
        Get number of matched objects
        """
        tab_ref_matched = Table.read(self.tab_path, format='ascii.csv')
        return len(tab_ref_matched)