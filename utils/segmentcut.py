import numpy as np
from dataclasses import dataclass


def cutout(img: np.ndarray, x: int, y: int, size: int):
    """
    Cutout image

    Parameters
    ----------
    img: np.ndarray
        Image array
    x: int
        x coordinate of the center of the cutout    
    y: int
        y coordinate of the center of the cutout
    size: int
        size of the cutout

    Returns
    -------
    cutout: np.ndarray
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
    from astropy.stats import sigma_clipped_stats
    return sigma_clipped_stats(data1d, sigma=3.0)

@dataclass
class SegmCut:
    segm_cut: np.ndarray
    size: int
    label_areas: np.ndarray
    labels_unmask: np.ndarray
    verbose: bool = False

    def __post_init__(self):
        self.labels = self._labels()
        self.mask_all = self._mask_all()
        self.mask_target = self._mask_target()

    def _labels(self):
        """
        Generate a list of labels for all sources.
        """
        labels = np.unique(self.segm_cut)
        return labels[labels!=0]

    def _mask_all(self):
        """
        Generate a mask for all sources.
        The mask is only determined by the input segmentation image.
        """
        mask_all = np.zeros_like(self.segm_cut, dtype=bool)
        for label_iter in self.labels:
            mask_iter = (self.segm_cut==label_iter)
            dilation_size = int(min(8,self.label_areas[label_iter-1]/np.pi))
            if self.verbose:
                print(f'dilation size for label {label_iter}: ', dilation_size)

            mask_iter_dilated = dilate_mask(mask_iter, size=dilation_size)
            mask_all = mask_all | mask_iter_dilated

        if self.verbose:
            print('total iteration: ', len(self.labels))
        return mask_all
    
    def _mask_target(self):

        mask_target = np.zeros_like(self.segm_cut, dtype=bool)
        for label_iter in self.labels_unmask:
            mask_iter = (self.segm_cut==label_iter)
            dilation_size = int(min(8,self.label_areas[label_iter-1]/np.pi))
            if self.verbose:
                print(f'dilation size for targeted label {label_iter}: ', dilation_size)

            mask_iter_dilated = dilate_mask(mask_iter, size=dilation_size)
            mask_target = mask_iter_dilated | mask_target

        return mask_target


    def _bpmask(self, mask_init: np.ndarray):
        """
        Generate a bad pixel mask.
        mask_init: initial mask (specific to every cutout) which exclude the artifacts in science cutout that are not recognized by segmentation,
        np.nan in the science cutout and np.inf in the error cutout.
        """
        return (self.mask_all & ~self.mask_target) | mask_init
        # bpmask = ((mask_all | (sci_cut> 7.*bkg_std)) & ~mask_source) | mask_nan # (sci_cut> 7.*bkg_std) is for masking contaminations that are not recognized by segmentation


    def _errmap(self, err_cut: np.ndarray, bkg_std: float):
        """
        Generate error map
        """
        if bkg_std!=bkg_std:
            raise ValueError('bkg_std is nan')
        
        mask_err = np.isinf(err_cut) | np.isnan(err_cut)
        err_mean = np.nanmean(err_cut[~self.mask_all & ~mask_err])
        scale_factor = bkg_std/err_mean
        if self.verbose:
            print('scaling error map to level of background std: ', bkg_std)
            print('scaling factor: ', scale_factor)

        err_scaled = np.choose(self.mask_all, [err_cut*scale_factor, err_cut])
        return err_scaled
    
    
    def gen_cutout(self, sci_cut: np.ndarray, err_cut: np.ndarray, hist: bool = False):
        """
        Generate cutout
        """
        mask_err = np.isinf(err_cut) | np.isnan(err_cut)
        bkg2d = sci_cut[~self.mask_all & ~mask_err]
        bkg_mean, _, bkg_std = estimate_local_background(bkg2d)

        mask_init = ((sci_cut> 7.*bkg_std) & ~self.mask_all) | np.isnan(sci_cut) | mask_err
        # (sci_cut> 7.*bkg_std) is for masking contaminations that are not recognized by segmentation
        # np.isinf(err_cut) is for masking pixels with infinite error
        # np.isnan(sci_cut) is for masking pixels with np.nan value

        bpmask = self._bpmask(mask_init)
        errmap = self._errmap(err_cut, bkg_std)
        scimap = sci_cut - bkg_mean

        if not hist:
            return scimap, errmap, bpmask
        else:
            if self.verbose:
                print('returning 1d data of masked science cutout on the rear')
            return scimap, errmap, bpmask, bkg2d.flatten()-bkg_mean