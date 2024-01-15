from astropy.io import fits
import numpy as np

def merge_image(img_list: list[str]):
    """
    Merge images
    
    Parameters
    ----------
    img_list : list
        List of image paths str
    """
    print('merging with:')
    img_cube = []
    wht_cube = []
    for img_path in img_list:
        print(img_path)
        with fits.open(img_path) as hdul:
            img = hdul['SCI_BKSUB'].data
            err = hdul['ERR'].data
        img_cube.append(img)
        wht_cube.append(np.reciprocal(err))
    img_cube = np.array(img_cube)
    wht_cube = np.array(wht_cube)
    img_merged_sum = np.nansum(img_cube*wht_cube, axis=0)
    wht_merged_sum = np.nansum(wht_cube, axis=0)
    mask_final = (wht_merged_sum == 0)
    wht_merged_sum_nan = np.choose(mask_final,(wht_merged_sum, np.nan))
    img_merged = img_merged_sum/wht_merged_sum_nan

    return img_merged