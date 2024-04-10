import glob
import argparse

from tqdm import tqdm
from astropy.io import fits
from astropy.table import Table

from galfitclass import GalfitClass, SersicComponent
from utils.zeropoint import calc_zpt

prep_base = './io/prep/'

band_labels = ['F115W','F150W','F200W','F277W','F356W','F410M','F444W']
band_wavelengths = [1.1540, 1.5000, 1.9880, 2.7610, 3.5680, 4.0820, 4.4040]
psf_list = [f'./io/psf/{band.lower()}_psf.fits' for band in band_labels]
constraint_path = './constraint.txt'

def genfeedme_sample(arguments :tuple[str, Table]) -> None:
    sample_dir, tab_ini = arguments
    id = sample_dir.split('/')[-2]

    img_list = [sample_dir + 'sci_' + band + '.fits' for band in band_labels]
    err_list = [sample_dir + 'err_' + band + '.fits' for band in band_labels]
    bpmask_list = [sample_dir + 'bpmask_' + band + '.fits' for band in band_labels]

    x_c = fits.getheader(img_list[0])['xc']
    y_c = fits.getheader(img_list[0])['yc']
    zpts = calc_zpt(x_c, y_c, corr_dict='CR', verbose=False)
    dict_band = dict(zip(band_labels, zip(band_wavelengths, zpts)))

    # create a galfit object
    gal_save_path = sample_dir + 'output.fits'
    gal_obj = GalfitClass(dict_band, gal_save_path, plate_scale=(0.03,0.03))

    gal_obj.load_data(img_list, err_list, bpmask_list, verbose=False)
    gal_obj.load_psf(psf_list, verbose=False)
    gal_obj.load_constraint(constraint_path, verbose=False)

    # add a sersic component
    sersic_comp = SersicComponent(gal_obj.nbands, skip=False)
    sersic_comp.config_x(7, [gal_obj.image_size[0]/2+1]*7) # 7 but restrain its offset
    sersic_comp.config_y(7, [gal_obj.image_size[1]/2+1]*7)

    row_ini = tab_ini[tab_ini['ID']==int(id)][0]
    mag_ini_list = [row_ini[f'KronPhot_{band}_mag'] for band in band_labels]
    sersic_comp.config_mag(7, mag_ini_list)
    sersic_comp.config_n(2, [2.0]*2) # 2 can be wrong
    sersic_comp.config_re(2, [5.0]*2)  # 2 can be wrong
    sersic_comp.config_q(1, [row_ini['eccentricity']])
    sersic_comp.config_pa(1, [row_ini['orientation']])

    gal_obj.add_component(sersic_comp)

    # generate feedme file
    feedme = gal_obj.genstr_feedme()
    feedme_path = sample_dir + f'{id}.galfit'
    with open(feedme_path, 'w') as f:
        f.write(feedme)


def main():
    parser = argparse.ArgumentParser(description='Cutout images')
    parser.add_argument('-i','--imgname',  help='image name to be cut', default='None')
    parser.add_argument('-o','--overwrite',  help='overwrite existing files', action='store_true')

    args = parser.parse_args()

    img_name = input("Enter image name: ") if args.imgname == 'None' else args.imgname

    pregalfit_path = prep_base + img_name + '/pregalfit.fits'
    with fits.open(pregalfit_path) as hdul:
        tab_ini = Table(hdul['INIPARAM'].data)

    sample_base = './io/sample/' + img_name + '/'
    sample_dir_list = glob.glob(sample_base + '*/')
    sample_dir_list.sort()
    
    for sample_dir in tqdm(sample_dir_list):
        genfeedme_sample((sample_dir, tab_ini))

if __name__ == '__main__':
    main()