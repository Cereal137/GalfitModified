import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ImageNormalize, stretch
from photutils.aperture import EllipticalAnnulus, EllipticalAperture, aperture_photometry

band_labels = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']
psf_path = ['~/Desktop/GALFIT/GalfitModified/io/psf/' + band.lower() + '_psf.fits' for band in band_labels]
psf_imglist = [fits.getdata(path) for path in psf_path]
pix_scale = 0.03 # arcsec/pixel

def _flux_profile(img, aperlist) -> np.ndarray:
    return _phot_profile(img, aperlist) # in MJy/sr

def _phot_profile(img, aperlist, err = None):
    _n_aper = len(aperlist)
    areas = np.array([aper.area for aper in aperlist])
    phot = aperture_photometry(img, aperlist, error=err)
    _flux = [phot[f'aperture_sum_{i}'].data[0] for i in range(_n_aper)]
    if err is not None:
        _fluxerr = [phot[f'aperture_sum_err_{i}'].data[0] for i in range(_n_aper)]
        return np.array(_flux) / areas, np.array(_fluxerr) / areas
    else:
        return np.array(_flux) / areas
    
def phot_profile(img, aperlist, err = None):
    return _phot_profile(img, aperlist, err)
    
def _sci_profile(img, err, aperlist, zero_point: float = 28.086519392283982):
    _flux, _fluxerr = _phot_profile(img, aperlist, err) 
    flux = _flux / (pix_scale**2)

    return flux_to_mag(flux, zero_point), 2.5/np.log(10) * _fluxerr / _flux

def flux_to_mag(flux, zero_point: float = 28.086519392283982):
    return -2.5 * np.log10(flux) + zero_point # in AB mag

def _mag_profile(img, aperlist, zero_point: float = 28.086519392283982) -> np.ndarray:
    flux = _flux_profile(img, aperlist) / (pix_scale**2)
    return flux_to_mag(flux, zero_point)

def _aperlist(xc: float, yc: float, pa: float, ar: float, radii: np.ndarray):
    theta = np.radians(pa-90)
    return [EllipticalAperture((xc, yc), radii[0], radii[0]*ar)] + [EllipticalAnnulus((xc, yc), a_in, a_out, a_out*ar, a_in*ar, theta) for a_in, a_out in zip(radii[:-1], radii[1:])]

def gen_aperlist(xc: float, yc: float, pa: float, ar: float, radii: np.ndarray):
    return _aperlist(xc, yc, pa, ar, radii)

def sci_to_profile(sci, err, xc: float, yc: float, pa: float, ar: float, radii: np.ndarray, zpt: float):
    """
    make sure xc, yc are 0-indexed, not directly from GALFIT output
    """
    aperlist = _aperlist(xc, yc, pa, ar, radii)
    return _sci_profile(sci, err, aperlist, zpt)

def image_to_profile(img, xc: float, yc: float, pa: float, ar: float, radii: np.ndarray, zpt: float):
    """
    make sure xc, yc are 0-indexed, not directly from GALFIT output
    """
    theta = np.radians(pa-90)
    aperlist = [EllipticalAperture((xc, yc), radii[0], radii[0]*ar)] + [EllipticalAnnulus((xc, yc), a_in, a_out, a_out*ar, a_in*ar, theta) for a_in, a_out in zip(radii[:-1], radii[1:])]
    return _mag_profile(img, aperlist, zpt)

def psfimage_to_profile(psf, xc: float, yc: float, pa: float, ar: float, radii: np.ndarray, zpt: float = 28.086519392283982):
    """
    make sure xc, yc are 0-indexed, not directly from GALFIT output
    """
    return image_to_profile(psf, xc, yc, pa, ar, radii, zpt)

def mpl_setup(func: callable):
    def wrapper(*args, **kwargs):
        mpl.rcParams['xtick.color'] = 'w'
        mpl.rcParams['ytick.color'] = 'w'
        try:
            fig, ax = func(*args, **kwargs)
        finally:
            mpl.rcParams['xtick.color'] = 'black'
            mpl.rcParams['ytick.color'] = 'black'

        return fig, ax
    return wrapper

@mpl_setup
def plot_result_singleband(sample_base: str, sample_id: str, n_aper: int = 20, max_radius: float = 35):
    sample_base = sample_base.rstrip('/') + '/'
    sample_dir = sample_base + sample_id + '/'

    # determine a constant list of radius for each image
    rlist = []
    for band in band_labels:
        output_path = sample_dir + f'output_{band}.fits'
        row_finalband = Table.read(output_path, hdu='FINAL_BAND')[0]
        rlist.append(row_finalband['COMP1_Re'])
    r_out = max(min(3 * np.median(rlist), max_radius), 7.5)
    S_sqr = stretch.SquaredStretch()
    lin = np.linspace(0,1,n_aper)
    radii = S_sqr(lin[1:]) * r_out

    fig, ax = plt.subplots(7,4,figsize=(17,30), facecolor='k')
    for i in range(7):
        band = band_labels[i]
        output_path = sample_dir + f'output_{band}.fits'
        with fits.open(output_path) as hdul:
            sci = hdul['INPUT'].data
            model = hdul['MODEL'].data
            res = hdul['RESIDUAL'].data
            zpt = hdul['MODEL'].header[f'MAGZPT_{band}']

            row_finalband = Table(hdul['FINAL_BAND'].data)[0]
            x, y, pa, ar, n = row_finalband['COMP1_XC']-1, row_finalband['COMP1_YC']-1, row_finalband['COMP1_PA'], row_finalband['COMP1_AR'], row_finalband['COMP1_n']

            # profile part
            sci_profile = image_to_profile(sci, x, y, pa, ar, radii, zpt)
            model_profile = image_to_profile(model, x, y, pa, ar, radii, zpt)
            psf_profile = image_to_profile(psf_imglist[i], 40, 40, pa, ar, radii, zpt)
            psf_scaled = psf_profile * model_profile[0] / psf_profile[0]
            
            ax[i,0].plot(radii*0.03, psf_scaled,'b--', label='PSF')
            ax[i,0].plot(radii*0.03, model_profile,'r--', label='Model')
            ax[i,0].plot(radii*0.03, sci_profile,'o', color='k', fillstyle='none', label='Data')
            ax[i,0].set_ylabel(r'$\mu$ [mag/arcsec$^2$]', color='w',)
            ax[i,0].xaxis.set_tick_params(labelbottom=False)
            ymin, ymax = np.min(model_profile), np.max(model_profile)
            ax[i,0].set_ylim(ymax+.6, ymin-1.4)
            ax[i,0].set_xlim(0, (r_out-.5)*0.03)

            ax[i,0].text(.05, .85, f'{band}', color='k', fontsize=14, transform = ax[i,0].transAxes)

            divider = make_axes_locatable(ax[i,0])
            ax_offset = divider.append_axes("bottom", 0.8, pad=0, sharex=ax[i,0], transform = ax[i,0].transData)
            ax_offset.plot(radii*0.03, sci_profile-model_profile, 'o', color='k', fillstyle='none')
            ax_offset.set_ylim(-0.6, 0.6)
            ax_offset.plot(ax_offset.get_xlim(),[0,0], 'r--')
            ax_offset.set_ylabel(r'$\Delta\mu$', color='w')


            # image part
            model_mean, model_std = np.mean(model), np.std(model)
            vmin = max(model_mean - 10*model_std, 0)
            vmax = np.max(model) * 3.
            norm = ImageNormalize(model, stretch=stretch.LogStretch(), vmin=vmin, vmax=vmax)
            ax[i,2].text(.1, .1, f'$n={n:.2f}$',color='w', fontsize=12, transform = ax[i,2].transAxes)
            ax[i,2].text(.1, .2, f'$R_e={rlist[i]*.03:.2f}\,$'+'\u2033', color='w', fontsize=12, transform = ax[i,2].transAxes)
            im1 = ax[i,1].imshow(sci, origin='lower', cmap='viridis', norm=norm)
            im2 = ax[i,2].imshow(model, origin='lower', cmap='viridis', norm=norm)
            im3 = ax[i,3].imshow(res, origin='lower', cmap='viridis', norm=norm)

            # chi2  
            chi2 = hdul['MODEL'].header['CHI2NU']
            t = ax[i,1].text(.1 ,.1 ,r'$\chi^2_\nu=$'+f'{chi2:.3f}',color='k', fontsize=12, ha='left', transform = ax[i,1].transAxes)
            t.set_bbox(dict(facecolor='white', alpha=.7))

            # scalebar
            scalebar = AnchoredSizeBar(ax[i,3].transData,
                            0.5/0.03, ' 0.5\u2033 ', 'lower right', 
                            pad=1,
                            color='w',
                            fontproperties={'size': 12},
                            frameon=False,
                            size_vertical=1,
                            )
            ax[i,1].add_artist(scalebar)
            # colorbar
            cax = ax[i,3].inset_axes([0.,0.,1.,0.03])
            fig.colorbar(im3, cax=cax, orientation='horizontal',)

            for j in [1,2,3]:
                ax[i,j].set_xticks([])
                ax[i,j].set_yticks([])
                for spine in ax[i,j].spines.values():
                    spine.set_edgecolor('w')
        
    ax[0,0].set_title('Profile', color='w', fontsize=20)
    ax[0,0].legend(loc='center right', fontsize=12)

    ax[0,1].set_title('Science', color='w', fontsize=20)
    ax[0,2].set_title('Model', color='w', fontsize=20)
    ax[0,3].set_title('Residual', color='w', fontsize=20)
    plt.subplots_adjust(wspace=0.04, hspace=0.08, left=0.05, right=0.95, top=0.95, bottom=0.05)

    return fig, ax

@mpl_setup
def plot_result_multiband(sample_base: str, sample_id: str, output_name: str = 'output', n_aper: int = 20, max_radius: float = 35):
    sample_base = sample_base.rstrip('/') + '/'
    sample_dir = sample_base + sample_id + '/'

    fig, ax = plt.subplots(7,4,figsize=(17,30), facecolor='k')
    #fig.suptitle(f'EGS-{sample_num}', fontsize=26, color='white')
    with fits.open(sample_dir + output_name + '.fits') as hdul:
        tab_finalband = Table(hdul['FINAL_BAND'].data)
        n_out, rlist = tab_finalband['COMP1_n'], tab_finalband['COMP1_Re']
        ar, pa, re = tab_finalband['COMP1_AR'][0], tab_finalband['COMP1_PA'][0], tab_finalband['COMP1_Re'][3]

        # determine a constant list of radius for each image
        r_out = max(min(3 * np.median(rlist), max_radius), 7.5)
        S_sqr = stretch.SquaredStretch()
        lin = np.linspace(0,1,n_aper)
        radii = S_sqr(lin[1:]) * r_out

        for i in range(7):
            sci = hdul[i+1].data
            model = hdul[i+8].data
            res = hdul[i+15].data
            zpt = hdul[i+8].header[f'MAGZPT_{band_labels[i]}']

            x, y= tab_finalband['COMP1_XC'][i]-1, tab_finalband['COMP1_YC'][i]-1
            re = rlist[i]

            # profile part
            sci_profile = image_to_profile(sci, x, y, pa, ar, radii, zpt)
            model_profile = image_to_profile(model, x, y, pa, ar, radii, zpt)
            psf_profile = image_to_profile(psf_imglist[i], 40, 40, pa, ar, radii, zpt)
            psf_scaled = psf_profile * model_profile[0] / psf_profile[0]

            ax[i,0].plot(radii*0.03, psf_scaled,'b--', label='PSF')
            ax[i,0].plot(radii*0.03, model_profile,'r--', label='Model')
            ax[i,0].plot(radii*0.03, sci_profile,'o', color='k', fillstyle='none', label='Data')
            ax[i,0].set_ylabel(r'$\mu$ [mag/arcsec$^2$]', color='w',)
            ax[i,0].xaxis.set_tick_params(labelbottom=False)
            ymin, ymax = np.min(model_profile), np.max(model_profile)
            ax[i,0].set_ylim(ymax+.6, ymin-1.4)
            ax[i,0].set_xlim(0, (r_out+.5)*0.03)

            ax[i,0].text(.05, .85, f'{band_labels[i]}', color='k', fontsize=14, transform = ax[i,0].transAxes)

            divider = make_axes_locatable(ax[i,0])
            ax_offset = divider.append_axes("bottom", 0.8, pad=0, sharex=ax[i,0], transform = ax[i,0].transData)
            ax_offset.plot(radii*0.03, sci_profile-model_profile, 'o', color='k', fillstyle='none')
            ax_offset.set_ylim(-0.6, 0.6)
            ax_offset.plot(ax_offset.get_xlim(),[0,0], 'r--')
            ax_offset.set_ylabel(r'$\Delta\mu$', color='w')


            # image part
            model_mean, model_std = np.mean(model), np.std(model)
            vmin = max(model_mean - 10*model_std, 0)
            vmax = np.max(model) * 3.
            norm = ImageNormalize(model, stretch=stretch.LogStretch(), vmin=vmin, vmax=vmax)
            ax[i,2].text(.1, .1, f'$n={n_out[i]:.2f}$',color='w', fontsize=12, transform = ax[i,2].transAxes)
            ax[i,2].text(.1, .2, f'$R_e={re*.03:.2f}\,$'+'\u2033', color='w', fontsize=12, transform = ax[i,2].transAxes)
            _ = ax[i,1].imshow(sci, origin='lower', cmap='viridis', norm=norm)
            _ = ax[i,2].imshow(model, origin='lower', cmap='viridis', norm=norm)
            im3 = ax[i,3].imshow(res, origin='lower', cmap='viridis', norm=norm)

            # chi2  
            chi2 = hdul[8+i].header['CHI2NU']
            t = ax[i,1].text(.1 ,.1 ,r'$\chi^2_\nu=$'+f'{chi2:.3f}',color='k', fontsize=12, ha='left', transform = ax[i,1].transAxes)
            t.set_bbox(dict(facecolor='white', alpha=.7))

            # scalebar
            scalebar = AnchoredSizeBar(ax[i,3].transData,
                            0.5/0.03, ' 0.5\u2033 ', 'lower right', 
                            pad=1,
                            color='w',
                            fontproperties={'size': 12},
                            frameon=False,
                            size_vertical=1,
                            )
            ax[i,1].add_artist(scalebar)
            
            # colorbar
            cax = ax[i,3].inset_axes([0.,0.,1.,0.03])
            fig.colorbar(im3, cax=cax, orientation='horizontal',)

            for j in [1,2,3]:
                ax[i,j].set_xticks([])
                ax[i,j].set_yticks([])
                for spine in ax[i,j].spines.values():
                    spine.set_edgecolor('w')
        
        chi2_tot = Table(hdul['FIT_INFO'].data)['CHI2NU'].data[0]
        ax[6,0].text(.8, .8, r'$\chi^2_\nu=$'+f'{chi2_tot:.3f}\n (All Bands)', color='k', fontsize=12, ha='center', transform = ax[6,0].transAxes)
            
    ax[0,0].set_title('Profile', color='w', fontsize=20)
    ax[0,0].legend(loc='center right', fontsize=12)

    ax[0,1].set_title('Science', color='w', fontsize=20)
    ax[0,2].set_title('Model', color='w', fontsize=20)
    ax[0,3].set_title('Residual', color='w', fontsize=20)
    plt.subplots_adjust(wspace=0.04, hspace=0.08, left=0.05, right=0.95, top=0.95, bottom=0.05)

    return fig, ax