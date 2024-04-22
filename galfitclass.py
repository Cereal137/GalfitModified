import subprocess
import numpy as np
from astropy.io import fits


class GalfitClass():
    def __init__(self, dict_band: dict, output_path: str, plate_scale: tuple =(0.030,0.030)):
        """
        Initialize GalfitClass
        args:
            dict_band: dict
                Dictionary of band labels, pivot wavelengths and magnitude zeropoints

            plate_scale: tuple
                Plate scale in arcsec/pixel
        """
        self.band_labels = dict_band.keys()
        self.band_wavelengths = [item[1][0] for item in dict_band.items()]
        self.nbands = len(dict_band)
        self.img_paths = []
        self.output_path = output_path
        self.err_paths = []
        self.bpmask_paths = []
        self.data_loaded = False
        self.psf_paths = []
        self.psf_loaded = False
        self.psf_sample_factor = 1.0
        self.constraint_path = ''
        self.image_size = (0,0)
        self.conv_box_size = (100,100)
        self.mag_zerpoints = [item[1][1] for item in dict_band.items()]
        self.plate_scale = plate_scale
        self.display_type = 0

        self.component_list = []

    def reload_data(self):
        """
        Reload data
        """
        self.img_paths = []
        self.err_paths = []
        self.bpmask_paths = []
        self.data_loaded = False

    def load_data(self, img_paths: list, err_paths: list, bpmask_paths: list, reload=False, verbose=True):
        """
        Load data
        args:
            img_paths: list
                List of paths to SCI images

            err_paths: list
                List of paths to ERR images

            bpmask_paths: list
                List of paths to bad pixel masks

            reload: bool
                If True, reload data
        """
        if reload:
            self.reload_data()
        
        if len(img_paths) != self.nbands:
            self.reload_data()
            raise ValueError('Number of images does not match the number of bands')
        self.img_paths = img_paths
        if verbose:
            print('SCI images loaded')
        
        if len(err_paths) != self.nbands:
            self.reload_data()
            raise ValueError('Number of error images does not match the number of bands')
        self.err_paths = err_paths
        if verbose:
            print('ERR images loaded')
        
        if len(bpmask_paths) != self.nbands:
            self.reload_data()
            raise ValueError('Number of bad pixel masks does not match the number of bands')
        self.bpmask_paths = bpmask_paths
        if verbose:
            print('Bad pixel masks loaded')

        with fits.open(self.img_paths[0]) as hdul:
            self.image_size = hdul['SCIENCE'].data.shape
        
        self.data_loaded = True
        pass

    def load_psf(self, psf_paths: list, verbose=True):
        """
        Load psf
        args:
            psf_paths: list
                List of paths to PSF images
        """
        if len(psf_paths) != self.nbands:
            raise ValueError('Number of psf images does not match the number of bands')
        self.psf_paths = psf_paths
        if verbose:
            print('PSF images loaded')
        self.psf_loaded = True

    def load_constraint(self,constraint_path: str, verbose=True):
        """ 
        Load constraint file
        args:
            constraint_path: str
                Path to constraint file
        """
        self.constraint_path = constraint_path
        if verbose:
            print('Constraint file loaded')

    def genstr_feedme(self):
        """
        Generate feedme file
        """
        feedme_A = 'A) ' + ','.join(self.img_paths) + ' # Input data image(s) (FITS file)\n'
        feedme_A1 = 'A1) ' + ','.join(self.band_labels) + ' # Band labels\n'
        feedme_A2 = 'A2) ' + ','.join([str(w) for w in self.band_wavelengths]) + ' # Band pivot wavelengths (microns)\n'
        feedme_B = f'B) {self.output_path} # Output data image block\n'
        feedme_C = 'C) ' + ','.join(self.err_paths) + ' # Input sigma image(s) (FITS file or NONE)\n'
        feedme_D = 'D) ' + ','.join(self.psf_paths) + ' # Input PSF image(s) (FITS file or NONE)\n'
        feedme_E = 'E) ' + str(self.psf_sample_factor) + ' # PSF oversampling factor relative to data\n'
        feedme_F = 'F) ' + ','.join(self.bpmask_paths) + ' # Bad pixel mask (FITS image or ASCII coord list)\n'
        feedme_G = 'G) ' + self.constraint_path + ' # File with parameter constraints (ASCII file)\n'
        feedme_H = 'H) ' + '1 ' + str(self.image_size[0]) + ' 1 ' + str(self.image_size[1]) + ' # Image region to fit (xmin xmax ymin ymax)\n'
        feedme_I = 'I) ' + str(self.conv_box_size[0]) + ' ' + str(self.conv_box_size[1]) + ' # Size of the convolution box (x y)\n'
        feedme_J = 'J) ' + ','.join([str(zp) for zp in self.mag_zerpoints]) + ' # Magnitude photometric zeropoint(s) (mag)\n'
        feedme_K = 'K) ' + str(self.plate_scale[0]) + ' ' + str(self.plate_scale[1]) + ' # Plate scale (dx dy)   [arcsec per pixel]\n'
        feedme_O = 'O) ' + str(self.display_type) + ' # Display type\n'

        feedme = feedme_A + feedme_A1 + feedme_A2 + feedme_B + feedme_C + feedme_D + feedme_E + feedme_F + feedme_G + feedme_H + feedme_I + feedme_J + feedme_K + feedme_O + '\n'

        if self.component_list:
            for component in self.component_list:
                feedme += component.genstr_feedme() + '\n'
        else:
            raise ValueError('No component is added')
        
        return feedme
    
    def add_component(self, component) -> None:
        """
        Add component to feedme file
        args:
            component: Component
                Component to be added
        """
        self.component_list.append(component)

    @staticmethod
    def run(feedme_path: str) -> None:
        """
        Run Galfit
        args:
            feedme_path: str
                Path to feedme file
        
        """
        cmd = 'galfitm '+ feedme_path
        child = subprocess.call(['/bin/zsh',"-i", "-c", cmd])

class SersicComponent():
    def __init__(self, nbands, skip=False):
        """
        Initialize SersicComponent
        args:
            nbands: int
                Number of bands
            skip: bool
                Skip this component
        """
        self.nbands = nbands
        self.x_configured = False
        self.y_configured = False
        self.mag_configured = False
        self.sersic_configured = False
        self.re_configured = False
        self.n_configured = False
        self.q_configured = False
        self.pa_configured = False
        self.skip = skip

    def config_x(self, dof_x: int, x_ini_list: list, cheb_mode: bool = False):
        """
        Configure x initial guess
        """
        if dof_x > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_x = dof_x
        if len(x_ini_list) != dof_x:
            raise ValueError('Number of x initial guesses does not match the number of degrees of freedom')
        
        self.x_ini_list = x_ini_list
        self.x_cheb_mode = cheb_mode
        self.x_configured = True
    
    def config_y(self, dof_y: int, y_ini_list: list, cheb_mode: bool = False):
        """
        Configure y initial guess
        """
        if dof_y > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_y = dof_y
        if len(y_ini_list) != dof_y:
            raise ValueError('Number of y initial guesses does not match the number of degrees of freedom')
        
        self.y_ini_list = y_ini_list
        self.y_cheb_mode = cheb_mode
        self.y_configured = True

    def config_mag(self, dof_mag: int, mag_ini_list: list, cheb_mode: bool = False):
        """
        Configure magnitude initial guess
        """
        if dof_mag > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_mag = dof_mag
        if len(mag_ini_list) != dof_mag:
            raise ValueError('Number of magnitude initial guesses does not match the number of degrees of freedom')
        
        mag_nan = np.isnan(mag_ini_list)
        if np.any(mag_nan):
            print('Warning: magnitude initial guess for band ' + str(np.where(mag_nan)[0]) + ' is NaN')
            print('Setting NaN to 25.0')
        mag_ini_list = np.nan_to_num(mag_ini_list, nan=25.0)

        self.mag_ini_list = mag_ini_list
        self.mag_cheb_mode = cheb_mode
        self.mag_configured = True

    def config_sersic(self, dof_sersic: int, sersic_ini_list: list, cheb_mode: bool = False):
        """
        Configure sersic initial guess
        """
        if dof_sersic > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_sersic = dof_sersic
        if len(sersic_ini_list) != dof_sersic:
            raise ValueError('Number of sersic initial guesses does not match the number of degrees of freedom')
        
        self.sersic_ini_list = sersic_ini_list
        self.sersic_cheb_mode = cheb_mode
        self.sersic_configured = True

    def config_re(self, dof_re: int, re_ini_list: list, cheb_mode: bool = False):
        """
        Configure re initial guess
        """
        if dof_re > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_re = dof_re
        if len(re_ini_list) != dof_re:
            raise ValueError('Number of re initial guesses does not match the number of degrees of freedom')
        
        self.re_ini_list = re_ini_list
        self.re_cheb_mode = cheb_mode
        self.re_configured = True

    def config_n(self, dof_n: int, n_ini_list: list, cheb_mode: bool = False):
        """
        Configure n initial guess
        """
        if dof_n > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_n = dof_n
        if len(n_ini_list) != dof_n:
            raise ValueError('Number of n initial guesses does not match the number of degrees of freedom')
        
        self.n_ini_list = n_ini_list
        self.n_cheb_mode = cheb_mode
        self.n_configured = True

    def config_q(self, dof_q: int, q_ini_list: list, cheb_mode: bool = False):
        """
        Configure q initial guess
        """
        if dof_q > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_q = dof_q
        if len(q_ini_list) != dof_q:
            raise ValueError('Number of q initial guesses does not match the number of degrees of freedom')
        
        self.q_ini_list = q_ini_list
        self.q_cheb_mode = cheb_mode
        self.q_configured = True

    def config_pa(self, dof_pa: int, pa_ini_list: list, cheb_mode: bool = False):
        """
        Configure pa initial guess
        """
        if dof_pa > self.nbands:
            raise ValueError('Number of degrees of freedom cannot be larger than the number of bands')
        self.dof_pa = dof_pa
        if len(pa_ini_list) != dof_pa:
            raise ValueError('Number of pa initial guesses does not match the number of degrees of freedom')
        
        self.pa_ini_list = pa_ini_list
        self.pa_cheb_mode = cheb_mode
        self.pa_configured = True


    def genstr_feedme(self):
        """
        Generate string for feedme file
        """
        feedme_sersic = ' 0) sersic\n'
        if self.x_configured:
            feedme_x = ' 1) '
            str_x_list = ','.join([str(x) for x in self.x_ini_list])
            feedme_x += str_x_list + ' ' + str(self.dof_x) + '\n'
        else:
            raise ValueError('x initial guess not configured')
        
        if self.y_configured:
            feedme_y = ' 2) '
            str_y_list = ','.join([str(y) for y in self.y_ini_list])
            feedme_y += str_y_list + ' ' + str(self.dof_y) + '\n'
        else:
            raise ValueError('y initial guess not configured')
        
        if self.mag_configured:
            feedme_mag = ' 3) '
            str_mag_list = ','.join([str(mag) for mag in self.mag_ini_list])
            feedme_mag += str_mag_list + ' ' + str(self.dof_mag) + '\n'
        else:
            feedme_mag = ' 3) 20.0 1 \n'

        if self.re_configured:
            feedme_re = ' 4) '
            str_re_list = ','.join([str(re) for re in self.re_ini_list])
            feedme_re += str_re_list + ' ' + str(self.dof_re) + '\n'
        else:
            feedme_re = ' 4) 10.0 3 \n'

        if self.n_configured:
            feedme_n = ' 5) '
            str_n_list = ','.join([str(n) for n in self.n_ini_list])
            feedme_n += str_n_list + ' ' + str(self.dof_n) + '\n'
        else:
            feedme_n = ' 5) 4.0 3 \n'

        if self.q_configured:
            feedme_q = ' 9) '
            str_q_list = ','.join([str(q) for q in self.q_ini_list])
            feedme_q += str_q_list + ' ' + str(self.dof_q) + '\n'
        else:
            feedme_q = ' 9) 0.5 1 \n'
        
        if self.pa_configured:
            feedme_pa = ' 10) '
            str_pa_list = ','.join([str(pa) for pa in self.pa_ini_list])
            feedme_pa += str_pa_list + ' ' + str(self.dof_pa) + '\n'
        else:
            feedme_pa = ' 10) 0.0 1 \n'

        feedme_skip = 'Z) 1 \n' if self.skip else 'Z) 0 \n'

        return feedme_sersic + feedme_x + feedme_y + feedme_mag + feedme_re + feedme_n + feedme_q + feedme_pa + feedme_skip