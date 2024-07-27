import numpy as np

# custom
def calc_zpt(x: int, y: int, corr_dict: str= None, **kwargs)-> list:
    """
    Calculate zeropoint
    args:
        x: int
            x position

        y: int
            y position

        corr_dict: str
            correction dictionary name
    
    return:
        list
            list of zeropoints
    """
    def _CR_zpt(**kwargs):
        if 'verbose' in kwargs.keys() and kwargs['verbose']:
            print('Using magnitude zeropoints for F150W from Boyer et al. (2022), and for the other 6 filters, we use the zero points from Brammer (2022).')
        
        CR = {'F115W':[0.914,0.908,0.862,0.813,0.959,0.897,0.995,0.831],
                'F150W':[1.005,0.966,0.932,0.855,1.000,0.968,0.988,0.865],
                'F200W':[0.865,0.889,0.800,0.791,0.873,0.829,0.901,0.812],
                'F277W':[1.107,1.107,1.107,1.107,1.000,1.000,1.000,1.000],
                'F356W':[1.088,1.088,1.088,1.088,1.000,1.000,1.000,1.000],
                'F410M':[1.031,1.031,1.031,1.031,1.042,1.042,1.042,1.042],
                'F444W':[1.036,1.036,1.036,1.036,1.096,1.096,1.096,1.096]
                }
        
        let = x<5400 #letter B for left half, A for right half
        index = let << 2
        if let:
            index = index | ((x>2350) << 1) | (y>2300) 
        else:
            index = index | ((x<8140) << 1) | (y<2400)
        #index = let * (index | ((x>2350) << 1) | (y>2300) ) + (1-let) * (index | ((x<8140) << 1) | (y<2400))

        return [28.086519392283982 - 2.5*np.log10(CR[band][index]) for band in CR.keys()]

    def _CR0_zpt(**kwargs): 
        if 'verbose' in kwargs.keys() and kwargs['verbose']:
            print('zeropoints are adopted from jwst_0995.pmap\n')

        CR0 = {'F115W':[1,1,1,1,1,1,1,1],
                'F150W':[1,1,1,1,1,1,1,1],
                'F200W':[1,1,1,1,1,1,1,1],
                'F277W':[1,1,1,1,1,1,1,1],
                'F356W':[1,1,1,1,1,1,1,1],
                'F410M':[1,1,1,1,1,1,1,1],
                'F444W':[1,1,1,1,1,1,1,1]
                }
        return [28.086519392283982]*len(CR0.keys())
    
    if corr_dict is None: # default
        kwargs['verbose']=True
        print('No correction dictionary name provided. Using default correction dictionary: CR. Automatically set verbose=True.')
        corr_dict = 'CR'
        
    corrfunc_name = '_' + corr_dict + '_zpt' # crucial

    if corrfunc_name not in dir():
        raise ValueError('Invalid correction dictionary name:', corr_dict)
    
    return eval(corrfunc_name)(**kwargs)

# test
if __name__ == '__main__':
    print(calc_zpt(0,0))
    print(calc_zpt(0,0, 'CR0'))