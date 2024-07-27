import glob
import argparse
from multiprocessing import Pool

import numpy as np

from tqdm import tqdm

from astropy.io import fits
from astropy.table import Table

from galfitclass import GalfitClass, SersicComponent
from utils.zeropoint import calc_zpt

def run_sample(sample_dir: str) -> None:
    id = sample_dir.split('/')[-2]

    feedme_path = sample_dir + f'{id}_linear.galfit'
    #print(f'Running {feedme_path}')
    GalfitClass.run(feedme_path)


def main():
    parser = argparse.ArgumentParser(description='Cutout images')
    parser.add_argument('-i','--imgname',  help='image name to be cut', default='none')
    parser.add_argument('-o','--overwrite',  help='overwrite existing files', action='store_true')
    parser.add_argument('-n','--nproc',  help='number of processes', type=int, default=1)

    args = parser.parse_args()

    img_name = input("Enter image name: ") if args.imgname == 'none' else args.imgname

    sample_base = './io/sample/' + img_name + '/'
    sample_dir_list = glob.glob(sample_base + '*/')
    sample_dir_list.sort()

    if args.nproc > 1:
        with Pool(args.nproc) as p:
            p.map(run_sample, sample_dir_list)
    else:
        for sample_dir in tqdm(sample_dir_list):
            run_sample(sample_dir)

if __name__ == '__main__':
    main()