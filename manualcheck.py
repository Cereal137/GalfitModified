import glob
import argparse

import warnings
warnings.filterwarnings('ignore')

from utils.matchcat import MatchCat

match_base = './io/match/'
image_base = './io/image/'
hst_base = match_base


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--imgname', help='image name to be checked', default='none')
    args = parser.parse_args()

    img_name = input("Enter image name: ") if args.imgname == 'none' else args.imgname
    
    matc = MatchCat(img_name)

    hst_dir = hst_base + img_name +'/'
    hst_path = glob.glob(hst_dir + 'egs*.fits')[0]
    
    matc.interactive_match(hst_path)


if __name__ == "__main__":
    main()