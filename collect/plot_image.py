import argparse
import sys
import glob
sys.path.append('../')

import matplotlib.pyplot as plt

from plot.plot_result import plot_result_multiband

from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='Plot image')
    parser.add_argument('-i', '--img_name', type=str, help='Image name')
    parser.add_argument('-s', '--save_dir', type=str, help='Save directory')
    args = parser.parse_args()

    img_name = args.img_name
    save_dir = args.save_dir
    save_dir = save_dir.rstrip('/') + '/'
    sample_base = '../io/sample/'+img_name+'/'

    img_dirs = glob.glob(sample_base+'*/')
    for img_dir in tqdm(img_dirs):
        id = img_dir.split('/')[-2]
        try:
            fig, _ = plot_result_multiband(sample_base, id, output_name='output_quad')
            fig.savefig(save_dir +id+'.png', dpi=100)
        except:
            print('Error:', id)
        plt.close(fig)

if __name__ == '__main__':
    main()