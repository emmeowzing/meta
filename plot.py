#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Plot processed image data """

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from argparse import ArgumentParser


def plot(fname: str) -> None:
    with Image.open(fname) as data:
        data = np.array(data)
        plt.figure(figsize=(9, 6), dpi=200)
        plt.imshow(data, interpolation='none', cmap='gist_heat')
        
        fig_name = fname.split('/')[-1] if '/' in fname else fname
        print('fig_{0}'.format(fig_name))
        plt.savefig('fig_{0}'.format(fig_name))


if __name__ == '__main__':
    parser = ArgumentParser(description='** ' + __doc__)
    parser.add_argument('inputfile', type=str, help='Input file name')
    args = parser.parse_args()
    del parser

    plot(args.inputfile)
