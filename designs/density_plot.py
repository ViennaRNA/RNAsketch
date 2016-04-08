from __future__ import print_function

try:
    from PyDesign import *
    import PyDesign
except ImportError, e:
    print(e.message)
    exit(1)

import RNAdesign as rd
import argparse
import sys
import time
import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import glob
import os
from decimal import Decimal

def plot_sequence_objective(args):

    x_obj = []
    y_obj = []
    path = '*.csv'
    if args.directory is not None:
        if not os.path.exists(args.directory):
            print("Path does not exist!", file=sys.stderr)
            exit()
        path = os.path.join(args.directory, path)

    filenames = []
    if args.file is not None:
        if ";" in args.file:
            filenames = args.file.split(";")
        else:
            filenames.append(args.file)
    else:
        filenames = glob.glob(path)

    for filename in filenames:
        with open(filename) as sampled_seq:
            reader = csv.reader(sampled_seq, delimiter=";")
            header = reader.next()

            # get information of intial sequence
            initial_sequence = reader.next()
            x_center = float(initial_sequence[0])
            y_center = float(initial_sequence[1])


            x_origin = 0 # = x_center-x_center
            y_origin = 0

            x_obj.append(x_origin)
            y_obj.append(y_origin)


            for line in reader:
                x_new = float(line[0])
                y_new = float(line[1])

                # calculate relative positions to initial sequence
                x = x_new - x_center
                y = y_new - y_center

                x_obj.append(x)
                y_obj.append(y)

    # find largest x/y and use as limit for both axes
    if max(x_obj) > max(y_obj):
        ax_lim = abs(max(x_obj))
    else:
        ax_lim = abs(max(y_obj))

    # density plot https://oceanpython.org/2013/02/25/2d-histogram/
    # Estimate the 2D histogram
    nbins = 50 # TODO: how to choose?
    H, xedges, yedges = np.histogram2d(x_obj, y_obj, bins=nbins)

    H = np.rot90(H)
    H = np.flipud(H)

    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

    # Plot 2D histogram
    plt.figure()
    ax_lim += 1  # +1 to avoid points on axes
    plt.xlim([-ax_lim, ax_lim])
    plt.ylim([-ax_lim, ax_lim])
    plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight, -ax_lim * args.weight], 'm')

    # counts per quadrants
    quad1 = 0
    quad2 = 0
    quad3 = 0
    quad4 = 0

    weighted1 = 0
    weighted2 = 0
    for i in range(0, len(x_obj)):
        if x_obj[i] > 0 and y_obj[i] > 0:
            quad1 += 1
        if x_obj[i] < 0 and y_obj[i] > 0:
            quad2 += 1
            if y_obj[i] < (-args.weight * x_obj[i]):
                weighted1 += 1
        if x_obj[i] < 0 and y_obj[i] < 0:
            quad3 += 1
        if x_obj[i] > 0 and y_obj[i] < 0:
            quad4 += 1
            if y_obj[i] < (-args.weight * x_obj[i]):
                weighted2 += 1

    sum_weighted = weighted1 + weighted2 + quad3

    plt.text(ax_lim - 4, ax_lim + 2, quad1, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), ax_lim + 2, quad2, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), -(ax_lim - 2), quad3, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(ax_lim - 4, -(ax_lim - 2), quad4, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), -(ax_lim - 5), sum_weighted, style='italic',
    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

    plt.pcolormesh(xedges,yedges,Hmasked)
    plt.xlabel('objective 1')
    plt.ylabel('objective 2')
    plt.axhline(color='k')
    plt.axvline(color='k')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig(args.out_file + '.png')
    plt.close()
    print('# Output file: %s'% args.out_file + '.png')
    print('# Quadrant 1: %s'% quad1)
    print('# Quadrant 2: %s'% quad2)
    print('# Quadrant 3: %s'% quad3)
    print('# Quadrant 4: %s'% quad4)
    print('# Quadrant 3++: %s'% sum_weighted)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot...')
    parser.add_argument("-w", "--weight", type=float, default=0.5, help='Define weighting-factor')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat(), help='Name graphic')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.csv format. Separate multiple filenames with ";". If not set reads all *.csv in current directory')
    parser.add_argument("-d", "--directory", type = str, default=None, help='Set directory')
    args = parser.parse_args()

    plot_sequence_objective(args)



