#!/usr/bin/env python

################################################################################
# SETUP
################################################################################
# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, numpy as np, pandas as pd, logging
import seaborn as sns, matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

sns.set_style('whitegrid')

from i_o import get_logger
from constants import SUBS, SUB_COLOR

import nimfa


# Load mutation signatures visualizations
this_dir = os.path.dirname(__file__)
viz_dir = os.path.join(this_dir, '../viz/src')
sys.path.append( viz_dir )
from mutation_signatures_visualization import sbs_signature_plot, BROAD


################################################################################
# MAIN
################################################################################
def run( args ):
    ############################################################################
    # I/O
    ############################################################################
    # Set up logger
    logger = get_logger(args.verbosity)
	
    ############################################################################
    # PLOT RESULTS AND OUTPUT
    ############################################################################
    logger.info('- Plotting')

    ind = range(30)
    PFS = [73, 157, 89, 327, 603, 210, 637, 171, 47, 342, 466, 1287, 166, 158, 147, 351, 698, 428, 98, 1776, 1136, 259, 131, 266, 18, 471, 104, 493, 438, 1028]
    Age = [53, 43, 57, 40, 67, 56, 81, 77, 60, 66, 64, 44, 63, 72, 66, 73, 66, 59, 60, 64, 56, 75, 66, 79, 61, 74, 69, 53, 74, 71]
    Sex = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]
	
    young = []
    old = []
    male = []
    female = []
    
    for i in range(30):
        if Age[i] <= 65:
            young.append(PFS[i])
        else:
            old.append(PFS[i])	
        if Sex[i] == 0:
            male.append(PFS[i])
        else:
            female.append(PFS[i])		
	
    data = [young, old]
    fig = plt.figure(1)
    f1 = fig.add_subplot(211)
    f1.set_title('PFS for <=65/>65 yrs')
    f1.boxplot(data, widths = 0.2)
    data = [male, female]
    f2 = fig.add_subplot(212)
    f2.set_title('PFS for male/female')
    f2.boxplot(data, widths = 0.2)
    fig.savefig(args.output_prefix)
	
	
if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-op', '--output_prefix', type=str, required=True)
    parser.add_argument('-v', '--verbosity', type=int, required=False, default=logging.INFO)

    # Run
    run(parser.parse_args(sys.argv[1:]))