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

    # Load the pre-treatment file
    logger.info('[Loading mutation counts]')
    sbs96_df = pd.read_csv(args.input_file, sep='\t', index_col=0)

    X = sbs96_df.values

    samples = sbs96_df.index
 
    N = len(samples)
    categories = list(sbs96_df.columns)
    cat_index  = dict(zip(categories, range(len(categories))))
    
    logger.info('- Loaded mutations in %s samples' % N)

    ############################################################################
    # Calculate total mutation load per patients
    ############################################################################
    load = list(range(N))
    for i in range(N):
        sum = 0
        for j in range(96):
            sum += X[i,j]
        load[i] = sum
		   
    ############################################################################
    # PLOT RESULTS AND OUTPUT
    ############################################################################
    logger.info('- Plotting')

    ind = range(30)
    PFS = [73, 157, 89, 327, 603, 210, 637, 171, 47, 342, 466, 1287, 166, 158, 147, 351, 698, 428, 98, 1776, 1136, 259, 131, 266, 18, 471, 104, 493, 438, 1028]

    if args.del_bool == 0:
        del PFS[12]
        ind = range(29)

    short = []
    long = []
	
    for i in range(len(PFS)):
       if load[i] <= np.median(load):
           short.append(PFS[i])
       else:
           long.append(PFS[i])
		   
    data = [short, long]
    fig = plt.figure(1)
    f1 = fig.add_subplot(211)
    f1.set_title('PFS for load <=median / > median')
    f1.boxplot(data, widths = 0.2)
    fig.savefig(args.output_prefix)
	
	
if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--input_file', type=str, required=True, 
                        help='96-category files')
    parser.add_argument('-d', '--del_bool', type=int, required=True) 
    parser.add_argument('-op', '--output_prefix', type=str, required=True)
    parser.add_argument('-v', '--verbosity', type=int, required=False, default=logging.INFO)

    # Run
    run(parser.parse_args(sys.argv[1:]))