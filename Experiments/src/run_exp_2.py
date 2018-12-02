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
	
    # Load the post-treatment file
    logger.info('[Loading mutation counts]')
    sbs96_df_post = pd.read_csv(args.input_post, sep='\t', index_col=0)

    X_post = sbs96_df_post.values
    samples_post = sbs96_df_post.index
 
    N_post = len(samples_post)
    
    logger.info('- Loaded mutations in %s samples for the pre-treatment case' % N_post)
	
    # Load Cosmic Signatures
    logger.info('[Loading mutation counts]')
    sbs96_cs = pd.read_csv(args.cosmic_file, sep='\t', index_col=0)

    Cosmic = np.array(sbs96_cs.values)    
    logger.info('- Loaded Cosmic Signatures')

    ############################################################################
    # INFER SIGNATURES
    ############################################################################
    # Run NMF for n_run iterations
    logger.info('[Running NMF]')
    np.random.seed(args.random_seed)
		   
	# Simulate the described procedure for the post-treatment case. Keep the factorization with the minimum residual error
    min_res_error = sys.maxsize
    for i in range(args.n_run):	
        nmf = nimfa.Nmf(X_post, rank=args.k_post, seed="random_vcol", max_iter=500, update='divergence', objective='div', track_error=True)
        nmf_fit = nmf()
        res_err = nmf_fit.fit.tracker.get_error()
        avg_res_err = np.sum(res_err)/len(res_err)
        if (avg_res_err < min_res_error):
           H_post = nmf_fit.coef()
           E_post = nmf_fit.basis()
           min_res_error = avg_res_err
    
	# Normalize the signature matrix 
    for i in range(args.k_post):
        sum_h = 0
        for j in range(96):
           sum_h += H_post[i,j]
        for j in range(96):
           H_post[i,j] = H_post[i,j] / sum_h	      

    for i in range(N_post):
        sum_e = 0
        for j in range(args.k_post):
            sum_e += E_post[i,j]
        for j in range(args.k_post):   
           E_post[i,j] = E_post[i,j] / sum_e

    		   
    # For each discovered signature find the corresponding cosmic one via highest cosine similarity			
    cosm_cos_post = list(range(args.k_post))
    cosm_ind_post = list(range(args.k_post))	
	
    cos_sim = 1.-cdist(H_post, Cosmic, metric='cosine')
    agr_cos_sim = 0
    for k in range(args.k_post):
        cosm_cos_post[k] = max(cos_sim[k])   
        cosm_ind_post[k] = np.argmax(cos_sim[k])	

    # Disambiguate mutational signatures 			
    apobec_post = age = cys = 0
    for k in range(args.k_post):
        if cosm_ind_post[k] == 0:
            age = k
        else:
            if cosm_ind_post[k] == 12:
                apobec_post = k
            else:
                cys = k			
		   
    ############################################################################
    # PLOT RESULTS AND OUTPUT
    ############################################################################
    logger.info('- Plotting')

    ind = range(29)
    PFS = [73, 157, 89, 327, 603, 210, 637, 171, 47, 342, 466, 1287, 158, 147, 351, 698, 428, 98, 1776, 1136, 259, 131, 266, 18, 471, 104, 493, 438, 1028]
	
    dom_apobec = []
    dom_age = []
    dom_cys = []
    
    for i in range(N_post):
        if (E_post[i,apobec_post] > E_post[i, age]) and (E_post[i,apobec_post] > E_post[i, cys]):
            dom_apobec.append(PFS[i])
        if (E_post[i,age] > E_post[i, apobec_post]) and (E_post[i,age] > E_post[i, cys]):
            dom_age.append(PFS[i])
        if (E_post[i,cys] >= E_post[i, age]) and (E_post[i,apobec_post] <= E_post[i, cys]):
            dom_cys.append(PFS[i])

	
    data = [dom_apobec, dom_age, dom_cys]
    fig = plt.figure(1)
    f1 = fig.add_subplot(111)
    f1.set_title('PFS for Apobec/Aging/Cisplatin')
    f1.boxplot(data)
    fig.savefig(args.output_prefix)
	
	
if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ia', '--input_post', type=str, required=True, 
                        help='96-category files')
    parser.add_argument('-sf', '--cosmic_file', type=str, required=True)
    parser.add_argument('-op', '--output_prefix', type=str, required=True)
    parser.add_argument('-nr', '--n_run', type=int, required=True)
    parser.add_argument('-ka', '--k_post', type=int, required=True)
    parser.add_argument('-rs', '--random_seed', type=int, required=False, default=81047)
    parser.add_argument('-v', '--verbosity', type=int, required=False, default=logging.INFO)

    # Run
    run(parser.parse_args(sys.argv[1:]))