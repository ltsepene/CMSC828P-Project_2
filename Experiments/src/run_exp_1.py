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
    sbs96_df_pre = pd.read_csv(args.input_pre, sep='\t', index_col=0)

    X_pre = sbs96_df_pre.values

    samples_pre = sbs96_df_pre.index
 
    N_pre = len(samples_pre)
    categories = list(sbs96_df_pre.columns)
    cat_index  = dict(zip(categories, range(len(categories))))
    
    logger.info('- Loaded mutations in %s samples for the pre-treatment case' % N_pre)
	
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

	# Simulate the described procedure for the pre-treatment case. Keep the factorization with the minimum residual error
    min_res_error = sys.maxsize
    for i in range(args.n_run):	
        nmf = nimfa.Nmf(X_pre, rank=args.k_pre, seed="random_vcol", max_iter=500, update='divergence', objective='div', track_error=True)
        nmf_fit = nmf()
        res_err = nmf_fit.fit.tracker.get_error()
        avg_res_err = np.sum(res_err)/len(res_err)
        if (avg_res_err < min_res_error):
           H_pre = nmf_fit.coef()
           E_pre = nmf_fit.basis()
           min_res_error = avg_res_err
    
	# Normalize the signature and exposures matrix 
    for i in range(args.k_pre):
        sum_h = 0
        for j in range(96):
           sum_h += H_pre[i,j]
        for j in range(96):
           H_pre[i,j] = H_pre[i,j] / sum_h	   

    for i in range(N_pre):
        sum_e = 0
        for j in range(args.k_pre):
            sum_e += E_pre[i,j]
        for j in range(args.k_pre):   
           E_pre[i,j] = E_pre[i,j] / sum_e	
		   
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
    cosm_cos_pre = list(range(args.k_pre))
    cosm_ind_pre = list(range(args.k_pre))	
	
    cos_sim = 1.-cdist(H_pre, Cosmic, metric='cosine')
    agr_cos_sim = 0
    for k in range(args.k_pre):
        cosm_cos_pre[k] = max(cos_sim[k])   
        cosm_ind_pre[k] = np.argmax(cos_sim[k])		
		
    cosm_cos_post = list(range(args.k_post))
    cosm_ind_post = list(range(args.k_post))	
	
    cos_sim = 1.-cdist(H_post, Cosmic, metric='cosine')
    agr_cos_sim = 0
    for k in range(args.k_post):
        cosm_cos_post[k] = max(cos_sim[k])   
        cosm_ind_post[k] = np.argmax(cos_sim[k])	

    # Disambiguate mutational signatures 
    apobec_pre = ner = 0
    for k in range(args.k_pre):
        if cosm_ind_pre[k] == 12:
            apobec_pre = k
        else:
            ner = k		
			
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
    E_pre = np.delete(E_pre, 12, 0)
    E_pre_tmp = np.transpose(E_pre)
    E_pre_tmp = E_pre_tmp.tolist()
	
    E_post_tmp = np.transpose(E_post)
    E_post_tmp = E_post_tmp.tolist()	

    bars = [E_post_tmp[apobec_post][j] + E_post_tmp[age][j] for j in range(N_post)]	
	
    fig = plt.figure(1)
    f1 = fig.add_subplot(211)	
    p1 = f1.bar(ind, E_pre_tmp[ner], width=0.4, label='Signature 1', color='blue', bottom=E_pre_tmp[apobec_pre])
    p2 = f1.bar(ind, E_pre_tmp[apobec_pre], width=0.4, label='Signature 2', color='red')
    plt.ylim([0,1.25])
    plt.legend((p1[0], p2[0]), ('NER', 'APOBEC'), fontsize=8, ncol=4, framealpha=0, fancybox=True)		
    plt.title('Pre-treatment')
    f2 = fig.add_subplot(212)
    p1 = f2.bar(ind, E_post_tmp[cys], width=0.4, label='Signature 1', color='green', bottom=bars)
    p2 = f2.bar(ind, E_post_tmp[apobec_post], width=0.4, label='Signature 2', color='silver', bottom=E_post_tmp[age])
    p3 = f2.bar(ind, E_post_tmp[age], width=0.4, label='Signature 3', color='yellow')
    plt.ylim([0,1.25])
    plt.legend((p1[0], p2[0], p3[0]), ('CYS', 'APOBEC', 'AGE'), fontsize=8, ncol=4, framealpha=0, fancybox=True)
    plt.title('Post-treatment')
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
    fig.savefig(args.output_prefix)
	
	
if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ib', '--input_pre', type=str, required=True, 
                        help='96-category files')
    parser.add_argument('-ia', '--input_post', type=str, required=True, 
                        help='96-category files')
    parser.add_argument('-sf', '--cosmic_file', type=str, required=True)
    parser.add_argument('-op', '--output_prefix', type=str, required=True)
    parser.add_argument('-nr', '--n_run', type=int, required=True)
    parser.add_argument('-kb', '--k_pre', type=int, required=True)
    parser.add_argument('-ka', '--k_post', type=int, required=True)
    parser.add_argument('-rs', '--random_seed', type=int, required=False, default=81047)
    parser.add_argument('-v', '--verbosity', type=int, required=False, default=logging.INFO)

    # Run
    run(parser.parse_args(sys.argv[1:]))