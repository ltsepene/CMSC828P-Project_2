#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(this_dir, '../viz/src'))
from mutation_signatures_visualization import BROAD, sbs_signature_plot

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-if', '--input_files', type=str, required=True, nargs='*')
parser.add_argument('-l', '--labels', type=str, required=True, nargs='*')
parser.add_argument('-of', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

assert( len(args.input_files) == len(args.labels))

# Load the input data
data = []
row_labels = args.labels
categories = None
for input_file in args.input_files:
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    data.append( df.values.sum(axis=0)  )
    if categories is None:
        categories = list(df.columns)
    else:
        assert( categories == list(df.columns) )

# Plot counts
sbs96_counts_df = pd.DataFrame(data=data, index=row_labels, columns=categories)
sbs_signature_plot(sbs96_counts_df, palette=BROAD, ylabel='Count')
plt.tight_layout(pad=2)
plt.savefig(args.output_file)
