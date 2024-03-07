#!/usr/bin/env python

'''
Scoring SourceApp tune results against a gold standard
Note, SourceApp's default installation instructions do
not include sklearn -- you'll want to install it 
e.g., pip install scikit-learn 
if you want to reproduce this workflow.
------------------------------------------------------
'''

import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_absolute_error

def score_att(est, key):
    # make est and key into binaries and then score as classification  
    scores = []
    key_abs = key.astype(bool).astype(int)
    est_abs = est.astype(bool).astype(int)
    for iteration in est_abs.columns:
        tn, fp, fn, tp = confusion_matrix(key_abs, est_abs[iteration]).ravel()
        sens = (tp / (tp+tn))
        spec = (tn / (tn+fp))
        scores.append(str(sens) + "|" + str(spec))  
    return scores

def score_app(est, key):
    # make est and key into composites of total fecal load, e.g., 1A + 2B
    scores = []
    key_app = key / key.sum()
    for iteration in est.columns:
        est_app = est[iteration] / est[iteration].sum()
        scores.append(mean_absolute_error(key_app, est_app))         
    return scores

def score_frac(est, key):
    # score as-is -- stretch goal of the work 
    scores = []
    for iteration in est.columns:
        scores.append(mean_absolute_error(key, est[iteration]))
    return scores

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input',
        help='File with cell fractions estimated by sourceapp_tune',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--key',
        help='File with true cell fractions',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output',
        help='Output ',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--sample',
        help='Sample name (e.g., MST124)',
        metavar='',
        type=str,
        required=True,
        default="Class"
        )
    args=vars(parser.parse_args())


    print('\n\Scoring SourceApp ...\n')
    est = pd.read_csv(args["input"])
    key = pd.read_csv(args["key"],header=0,index_col=0).loc[args["sample"]]
    sources = est.index[0:-3][0::2].tolist()
    est_nocrx=pd.DataFrame()
    est_crx = pd.DataFrame()
    for source in sources:
        est_nocrx[source]=est.loc[source]
        est_crx[source]=est.loc[source] + est.loc[source + "_crx"]
        
    scores = score_att(est_nocrx, key)
    
    est["score_att"] = scores
    
    scores = score_att(est_crx, key)
    
    est["score_att_crx"] = scores
    
    scores = score_app(est_nocrx, key)
    
    est["score_app"] = scores
    
    scores = score_app(est_crx, key)
    
    est["score_app_crx"] = scores
    
    scores = score_frac(est_nocrx, key)
    
    est["score_frac"] = scores
    
    scores = score_frac(est_crx, key)
    
    est["score_frac_crx"] = scores
    
    est.to_csv(args["output"]+".scores.csv",sep=",")
    
    ####
    #           iteration0, ..., iterationN
    # source1
    # ...
    # source N
    # param1
    # param...
    # paramN
    
    est_score_nocrx
    est_score_crx


if __name__ == "__main__":
    main()
