#!/usr/bin/env python

'''
Scoring SourceApp tune results against a gold standard. Note, SourceApp's default
installation instructions do not include sklearn -- you'll want to install it if 
you want to run this workflow
e.g., pip install scikit-learn 
----------------------------------------------------------------------------------
'''
import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import root_mean_squared_error

def score_att(est, key):
    # make est and key into binaries and then score as classification  
    sensitivity = []
    specificity = []
    key_abs = key.astype(bool).astype(int)
    est_abs = est.astype(bool).astype(int)
    for iteration in est_abs.index:
        tn, fp, fn, tp = confusion_matrix(key_abs, est_abs.loc[iteration]).ravel()
        if key.sum() == 0: # i.e., a negative control 
            sensitivity.append(0)  # how many positives were selected as positive?
            specificity.append(tn / (tn+fp))  # how many negatives are truly negative?
        else:
            sensitivity.append(tp / (tp+fn))  # how many positives were selected as positive?
            specificity.append(tn / (tn+fp))  # how many negatives are truly negative?
    return sensitivity, specificity

def score_app(est, key):
    # make est and key into composites of total fecal load, e.g., 1A + 2B
    mae = []
    mse = []
    rmse = []
    
    if key.sum() == 0: # i.e., a negative control
        key_app = key
    else:
        key_app = key / key.sum()
    for iteration in est.index:
        est_app = est.loc[iteration] / est.loc[iteration].sum()
        mae.append(mean_absolute_error(key_app, est_app))     
        mse.append(mean_squared_error(key_app, est_app))      
        rmse.append(root_mean_squared_error(key_app, est_app))      
    return mae, mse, rmse

def score_frac(est, key):
    # score as-is -- stretch goal of the work 
    mae = []
    mse = []
    rmse = []
    for iteration in est.index:
        mae.append(mean_absolute_error(key, est.loc[iteration]))
        mse.append(mean_squared_error(key, est.loc[iteration]))
        rmse.append(root_mean_squared_error(key, est.loc[iteration]))
    return mae, mse, rmse

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
    parser.add_argument(
        '-r', '--relabd',
        help='Are the input values percentages?',
        action="store_true",
        required=False
        )
    args=vars(parser.parse_args())

    print('Scoring SourceApp for sample: ' + args['sample'])
    est = pd.read_csv(args["input"],header=0,index_col=0).T
    key = pd.read_csv(args["key"],header=0,index_col=0).loc[args["sample"]]
    sources = est.columns[0:-3][0::2].tolist()
    est_nocrx=pd.DataFrame()
    est_crx = pd.DataFrame()
    for source in sources:
        est_nocrx[source]=est[source]
        est_crx[source]=est[source] + est[source + "_crx"]

    if args["relabd"]: 
        est_nocrx = est_nocrx / 100
        est_crx =  est_crx / 100

    print('Scoring attribution ...\n')
    sensitivity, specificity = score_att(est_nocrx, key)
    est["att_nocrx_sens"] = sensitivity
    est["att_nocrx_spec"] = specificity
    
    sensitivity, specificity = score_att(est_crx, key)
    est["att_wcrx_sens"] = sensitivity
    est["att_wcrx_spec"] = specificity

    print('Scoring apportioning ...\n')
    mae, mse, rmse = score_app(est_nocrx, key)
    est["app_nocrx_mae"] = mae
    est["app_nocrx_mse"] = mse
    est["app_nocrx_rmse"] = rmse
    
    mae, mse, rmse = score_app(est_crx, key)
    est["app_wcrx_mae"] = mae
    est["app_wcrx_mse"] = mse
    est["app_wcrx_rmse"] = rmse

    print('Scoring cell fractioning ...\n')
    mae, mse, rmse = score_frac(est_nocrx, key)
    est["frac_nocrx_mae"] = mae
    est["frac_nocrx_mse"] = mse
    est["frac_nocrx_rmse"] = rmse
    
    mae, mse, rmse = score_frac(est_crx, key)
    est["frac_wcrx_mae"] = mae
    est["frac_wcrx_mse"] = mse
    est["frac_wcrx_rmse"] = rmse

    est["min_frac"] = key[key>0].min()

    est["sample"] = args["sample"]

    est.to_csv(args["output"]+".scores.csv",sep=",",index=False)

    print("Parameter sets found for the following...", flush=True)
    
    row_ind = est["app_wcrx_rmse"].idxmin() 
    print("Apportioning with crx:", flush=True)
    print(est.loc[row_ind]["percent_identity"], flush=True)
    print(est.loc[row_ind]["query_coverage"], flush=True)
    print(est.loc[row_ind]["limit_threshold"], flush=True)

    row_ind = est["app_nocrx_rmse"].idxmin()
    print("Apportioning without crx:", flush=True)
    print(est.loc[row_ind]["percent_identity"], flush=True)
    print(est.loc[row_ind]["query_coverage"], flush=True)
    print(est.loc[row_ind]["limit_threshold"], flush=True)

    row_ind = est["frac_wcrx_rmse"].idxmin() 
    print("Fractioning with crx:", flush=True)
    print(est.loc[row_ind]["percent_identity"], flush=True)
    print(est.loc[row_ind]["query_coverage"], flush=True)
    print(est.loc[row_ind]["limit_threshold"], flush=True)

    row_ind = est["frac_nocrx_rmse"].idxmin() 
    print("Fractioning without crx:", flush=True)
    print(est.loc[row_ind]["percent_identity"], flush=True)
    print(est.loc[row_ind]["query_coverage"], flush=True)
    print(est.loc[row_ind]["limit_threshold"], flush=True)
    
    est.min

if __name__ == "__main__":
    main()
