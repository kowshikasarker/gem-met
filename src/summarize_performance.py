import argparse
import sys
import os
import json
import dataframe_image as dfi
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--working_dir", type=str,
                        help="path to the directory containing features in .tsv format computed by the script preprocess_and_compute_features.py",
                        required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="path to the output directory",
                        required=True, default=None)
    parser.add_argument("-l", "--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    args = parser.parse_args()
    return args

def color_cells(x):
    if(x >= 0.1):
        return 'background-color: darkturquoise'
    if(x > 0):
        return 'background-color: paleturquoise'
    if(x == 0):
        return 'background-color: khaki'
    if(x <= -0.1):
        return 'background-color: coral'
    if(x < 0):
        return 'background-color: mistyrose'


def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    sys.stderr = log_file
    
    treatments = ['Almond', 'Avocado', 'Barley', 'Broccoli', 'Oats', 'Walnut']
    metrics = ['Accuracy', 'AUROC', 'AUPRC']
    react_feat_map = {
        1: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        2: ['Change', 'Prob/Reaction', 'Prob/Subsystem', 'Ratio'],
        3: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        4: ['Change', 'Prob/Reaction', 'Prob/Subsystem', 'Ratio'],
        5: ['Prob/Reaction', 'Prob/Subsystem'],
        6: ['Prob/Reaction', 'Prob/Subsystem'],
        7: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        8: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        9: ['Prob/Reaction', 'Prob/Subsystem']
    }
    
    summary_dict = {'index': [],
        'columns': [],
        'data': [],
        'index_names': ['Reaction Set', 'Input Features'],
        'column_names': ['Treatment', 'Metric']}
    
    for i in range(len(treatments)):
        for j in range(len(metrics)):
            summary_dict['columns'].append((treatments[i], metrics[j]))
    print(summary_dict['columns'])
    metrics = ['accuracy', 'auroc', 'auprc']
    for react_set_no in range(1, 10):
        react_set_dir = args.working_dir + '/Reaction-Set-' + str(react_set_no)
        print()
        
        for feature_path in react_feat_map[react_set_no]:
            feature_dir = react_set_dir + '/' + feature_path
            feature_type = (feature_path.split('/')[0]).lower()
            print(feature_type)
            summary_dict['index'].append(('Reaction-Set-' + str(react_set_no), feature_path.replace('/', '-')))
            data = []
            for i in range(len(treatments)):
                treatment = treatments[i]
                control = 'No' + treatment
                treatment_dir = feature_dir + '/' + treatment
                json_file = treatment_dir + '/' + treatment + '.' + control + '.classifier_performance.json'
                metric_dict = json.load(open(json_file))
                for j in range(len(metrics)):
                    data.append(round(metric_dict[metrics[j]], 2))
            summary_dict['data'].append(data)
    
    summary_df = pd.DataFrame.from_dict(summary_dict, orient='tight')
    
    summary_json = summary_df.to_json(orient="split")
    json_file = open(args.out_dir + '/summary.json', 'w')
    json.dump(json.loads(summary_json), json_file, indent=4)
    json_file.close()
    
    baseline = [0.79, 0.81, 0.82, 0.61, 0.56, 0.69, 0.54, 0.57, 0.59, 0.43, 0.45, 0.49, 0.45, 0.45, 0.58, 0.94, 0.98, 0.99]
    diff = summary_df.subtract(baseline)
    summary_df.loc[('Baseline (Met-75)', 'ΔConcentration'), :] = baseline
    
    row_idx = summary_df.index[:-1]
    print(row_idx)
    col_idx = summary_df.columns
    print(col_idx)
    subset = pd.IndexSlice[row_idx, col_idx]

    df_styled = summary_df.style.format(precision=2).set_table_styles(
    [{"selector": "", "props": [("border", "1px solid grey"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap')]},
      {"selector": "tbody td", "props": [("border", "1px solid grey"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap')]},
     {"selector": "th", "props": [("border", "1px solid grey"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap')]},
    ]).set_properties(**{'max-width': '100%', 'white-space': 'nowrap'}).apply(lambda x: diff.applymap(color_cells), axis=None).highlight_max(axis=0, props='font-weight: bold; color: blue', subset=subset)

    dfi.export(df_styled, args.out_dir + '/summary.png')

    sys.stdout = orig_stdout
    sys.stderr = orig_stderr 
    log_file.close()
                
if __name__ == "__main__":
    main(parse_args())