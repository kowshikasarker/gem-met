import argparse
import sys
import os
import json
import dataframe_image as dfi
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--case", type=str,
                        required=True, default=None)
    parser.add_argument("--control", type=str,
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--log_path", type=str,
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
    
    json_filename = args.case + "." + args.control + ".classifier_performance.json"
    
    baseline_dict = json.load(open(args.in_dir + '/metabolite/' + json_filename))
    baseline = [round(baseline_dict['accuracy'], 2),
                round(baseline_dict['auroc'], 2),
                round(baseline_dict['auprc'], 2)]
    
    print("baseline", baseline)
    
    react_feat_map = {
        1: ['Change', 'Prob'],
        2: ['Change', 'Ratio', 'Prob'],
        3: ['Change', 'Prob'],
        4: ['Change', 'Ratio', 'Prob'],
        5: ['Prob'],
        6: ['Prob'],
        7: ['Change', 'Ratio', 'Prob'],
        8: ['Change', 'Ratio', 'Prob'],
        9: ['Prob']
    }
    
    #reaction
    
    summary_dict = {'index': ['Accuracy', 'AUROC', 'AUPRC'],
        'columns': [],
        'data': [],
        'index_names': [''],
        'column_names': ['Reaction set', 'Feature']}
    
    for r in range (1, 10):
        feat_map = react_feat_map[r]
        for f in range(len(feat_map)):
            summary_dict['columns'].append((r, feat_map[f]))
    summary_dict['columns'].append(('Baseline', 'Metabolite'))
            
    metrics = ['accuracy', 'auroc', 'auprc']
    
    for m in range(len(metrics)):
        metric = metrics[m]
        print(m, metric)
        data = []
        for react_set_no in range(1, 10):
            react_set_dir = args.in_dir + '/reaction/reaction-set-' + str(react_set_no)
            print('react_set_no', react_set_no)
            feat_map = react_feat_map[react_set_no]
            for f in range(len(feat_map)):
                feature_type = feat_map[f].lower()
                print('feature_type', feature_type)
                feature_dir = react_set_dir + '/' + feature_type
                
                json_file = feature_dir + '/' + json_filename
                metric_dict = json.load(open(json_file))
                data.append(round(metric_dict[metric], 2))
        data.append(baseline[m])
        summary_dict['data'].append(data)
    print('columns', summary_dict['columns'])
    print('summary dict created')
    print(summary_dict)
    
    summary_df = pd.DataFrame.from_dict(summary_dict, orient='tight')
    summary_df.index.name = None
    summary_json = summary_df.to_json(orient="split")
    json_file = open(args.out_dir + '/summary.json', 'w')
    json.dump(json.loads(summary_json), json_file, indent=4)
    json_file.close()
    

    diff = summary_df.subtract(baseline, axis=0)
    
    row_idx = summary_df.index
    print(row_idx)
    col_idx = summary_df.columns[:-1]
    print(col_idx)
    subset = pd.IndexSlice[row_idx, col_idx]

    df_styled = summary_df.style.format(precision=2).set_table_styles(
    [{"selector": "", "props": [("border", "1px solid black"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap')]},
      {"selector": "tbody td", "props": [("border", "1px solid black"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap'), ('margin', '5'), ('padding', '5'), ('font-size', 'large')]},
     {"selector": "th", "props": [("border", "1px solid black"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap'), ('margin', '5'), ('padding', '5'), ('font-size', 'large'), ('font-weight', 'normal'), ('white-space', 'pre')]},
    ]).set_properties(**{'max-width': '100%', 'white-space': 'nowrap'}).apply(lambda x: diff.applymap(color_cells), axis=None).highlight_max(axis=1, props='font-weight: bold; color: blue', subset=subset)

    dfi.export(df_styled, args.out_dir + '/reaction.png', table_conversion="selenium")
    
    
    
    
    #metabolite+reaction
    summary_dict = {'index': ['Accuracy', 'AUROC', 'AUPRC'],
        'columns': [],
        'data': [],
        'index_names': [''],
        'column_names': ['Reaction set', 'Feature']}
    
    for r in range (1, 10):
        feat_map = react_feat_map[r]
        for f in range(len(feat_map)):
            summary_dict['columns'].append((r, feat_map[f]))
    summary_dict['columns'].append(('Baseline', 'Metabolite'))
            
    metrics = ['accuracy', 'auroc', 'auprc']
    
    for m in range(len(metrics)):
        metric = metrics[m]
        print(m, metric)
        data = []
        for react_set_no in range(1, 10):
            react_set_dir = args.in_dir + '/metabolite+reaction/reaction-set-' + str(react_set_no)
            print('react_set_no', react_set_no)
            feat_map = react_feat_map[react_set_no]
            for f in range(len(feat_map)):
                feature_type = feat_map[f].lower()
                print('feature_type', feature_type)
                feature_dir = react_set_dir + '/' + feature_type
                
                json_file = feature_dir + '/' + json_filename
                metric_dict = json.load(open(json_file))
                data.append(round(metric_dict[metric], 2))
        data.append(baseline[m])
        summary_dict['data'].append(data)
    print('columns', summary_dict['columns'])
    print('summary dict created')
    print(summary_dict)
    
    summary_df = pd.DataFrame.from_dict(summary_dict, orient='tight')
    summary_df.index.name = None
    summary_json = summary_df.to_json(orient="split")
    json_file = open(args.out_dir + '/summary.json', 'w')
    json.dump(json.loads(summary_json), json_file, indent=4)
    json_file.close()
    

    diff = summary_df.subtract(baseline, axis=0)
    
    row_idx = summary_df.index
    print(row_idx)
    col_idx = summary_df.columns[:-1]
    print(col_idx)
    subset = pd.IndexSlice[row_idx, col_idx]

    df_styled = summary_df.style.format(precision=2).set_table_styles(
    [{"selector": "", "props": [("border", "1px solid black"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap')]},
      {"selector": "tbody td", "props": [("border", "1px solid black"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap'), ('margin', '5'), ('padding', '5'), ('font-size', 'large')]},
     {"selector": "th", "props": [("border", "1px solid black"), ('text-align', 'center'), ('max-width', '100%'), ('white-space', 'nowrap'), ('margin', '5'), ('padding', '5'), ('font-size', 'large'), ('font-weight', 'normal'), ('white-space', 'pre')]},
    ]).set_properties(**{'max-width': '100%', 'white-space': 'nowrap'}).apply(lambda x: diff.applymap(color_cells), axis=None).highlight_max(axis=1, props='font-weight: bold; color: blue', subset=subset)

    dfi.export(df_styled, args.out_dir + '/metabolite+reaction.png', table_conversion="selenium")

    sys.stdout = orig_stdout
    sys.stderr = orig_stderr 
    log_file.close()
                
if __name__ == "__main__":
    main(parse_args())
