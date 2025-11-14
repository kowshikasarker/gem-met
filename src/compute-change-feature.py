import argparse
import pandas as pd
from pathlib import Path
import sys
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--react_set_path", type=str,
                        help="path to a reaction set in .tsv format",
                        required=True, default=None)
    parser.add_argument("--met_change_path", type=str,
                        help="path to metabolite concentration change in .tsv format", required=True, default=None)
    parser.add_argument("--valid_met_path", type=str,
                        help="path to list of human-gem overlapped metabolites in .tsv format", required=True, default=None)
    
    parser.add_argument("--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    args = parser.parse_args()
    return args

def str_to_set(cell):
    cell = ''.join(c for c in cell if c not in "'{}")
    if(len(cell) == 0):
        return {}
    cell = set(cell.split(', '))
    if('set()' in cell):
        cell.remove('set()')
    return cell

def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    
    print('react_set_path', args.react_set_path)
    print('met_change_path', args.met_change_path)
    print('valid_met_path', args.valid_met_path)
    print('log_path', args.log_path)
    print('out_dir', args.out_dir)
    
    change_df = pd.read_csv(args.met_change_path, sep='\t', index_col='key')
       
    change_df = change_df.sort_index().sort_index(axis=1)

    id_df = pd.read_csv(args.valid_met_path, sep='\t')
    met_to_id = dict(zip(id_df.MET_ID, id_df.ID))
    
    

    react_set_df = pd.read_csv(args.react_set_path, sep='\t')

    react_set_df['Measured_Substrate'] = react_set_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_set_df['Measured_Product'] = react_set_df['Measured_Product'].apply(lambda x: str_to_set(x))

    rc1_df_dict = {}
    for idx, rxn in react_set_df.iterrows():
        rc1_col = 0
        print(rxn)
        for s in rxn['Measured_Substrate']:
            if(not s in met_to_id):
                print('s', s, idx, rxn['Measured_Substrate'], rxn)
            rc1_col = rc1_col - change_df[met_to_id[s]]
        for p in rxn['Measured_Product']:
            if(not p in met_to_id):
                print('p', p, idx, rxn['Measured_Product'], rxn)
            rc1_col = rc1_col + change_df[met_to_id[p]]
        rc1_df_dict[rxn['RXN_ID']] = rc1_col
    rc1_df = pd.DataFrame(rc1_df_dict)
    rc1_df.index = change_df.index
    
    rc1_df = rc1_df.iloc[:, :2]
    
    df1 = rc1_df.reset_index(names='index')
    sample_split = df1['index'].str.rsplit(":", n=1, expand=True)
    if sample_split.shape[1] != 2:
        raise ValueError(
            "Unexpected sample index format encountered while splitting into sample_id and sample_group"
        )
    df1[['sample_id', 'sample_group']] = sample_split
    df1 = df1.set_index(['sample_id', 'sample_group'])
    df1 = df1.drop(columns='index')
    df1.to_csv(args.out_dir + '/reaction.change.tsv', sep='\t', index=True)
    
    df2 = pd.concat([change_df, rc1_df], axis=1)
    df2 = df2.reset_index(names='index')
    sample_split = df2['index'].str.rsplit(":", n=1, expand=True)
    if sample_split.shape[1] != 2:
        raise ValueError(
            "Unexpected sample index format encountered while splitting into sample_id and sample_group"
        )
    df2[['sample_id', 'sample_group']] = sample_split
    df2 = df2.set_index(['sample_id', 'sample_group'])
    df2 = df2.drop(columns='index')
    df2.to_csv(args.out_dir + '/metabolite.reaction.change.tsv', sep='\t', index=True)
        
    sys.stdout = orig_stdout
    log_file.close()
    
if __name__ == "__main__":
    main(parse_args())
