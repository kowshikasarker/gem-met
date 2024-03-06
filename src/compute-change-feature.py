import argparse
import pandas as pd
from pathlib import Path
import sys
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--react_set_path", type=str,
                        help="path to reaction set in .tsv format",
                        required=True, default=None)
    parser.add_argument("-c", "--met_change_path", type=str,
                        help="path to metabolite concentration change in .tsv format", required=True, default=None)
    parser.add_argument("-m", "--met_path", type=str,
                        help="path to list of human-gem overlapped metabolites in .tsv format", required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    parser.add_argument("-l", "--log_path", type=str,
                        help="path to log file",
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
    print('met_path', args.met_path)
    print('out_dir', args.out_dir)
    print('log_path', args.log_path)
    
    
    change_df = pd.read_csv(args.met_change_path, sep='\t')
    
    change_df[['person_id', 'diet']] = change_df['key'].str.split('.', expand=True)
    
    treatments = set(change_df.diet)
    treatments = [treatment for treatment in treatments if not treatment.startswith('No')]
    print('treatments', treatments)
    
    
    
    change_df = change_df.drop(columns=['person_id', 'diet']).set_index('key')
    
    
    change_df = change_df.sort_index().sort_index(axis=1)

    id_df = pd.read_csv(args.met_path, sep='\t')
    met_to_hmdb = dict(zip(id_df.MET_ID, id_df.HMDB_ID))
    
    

    react_set_df = pd.read_csv(args.react_set_path, sep='\t')

    react_set_df['Measured_Substrate'] = react_set_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_set_df['Measured_Product'] = react_set_df['Measured_Product'].apply(lambda x: str_to_set(x))

    #change_dir = args.out_dir +
    #pathlib.Path(change_dir).mkdir(parents=True, exist_ok=True)
    
    

    for treatment in treatments:
        control = 'No' + treatment
        study_change_df = change_df[change_df.index.str.contains(treatment)]
        print('study_change_df', study_change_df.shape)
        # rc1
        #rc1_df = study_change_df.copy().drop(columns=study_change_df.columns)
        rc1_df_dict = {}
        for idx, rxn in react_set_df.iterrows():
            rc1_col = 0
            for s in rxn['Measured_Substrate']:
                if(not s in met_to_hmdb):
                    print('s', s, idx, rxn['Measured_Substrate'], rxn)
                rc1_col = rc1_col - study_change_df[met_to_hmdb[s]]
            for p in rxn['Measured_Product']:
                if(not p in met_to_hmdb):
                    print('p', p, idx, rxn['Measured_Product'], rxn)
                rc1_col = rc1_col + study_change_df[met_to_hmdb[p]]
            rc1_df_dict[rxn['RXN_ID']] = rc1_col
            #rc1_df[rxn['RXN_ID']] = rc1_col
        rc1_df = pd.DataFrame(rc1_df_dict)
        rc1_df.index = study_change_df.index

        rc1_path = args.out_dir + '/' + treatment + '.' + control + '.change.tsv'
        rc1_df.to_csv(rc1_path, sep='\t', index=True)
        
    sys.stdout = orig_stdout
    log_file.close()
    
if __name__ == "__main__":
    main(parse_args())
