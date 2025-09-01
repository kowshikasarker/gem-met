import argparse
import pandas as pd
from pathlib import Path
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--react_set_path", type=str,
                        help="path to reaction set in .tsv format",
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
    cell = set(cell.split(', '))
    if('set()' in cell):
        cell.remove('set()')
    return cell

def cross_product(set1, set2):
    product = []
    for s1 in set1:
        for s2 in set2:
            product.append((s1, s2))
    return product

def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    sys.stderr = log_file
    
    change_df = pd.read_csv(args.met_change_path, sep='\t')
    change_df[['person_id', 'diet']] = change_df['key'].str.split('.', expand=True)
    
    treatments = set(change_df.diet)
    treatments = [treatment for treatment in treatments if not treatment.startswith('No')]
    print('treatments', treatments)
    
    change_df = change_df.drop(columns=['person_id', 'diet']).set_index('key')
    
    change_df = change_df.sort_index().sort_index(axis=1)

    id_df = pd.read_csv(args.valid_met_path, sep='\t')
    met_to_hmdb = dict(zip(id_df.MET_ID, id_df.ID))

    react_set_df = pd.read_csv(args.react_set_path, sep='\t')
    react_set_df['Measured_Substrate'] = react_set_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_set_df['Measured_Product'] = react_set_df['Measured_Product'].apply(lambda x: str_to_set(x))

    react_set_df['Measured_Product_Substrate_Pair'] = react_set_df.apply(lambda x: cross_product(x['Measured_Product'], x['Measured_Substrate']), axis=1)
    
    react_col = 0
        
    for treatment in treatments:
        control = 'No' + treatment
        study_change_df = change_df[change_df.index.str.contains(treatment)]
        print('study_change_df', study_change_df.shape)
        er1_df_dict = {}
        for idx, rxn in react_set_df.iterrows():
            er1_col = 0
            for pair in rxn['Measured_Product_Substrate_Pair']:
                if not (pair[0] in met_to_hmdb):
                    print('pair[0]', pair[0], rxn['Measured_Product_Substrate_Pair'])
                hmdb_product = met_to_hmdb[pair[0]]
                if not (pair[1] in met_to_hmdb):
                    print('pair[1]', pair[1], rxn['Measured_Product_Substrate_Pair'])
                hmdb_substrate = met_to_hmdb[pair[1]]
                er1_col = er1_col + (study_change_df[hmdb_product] / study_change_df[hmdb_substrate])
            er1_df_dict[rxn['RXN_ID']] = er1_col

        er1_df = pd.DataFrame(er1_df_dict)
        er1_df.to_csv(args.out_dir + '/' + treatment + '.' + control + '.ratio.tsv', sep='\t', index=True)
        
    sys.stdout = orig_stdout
    sys.stderr = orig_stderr
    log_file.close()

if __name__ == "__main__":
    main(parse_args())
