import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--react_set_path", type=str,
                        help="path to reaction set in tsv format",
                        required=True, default=None)
    parser.add_argument("-c", "--met_change_path", type=str,
                        help="path to metabolite concentration change in tsv format", required=True, default=None)
    parser.add_argument("-m", "--met_75_path", type=str,
                        help="path to list of 75 metabolites in tsv format", required=True, default=None)
    parser.add_argument("-i", "--id_path", type=str,
                        help="path to metabolite identifiers in tsv format", required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
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
    '''print(set1)
    print(set2)
    print(product)'''
    return product

def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    change_df = pd.read_csv(args.met_change_path, sep='\t', index_col='Key')
    
    change_df = change_df.sort_index().sort_index(axis=1)

    id_df = pd.read_csv(args.id_path, sep='\t')
    
    met_to_hmdb = dict(zip(id_df.MET_ID, id_df.HMDB_ID))

    react_set_df = pd.read_csv(args.react_set_path, sep='\t')
    react_set_df['Measured_Substrate'] = react_set_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_set_df['Measured_Product'] = react_set_df['Measured_Product'].apply(lambda x: str_to_set(x))

    react_set_df['Measured_Product_Substrate_Pair'] = react_set_df.apply(lambda x: cross_product(x['Measured_Product'], x['Measured_Substrate']), axis=1)
    react_set_df.to_csv(react_set_dir + '/' + react_set.lower() + '.er.tsv', sep='\t', index=False)
    react_col = 0
    
    ratio_dir = react_set_dir + '/Ratio'
    pathlib.Path(er1_dir).mkdir(parents=True, exist_ok=True)
    print('ratio_dir', er1_dir)
    
    
    for treatment in treatments:
        control = 'No' + treatment
        study_change_df = change_df[change_df.index.str.contains(treatment)]
        print('study_change_df', study_change_df.shape)

    
        er1_df = study_change_df.copy().drop(columns=study_change_df.columns)
    
        for idx, rxn in react_set_df.iterrows():
            er1_col = 0
            for pair in rxn['Measured_Product_Substrate_Pair']:
                hmdb_product = met_to_hmdb[pair[0]]
                hmdb_substrate = met_to_hmdb[pair[1]]
                er1_col = er1_col + (study_change_df[hmdb_product] / study_change_df[hmdb_substrate])
            er1_df[rxn['RXN_ID']] = er1_col

        er1_df.to_csv(ratio_dir + '/' + treatment + '.' + control + '.ratio.tsv', sep='\t', index=True)

if __name__ == "__main__":
    main(parse_args())