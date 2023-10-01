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

def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    change_df = pd.read_csv(args.met_change_path, sep='\t', index_col='Key')
    
    change_df = change_df.sort_index().sort_index(axis=1)

    id_df = pd.read_csv(args.met_75_path, sep='\t')
    met_to_hmdb = dict(zip(id_df.MET_ID, id_df.HMDB_ID))

    treatments = ['Walnut', 'Broccoli', 'Barley', 'Oats', 'Avocado', 'Almond']

    react_set_df = pd.read_csv(args.react_set_path, sep='\t')

    react_set_df['Measured_Substrate'] = react_set_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_set_df['Measured_Product'] = react_set_df['Measured_Product'].apply(lambda x: str_to_set(x))

    change_dir = args.out_dir + '/Change'
    pathlib.Path(change_dir).mkdir(parents=True, exist_ok=True)

    for treatment in treatments:
        control = 'No' + treatment
        study_change_df = change_df[change_df.index.str.contains(treatment)]
        print('study_change_df', study_change_df.shape)
        # rc1
        rc1_df = study_change_df.copy().drop(columns=study_change_df.columns)
        for idx, rxn in react_set_df.iterrows():
            rc1_col = 0
            for s in rxn['Measured_Substrate']:

                rc1_col = rc1_col - study_change_df[met_to_hmdb[s]]
            for p in rxn['Measured_Product']:
                rc1_col = rc1_col + study_change_df[met_to_hmdb[p]]

            rc1_df[rxn['RXN_ID']] = rc1_col

        rc1_path = change_dir + '/' + treatment + '.' + control + '.change.tsv'
        rc1_df.to_csv(rc1_path, sep='\t', index=True)
    
if __name__ == "__main__":
    main(parse_args())