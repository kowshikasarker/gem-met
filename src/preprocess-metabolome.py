import argparse
import sys
import pathlib
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocessing of user-provided baseline (before diet) and end (after diet) metabolomic profiles (expected file format: .tsv).' + \
        '\n' + 'Two columns sample_id and sample_group are expected to have sample id and sample groups (case/control).' + \
        '\n' + 'Other columns are expected to have metabolite abundances. These columns may have arbitrary names.' + \
        '\n' + 'Rows or columns absent in any profile are discarded from the other profile.' + \
        '\n' + 'Columns are discarded based on missing value proportion too.' + \
        '\n' + 'In the output directory, 10 .tsv files are generated.')
    
    parser.add_argument("--base_path", type=str,
                        help="path to baseline metabolomics in tsv format",
                        required=True, default=None)
    parser.add_argument("--end_path", type=str,
                        help="path to end metabolomics in tsv format",
                        required=True, default=None)
    
    parser.add_argument("--missing_pct", type=float,
                        help="threshold of missing value percentage for dropping columns",
                        required=True, default=None)
    
    parser.add_argument("--user_met_id_path", type=str,
                        help="path to list of user metabolites (columns in in base_path and end_path) with standard identifiers in tsv format",
                        required=True, default=None)
    parser.add_argument("--user_met_name_col", type=str,
                        help="which column in user_met_id_path contains the metabolite names",
                        required=True, default=None)
    parser.add_argument("--user_met_id_col", type=str,
                        help="which column in user_met_id_path contains the metabolite standard identifiers",
                        required=True, default=None)
    
    parser.add_argument("--gem_path", type=str,
                        help="path to a gem model in .xlsx format",
                        required=True, default=None)
    parser.add_argument("--gem_met_id_path", type=str,
                        help="path to list of human gem metabolites with standard identifiers list in tsv format",
                        required=True, default=None)
    parser.add_argument("--gem_met_id_col", type=str,
                        help="which column in gem_met_id_path contains the metabolite standard identifiers",
                        required=True, default=None)
    
    parser.add_argument("--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    
    args = parser.parse_args()
    return args

# Imputes missing values to uniform random values between [0, mm * minimum observed] for every feature
def impute_missing_values(base_df, end_df, coeff):

    # Compute per-feature minimums for dataset
    base_df_feature_mins = np.min(base_df, axis=0)
    base_df_nan_dict = {}
    end_df_feature_mins = np.min(end_df, axis=0)
    end_df_nan_dict = {}
    
    print('base_df_feature_mins', base_df_feature_mins)
    print('end_df_feature_mins', end_df_feature_mins)

    # Create new datasets that contains random values for each subject for each feature,
    # between 0 and mm * the minimum for that feature
    for feature, minimum in base_df_feature_mins.items():
        base_df_nan_dict[feature] = np.random.uniform(
            low=0, high=coeff*minimum, size=len(base_df)
        )

    for feature, minimum in end_df_feature_mins.items():
        end_df_nan_dict[feature] = np.random.uniform(
            low=0, high=coeff*minimum, size=len(end_df)
        )

    # Update original dataset with new values for any missing entries
    # Original values should be preserved
    base_df_nan = pd.DataFrame(base_df_nan_dict)
    base_df_nan.index = base_df.index

    end_df_nan = pd.DataFrame(end_df_nan_dict)
    end_df_nan.index = end_df.index

    base_df.update(base_df_nan, overwrite=False)
    end_df.update(end_df_nan, overwrite=False)

    return base_df, end_df


def main(args):
    pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    
    print('base_path', args.base_path)
    print('end_path', args.end_path)
    print('missing_pct', args.missing_pct)
    print('user_met_id_path', args.user_met_id_path)
    print('user_met_name_col', args.user_met_name_col)
    print('user_met_id_col', args.user_met_id_col)
    print('gem_path', args.gem_path)
    print('gem_met_id_path', args.gem_met_id_path)
    print('gem_met_id_col', args.gem_met_id_col)
    print('log_path', args.log_path)
    print('out_dir', args.out_dir)
    
    base_df = pd.read_csv(args.base_path, sep='\t')
    print(base_df.shape)
    end_df = pd.read_csv(args.end_path, sep='\t')
    print(end_df.shape)
    
    base_df.columns = map(str.lower, base_df.columns)
    end_df.columns = map(str.lower, end_df.columns)
    
    base_df['key'] = base_df['sample_id'].astype(str) + ':' + base_df['sample_group'].astype(str)
    end_df['key'] = end_df['sample_id'].astype(str) + ':' + end_df['sample_group'].astype(str)
    
    base_df = base_df.set_index('key')
    end_df = end_df.set_index('key')
    
    common_rows = set(base_df.index).intersection(set(end_df.index))
    print('common_rows', len(common_rows), common_rows)
    
    base_df = base_df[base_df.index.isin(common_rows)]
    end_df = end_df[end_df.index.isin(common_rows)]
     
    common_cols = list(set(base_df.columns).intersection(set(end_df.columns)))
    print('common_cols', len(common_cols))
    
    base_df = base_df[common_cols]
    end_df = end_df[common_cols]
    
    print('base_df', base_df.shape, 'end_df', end_df.shape)
    
    base_missing = base_df.isnull().mean()
    end_missing = end_df.isnull().mean()
    
    print('base_missing', base_missing)
    print('end_missing', end_missing)
    
    base_drop_cols = base_missing[base_missing >= args.missing_pct].index
    end_drop_cols = end_missing[end_missing >= args.missing_pct].index
    
    print('base_drop_cols', len(base_drop_cols), base_drop_cols)
    print('end_drop_cols', len(end_drop_cols), end_drop_cols)
    
    drop_cols = list(set(base_drop_cols).union(set(end_drop_cols))) + ['sample_id', 'sample_group']
    print('drop_cols', len(drop_cols), drop_cols)
    
    base_df = base_df.drop(columns=drop_cols)
    end_df = end_df.drop(columns=drop_cols)
    
    print('base_df', base_df.shape, 'end_df', end_df.shape)
    
    base_df, end_df = impute_missing_values(base_df, end_df, coeff=0.25)
    print(base_df.shape, end_df.shape)    
    
    common_cols = list(set(base_df.columns).intersection(set(end_df.columns)))
    print('common_cols', len(common_cols))
    base_df = base_df[common_cols]
    end_df = end_df[common_cols]
    
    print('base_df',base_df.shape, 'end_df', end_df.shape)
    
    base_df.to_csv(args.out_dir + '/preprocessed_base_name.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/preprocessed_end_name.tsv', sep='\t')
    
    base_df = base_df.sort_index(axis=0).sort_index(axis=1)
    end_df = end_df.sort_index(axis=0).sort_index(axis=1)
    
    print('base_df',base_df.shape, 'end_df', end_df.shape)
    
    change_df = end_df.sub(base_df)
    print('change_df', change_df.shape)
    change_df.to_csv(args.out_dir + '/preprocessed_change_name.tsv', sep='\t')
    
    id_df = pd.read_csv(args.user_met_id_path, sep='\t', usecols=[args.user_met_name_col, args.user_met_id_col])
    id_df[args.user_met_name_col] = id_df[args.user_met_name_col].apply(lambda x: x.lower())
    id_df = id_df[id_df[args.user_met_name_col].isin(common_cols)]
    id_df = id_df.dropna(subset=[args.user_met_id_col])
    
    print('id_df', id_df.shape)
        
    name_to_id = dict(zip(id_df[args.user_met_name_col], id_df[args.user_met_id_col]))
    print('name_to_id', name_to_id)
    
    id_to_name = dict(zip(id_df[args.user_met_id_col], id_df[args.user_met_name_col]))
    print('id_to_name', id_to_name)
    
    common_cols_id = set([name_to_id[col] for col in common_cols if col in name_to_id.keys()])
    print(len(common_cols_id))
    
    gem_met = pd.read_csv(args.gem_met_id_path, sep='\t')
    gem_met = gem_met.dropna(subset=[args.gem_met_id_col])
    print(gem_met.shape)
    
    common_mets = common_cols_id.intersection(set(gem_met[args.gem_met_id_col]))
    print('common_mets', len(common_mets))
    
    common_mets_names = [id_to_name[met] for met in common_mets]
    print(len(common_mets_names))
    print(common_mets_names)
    
    base_df = base_df[common_mets_names]
    end_df = end_df[common_mets_names]
    print(base_df.shape, end_df.shape)
    
    base_df.to_csv(args.out_dir + '/gem_overlapped_base_name.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/gem_overlapped_end_name.tsv', sep='\t')
    
    base_df = base_df.rename(columns=name_to_id)
    end_df = end_df.rename(columns=name_to_id)
    
    base_df.to_csv(args.out_dir + '/gem_overlapped_base_id.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/gem_overlapped_end_id.tsv', sep='\t')
     
    base_df = base_df.sort_index(axis=0).sort_index(axis=1)
    end_df = end_df.sort_index(axis=0).sort_index(axis=1)
    
    change_df = end_df.sub(base_df)
    
    change_df.to_csv(args.out_dir + '/gem_overlapped_change_id.tsv', sep='\t')
    change_df = change_df.rename(columns=id_to_name)
    change_df.to_csv(args.out_dir + '/gem_overlapped_change_name.tsv', sep='\t')
    print('change_df', change_df.shape)
    
    gem_overlapeed_met_names = list(change_df.columns)
    gem_met = gem_met[gem_met[args.gem_met_id_col].isin(id_to_name.keys())]
    id_to_mam = dict(zip(gem_met[args.gem_met_id_col], gem_met.metsNoComp))
    print('id_to_mam', id_to_mam)
    
    gem_model = pd.read_excel(args.gem_path, sheet_name='METS', usecols=['NAME', 'REPLACEMENT ID'])
    print('gem_path', args.gem_path)
    #pattern = '|'.join(['e', 'x', 'm', 'c', 'l', 'r', 'g', 'n', 'i'])
    #gem_model['MAM_ID'] = gem_model['REPLACEMENT ID'].str.replace(pattern, '')
    gem_model['MAM_ID'] = gem_model['REPLACEMENT ID'].str[:-1]

    print('gem_model', gem_model.shape)
    print(list(id_to_mam.values()))
    print(gem_model['MAM_ID'])
    gem_model = gem_model[gem_model['MAM_ID'].isin(id_to_mam.values())]
    print('gem_model', gem_model.shape)
    mam_to_met = dict(zip(gem_model['MAM_ID'], gem_model['NAME']))
    print('mam_to_met', mam_to_met)
    
    
    file = open(args.out_dir + '/gem_overlapped_metabolites.tsv', 'w')
    file.write('Name' + '\t' + 'ID' + '\t' + 'MAM_ID' + '\t' + 'MET_ID' + '\n')
    for name in gem_overlapeed_met_names:
        id = name_to_id[name]
        mam = id_to_mam[id]
        met = mam_to_met[mam]
        file.write(name + '\t' + id + '\t' + mam + '\t' + met + '\n')
    file.flush()
    file.close()
    
    df = pd.read_csv(args.out_dir + '/gem_overlapped_change_id.tsv', sep='\t')
    df[['sample_id', 'sample_group']] = df['key'].str.split(":", expand=True)
    df = df.set_index(['sample_id', 'sample_group'])
    df = df.drop(columns='key')
    df.to_csv(args.out_dir + '/metabolite.tsv', sep='\t', index=True)
    
    sys.stdout = orig_stdout
    log_file.close()
        
if __name__ == "__main__":
    main(parse_args())
