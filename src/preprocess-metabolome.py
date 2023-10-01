import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--base_path", type=str,
                        help="path to baseline metabolomics profile in csv format",
                        required=True, default=None)
    parser.add_argument("-e", "--end_path", type=str,
                        help="path to end metabolomics profile in csv format",
                        required=True, default=None)
    parser.add_argument("-m", "--meta_path", type=str,
                        help="path to metadata for metabolomics profiles in csv format", required=True, default=None)
    parser.add_argument("-i", "--id_path", type=str,
                        help="path to metabolite identifiers in tsv format", required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    args = parser.parse_args()
    return args

def split_grains(df, barley_key, oats_key):
    no_barley = df[df.index.str.contains('NoGrains')]
    print('no_barley', no_barley.shape)
    #print(no_barley)
    
    no_barley.index = no_barley.index.str.replace('Grains', 'Barley')
    print('no_barley', no_barley.shape)
    
    no_oats = no_barley.copy()
    no_oats.index = no_oats.index.str.replace('Barley', 'Oats')
    print('no_oats', no_oats.shape)
    
    barley = df[df.index.isin(barley_key)]
    barley.index = barley.index.str.replace('Grains', 'Barley')
    print('barley', barley.shape)
    
    oats = df[df.index.isin(oats_key)]
    oats.index = oats.index.str.replace('Grains', 'Oats')
    print('oats', oats.shape)
    
    df = df[~df.index.str.contains('Grains')]
    print('df', df.shape)
    
    df = pd.concat([df, barley, no_barley, oats, no_oats], axis=0)
    print('df', df.shape)
    return df

# Keep features that have proportion of missing values <  p for any food
def drop_missing_values(
    met_baseline, met_end, p=0.2
):

    columns_to_keep_baseline = set()
    columns_to_keep_end = set()

    foods = ['Almond', 'Walnut', 'Avocado', 'Barley', 'Oats', 'Broccoli']

    for study in foods:
        # Select dataset for this study
        met_baseline_study = met_baseline.filter(like=study, axis=0)
        met_end_study = met_end.filter(like=study, axis=0)

        # Compute percent of missing values for the datasets
        p_baseline = met_baseline_study.isnull().sum() / len(met_baseline_study)
        p_end = met_end_study.isnull().sum() / len(met_end_study)

        # Keep all features that have < p percent missing
        # i.e. have > p percent features present
        p_baseline = p_baseline < p
        p_end = p_end < p
        
        # Subset feature list to only include those features
        p_baseline = p_baseline.where(lambda a: a).dropna().index
        p_end = p_end.where(lambda a: a).dropna().index

        # Add column to keep list
        columns_to_keep_baseline.update(list(p_baseline))
        columns_to_keep_end.update(list(p_end))

    # Subset columns
    met_baseline = met_baseline[list(columns_to_keep_baseline)]
    met_end = met_end[list(columns_to_keep_end)]

    return met_baseline, met_end

# Imputes missing values to uniform random values between [0, mm * minimum observed] for every feature
def impute_missing_values(met_baseline, met_end, mm=0.25, verbose=False):

    # Compute per-feature minimums for dataset
    met_baseline_feature_mins = np.min(met_baseline, axis=0)
    met_baseline_nan_dict = {}
    met_end_feature_mins = np.min(met_end, axis=0)
    met_end_nan_dict = {}

    # Create new datasets that contains random values for each subject for each feature,
    # between 0 and mm * the minimum for that feature
    for feature, minimum in met_baseline_feature_mins.iteritems():
        met_baseline_nan_dict[feature] = np.random.uniform(
            low=0, high=mm * minimum, size=len(met_baseline)
        )

    for feature, minimum in met_end_feature_mins.iteritems():
        met_end_nan_dict[feature] = np.random.uniform(
            low=0, high=mm * minimum, size=len(met_end)
        )

    # Update original dataset with new values for any missing entries
    # Original values should be preserved
    met_baseline_nan = pd.DataFrame(met_baseline_nan_dict)
    met_baseline_nan.index = met_baseline.index

    met_end_nan = pd.DataFrame(met_end_nan_dict)
    met_end_nan.index = met_end.index

    met_baseline.update(met_baseline_nan, overwrite=False)
    met_end.update(met_end_nan, overwrite=False)

    return met_baseline, met_end

def main(args):
    base_df = pd.read_csv(args.base_path, sep=',', index_col='Key')
    print(base_df.shape)
    end_df = pd.read_csv(args.end_path, sep=',', index_col='Key')
    print(end_df.shape)
    
    base_df.columns = map(str.lower, base_df.columns)
    end_df.columns = map(str.lower, end_df.columns)
    
    base_df.replace(0, np.nan, inplace=True)

    base_df['3-methyl-2-oxobutanoic acid'] = base_df['2-oxoisovaleric  acid'].fillna(0) + base_df['3-methyl-2-oxobutanoic acid'].fillna(0)# + base_df['butanoic acid, 3-methyl-2-oxo-']
    base_df['indole-3-acetic acid'] = base_df['indole-3-acetic acid'].fillna(0) + base_df['3-indole lactic acid'].fillna(0)
    base_df['c16:1 (9)'] = base_df['c16:1'].fillna(0) + base_df['c16:1 (9)'].fillna(0)
    base_df = base_df.drop(columns=['2-oxoisovaleric  acid', '3-indole lactic acid', 'c16:1'])
    base_df.replace(0, np.nan, inplace=True)
    print(base_df.shape)
    
    
    end_df.replace(0, np.nan, inplace=True)
    end_df['3-methyl-2-oxobutanoic acid'] = end_df['2-oxoisovaleric  acid'].fillna(0) + end_df['3-methyl-2-oxobutanoic acid'].fillna(0) + end_df['butanoic acid, 3-methyl-2-oxo-'].fillna(0)
    end_df['indole-3-acetic acid'] = end_df['indole-3-acetic acid'].fillna(0) + end_df['3-indole lactic acid'].fillna(0)
    end_df['c16:1 (9)'] = end_df['c16:1'].fillna(0) + end_df['c16:1 (9)'].fillna(0)
    end_df = end_df.drop(columns=['2-oxoisovaleric  acid', 'butanoic acid, 3-methyl-2-oxo-', '3-indole lactic acid', 'c16:1'])
    end_df.replace(0, np.nan, inplace=True)
    print(end_df.shape)
    
    
    metadata = pd.read_csv(args.meta_path, sep=',')
    barley_key = metadata[metadata.Treatment2 == 'Barley']['Key']
    oats_key = metadata[metadata.Treatment2 == 'Oats']['Key']
    
    base_df = split_grains(base_df, barley_key, oats_key)
    end_df = split_grains(end_df, barley_key, oats_key)
    
    base_df.index = base_df.index.str.replace('.P', ':P')
    base_df['person_id'] = base_df.index.str.replace('.Baseline', '')

    end_df.index = end_df.index.str.replace('.P', ':P')
    end_df['person_id'] = end_df.index.str.replace('.End', '')
    
    common_rows = set(base_df.person_id).intersection(set(end_df.person_id))
    print('common_rows', len(common_rows))
    base_df = base_df[base_df.person_id.isin(common_rows)]
    end_df = end_df[end_df.person_id.isin(common_rows)]
    base_df = base_df.drop(columns=['person_id'])
    end_df = end_df.drop(columns=['person_id'])
    
    common_cols = set(base_df.columns).intersection(set(end_df.columns))
    print('common_cols', len(common_cols))
    base_df = base_df[common_cols]
    end_df = end_df[common_cols]
    
    print(base_df.shape)
    print(end_df.shape)
    
    #base_df.to_csv(args.out_dir + '/base_name_109.tsv', sep='\t')
    #end_df.to_csv(args.out_dir + '/end_name_109.tsv', sep='\t')
    
    base_df, end_df = drop_missing_values(base_df, end_df)
    print(base_df.shape, end_df.shape)
    base_df, end_df = impute_missing_values(base_df, end_df)
    print(base_df.shape, end_df.shape)
    
    
    common_cols = set(base_df.columns).intersection(set(end_df.columns))
    print('common_cols', len(common_cols))
    base_df = base_df[common_cols]
    end_df = end_df[common_cols]
    
    id_df = pd.read_csv(args.id_path, sep='\t',
                   usecols=['Name', 'Name_Lower', 'HMDB_ID'])
    id_df = id_df[id_df.Name_Lower.isin(common_cols)]
    id_df = id_df.dropna(subset=['HMDB_ID'])
    
    name_to_hmdb = dict(zip(id_df.Name_Lower, id_df.HMDB_ID))
    hmdb_to_name = dict(zip(id_df.HMDB_ID, id_df.Name_Lower))
    
    common_cols_hmdb = set([name_to_hmdb[col] for col in common_cols if col in name_to_hmdb.keys()])
    print(len(common_cols_hmdb))
    
    gem_met = pd.read_csv('/home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Human-GEM/model/metabolites.tsv', sep='\t')
    gem_met = gem_met.dropna(subset=['metHMDBID'])
    print(gem_met.shape)
    
    common_mets = common_cols_hmdb.intersection(set(gem_met['metHMDBID']))
    print(len(common_mets))
    
    common_mets_names = [hmdb_to_name[met] for met in common_mets]
    print(len(common_mets_names))
    print(common_mets_names)
    
    base_df = base_df[common_mets_names]
    end_df = end_df[common_mets_names]
    print(base_df.shape, end_df.shape)
    
    base_df.to_csv(args.out_dir + '/base_name_75.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/end_name_75.tsv', sep='\t')
    
    base_df = base_df.rename(columns=name_to_hmdb)
    end_df = end_df.rename(columns=name_to_hmdb)
    
    base_df.to_csv(args.out_dir + '/base_hmdb_75.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/end_hmdb_75.tsv', sep='\t')
    
    base_df.index = base_df.index.str.replace('.Baseline', '')
    end_df.index = end_df.index.str.replace('.End', '')
    
    base_df = base_df.sort_index(axis=0).sort_index(axis=1)
    end_df = end_df.sort_index(axis=0).sort_index(axis=1)
    
    change_df = end_df.sub(base_df)
    
    change_df.to_csv(args.out_dir + '/change_hmdb_75.tsv', sep='\t')
    
if __name__ == "__main__":
    main(parse_args())