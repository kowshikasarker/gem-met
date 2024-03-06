import argparse
import pandas as pd
import sys
import pathlib
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocessing of basleine and end metabolomics profiles.\nBoth profiles are expected to be .tsv files with two columns \person_id\' and \'diet\' indicating sample identifier and assigned diet.\n Other columns are expected to be concentrations of metabolic compounds.\n\n Only rows and coulmns common in both prifles are kept and columns are also filtered based on missing value proportion.\nIn the output directory, 10 .tsv files are generated.')
    
    parser.add_argument("-b", "--base_path", type=str,
                        help="path to baseline metabolomics profile in tsv format",
                        required=True, default=None)
    parser.add_argument("-e", "--end_path", type=str,
                        help="path to end metabolomics profile in tsv format",
                        required=True, default=None)
    parser.add_argument("-g", "--gem_met_path", type=str,
                        help="path to human gem metabolites list in tsv format", required=True, default=None)
    parser.add_argument("-x", "--human_gem_path", type=str,
                        help="path to human gem model in .xlsx format", required=True, default=None)
    parser.add_argument("-i", "--id_path", type=str,
                        help="path to metabolite identifiers in tsv format", required=True, default=None)
    parser.add_argument("-l", "--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    parser.add_argument("-m", "--missing_pct", type=float,
                        help="threshold of missing value percentage for dropping columns", required=True, default=None)
    
    args = parser.parse_args()
    return args


# Keep features that have proportion of missing values <  p for any food
def drop_missing_values(base_df, end_df, missing_pct):
    print('drop_missing_values')
    diets = set(base_df['diet']).union(set(end_df['diet']))
    diets = [diet for diet in diets if not diet.startswith('No')]
    print(len(diets), ' diets', diets)
    
    columns_to_keep_baseline = set()
    columns_to_keep_end = set()

    for diet in diets:
        # Select dataset for this study
        base_df_diet = base_df[base_df['diet'].str.contains(diet)]
        end_df_diet = end_df[end_df['diet'].str.contains(diet)]
        
        # Compute percent of missing values for the datasets
        missing_pct_base = base_df_diet.isnull().sum() / len(base_df_diet)
        missing_pct_end = end_df_diet.isnull().sum() / len(end_df_diet)

        # Keep all features that have < p percent missing
        # i.e. have > p percent features present
        missing_pct_base = missing_pct_base < missing_pct
        missing_pct_end = missing_pct_end < missing_pct

        # Subset feature list to only include those features
        missing_pct_base = missing_pct_base.where(lambda a: a).dropna().index
        missing_pct_end = missing_pct_end.where(lambda a: a).dropna().index
        
        columns_to_keep_baseline.update(list(missing_pct_base))
        columns_to_keep_end.update(list(missing_pct_end))
                                   
    # Subset columns
    base_df = base_df[list(columns_to_keep_baseline)]
    end_df = end_df[list(columns_to_keep_end)]
   
    return base_df, end_df

# Imputes missing values to uniform random values between [0, mm * minimum observed] for every feature
def impute_missing_values(base_df, end_df, coeff):

    # Compute per-feature minimums for dataset
    base_df_feature_mins = np.min(base_df, axis=0)
    base_df_nan_dict = {}
    end_df_feature_mins = np.min(end_df, axis=0)
    end_df_nan_dict = {}

    # Create new datasets that contains random values for each subject for each feature,
    # between 0 and mm * the minimum for that feature
    for feature, minimum in base_df_feature_mins.iteritems():
        base_df_nan_dict[feature] = np.random.uniform(
            low=0, high=coeff*minimum, size=len(base_df)
        )

    for feature, minimum in end_df_feature_mins.iteritems():
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
    print('id_path', args.id_path)
    print('log_path', args.log_path)
    print('out_dir', args.out_dir)
    
    base_df = pd.read_csv(args.base_path, sep='\t')
    print(base_df.shape)
    end_df = pd.read_csv(args.end_path, sep='\t')
    print(end_df.shape)
    
    base_df.columns = map(str.lower, base_df.columns)
    end_df.columns = map(str.lower, end_df.columns)
    
    base_df['key'] = base_df['person_id'] + '.' + base_df['diet']
    end_df['key'] = end_df['person_id'] + '.' + end_df['diet']
    
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
    
    print(base_df.shape)
    print(end_df.shape)
    
    base_df, end_df = drop_missing_values(base_df, end_df, args.missing_pct)
    
    base_df = base_df.drop(['person_id', 'diet'], axis=1)
    end_df = end_df.drop(['person_id', 'diet'], axis=1)
   
    
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
    print()
    
    change_df = end_df.sub(base_df)
    print('change_df', change_df.shape)
    change_df.to_csv(args.out_dir + '/preprocessed_change_name.tsv', sep='\t')
    
    
    id_df = pd.read_csv(args.id_path, sep='\t',
                   usecols=['Name', 'HMDB_ID'])
    id_df['Name'] = id_df['Name'].apply(lambda x: x.lower())
    id_df = id_df[id_df.Name.isin(common_cols)]
    id_df = id_df.dropna(subset=['HMDB_ID'])
    
    print('id_df', id_df.shape)
        
    name_to_hmdb = dict(zip(id_df.Name, id_df.HMDB_ID))
    hmdb_to_name = dict(zip(id_df.HMDB_ID, id_df.Name))
    
    common_cols_hmdb = set([name_to_hmdb[col] for col in common_cols if col in name_to_hmdb.keys()])
    print(len(common_cols_hmdb))
    
    #gem_met = pd.read_csv('/home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Human-GEM/model/metabolites.tsv', sep='\t')
    gem_met = pd.read_csv(args.gem_met_path, sep='\t')
    gem_met = gem_met.dropna(subset=['metHMDBID'])
    print(gem_met.shape)
    
    common_mets = common_cols_hmdb.intersection(set(gem_met['metHMDBID']))
    print('common_mets', len(common_mets))
    
    common_mets_names = [hmdb_to_name[met] for met in common_mets]
    print(len(common_mets_names))
    print(common_mets_names)
    
    base_df = base_df[common_mets_names]
    end_df = end_df[common_mets_names]
    print(base_df.shape, end_df.shape)
    
    base_df.to_csv(args.out_dir + '/gem_overlapped_base_name.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/gem_overlapped_end_name.tsv', sep='\t')
    
    base_df = base_df.rename(columns=name_to_hmdb)
    end_df = end_df.rename(columns=name_to_hmdb)
    
    base_df.to_csv(args.out_dir + '/gem_overlapped_base_hmdb.tsv', sep='\t')
    end_df.to_csv(args.out_dir + '/gem_overlapped_end_hmdb.tsv', sep='\t')
     
    base_df = base_df.sort_index(axis=0).sort_index(axis=1)
    end_df = end_df.sort_index(axis=0).sort_index(axis=1)
    
    change_df = end_df.sub(base_df)
    
    change_df.to_csv(args.out_dir + '/gem_overlapped_change_hmdb.tsv', sep='\t')
    change_df = change_df.rename(columns=hmdb_to_name)
    change_df.to_csv(args.out_dir + '/gem_overlapped_change_name.tsv', sep='\t')
    print('change_df', change_df.shape)
    
    gem_overlapeed_met_names = list(change_df.columns)
    #print('gem_overlapeed_met_names', gem_overlapeed_met_names)
    gem_met = gem_met[gem_met['metHMDBID'].isin(hmdb_to_name.keys())]
    hmdb_to_mam = dict(zip(gem_met.metHMDBID, gem_met.metsNoComp))
    #print(len(hmdb_to_mam))
    #print(hmdb_to_mam)
    
    gem_model = pd.read_excel(args.human_gem_path, sheet_name='METS', usecols=['NAME', 'REPLACEMENT ID'])
    print('human_gem_path', args.human_gem_path)
    #gem_model = pd.read_excel('/home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Human-GEM/model/Human-GEM.xlsx', sheet_name='METS', usecols=['NAME', 'REPLACEMENT ID'])
    pattern = '|'.join(['e', 'x', 'm', 'c', 'l', 'r', 'g', 'n', 'i'])
    gem_model['MAM_ID'] = gem_model['REPLACEMENT ID'].str.replace(pattern, '')
    gem_model = gem_model[gem_model['MAM_ID'].isin(hmdb_to_mam.values())]
    #print('gem_model', gem_model.shape)
    mam_to_met = dict(zip(gem_model['MAM_ID'], gem_model['NAME']))
    #print(len(mam_to_met))
    #print(mam_to_met)
    
    file = open(args.out_dir + '/gem_overlapped_metabolites.tsv', 'w')
    file.write('Name' + '\t' + 'HMDB_ID' + '\t' + 'MAM_ID' + '\t' + 'MET_ID' + '\n')
    for name in gem_overlapeed_met_names:
        hmdb = name_to_hmdb[name]
        mam = hmdb_to_mam[hmdb]
        met = mam_to_met[mam]
        file.write(name + '\t' + hmdb + '\t' + mam + '\t' + met + '\n')
    file.flush()
    file.close()
    
    sys.stdout = orig_stdout
    log_file.close()
        
if __name__ == "__main__":
    main(parse_args())
