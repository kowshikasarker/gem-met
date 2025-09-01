import argparse
import sys
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Discards and prerpocesses GEM reactions based on reaction reversibility and metabolite overlap filters.')
    
    parser.add_argument("--gem_path", type=str,
                        help="path to a gem model in .xlsx format",
                        required=True, default=None)
    parser.add_argument("--valid_met_path", type=str,
                        help="path to list of valid (common between user data and gem) metabolites in tsv format",
                        required=True, default=None)
    
    parser.add_argument("--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    args = parser.parse_args()
    return args

def extract_metabolites(eqn):
    print('eqn', eqn)
    old = eqn.split(' + ')
    print('old', old)
    new = set()
    flag = False
    for i_old in old:
        print(i_old)
        i_old_split = i_old.split(' ')
        print(i_old)
        if(len(i_old_split)>1):
            print('if')
            if(i_old_split[0].replace(".", "").isnumeric()): # If a metabolite   a coefficient in reaction equation
                i_new = i_old[len(i_old_split[0])+1:]
                print(i_new, i_old)
                
                new.add(i_new)
                flag = True
            else:
                new.add(i_old)
        else:
            
            new.add(i_old)
    print('new', new)
    return new
            

def all_reactions(valid_met_path, gem_path, out_dir):
    valid_mets = pd.read_csv(valid_met_path, sep='\t')
    valid_mets = set(valid_mets['MET_ID'])
    react_gem = pd.read_excel(gem_path, sheet_name='RXNS', usecols=['ID', 'EQUATION', 'SUBSYSTEM'])
    react_gem = react_gem.rename(columns={'ID': 'RXN_ID'})
    
    print(react_gem.head())
    pattern = '|'.join(['\[e\]', '\[x\]', '\[m\]', '\[c\]', '\[l\]', '\[r\]', '\[g\]', '\[n\]', '\[i\]'])
    react_gem['EQUATION'] = react_gem['EQUATION'].str.replace(pattern, '', regex=True)
    print(react_gem.head())
    react_gem['Direction'] = react_gem['EQUATION'].apply(lambda x: 2 if '<=>' in x else 1)
    react_gem['EQUATION'] = react_gem['EQUATION'].str.replace('<=>', '=>')
    react_gem[['EQUATION_LHS', 'EQUATION_RHS']] = react_gem['EQUATION'].str.split(' => ', expand=True)
    react_gem['Substrate_Set'] = react_gem['EQUATION_LHS'].apply(lambda x: extract_metabolites(x))
    react_gem['Product_Set'] = react_gem['EQUATION_RHS'].apply(lambda x: extract_metabolites(x))
    react_gem['Metabolite_Set'] = react_gem.apply(lambda x: x['Product_Set'].union(x['Substrate_Set']), axis=1)
    
    print(react_gem)

    react_gem['Substrate_Count'] = react_gem['Substrate_Set'].apply(lambda x: len(x))
    react_gem['Product_Count'] = react_gem['Product_Set'].apply(lambda x: len(x))
    react_gem['Metabolite_Count'] = react_gem['Metabolite_Set'].apply(lambda x: len(x))
    
    react_gem['Measured_Substrate'] = react_gem['Substrate_Set'].apply(lambda x: x.intersection(valid_mets))
    react_gem['Measured_Product'] = react_gem['Product_Set'].apply(lambda x: x.intersection(valid_mets))
    react_gem['Measured_Metabolite'] = react_gem['Metabolite_Set'].apply(lambda x: x.intersection(valid_mets))

    react_gem['Measured_Substrate_Count'] = react_gem['Measured_Substrate'].apply(lambda x: len(x))
    react_gem['Measured_Product_Count'] = react_gem['Measured_Product'].apply(lambda x: len(x))
    react_gem['Measured_Metabolite_Count'] = react_gem['Measured_Metabolite'].apply(lambda x: len(x))
    #react_gem.head()
    react_gem.to_csv(out_dir + '/all-reactions.tsv', sep='\t', index=False)
    return react_gem

def react_set_1(react_gem, out_dir):
    react_set_1_dir = out_dir + '/reaction-set-1'
    Path(react_set_1_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_gem = react_gem[react_gem['Direction'] == 1]
    one_more_mets_in = react_gem[react_gem['Measured_Metabolite_Count'] > 0]
    print(one_more_mets_in.shape)
    one_more_mets_in.to_csv(react_set_1_dir + '/reaction-set-1.tsv', sep='\t', index=False)
    print(one_more_mets_in.shape)
    
def react_set_2(react_gem, out_dir):
    react_set_2_dir = out_dir + '/reaction-set-2'
    Path(react_set_2_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_gem = react_gem[react_gem['Direction'] == 1]
    one_more_mets_in = react_gem[(react_gem['Measured_Substrate_Count'] > 0) & (react_gem['Measured_Product_Count'] > 0)]
    print(one_more_mets_in.shape)
    one_more_mets_in.to_csv(react_set_2_dir + '/reaction-set-2.tsv', sep='\t', index=False)
    print(one_more_mets_in.shape)
    
def reverse_equation(eqn):
    eqn_split = eqn.split(' => ')
    if(len(eqn_split) != 2):
        print('Something wrong!')
    new_eqn = eqn_split[1] + ' => ' + eqn_split[0]
    #print(eqn, eqn_split, new_eqn)
    return new_eqn

def react_set_3(react_gem, out_dir):
    react_set_3_dir = out_dir + '/reaction-set-3'
    Path(react_set_3_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_df = react_gem[react_gem['Measured_Metabolite_Count'] > 0]
    print(react_df.shape)
    react_df.to_csv(react_set_3_dir + '/reaction-set-3-basic.tsv', sep='\t', index=False)
    print(react_df.shape)
    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]
    react3_df = react2_df.copy()
    react2_df['RXN_ID'] = react2_df['RXN_ID'] + 'F'
    react2_df['EQUATION'] = react2_df['EQUATION'].str.replace('<=>','=>')
    print(react2_df.head())




    react3_df['EQUATION'] = react3_df['EQUATION'].apply(lambda x: reverse_equation(x))
    react3_df['RXN_ID'] = react3_df['RXN_ID'] + 'B'
    #swapping columns
    react3_df[['EQUATION_LHS', 'EQUATION_RHS']] = react3_df[['EQUATION_RHS', 'EQUATION_LHS']]
    react3_df[['Substrate_Set', 'Product_Set']] = react3_df[['Product_Set', 'Substrate_Set']]
    react3_df[['Substrate_Count', 'Product_Count']] = react3_df[['Product_Count', 'Substrate_Count']]
    react3_df[['Measured_Substrate', 'Measured_Product']] = react3_df[['Measured_Product', 'Measured_Substrate']]
    react3_df[['Measured_Substrate_Count', 'Measured_Product_Count']] = react3_df[['Measured_Product_Count', 'Measured_Substrate_Count']]

    print(react3_df.head())

    react_df = pd.concat([react1_df, react2_df, react3_df], axis=0)
    react_df.to_csv(react_set_3_dir + '/reaction-set-3.tsv', sep='\t', index=False)
    

def react_set_4(react_gem, out_dir):
    react_set_4_dir = out_dir + '/reaction-set-4'
    Path(react_set_4_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_df = react_gem[(react_gem['Measured_Substrate_Count'] > 0) & (react_gem['Measured_Product_Count'] > 0)]
    print(react_df.shape)
    react_df.to_csv(react_set_4_dir + '/reaction-set-4-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]
    react3_df = react2_df.copy()

    react2_df['RXN_ID'] = react2_df['RXN_ID'] + 'F'
    react2_df['EQUATION'] = react2_df['EQUATION'].str.replace('<=>','=>')
    print(react2_df.head())

    react3_df['EQUATION'] = react3_df['EQUATION'].apply(lambda x: reverse_equation(x))
    react3_df['RXN_ID'] = react3_df['RXN_ID'] + 'B'
    #swapping columns
    react3_df[['EQUATION_LHS', 'EQUATION_RHS']] = react3_df[['EQUATION_RHS', 'EQUATION_LHS']]
    react3_df[['Substrate_Set', 'Product_Set']] = react3_df[['Product_Set', 'Substrate_Set']]
    react3_df[['Substrate_Count', 'Product_Count']] = react3_df[['Product_Count', 'Substrate_Count']]
    react3_df[['Measured_Substrate', 'Measured_Product']] = react3_df[['Measured_Product', 'Measured_Substrate']]
    react3_df[['Measured_Substrate_Count', 'Measured_Product_Count']] = react3_df[['Measured_Product_Count', 'Measured_Substrate_Count']]

    print(react3_df.head())

    react_df = pd.concat([react1_df, react2_df, react3_df], axis=0)
    react_df.to_csv(react_set_4_dir + '/reaction-set-4.tsv', sep='\t', index=False)
    
def react_set_5(react_gem, out_dir):
    react_set_5_dir = out_dir + '/reaction-set-5'
    Path(react_set_5_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_df = react_gem[react_gem['Measured_Metabolite_Count'] > 0]
    print(react_df.shape)
    react_df.to_csv(react_set_5_dir + '/reaction-set-5-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]

    react2_df['Substrate_Set'] = react2_df['Metabolite_Set']
    react2_df['Product_Set'] = react2_df['Metabolite_Set']

    react2_df['Measured_Substrate'] = react2_df['Measured_Metabolite']
    react2_df['Measured_Product'] = react2_df['Measured_Metabolite']

    react_df = pd.concat([react1_df, react2_df], axis=0)

    react_df.to_csv(react_set_5_dir + '/reaction-set-5.tsv', sep='\t', index=False)
    
def react_set_6(react_gem, out_dir):
    react_set_6_dir = out_dir + '/reaction-set-6'
    Path(react_set_6_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_df = react_gem[(react_gem['Measured_Substrate_Count'] > 0) & (react_gem['Measured_Product_Count'] > 0)]
    print(react_df.shape)
    react_df.to_csv(react_set_6_dir + '/reaction-set-6-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]

    react2_df['Substrate_Set'] = react2_df['Metabolite_Set']
    react2_df['Product_Set'] = react2_df['Metabolite_Set']

    react2_df['Measured_Substrate'] = react2_df['Measured_Metabolite']
    react2_df['Measured_Product'] = react2_df['Measured_Metabolite']

    react_df = pd.concat([react1_df, react2_df], axis=0)
    react_df.to_csv(react_set_6_dir + '/reaction-set-6.tsv', sep='\t', index=False)

def react_set_7(react_gem, out_dir):
    react_set_7_dir = out_dir + '/reaction-set-7'
    Path(react_set_7_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_gem = react_gem[react_gem['Direction'] == 1]
    two_more_mets_in = react_gem[react_gem['Measured_Metabolite_Count'] > 1]
    print(two_more_mets_in.shape)
    two_more_mets_in.to_csv(react_set_7_dir + '/reaction-set-7.tsv', sep='\t', index=False)
    print(two_more_mets_in.shape)
    
def react_set_8(react_gem, out_dir):
    react_set_8_dir = out_dir + '/reaction-set-8'
    Path(react_set_8_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_df = react_gem[react_gem['Measured_Metabolite_Count'] > 1]
    print(react_df.shape)
    react_df.to_csv(react_set_8_dir + '/reaction-set-8-basic.tsv', sep='\t', index=False)
    print(react_df.shape)
    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]
    react3_df = react2_df.copy()
    react2_df['RXN_ID'] = react2_df['RXN_ID'] + 'F'
    react2_df['EQUATION'] = react2_df['EQUATION'].str.replace('<=>','=>')
    print(react2_df.head())




    react3_df['EQUATION'] = react3_df['EQUATION'].apply(lambda x: reverse_equation(x))
    react3_df['RXN_ID'] = react3_df['RXN_ID'] + 'B'
    #swapping columns
    react3_df[['EQUATION_LHS', 'EQUATION_RHS']] = react3_df[['EQUATION_RHS', 'EQUATION_LHS']]
    react3_df[['Substrate_Set', 'Product_Set']] = react3_df[['Product_Set', 'Substrate_Set']]
    react3_df[['Substrate_Count', 'Product_Count']] = react3_df[['Product_Count', 'Substrate_Count']]
    react3_df[['Measured_Substrate', 'Measured_Product']] = react3_df[['Measured_Product', 'Measured_Substrate']]
    react3_df[['Measured_Substrate_Count', 'Measured_Product_Count']] = react3_df[['Measured_Product_Count', 'Measured_Substrate_Count']]

    print(react3_df.head())

    react_df = pd.concat([react1_df, react2_df, react3_df], axis=0)
    react_df.to_csv(react_set_8_dir + '/reaction-set-8.tsv', sep='\t', index=False)
    
def react_set_9(react_gem, out_dir):
    react_set_9_dir = out_dir + '/reaction-set-9'
    Path(react_set_9_dir).mkdir(parents=True, exist_ok=True)
    
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    react_df = react_gem[react_gem['Measured_Metabolite_Count'] > 1]
    print(react_df.shape)
    react_df.to_csv(react_set_9_dir + '/reaction-set-9-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]

    react2_df['Substrate_Set'] = react2_df['Metabolite_Set']
    react2_df['Product_Set'] = react2_df['Metabolite_Set']

    react2_df['Measured_Substrate'] = react2_df['Measured_Metabolite']
    react2_df['Measured_Product'] = react2_df['Measured_Metabolite']

    react_df = pd.concat([react1_df, react2_df], axis=0)

    react_df.to_csv(react_set_9_dir + '/reaction-set-9.tsv', sep='\t', index=False)
    
def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    
    print('gem_path', args.gem_path)
    print('valid_met_path', args.valid_met_path)
    print('log_path', args.log_path)
    print('out_dir', args.out_dir)
    
    react_gem = all_reactions(args.valid_met_path, args.gem_path, args.out_dir)
    
    react_set_1(react_gem, args.out_dir)
    react_set_2(react_gem, args.out_dir)
    react_set_3(react_gem, args.out_dir)
    react_set_4(react_gem, args.out_dir)
    react_set_5(react_gem, args.out_dir)
    react_set_6(react_gem, args.out_dir)
    react_set_7(react_gem, args.out_dir)
    react_set_8(react_gem, args.out_dir)
    react_set_9(react_gem, args.out_dir)
    
    sys.stdout = orig_stdout
    log_file.close()
    
if __name__ == "__main__":
    main(parse_args())
