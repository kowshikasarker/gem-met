import argparse
import pandas as pd
from pathlib import Path
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocesses Human-GEM')
    parser.add_argument("-x", "--gem_path", type=str,
                        help="path to human gem model in .xlsx format",
                        required=True, default=None)
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

def extract_metabolites(eqn):
    old = eqn.split(' + ')
    new = set()
    flag = False
    for i_old in old:
        i_old_split = i_old.split(' ')
        if(len(i_old_split)>1):
            if(i_old_split[0].replace(".", "").isnumeric()): # If a metabolite has a coefficient in reaction equation
                i_new = i_old[len(i_old_split[0])+1:]
                #print(i_new, i_old)
                #i_new = i_old.replace(i_old_split[0] + ' ', '')
                #print(i_old, '->', i_new)
                new.add(i_new)
                flag = True
            else:
                new.add(i_old)
        else:
            new.add(i_old)
    return new
            

def prepare_all_reactions(met_path, gem_path, out_dir):
    gem_overlapped_mets = pd.read_csv(met_path, sep='\t')
    gem_overlapped_mets = set(gem_overlapped_mets['MET_ID'])
    react_gem = pd.read_excel(gem_path, sheet_name='RXNS', usecols=['ID', 'EQUATION', 'SUBSYSTEM'])

    react_gem = react_gem.rename(columns={'ID': 'RXN_ID'})
    pattern = '|'.join(['\[e\]', '\[x\]', '\[m\]', '\[c\]', '\[l\]', '\[r\]', '\[g\]', '\[n\]', '\[i\]'])
    react_gem['EQUATION'] = react_gem['EQUATION'].str.replace(pattern, '')
    react_gem['Direction'] = react_gem['EQUATION'].apply(lambda x: 2 if '<=>' in x else 1)
    react_gem['EQUATION'] = react_gem['EQUATION'].str.replace('<=>', '=>')
    react_gem[['EQUATION_LHS', 'EQUATION_RHS']] = react_gem['EQUATION'].str.split(' => ', expand=True)
    react_gem['Substrate_Set'] = react_gem['EQUATION_LHS'].apply(lambda x: extract_metabolites(x))
    react_gem['Product_Set'] = react_gem['EQUATION_RHS'].apply(lambda x: extract_metabolites(x))
    react_gem['Metabolite_Set'] = react_gem.apply(lambda x: x['Product_Set'].union(x['Substrate_Set']), axis=1)

    react_gem['Substrate_Count'] = react_gem['Substrate_Set'].apply(lambda x: len(x))
    react_gem['Product_Count'] = react_gem['Product_Set'].apply(lambda x: len(x))
    react_gem['Metabolite_Count'] = react_gem['Metabolite_Set'].apply(lambda x: len(x))
    react_gem['Measured_Substrate'] = react_gem['Substrate_Set'].apply(lambda x: x.intersection(gem_overlapped_mets))
    react_gem['Measured_Product'] = react_gem['Product_Set'].apply(lambda x: x.intersection(gem_overlapped_mets))
    react_gem['Measured_Metabolite'] = react_gem['Metabolite_Set'].apply(lambda x: x.intersection(gem_overlapped_mets))

    react_gem['Measured_Substrate_Count'] = react_gem['Measured_Substrate'].apply(lambda x: len(x))
    react_gem['Measured_Product_Count'] = react_gem['Measured_Product'].apply(lambda x: len(x))
    react_gem['Measured_Metabolite_Count'] = react_gem['Measured_Metabolite'].apply(lambda x: len(x))
    #react_gem.head()
    react_gem.to_csv(out_dir + '/all-reactions.tsv', sep='\t', index=False)
    react_gem = react_gem[~react_gem['SUBSYSTEM'].isin(['Transport reactions', 'Exchange/demand reactions'])]
    return react_gem

def preapre_reaction_set_1(react_gem, react_set_1_dir):
    print('preapre_reaction_set_1')
    react_gem = react_gem[react_gem['Direction'] == 1]
    one_more_mets_in = react_gem[react_gem['Measured_Metabolite_Count'] > 0]
    print(one_more_mets_in.shape)
    one_more_mets_in.to_csv(react_set_1_dir + '/reaction-set-1.tsv', sep='\t', index=False)
    print(one_more_mets_in.shape)
    
def preapre_reaction_set_2(react_gem, react_set_2_dir):
    print('preapre_reaction_set_2')
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
    return new_eqn

def preapre_reaction_set_3(react_gem, react_set_3_dir):
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
    

def preapre_reaction_set_4(react_gem, react_set_4_dir):
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
    
def str_to_set(cell):
    if isinstance(cell, set):
        return cell
    cell = ''.join(c for c in cell if c not in "'{}")
    print(cell)
    cell = set(cell.split(', '))
    print(cell)
    return cell
    
def preapre_reaction_set_5(react_gem, react_set_5_dir):
    react_df = react_gem[react_gem['Measured_Metabolite_Count'] > 0]
    print(react_df.shape)
    react_df.to_csv(react_set_5_dir + '/reaction-set-5-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react_df['Substrate_Set'] = react_df['Substrate_Set'].apply(lambda x: str_to_set(x))
    react_df['Product_Set'] = react_df['Product_Set'].apply(lambda x: str_to_set(x))
    react_df['Measured_Substrate'] = react_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_df['Measured_Product'] = react_df['Measured_Product'].apply(lambda x: str_to_set(x))

    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]

    react2_df['Substrate_Set'] = react2_df['Metabolite_Set']
    react2_df['Product_Set'] = react2_df['Metabolite_Set']

    react2_df['Measured_Substrate'] = react2_df['Measured_Metabolite']
    react2_df['Measured_Product'] = react2_df['Measured_Metabolite']

    react_df = pd.concat([react1_df, react2_df], axis=0)

    react_df.to_csv(react_set_5_dir + '/reaction-set-5.tsv', sep='\t', index=False)
    
def preapre_reaction_set_6(react_gem, react_set_6_dir):
    react_df = react_gem[(react_gem['Measured_Substrate_Count'] > 0) & (react_gem['Measured_Product_Count'] > 0)]
    print(react_df.shape)
    react_df.to_csv(react_set_6_dir + '/reaction-set-6-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react_df['Substrate_Set'] = react_df['Substrate_Set'].apply(lambda x: str_to_set(x))
    react_df['Product_Set'] = react_df['Product_Set'].apply(lambda x: str_to_set(x))
    react_df['Measured_Substrate'] = react_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_df['Measured_Product'] = react_df['Measured_Product'].apply(lambda x: str_to_set(x))

    react1_df = react_df[react_df.Direction == 1]
    react2_df = react_df[react_df.Direction == 2]

    react2_df['Substrate_Set'] = react2_df['Metabolite_Set']
    react2_df['Product_Set'] = react2_df['Metabolite_Set']

    react2_df['Measured_Substrate'] = react2_df['Measured_Metabolite']
    react2_df['Measured_Product'] = react2_df['Measured_Metabolite']

    react_df = pd.concat([react1_df, react2_df], axis=0)
    react_df.to_csv(react_set_6_dir + '/reaction-set-6.tsv', sep='\t', index=False)

def preapre_reaction_set_7(react_gem, react_set_7_dir):
    react_gem = react_gem[react_gem['Direction'] == 1]
    two_more_mets_in = react_gem[react_gem['Measured_Metabolite_Count'] > 1]
    print(two_more_mets_in.shape)
    two_more_mets_in.to_csv(react_set_7_dir + '/reaction-set-7.tsv', sep='\t', index=False)
    print(two_more_mets_in.shape)
    
def preapre_reaction_set_8(react_gem, react_set_8_dir):
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
    
def preapre_reaction_set_9(react_gem, react_set_9_dir):
    react_df = react_gem[react_gem['Measured_Metabolite_Count'] > 1]
    print(react_df.shape)
    react_df.to_csv(react_set_9_dir + '/reaction-set-9-basic.tsv', sep='\t', index=False)
    print(react_df.shape)

    react_df['Substrate_Set'] = react_df['Substrate_Set'].apply(lambda x: str_to_set(x))
    react_df['Product_Set'] = react_df['Product_Set'].apply(lambda x: str_to_set(x))
    react_df['Measured_Substrate'] = react_df['Measured_Substrate'].apply(lambda x: str_to_set(x))
    react_df['Measured_Product'] = react_df['Measured_Product'].apply(lambda x: str_to_set(x))

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
    print('met_path', args.met_path)
    print('log_path', args.log_path)
    print('out_dir', args.out_dir)
    
    react_gem = prepare_all_reactions(args.met_path, args.gem_path, args.out_dir)
    
    react_set_1_dir = args.out_dir + '/Reaction-Set-1'
    Path(react_set_1_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_1(react_gem, react_set_1_dir)
    
    react_set_2_dir = args.out_dir + '/Reaction-Set-2'
    Path(react_set_2_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_2(react_gem, react_set_2_dir)
    
    react_set_3_dir = args.out_dir + '/Reaction-Set-3'
    Path(react_set_3_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_3(react_gem, react_set_3_dir)
    
    react_set_4_dir = args.out_dir + '/Reaction-Set-4'
    Path(react_set_4_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_4(react_gem, react_set_4_dir)
    
    react_set_5_dir = args.out_dir + '/Reaction-Set-5'
    Path(react_set_5_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_5(react_gem, react_set_5_dir)
    
    react_set_6_dir = args.out_dir + '/Reaction-Set-6'
    Path(react_set_6_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_6(react_gem, react_set_6_dir)
    
    react_set_7_dir = args.out_dir + '/Reaction-Set-7'
    Path(react_set_7_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_7(react_gem, react_set_7_dir)
    
    react_set_8_dir = args.out_dir + '/Reaction-Set-8'
    Path(react_set_8_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_8(react_gem, react_set_8_dir)
    
    react_set_9_dir = args.out_dir + '/Reaction-Set-9'
    Path(react_set_9_dir).mkdir(parents=True, exist_ok=True)
    preapre_reaction_set_9(react_gem, react_set_9_dir)
    
    sys.stdout = orig_stdout
    log_file.close()
    
if __name__ == "__main__":
    main(parse_args())
