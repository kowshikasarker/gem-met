import argparse
import pandas as pd
import networkx as nx
import pickle
import numpy as np
from pathlib import Path
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_path", type=str,
                        help="path to baseline metabolomics profile in .tsv format",
                        required=True, default=None)
    parser.add_argument("--end_path", type=str,
                        help="path to end metabolomics profile in .tsv format",
                        required=True, default=None)
    parser.add_argument("--met_change_path", type=str,
                        help="path to metabolite concentration change in .tsv format", required=True, default=None)
    
    parser.add_argument("--react_set_path", type=str,
                        help="path to reaction set in .tsv format",
                        required=True, default=None)
    
    parser.add_argument("--valid_met_path", type=str,
                        help="path to list of human-gem overlapped metabolites in .tsv format", required=True, default=None)
    
    parser.add_argument("--alpha", type=float,
                        help="Parameter alpha for RWR algorithm", required=True, default=None)
    
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
    return cell

def analyze_single_study_per_sample(G, treatment, alpha, out_dir):
    control = 'No' + treatment
    print('***', treatment, control, '***')
    
    G = G.copy()
    
    node_label = nx.get_node_attributes(G, 'label')
    
    react_nodes = set([node for node, label in node_label.items() if label=='reaction'])
    #print(len(react_nodes), 'react_nodes')
    
    met_nodes = set([node for node, label in node_label.items() if label=='metabolite'])
    #print(len(met_nodes), 'met_nodes')
    
    sample_nodes = set([node for node, label in node_label.items() if label=='sample'])
    #print(len(sample_nodes), 'sample_nodes')
    #print(sample_nodes)
    
    study_samples = set([s for s in sample_nodes if (treatment in s)])
    #print('study_samples', len(study_samples))
    
    control_samples = set([s for s in study_samples if (control in s)])
    #print('control_samples', len(control_samples))
    
    treatment_samples = study_samples.difference(control_samples)
    #print('treatment_samples', len(treatment_samples))
    
    invalid_samples = sample_nodes.difference(study_samples)
    #print('invalid_samples', len(invalid_samples))
    #print('invalid_samples', len(invalid_samples), invalid_samples)
    
    #print('Before removing invalid_samples and isolates', G.number_of_nodes(), 'nodes', G.number_of_edges(), 'edges')
    
    G.remove_nodes_from(invalid_samples)
    G.remove_nodes_from(list(nx.isolates(G)))
    
    #print('After removing invalid_samples and isolaets', G.number_of_nodes(), 'nodes', G.number_of_edges(), 'edges')
    
    eq_prob = []
    
    #print('study_samples', len(study_samples))
    
    for key in study_samples:
        prob = nx.pagerank(G, alpha=alpha, personalization={key: 1}, weight="weight", max_iter=1000)
        prob['key'] = key
        eq_prob.append(prob)
        
    #print('eq_prob', eq_prob)
    df = pd.DataFrame.from_records(eq_prob, index='key')
    df.to_csv(out_dir + '/' + treatment + '.' + control + '.prob.tsv', sep='\t')
    

def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    
    base_df = pd.read_csv(args.base_path, sep='\t', index_col='key')
    base_df = base_df.sort_index().sort_index(axis=1)

    end_df = pd.read_csv(args.end_path, sep='\t', index_col='key')
    end_df = end_df.sort_index().sort_index(axis=1)
    
    change_df = pd.read_csv(args.met_change_path, sep='\t')
    change_df[['person_id', 'diet']] = change_df['key'].str.split('.', expand=True)
    
    
    treatments = set(change_df.diet)
    treatments = [treatment for treatment in treatments if not treatment.startswith('No')]
    print('treatments', treatments)
    change_df = change_df.drop(columns=['person_id', 'diet']).set_index('key')

    change_df = change_df.sort_index().sort_index(axis=1)

    change_direction = change_df > 0 # True(1) -> increase (M+ node), False(0) -> decrease (M- node)
    change_magnitude = change_df.abs()

    sample_met_inc = change_magnitude.div(end_df).fillna(0) # if end_df concentration is zero
    sample_met_inc = sample_met_inc.replace([np.inf, -np.inf], 0)

    sample_met_dec = change_magnitude.div(base_df).fillna(0)
    sample_met_dec = sample_met_dec.replace([np.inf, -np.inf], 0)

    sample_met_weight = sample_met_inc * change_direction + sample_met_dec * (1 - change_direction)
    sample_met_weight = np.exp(sample_met_weight)
    nomalizer = sample_met_weight.sum(axis=1)
    sample_met_prob = sample_met_weight.div(nomalizer, axis=0)

    hmdb_list = sample_met_prob.columns
    print('hmdb_list', len(hmdb_list))

    sample_met_edges = []
    for key, metabolites in sample_met_prob.iterrows():
        #print(person_id)
        for met in hmdb_list:
            #print(type(met), met)
            if(change_direction.loc[key][met]): # true -> increase
                sample_met_edges.append((key, met + '+', metabolites[met]))
            else:
                sample_met_edges.append((key, met + '-', metabolites[met]))
    sample_met_edges[:5]

    met_sample_edges = []
    for met in hmdb_list:
        #print(change_direction[met])
        met_inc = sample_met_inc[change_direction[met]][met]
        #print(met_inc)
        met_inc = np.exp(met_inc)
        met_inc = met_inc / met_inc.sum()
        #print('met_inc.sum', met_inc.sum())

        for key, inc_val in met_inc.items():
            met_sample_edges.append((met + '+', key, inc_val))

        met_dec = sample_met_dec[~change_direction[met]][met]
        met_dec = np.exp(met_dec)
        met_dec = met_dec / met_dec.sum()
        #print('met_dec.sum', met_dec.sum())

        for key, dec_val in met_dec.items():
            met_sample_edges.append((met + '-', key, dec_val))

    react_df = pd.read_csv(args.react_set_path, sep='\t', index_col='RXN_ID')

    react_df['Measured_Substrate'] = react_df['Measured_Substrate'].apply(str_to_set)
    react_df['Measured_Product'] = react_df['Measured_Product'].apply(str_to_set)


    id_df = pd.read_csv(args.valid_met_path, sep='\t')

    met_to_hmdb = dict(zip(id_df.MET_ID, id_df.ID))

    met_react_map = {}
    for hmdb in hmdb_list:
        met_react_map[hmdb+'+'] = []
        met_react_map[hmdb+'-'] = []

    react_met_edges = []

    met_to_hmdb_keys = met_to_hmdb.keys() # measured mets # for discarding

    for rxn_id, rxn_details in react_df.iterrows():
        w = 1 / rxn_details['Measured_Metabolite_Count']
        substrate_set = rxn_details['Measured_Substrate']
        product_set = rxn_details['Measured_Product']

        substrate_set = [met_to_hmdb[s] for s in substrate_set if (s in met_to_hmdb_keys)]

        for hmdb in substrate_set:
            react_met_edges.append((rxn_id, hmdb+'-', w))
            met_react_map[hmdb+'-'].append(rxn_id)


        product_set = [met_to_hmdb[p] for p in product_set if (p in met_to_hmdb_keys)]
        for hmdb in product_set:
            react_met_edges.append((rxn_id, hmdb+'+', w))
            met_react_map[hmdb+'+'].append(rxn_id)

        if(rxn_details['Measured_Metabolite_Count'] != (len(substrate_set) + len(product_set))):
            print('Count mismatch', rxn_id)

    met_react_edges = []
    no_react_count = 0
    no_react_mets = set()
    for met in met_react_map.keys():
        react_map = met_react_map[met]
        if(len(react_map) > 0):
            w = 1 / len(react_map)
            for rxn in react_map:
                met_react_edges.append((met, rxn, w))
        else:
            print('No reaction mapped to', met)
            no_react_count += 1
            no_react_mets.add(met)
    #print('no_react_count', no_react_count)

    all_edges = sample_met_edges + met_sample_edges + react_met_edges + met_react_edges


    G = nx.DiGraph()
    G.add_weighted_edges_from(all_edges)
    G.remove_nodes_from(list(nx.isolates(G)))
    nodes = set(G.nodes)
    react_nodes = set([n for n in nodes if n.startswith('MAR')])
    met_nodes = set([n for n in nodes if (n.endswith('+') or n.endswith('-'))])
    sample_nodes = nodes.difference(react_nodes).difference(met_nodes)
    node_label = {}
    for n in sample_nodes:
        node_label[n] = 'sample'
    for n in react_nodes:
        node_label[n] = 'reaction'
    for n in met_nodes:
        node_label[n] = 'metabolite'

    nx.set_node_attributes(G, node_label, name='label')
    sample_nodes = nodes.difference(react_nodes).difference(met_nodes)
    pickle.dump(G, open(args.out_dir + '/network.pickle', 'wb'))
    
    for treatment in treatments:
        analyze_single_study_per_sample(G, treatment, args.alpha, args.out_dir)
        
    sys.stdout = orig_stdout
    log_file.close()
        
if __name__ == "__main__":
    main(parse_args())
    
    

