import argparse
import os
import pathlib
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocesses user-provided metabolomic profiles and Human-GEM, computes novel features integerating these two processed sources and performs classification with the features.')
    
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
    
    parser.add_argument("--alpha", type=float,
                        help="Parameter alpha for RWR algorithm",
                        required=True, default=None)
    
    parser.add_argument("--script_path", type=str,
                        help="path to folder ith all scripts",
                        required=True, default=None)
    
    parser.add_argument("--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    
    args = parser.parse_args()
    return args

def main(args):
    pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    
    #preprocess metabolome
    met_out_dir = args.out_dir + '/preprocess/metabolome'
    met_log_path =  met_out_dir + '/preprocess-metabolome.log'
    command = 'python3 -W ignore ' + args.script_path + '/preprocess-metabolome.py' + \
        ' --base_path ' + args.base_path + \
        ' --end_path ' + args.end_path + \
        ' --missing_pct ' + str(args.missing_pct) + \
        ' --user_met_id_path ' + args.user_met_id_path + \
        ' --user_met_name_col ' + args.user_met_name_col + \
        ' --user_met_id_col ' + args.user_met_id_col + \
        ' --gem_path ' + args.gem_path + \
        ' --gem_met_id_path ' + args.gem_met_id_path + \
        ' --gem_met_id_col ' + args.gem_met_id_col + \
        ' --log_path ' + met_log_path + \
        ' --out_dir ' + met_out_dir
      
    os.system(command)
    print("Preprocessing metabolome complete.")
        
    #preprocessing human-gem
    valid_met_path = met_out_dir + '/gem_overlapped_metabolites.tsv'
    gem_out_dir = args.out_dir + '/preprocess/gem'
    gem_log_path =  gem_out_dir + '/preprocess-gem.log'
    
    command = 'python3 -W ignore ' + args.script_path + '/preprocess-gem.py' + \
        ' --gem_path ' + args.gem_path + \
        ' --valid_met_path ' + valid_met_path + \
        ' --log_path ' + gem_log_path + \
        ' --out_dir ' + gem_out_dir
    
    os.system(command)
    print("Preprocessing gem complete.")
    
    #compute features
    met_change_path = met_out_dir + '/gem_overlapped_change_id.tsv'
    feature_out_dir = args.out_dir + '/feature'
    
    react_set_paths = {}
    react_set_out_dirs = {}
    for react_set_no in range(1, 10):
        react_set_paths[react_set_no] = gem_out_dir + '/reaction-set-' + str(react_set_no) + '/reaction-set-' + str(react_set_no) + '.tsv'
        react_set_out_dirs[react_set_no] = feature_out_dir + '/reaction-set-' + str(react_set_no)
    #change
    for react_set_no in [1, 2, 3, 4, 7, 8, 9]:
        change_out_dir = react_set_out_dirs[react_set_no] + '/change'
        change_log_path = change_out_dir + '/compute-change-feature.log'
        
        command = 'python3 -W ignore ' + args.script_path + '/compute-change-feature.py' + \
            ' --react_set_path ' + react_set_paths[react_set_no] + \
            ' --met_change_path ' + met_change_path + \
            ' --valid_met_path ' + valid_met_path + \
            ' --log_path ' + change_log_path + \
            ' --out_dir ' + change_out_dir
        os.system(command)
        
    #ratio
    for react_set_no in [2, 4, 7, 8, 9]:
        ratio_out_dir = react_set_out_dirs[react_set_no] + '/ratio'
        ratio_log_path = ratio_out_dir + '/compute-ratio-feature.log'
        
        command = 'python3 -W ignore ' + args.script_path + '/compute-ratio-feature.py' + \
            ' --react_set_path ' + react_set_paths[react_set_no] + \
            ' --met_change_path ' + met_change_path + \
            ' --valid_met_path ' + valid_met_path + \
            ' --log_path ' + ratio_log_path + \
            ' --out_dir ' + ratio_out_dir
        os.system(command)
        
    #prob    
    for react_set_no in range(1, 10):
        prob_out_dir = react_set_out_dirs[react_set_no] + '/prob'
        prob_log_path = prob_out_dir + '/compute-prob-feature.log'
        command = 'python3 -W ignore ' + args.script_path + '/compute-prob-feature.py' + \
            ' --base_path ' + met_out_dir + '/gem_overlapped_base_id.tsv' + \
            ' --end_path ' + met_out_dir + '/gem_overlapped_end_id.tsv' + \
            ' --met_change_path ' + met_change_path + \
            ' --react_set_path ' + react_set_paths[react_set_no] + \
            ' --valid_met_path ' + valid_met_path + \
            ' --alpha ' + str(args.alpha) + \
            ' --log_path ' + prob_log_path + \
            ' --out_dir ' + prob_out_dir    
        os.system(command)
        
    sys.stdout = orig_stdout
    log_file.close()
        
if __name__ == "__main__":
    main(parse_args())
    
