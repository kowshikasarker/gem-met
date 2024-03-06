import argparse
import os
import pathlib
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocesses metablomic profiles and Human-GEM, computes features integerating these two processed sources and performs classification with the features.')
    
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
    parser.add_argument("-a", "--alpha", type=float,
                        help="Parameter alpha for RWR algorithm", required=True, default=None)
    args = parser.parse_args()
    return args

def main(args):
    pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    
    #preprocess metabolome
    b = args.base_path
    e = args.end_path
    g = args.gem_met_path
    x = args.human_gem_path
    i = args.id_path
    o = args.out_dir + '/Preprocessed-Metabolome'
    l = o + '/preprocess-metabolome.log'
    m = str(args.missing_pct)
    
    command = 'python3 -W ignore preprocess-metabolome.py -b ' + b + ' -e ' + e + ' -g ' + g + ' -x ' + x + ' -i ' + i + ' -l ' + l + ' -o ' + o + ' -m ' + m
    
    os.system(command)
    print("Preprocessing metabolome complete.")
    
    #preprocessing human-gem
    m = args.out_dir + '/Preprocessed-Metabolome/gem_overlapped_metabolites.tsv'
    o = args.out_dir + '/Preprocessed-Human-GEM'
    l = o + '/preprocess-human-gem.log'
    x = args.human_gem_path
    
    command = 'python3 -W ignore preprocess-human-gem.py -x ' + x + ' -m ' + m + ' -o ' + o + ' -l ' + l
    
    os.system(command)
    
    print("Preprocessing human-gem complete.")
    
    #compute features
    c = args.out_dir + '/Preprocessed-Metabolome/gem_overlapped_change_hmdb.tsv'
    m = args.out_dir + '/Preprocessed-Metabolome/gem_overlapped_metabolites.tsv'
    #change
    for react_set_no in [1, 2, 3, 4, 7, 8, 9]:
        r = args.out_dir + '/Preprocessed-Human-GEM/Reaction-Set-' + str(react_set_no) + '/reaction-set-' + str(react_set_no) + '.tsv'
        o = args.out_dir + '/Features/Reaction-Set-' + str(react_set_no) + '/Change'
        l = o + '/compute-change-feature.log'
        
        command = 'python3 -W ignore compute-change-feature.py -r ' + r + ' -c ' + c + ' -m' + m + ' -o ' + o + ' -l ' + l
        os.system(command)
        
    #ratio
    for react_set_no in [2, 4, 7, 8, 9]:
        r = args.out_dir + '/Preprocessed-Human-GEM/Reaction-Set-' + str(react_set_no) + '/reaction-set-' + str(react_set_no) + '.tsv'
        o = args.out_dir + '/Features/Reaction-Set-' + str(react_set_no) + '/Ratio'
        l = o + '/compute-ratio-feature.log'
        
        command = 'python3 -W ignore compute-ratio-feature.py -r ' + r + ' -c ' + c + ' -m' + m + ' -o ' + o + ' -l ' + l
        os.system(command)
        
    #prob    
    b = args.out_dir + '/Preprocessed-Metabolome/gem_overlapped_base_hmdb.tsv'
    e = args.out_dir + '/Preprocessed-Metabolome/gem_overlapped_end_hmdb.tsv'
    a = str(args.alpha)
    for react_set_no in range(1, 10):
        r = args.out_dir + '/Preprocessed-Human-GEM/Reaction-Set-' + str(react_set_no) + '/reaction-set-' + str(react_set_no) + '.tsv'
        o = args.out_dir + '/Features/Reaction-Set-' + str(react_set_no) + '/Prob'
        l = o + '/compute-prob-feature.log'
        command = 'python3 -W ignore compute-prob-feature.py -b ' + b + ' -e ' + e + ' -a ' + a + ' -r ' + r + ' -c ' + c + ' -m' + m + ' -o ' + o + ' -l ' + l
        os.system(command)
        
    sys.stdout = orig_stdout
    log_file.close()
        
if __name__ == "__main__":
    main(parse_args())
    
