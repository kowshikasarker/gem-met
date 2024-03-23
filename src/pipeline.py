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
    orig_stderr = sys.stderr
    
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    sys.stderr = log_file
    
    command = 'python3 preprocess-and-compute-features.py -b ' + args.base_path + ' -e ' + args.end_path + ' -g ' + args.gem_met_path + ' -x ' + args.human_gem_path + ' -i ' + args.id_path + ' -l ' + args.log_path + ' -o ' + args.out_dir + ' -m ' + str(args.missing_pct) + ' -a ' + str(args.alpha)
    
    os.system(command)
    
    command = 'python3 run_classification.py -d ' + args.out_dir + '/Features' + ' -l ' + args.out_dir + '/run_classification.log'
    os.system(command)
    
    command = 'python3 summarize_performance.py -d ' + args.out_dir + '/Features' + ' -o ' + args.out_dir + '/Summary' + ' -l ' + args.out_dir + '/summarize_performance.log'
    os.system(command)
    
    sys.stdout = orig_stdout
    sys.stderr = orig_stderr 
    log_file.close()
        
if __name__ == "__main__":
    main(parse_args())
    