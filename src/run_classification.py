import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--working_dir", type=str,
                        help="path to the directory containing features in .tsv format computed by the script preprocess_and_compute_features.py",
                        required=True, default=None)
    parser.add_argument("-l", "--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    args = parser.parse_args()
    return args

def main(args):
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    sys.stderr = log_file
    
    treatments = ['Almond', 'Walnut', 'Avocado', 'Barley', 'Broccoli', 'Oats']
    react_feat_map = {
        1: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        2: ['Change', 'Prob/Reaction', 'Prob/Subsystem', 'Ratio'],
        3: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        4: ['Change', 'Prob/Reaction', 'Prob/Subsystem', 'Ratio'],
        5: ['Prob/Reaction', 'Prob/Subsystem'],
        6: ['Prob/Reaction', 'Prob/Subsystem'],
        7: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        8: ['Change', 'Prob/Reaction', 'Prob/Subsystem'],
        9: ['Prob/Reaction', 'Prob/Subsystem']
    }
    for react_set_no in range(1, 10):
        react_set_dir = args.working_dir + '/Reaction-Set-' + str(react_set_no)
        print()
        for feature_path in react_feat_map[react_set_no]:
            feature_dir = react_set_dir + '/' + feature_path
            feature_type = (feature_path.split('/')[0]).lower()
            print(feature_type)
            for treatment in treatments:
                control = 'No' + treatment
                
                i = feature_dir + '/' + treatment + '.' + control + '.' + feature_type + '.tsv'
                o = feature_dir + '/' + treatment
                l = o + '/log.txt'
                
                command = 'python3 /home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Analysis/TCBB-New-Experiments/classification.py -i ' + i + ' -o ' + o + ' -l ' + l + ' -t ' + treatment
                os.system(command)
        
        
    sys.stdout = orig_stdout
    sys.stderr = orig_stderr 
    log_file.close()
                
if __name__ == "__main__":
    main(parse_args())
    