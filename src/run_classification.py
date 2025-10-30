import argparse
import sys
import os
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--met_path", type=str,
                        required=True, default=None)
    parser.add_argument("--feature_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--case", type=str,
                        help="name of the case group",
                        required=True, default=None)
    parser.add_argument("--control", type=str,
                        help="name of the control group",
                        required=True, default=None)
    parser.add_argument("--script_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--log_path", type=str,
                        required=True, default=None)
    args = parser.parse_args()
    return args

def main(args):
    Path(args.out_dir).mkdir(exist_ok=True, parents=True)
    
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    sys.stderr = log_file
    
    #metabolite
    print("Classification with metabolite")
    met_out_dir = args.out_dir + "/metabolite"
    met_log_path = met_out_dir + '/classification.log'
    Path(met_out_dir).mkdir(exist_ok=True, parents=True)
    
    command = "python -W ignore " + args.script_dir + "/classification.py" + \
        " --in_path " + args.met_path + \
        " --case " + args.case + \
        " --control " + args.control + \
        " --out_dir " + met_out_dir + \
        " --log_path " + met_log_path
    
    os.system(command)
    
    #reaction
    print("Classification with reaction")
    react_out_dir = args.out_dir + "/reaction"
    for entry in os.scandir(args.feature_dir):
        if (entry.is_dir()):
            for sub_entry in os.scandir(entry.path):
                if (sub_entry.is_dir()):
                    command = "python -W ignore " + args.script_dir + "/classification.py" + \
                        " --in_path " + sub_entry.path + '/reaction.' + sub_entry.name + '.tsv' + \
                        " --case " + args.case + \
                        " --control " + args.control + \
                        " --out_dir " + react_out_dir + "/" + entry.name + "/" + sub_entry.name + \
                        " --log_path " + react_out_dir + "/" + entry.name + "/" + sub_entry.name + "/classification.log"
                    os.system(command)
                    
    #metabolite+reaction
    print("Classification with metabolite and reaction")
    met_react_out_dir = args.out_dir + "/metabolite+reaction"
    print("Classification with reaction")
    react_out_dir = args.out_dir + "/reaction"
    for entry in os.scandir(args.feature_dir):
        if (entry.is_dir()):
            for sub_entry in os.scandir(entry.path):
                if (sub_entry.is_dir()):
                    command = "python -W ignore " + args.script_dir + "/classification.py" + \
                        " --in_path " + sub_entry.path + '/metabolite.reaction.' + sub_entry.name + '.tsv' + \
                        " --case " + args.case + \
                        " --control " + args.control + \
                        " --out_dir " + met_react_out_dir + "/" + entry.name + "/" + sub_entry.name + \
                        " --log_path " + met_react_out_dir + "/" + entry.name + "/" + sub_entry.name + "/classification.log"
                    os.system(command) 
    sys.stdout = orig_stdout
    sys.stderr = orig_stderr 
    log_file.close()
                
if __name__ == "__main__":
    main(parse_args())
    
