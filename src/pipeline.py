import argparse
import os
import pathlib
import sys
import subprocess


def execute_command(command, log_file):
    """Executes a command, logs its output, and prints to stderr on failure."""
    log_file.write(f"--- Running command: {command} ---\n")

    process = subprocess.Popen(
        command, shell=True, stdout=log_file, stderr=log_file, text=True
    )
    process.communicate()

    if process.returncode != 0:
        # On failure, also write to the main script's stderr
        sys.stderr.write(f"--- Command failed: {command} ---\n")
        sys.stderr.write(f"--- Return Code: {process.returncode} ---\n")
        # The output is already in the log file, so we don't write it to stderr here.
        # We could, however, read the last few lines of the log file and print them.
        # For now, we'll just point the user to the log file.
        sys.stderr.write(f"--- See log file for details: {log_file.name} ---\n")
        sys.stderr.flush()

    return process.returncode


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocesses user-provided metabolomic profiles and Human-GEM, computes novel features integerating these two processed sources, performs downstream classification with the features, and generates benchmark plots."
    )

    parser.add_argument(
        "--base_path",
        type=str,
        help="path to baseline metabolomics in tsv format",
        required=True,
        default=None,
    )
    parser.add_argument(
        "--end_path",
        type=str,
        help="path to end metabolomics in tsv format",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--case", type=str, help="name of the case group", required=True, default=None
    )
    parser.add_argument(
        "--control",
        type=str,
        help="name of the control group",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--missing_pct",
        type=float,
        help="threshold of missing value percentage for dropping columns",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--user_met_id_path",
        type=str,
        help="path to list of user metabolites (columns in in base_path and end_path) with standard identifiers in tsv format",
        required=True,
        default=None,
    )
    parser.add_argument(
        "--user_met_name_col",
        type=str,
        help="which column in user_met_id_path contains the metabolite names",
        required=True,
        default=None,
    )
    parser.add_argument(
        "--user_met_id_col",
        type=str,
        help="which column in user_met_id_path contains the metabolite standard identifiers",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--gem_path",
        type=str,
        help="path to a gem model in .xlsx format",
        required=True,
        default=None,
    )
    parser.add_argument(
        "--gem_met_id_path",
        type=str,
        help="path to list of human gem metabolites with standard identifiers list in tsv format",
        required=True,
        default=None,
    )
    parser.add_argument(
        "--gem_met_id_col",
        type=str,
        help="which column in gem_met_id_path contains the metabolite standard identifiers",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--alpha",
        type=float,
        help="Parameter alpha for RWR algorithm",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--script_dir",
        type=str,
        help="path to folder ith all scripts",
        required=True,
        default=None,
    )

    parser.add_argument(
        "--log_path", type=str, help="path to log file", required=True, default=None
    )
    parser.add_argument(
        "--out_dir", type=str, help="path to output dir", required=True, default=None
    )

    args = parser.parse_args()
    return args


def run_pipeline(args, log_file):
    """
    Executes the full pipeline, returning True on success and False on failure.
    """
    # preprocess metabolome
    met_out_dir = os.path.join(args.out_dir, "preprocess", "metabolome")
    met_log_path = os.path.join(met_out_dir, "preprocess-metabolome.log")
    command = (
        "python3 -W ignore "
        + os.path.join(args.script_dir, "preprocess-metabolome.py")
        + " --base_path "
        + args.base_path
        + " --end_path "
        + args.end_path
        + " --missing_pct "
        + str(args.missing_pct)
        + " --user_met_id_path "
        + args.user_met_id_path
        + " --user_met_name_col "
        + args.user_met_name_col
        + " --user_met_id_col "
        + args.user_met_id_col
        + " --gem_path "
        + args.gem_path
        + " --gem_met_id_path "
        + args.gem_met_id_path
        + " --gem_met_id_col "
        + args.gem_met_id_col
        + " --log_path "
        + met_log_path
        + " --out_dir "
        + met_out_dir
    )

    if execute_command(command, log_file) != 0:
        return False
    log_file.write("--- Preprocessing metabolome complete. ---\n")

    # preprocessing human-gem
    valid_met_path = os.path.join(met_out_dir, "gem_overlapped_metabolites.tsv")
    gem_out_dir = os.path.join(args.out_dir, "preprocess", "gem")
    gem_log_path = os.path.join(gem_out_dir, "preprocess-gem.log")

    command = (
        "python3 -W ignore "
        + os.path.join(args.script_dir, "preprocess-gem.py")
        + " --gem_path "
        + args.gem_path
        + " --valid_met_path "
        + valid_met_path
        + " --log_path "
        + gem_log_path
        + " --out_dir "
        + gem_out_dir
    )

    if execute_command(command, log_file) != 0:
        return False
    log_file.write("--- Preprocessing gem complete. ---\n")

    # compute features
    met_change_path = os.path.join(met_out_dir, "gem_overlapped_change_id.tsv")
    feature_out_dir = os.path.join(args.out_dir, "feature")

    react_set_paths = {}
    react_set_out_dirs = {}
    for react_set_no in range(1, 10):
        react_set_paths[react_set_no] = os.path.join(
            gem_out_dir,
            f"reaction-set-{react_set_no}",
            f"reaction-set-{react_set_no}.tsv",
        )
        react_set_out_dirs[react_set_no] = os.path.join(
            feature_out_dir, f"reaction-set-{react_set_no}"
        )

    # change
    for react_set_no in [1, 2, 3, 4, 7, 8, 9]:
        change_out_dir = os.path.join(react_set_out_dirs[react_set_no], "change")
        change_log_path = os.path.join(change_out_dir, "compute-change-feature.log")

        command = (
            "python3 -W ignore "
            + os.path.join(args.script_dir, "compute-change-feature.py")
            + " --react_set_path "
            + react_set_paths[react_set_no]
            + " --met_change_path "
            + met_change_path
            + " --valid_met_path "
            + valid_met_path
            + " --log_path "
            + change_log_path
            + " --out_dir "
            + change_out_dir
        )
        if execute_command(command, log_file) != 0:
            return False

    # ratio
    for react_set_no in [2, 4, 7, 8, 9]:
        ratio_out_dir = os.path.join(react_set_out_dirs[react_set_no], "ratio")
        ratio_log_path = os.path.join(ratio_out_dir, "compute-ratio-feature.log")

        command = (
            "python3 -W ignore "
            + os.path.join(args.script_dir, "compute-ratio-feature.py")
            + " --react_set_path "
            + react_set_paths[react_set_no]
            + " --met_change_path "
            + met_change_path
            + " --valid_met_path "
            + valid_met_path
            + " --log_path "
            + ratio_log_path
            + " --out_dir "
            + ratio_out_dir
        )
        if execute_command(command, log_file) != 0:
            return False

    # prob
    for react_set_no in range(1, 10):
        prob_out_dir = os.path.join(react_set_out_dirs[react_set_no], "prob")
        prob_log_path = os.path.join(prob_out_dir, "compute-prob-feature.log")
        command = (
            "python3 -W ignore "
            + os.path.join(args.script_dir, "compute-prob-feature.py")
            + " --base_path "
            + os.path.join(met_out_dir, "gem_overlapped_base_id.tsv")
            + " --end_path "
            + os.path.join(met_out_dir, "gem_overlapped_end_id.tsv")
            + " --met_change_path "
            + met_change_path
            + " --react_set_path "
            + react_set_paths[react_set_no]
            + " --valid_met_path "
            + valid_met_path
            + " --alpha "
            + str(args.alpha)
            + " --case "
            + args.case
            + " --control "
            + args.control
            + " --log_path "
            + prob_log_path
            + " --out_dir "
            + prob_out_dir
        )
        if execute_command(command, log_file) != 0:
            return False

    # classification
    classification_out_dir = os.path.join(args.out_dir, "classification")
    classification_log_path = os.path.join(classification_out_dir, "classification.log")
    command = (
        "python3 -W ignore "
        + os.path.join(args.script_dir, "run_classification.py")
        + " --met_path "
        + os.path.join(met_out_dir, "metabolite.tsv")
        + " --feature_dir "
        + feature_out_dir
        + " --case "
        + args.case
        + " --control "
        + args.control
        + " --script_dir "
        + args.script_dir
        + " --out_dir "
        + classification_out_dir
        + " --log_path "
        + classification_log_path
    )

    if execute_command(command, log_file) != 0:
        return False

    # summary
    summary_out_dir = os.path.join(args.out_dir, "summary")
    summary_log_path = os.path.join(summary_out_dir, "summary.log")
    command = (
        "python3 -W ignore "
        + os.path.join(args.script_dir, "summarize_performance.py")
        + " --in_dir "
        + classification_out_dir
        + " --case "
        + args.case
        + " --control "
        + args.control
        + " --out_dir "
        + summary_out_dir
        + " --log_path "
        + summary_log_path
    )
    if execute_command(command, log_file) != 0:
        return False

    return True


def main(args):
    pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    with open(args.log_path, "w") as log_file:
        success = run_pipeline(args, log_file)

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main(parse_args())

