import argparse
import sys
import os
import json
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--case", type=str,
                        required=True, default=None)
    parser.add_argument("--control", type=str,
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        required=True, default=None)
    parser.add_argument("--log_path", type=str,
                        required=True, default=None)
    args = parser.parse_args()
    return args


def build_summary_dataframe(parent_dir, feature_root, react_feat_map, metrics, json_filename, baseline_metrics):
    summary_dict = {
        "index": ["Accuracy", "AUROC", "AUPRC"],
        "columns": [],
        "data": [],
        "index_names": [""],
        "column_names": ["Reaction set", "Feature"],
    }

    for react_set_no in range(1, 10):
        for feature_name in react_feat_map[react_set_no]:
            summary_dict["columns"].append((react_set_no, feature_name))
    summary_dict["columns"].append(("Baseline", "Metabolite"))

    for metric in metrics:
        metric_values = []
        for react_set_no in range(1, 10):
            react_set_dir = os.path.join(parent_dir, feature_root, f"reaction-set-{react_set_no}")
            for feature_name in react_feat_map[react_set_no]:
                feature_dir = os.path.join(react_set_dir, feature_name.lower())
                json_path = os.path.join(feature_dir, json_filename)
                metric_dict = json.load(open(json_path))
                metric_values.append(round(metric_dict[metric], 2))
        metric_values.append(baseline_metrics[metric])
        summary_dict["data"].append(metric_values)

    summary_df = pd.DataFrame.from_dict(summary_dict, orient="tight")
    summary_df.index.name = None
    return summary_df


def save_summary_outputs(summary_df, out_dir, prefix, title):
    summary_payload = json.loads(summary_df.to_json(orient="split"))
    with open(os.path.join(out_dir, f"{prefix}.json"), "w") as json_file:
        json.dump(summary_payload, json_file, indent=4)

    plot_df = summary_df.copy()
    plot_df.columns = [
        f"{col[0]} | {col[1]}" if isinstance(col, tuple) else str(col)
        for col in plot_df.columns
    ]

    plt.figure(figsize=(max(8, len(plot_df.columns) * 0.6), 4.5))
    ax = sns.heatmap(
        plot_df,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        cbar=True,
        linewidths=0.5,
        linecolor="white",
    )
    ax.set_title(title)
    ax.set_xlabel("Feature configuration")
    ax.set_ylabel("Metric")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}.png"), dpi=300)
    plt.close()


def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    
    log_file = open(args.log_path, 'w')
    sys.stdout = log_file
    sys.stderr = log_file
    
    json_filename = args.case + "." + args.control + ".classifier_performance.json"
    
    baseline_dict = json.load(open(args.in_dir + '/metabolite/' + json_filename))
    baseline = [round(baseline_dict['accuracy'], 2),
                round(baseline_dict['auroc'], 2),
                round(baseline_dict['auprc'], 2)]
    
    print("baseline", baseline)
    
    react_feat_map = {
        1: ['Change', 'Prob'],
        2: ['Change', 'Ratio', 'Prob'],
        3: ['Change', 'Prob'],
        4: ['Change', 'Ratio', 'Prob'],
        5: ['Prob'],
        6: ['Prob'],
        7: ['Change', 'Ratio', 'Prob'],
        8: ['Change', 'Ratio', 'Prob'],
        9: ['Prob']
    }
    
    metrics = ['accuracy', 'auroc', 'auprc']

    baseline_metrics = dict(zip(metrics, baseline))

    reaction_summary_df = build_summary_dataframe(
        args.in_dir,
        "reaction",
        react_feat_map,
        metrics,
        json_filename,
        baseline_metrics,
    )
    save_summary_outputs(
        reaction_summary_df,
        args.out_dir,
        "reaction",
        "Benchmark: Reaction features only",
    )

    combo_summary_df = build_summary_dataframe(
        args.in_dir,
        "metabolite+reaction",
        react_feat_map,
        metrics,
        json_filename,
        baseline_metrics,
    )
    save_summary_outputs(
        combo_summary_df,
        args.out_dir,
        "metabolite+reaction",
        "Benchmark: Reaction and metabolite features together",
    )

    sys.stdout = orig_stdout
    sys.stderr = orig_stderr 
    log_file.close()
                
if __name__ == "__main__":
    main(parse_args())
