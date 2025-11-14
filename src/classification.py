# export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
# to resolve
# ImportError: /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by /home/ksarker2/miniconda3/envs/metabolite/lib/python3.10/site-packages/kiwisolver/_cext.cpython-310-x86_64-linux-gnu.so)

# imports
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv, os, warnings, sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict, cross_validate, ParameterGrid
from sklearn.metrics import classification_report, roc_auc_score, average_precision_score
import shutil
import json
from pathlib import Path

# turn off spammy warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def random_forest_classifier(in_path, case, control, out_dir, log_path):
    classes = {case: 1, control: 0}
    df = pd.read_csv(in_path, sep='\t')
        
    df['class'] = df['sample_group'].apply(lambda x: 1 if x == case else 0)
    
    y_true = np.array(df['class'].values)
    X = df.set_index('sample_id').drop(columns=['class', 'sample_group'])
    print(X.shape, y_true.shape)
    
    old_cwd = os.getcwd()
    try:
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.makedirs(out_dir)
        os.chdir(out_dir)

        # Param grid to search for each food
        # Note: keep this lightweight to avoid very long runtimes and
        # avoid joblib multiprocessing issues in some environments.
        param_grid = {
            "n_estimators": [200],
            "oob_score": [True],
            "n_jobs": [1],
            "random_state": [1],
            "max_features": ["sqrt"],
            "min_samples_leaf": [1, 3],
        }

        # Ensure the log directory exists before opening the file
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        log_file = open(log_path, 'w')

        original_stdout = sys.stdout
        original_stderr = sys.stderr
        try:
            sys.stdout = log_file
            sys.stderr = log_file

            # Grid search
            best_rf = None
            best_params = None
            for params in ParameterGrid(param_grid):
                rfc = RandomForestClassifier()
                rfc.set_params(**params)

                # Perform LOO evaluation for this parameter set
                cv_result = cross_validate(
                    rfc,
                    X.values,
                    y_true,
                    scoring=None,
                    cv=LeaveOneOut(),
                    n_jobs=1,
                    return_estimator=True,
                )

                # Update the best parameters
                estimators = cv_result["estimator"]
                for estimator in estimators:
                    if best_rf is None or estimator.oob_score_ > best_rf.oob_score_:
                        best_rf = estimator
                        best_params = params

            print(case, 'vs', control, "-> Best parameters from grid search:", best_params)

            # Cross-val predict probabilities using leave one out and our new best parameters
            rfc = RandomForestClassifier()
            rfc.set_params(**best_params)
            y_proba = cross_val_predict(
                rfc,
                X.values,
                y_true,
                cv=LeaveOneOut(),
                n_jobs=1,
                method="predict_proba",
            )

            np.savetxt(fname=case + '.' + control + '.classifier_prediction.tsv', header=case + '\t' + control, X=y_proba, delimiter='\t')
            
            y_pred = [score.argmax() for score in y_proba]
            
            metric = classification_report(y_true, y_pred, target_names=[control, case], output_dict=True)
            metric['auroc'] = round(roc_auc_score(y_true, y_proba[:, 1]), 2)
            metric['auprc'] = round(average_precision_score(y_true, y_proba[:, 1]), 2)
            
            metric_file = open(case + '.' + control + '.classifier_performance.json', 'w')
            # Write metrics and force flush/close immediately to avoid data loss
            json.dump(metric, metric_file)
            metric_file.flush()
            os.fsync(metric_file.fileno())
            metric_file.close()
            
            print('Accuracy', metric['accuracy'])
            print('AUROC', metric['auroc'])
            print('AUPRC', metric['auprc'])

            # Plot feature importance graph
            feature_idxs = np.argsort(best_rf.feature_importances_)[::-1]
            plt.figure(figsize=(5, 5))
            plt.title(case + 'Feature Importances')
            plt.xlabel("Feature #")
            plt.ylabel("Feature Importance")
            plt.plot(sorted(best_rf.feature_importances_, reverse=True))
            plt.savefig(case + '.' + control + '.feature-importance.png')
            plt.close()
            best_features = X.columns[feature_idxs[:10]]
            print('Top-10 features for', case, ':', best_features)

            # feature means per group write-out
            classes = {0: control, 1: case}
            X_gb = X.copy().iloc[:, feature_idxs]
            X_gb["group"] = list(map(lambda i: classes[i], y_true))
            print("group", X_gb["group"])
            X_gb.groupby("group").mean().to_csv(case + '.' + control + '.feature-mean.csv')
            X_gb.groupby("group").std().to_csv(case + '.' + control + '.feature-std.csv')

            # Feautre importances write out + figure generation
            best_features_list = list(
                zip(
                    [X.columns[idx] for idx in feature_idxs],
                    [best_rf.feature_importances_[idx] for idx in feature_idxs],
                )
            )
            #best_features_per_food[food] = best_features_list
            with open(case + '.' + control + '.feature-importance.csv', "w") as f:
                w = csv.writer(f)
                w.writerow(["Feature", "Importance"])
                for idx in feature_idxs:
                    w.writerow([X.columns[idx], best_rf.feature_importances_[idx]])

            # plot features
            for feature in best_features:
                fig, ax = plt.subplots(figsize=(5, 5))
                sns.boxplot(
                    data=X_gb, x="group", y=feature, linewidth=2.5, width=0.4, ax=ax
                )
                ax.set_xlabel("Treatment group")
                ax.set_ylabel("Relative concentration")
                fig.suptitle(feature, size=22)
                fig.savefig(case + '.' + control + '.' + feature + '.boxplot.png', bbox_inches="tight")
                plt.close(fig)
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            log_file.flush()
            log_file.close()
    finally:
        os.chdir(old_cwd)
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_path", type=str,
                        help="path to input features in .tsv format",
                        required=True, default=None)
    parser.add_argument("--case", type=str,
                        help="name of the case group",
                        required=True, default=None)
    parser.add_argument("--control", type=str,
                        help="name of the control group",
                        required=True, default=None)
    parser.add_argument("--out_dir", type=str,
                        help="path to output dir",
                        required=True, default=None)
    parser.add_argument("--log_path", type=str,
                        help="path to log file",
                        required=True, default=None)
    args = parser.parse_args()
    return args

def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    # Run classifier first (may recreate the out_dir). Then copy the input for provenance.
    random_forest_classifier(args.in_path, args.case, args.control, args.out_dir, args.log_path)
    try:
        shutil.copy(args.in_path, args.out_dir)
    except Exception as e:
        print(f"Warning: failed to copy input to output dir: {e}")
    
if __name__ == "__main__":
    main(parse_args())
