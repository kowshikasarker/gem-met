# imports
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import csv, scipy, os, warnings, sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_predict, cross_validate
from sklearn.model_selection import ParameterGrid
from sklearn.metrics import classification_report, roc_auc_score
import shutil

#np.random.seed(1)

# formatting
#sns.set()

'''params = {
    "legend.fontsize": "x-large",
    "figure.figsize": (15, 10),
    "axes.labelsize": "x-large",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "x-large",
    "ytick.labelsize": "x-large",
}
plt.rcParams.update(params)'''

# turn off spammy warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def random_forest_classifier(treatment):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import LeaveOneOut, cross_val_predict, cross_validate, ParameterGrid
    from sklearn.metrics import classification_report, roc_auc_score, average_precision_score
    import numpy as np
    import json
    import time
    
    
    start_time = time.time()
    
    control = 'No' + treatment
    classes = {treatment: 1, control: 0}
    
    df = pd.read_csv(treatment + '.' + control + '.tsv', sep='\t')
    if(df.shape[1] == 1):
        print('No features')
        return
    '''df['class'] = df['sample'].apply(lambda x: 0 if control in x else 1)
    print(list(zip(df['sample'], df['class'])))
    y_true = np.array(df['class'].values)
    X = df.set_index('sample').drop(columns=['class'])
    print(X.shape, y_true.shape)'''
    
    df['class'] = df['Key'].apply(lambda x: 0 if control in x else 1)
    print(list(zip(df['Key'], df['class'])))
    y_true = np.array(df['class'].values)
    X = df.set_index('Key').drop(columns=['class'])
    print(X.shape, y_true.shape)
    
    
    
    result_dir = treatment
    if os.path.exists(result_dir):
        shutil.rmtree(result_dir)
    os.makedirs(result_dir)
    os.chdir(result_dir)
    

    # Param grid to search for each food
    param_grid = {
        "n_estimators": [1000],
        "oob_score": [True],
        "n_jobs": [-1],
        "random_state": [1],
        "max_features": [None, "sqrt", "log2"],
        "min_samples_leaf": [1, 3, 5],
    }

    original_stdout = sys.stdout
    log_file = open(treatment + '.' + control + '.random_forest.log', 'w')
    sys.stdout = log_file
    print("-------", treatment, "-------")

    # make directory
    

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
            n_jobs=-1,
            return_estimator=True,
        )

        # Update the best parameters
        estimators = cv_result["estimator"]
        for estimator in estimators:
            if best_rf is None or estimator.oob_score_ > best_rf.oob_score_:
                best_rf = estimator
                best_params = params

    print(treatment, "-> Best parameters from grid search:", best_params)

    # Cross-val predict probabilities using leave one out and our new best parameters
    rfc = RandomForestClassifier()
    rfc.set_params(**best_params)
    y_proba = cross_val_predict(
        rfc,
        X.values,
        y_true,
        cv=LeaveOneOut(),
        n_jobs=-1,
        method="predict_proba",
    )

    np.savetxt(fname=treatment + '.' + control + '.classifier_prediction.tsv', header=treatment + '\t' + control, X=y_proba, delimiter='\t')
    
    y_pred = [score.argmax() for score in y_proba]
    
    metric = classification_report(y_true, y_pred, target_names=[control, treatment], output_dict=True)
    metric['auroc'] = round(roc_auc_score(y_true, y_proba[:, 1]), 2)
    metric['auprc'] = round(average_precision_score(y_true, y_proba[:, 1]), 2)
    
    metric_file = open(treatment + '.' + control + '.classifier_performance.json', 'w')
    
    json.dump(metric, metric_file)
    
    print('Accuracy', metric['accuracy'])
    print('AUROC', metric['auroc'])
    print('AUPRC', metric['auprc'])

    # Plot feature importance graph
    feature_idxs = np.argsort(best_rf.feature_importances_)[::-1]
    plt.figure(figsize=(5, 5))
    plt.title(treatment + 'Feature Importances')
    plt.xlabel("Feature #")
    plt.ylabel("Feature Importance")
    plt.plot(sorted(best_rf.feature_importances_, reverse=True))
    plt.savefig(treatment + '.' + control + '.feature-importance.png')
    plt.close()
    best_features = X.columns[feature_idxs[:10]]
    print('Top-10 features for', treatment, ':', best_features)

    # feature means per group write-out
    classes = {0: 'Control', 1: treatment}
    X_gb = X.copy().iloc[:, feature_idxs]
    X_gb["group"] = list(map(lambda i: classes[i], y_true))
    print("group", X_gb["group"])
    X_gb.groupby("group").mean().to_csv(treatment + '.' + control + '.feature-mean.csv')
    X_gb.groupby("group").std().to_csv(treatment + '.' + control + '.feature-std.csv')

    # Feautre importances write out + figure generation
    best_features_list = list(
        zip(
            [X.columns[idx] for idx in feature_idxs],
            [best_rf.feature_importances_[idx] for idx in feature_idxs],
        )
    )
    #best_features_per_food[food] = best_features_list
    with open(treatment + '.' + control + '.feature-importance.csv', "w") as f:
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
        fig.savefig(treatment + '.' + control + '.' + feature + '.boxplot.png', bbox_inches="tight")
        plt.close(fig)
    end_time = time.time()
    print(end_time - start_time, 'seconds')
    
    sys.stdout = original_stdout
    log_file.flush()
    log_file.close()
    

def main():
    #alpha_list = [0.99, 0.85, 0.75, 0.50, 0.25, 0.10]
    #treatments = ['Almond', 'Avocado', 'Barley', 'Broccoli', 'Oats', 'Walnut']
    #feature_types = ['reaction', 'metabolite']
    '''react_sets = ['Reaction-Set-' + str(i) for i in range(1, 5)]
    rwr_list = ['RWR-3.1', 'RWR-3.2']
    alpha_list = [0.85]
    treatments = ['Almond', 'Avocado', 'Barley', 'Broccoli', 'Oats', 'Walnut']
    feature_types = ['reaction']
    
    for react_set in react_sets:
        for rwr in rwr_list:
            for treatment in treatments:
                for feature_type in feature_types:
                    for alpha in alpha_list:

                        print('===', 'alpha', alpha, 'feature_type', feature_type, 'treatment', treatment, '===', flush=True)
                        print()
                        working_dir = '/home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Analysis/Reactions/' + react_set + '/' + rwr + '/alpha=' + str(alpha) + '/classification/' + feature_type 
                        os.chdir(working_dir)
                        random_forest_classifier(treatment)
                        print('\n')'''
    react_sets = ['Reaction-Set-' + str(i) for i in range(1, 5)]
    rc_list = ['RC-1/classification', 'RC-1/classification-met34']
    treatments = ['Almond', 'Avocado', 'Barley', 'Broccoli', 'Oats', 'Walnut']
    
    for react_set in react_sets:
        for rc in rc_list:
            for treatment in treatments:
                print('===', react_set, rc, 'treatment', treatment, '===', flush=True)
                print()
                working_dir = '/home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Analysis/Reactions/' + react_set + '/' + rc
                os.chdir(working_dir)
                random_forest_classifier(treatment)
                print('\n')
    
    react_sets = ['Reaction-Set-' + str(i) for i in range(1, 5, 2)]
    er_list = ['ER-1/classification', 'ER-1/classification-met34']
    treatments = ['Almond', 'Avocado', 'Barley', 'Broccoli', 'Oats', 'Walnut']
    
    for react_set in react_sets:
        for er in er_list:
            for treatment in treatments:
                print('===', react_set, er, 'treatment', treatment, '===', flush=True)
                print()
                working_dir = '/home/ksarker2/Nutrition/New-Metabolomic-Data-Analysis-Phase-2/Analysis/Reactions/' + react_set + '/' + er
                os.chdir(working_dir)
                random_forest_classifier(treatment)
                print('\n')
main()