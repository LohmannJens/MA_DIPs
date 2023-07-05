'''
    Tests different classifers on different data sets. Also tests which
    combination of features is the best.
'''
import os
import sys
import shap
import joblib
import logging
import argparse
import datetime
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix, make_scorer

from ml_utils import generate_features
from ml_utils import CHARS, CHARS_COUNT, MAX_LEN
from classifier import select_classifier

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_pelz_dataset


def timepoints_set_labels(df: pd.DataFrame)-> pd.Series:
    '''
        Set the labels for the classification by labeling them at each time
        point. Then the sum of all "high" is calculated and normalized.
        :param df: Dataframe including the NGS count for the single time points

        :return: Pandas Series with the labels
    '''
    old_col_n = df.shape[1]
    for i in range(3, old_col_n):
        col = df.iloc[:,i]
        median = col.replace(0, np.nan).median()
      
        df[f"{i}"] = np.where((col == 0), 0, np.where((col < median), "low", "high"))


        # count numbers of "high" and "low"
    count_df = df.iloc[:,old_col_n:df.shape[1]].apply(pd.Series.value_counts, axis=1)
    count_df.replace(np.nan, 0, inplace=True)

    score = count_df["high"] / (count_df["high"] + count_df["low"])
    labels = pd.Series(np.where((score < 0.5), "low", "high"))
    return labels


def test_classifiers(df: pd.DataFrame,
                     folder: str,
                     perform_grid_search: bool
                     )-> None:
    '''
        Tests different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param folder: folder where to save the plots
        :param perform_grid_search: True if a grid search should be performed

        :return: None
    '''
    y = timepoints_set_labels(df)

    df["Strain"] = "PR8"

    # add features
    features = ["Segment", "DI_Length", "Direct_repeat","Junction", "3_5_ratio", "length_proportion", "delta_G", "Peptide_Length"]
    features = ["Segment", "DI_Length", "Direct_repeat","Junction", "3_5_ratio", "length_proportion", "Peptide_Length"]
    df, feature_cols = generate_features(df, features, load_precalc=False)

    features_out = "\n\t".join(features)
    logging.info(f"\nUsed features:\n\t{features_out}\n#####\n")

    # Selecting train/test and validation data sets
    X, X_val, y, y_val = train_test_split(df[feature_cols], y, test_size=0.2, random_state=42)

    logging.info("Distribution of labels:")
    logging.info(y.value_counts())
    logging.info(y_val.value_counts())
    logging.info("#####\n")

    # Testing different classifiers
    clf_names = ["logistic_regression", "svc", "random_forest", "mlp", "ada_boost", "naive_bayes"]
    data_dict = dict()
    data_dict["param"] = ["accuracy"]
    for clf_name in clf_names:
        print(clf_name)
        logging.info(f"\n### {clf_name} ###")

        data_dict[clf_name] = list()
        # setting up classifier and k-fold validation
        clf, param_grid = select_classifier(clf_name, grid_search=perform_grid_search)
        skf = StratifiedKFold(n_splits=5)
        scorers = {"accuracy_score": make_scorer(accuracy_score)}

        # perform grid search for best parameters
        grid_search = GridSearchCV(clf, param_grid, scoring=scorers, cv=skf, return_train_score=True, refit="accuracy_score")
        grid_search.fit(X, y)

        print(f"training acc.:\t{grid_search.best_score_}")
        logging.info(f"training acc.:\t{grid_search.best_score_}")

        if perform_grid_search:
            print(grid_search.best_params_)
            logging.info(grid_search.best_params_)

        # save model as pickle file (for R shiny application)
        path = os.path.join(folder, f"{clf_name}.pkl")
        joblib.dump(grid_search, path, compress=3)

        # fit on overall model and create confusion matrix for validation set
        predicted_val = grid_search.predict(X_val)
        acc_score = accuracy_score(predicted_val, y_val)
        confusion_m = confusion_matrix(predicted_val, y_val)

        print(f"validation acc.:{acc_score}")
        print(confusion_m)
        logging.info(f"validation acc.:{acc_score}")
        logging.info(confusion_m)
        data_dict[clf_name].append(acc_score)

        # if two classes given create a ROC
        if len(y.unique()) == 2:
            plt.rc("font", size=14)
            fig, axs = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)

            y_shuffled = y_val.sample(frac=1, random_state=42, ignore_index=True).to_numpy()

            RocCurveDisplay.from_estimator(grid_search, X_val, y_val, name=clf_name, ax=axs)
            RocCurveDisplay.from_estimator(grid_search, X_val, y_shuffled, name="shuffled", ax=axs)
            plt.plot([0,1], [0,1])
            path = os.path.join(folder, f"{clf_name}_roc_curve.png")
            plt.savefig(path)
            plt.close()
        # compare to shuffled labels, when three labels are selected
        elif len(y.unique()) == 3:
            y_shuffled = y_val.sample(frac=1, random_state=42, ignore_index=True).to_numpy()
            shuf_acc_score = accuracy_score(predicted_val, y_shuffled)
            print(f"validation acc. shuffled:{shuf_acc_score}")
            logging.info(f"validation acc. shuffled:{shuf_acc_score}")

    o_df = pd.DataFrame(data_dict)
    o_df["mean"] = o_df.mean(axis=1)
    path = os.path.join(folder, f"means.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)

    # initialize logger
    folder_name = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    folder = os.path.join(RESULTSPATH, "ML", folder_name)
    os.makedirs(folder)
    path = os.path.join(folder, "results.log")

    logging.basicConfig(filename=path, filemode="w", format="%(message)s", level=logging.INFO)
    logging.info(f"Script started on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info("using this command:")
    a = sys.argv
    a.insert(0, "python")
    logging.info(f"\t{' '.join(a)}\n#####")

    # Loading the dataset
    df = load_pelz_dataset(by_time=True)["PR8"]
    
    test_classifiers(df, folder, True)

    logging.info(f"Script ended on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

