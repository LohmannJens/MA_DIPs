'''
    Trains and tests different classifiers using different labels as output.
    The labels are: gain, loss, de novo gain, and de novo loss. These are
    taken from Pelz et al. 2021
'''
import os
import sys
import warnings

import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import GridSearchCV, StratifiedKFold , train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, make_scorer, f1_score
from sklearn.cluster import KMeans
from sklearn.preprocessing import power_transform

from classifier import select_classifier

from ml_utils import select_datasets, generate_features

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from utils import load_pelz_dataset


def test_classifiers(df: pd.DataFrame,
                     dataset_name: str,
                     perform_grid_search: bool
                     )-> None:
    '''
        Tests six different classifiers on a given dataset.
        :param df: data frame containing pelz dataset
        :param dataset_name: string indicating which datasets to use as train/
                             test and validation data set
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from
        :param perform_grid_search: True if a grid search should be performed

        :return: None
    '''
    # add features
    features = ["DI_length", "Direct_repeat", "Segment", "Junction", "3_5_ratio", "length_proportion" ,"full_sequence",]# "delta_G"]
    df, feature_cols = generate_features(df, features, load_precalc=False)

    df["test"] = "low"

    file_path = os.path.join(DATAPATH, "Pelz2021", "top_candidates_ranked.xlsx")
    top_df = pd.read_excel(io=file_path,
                            sheet_name=None,
                            header=0,
                            na_values=["", "None"],
                            keep_default_na=False)["PR8"]

    all_candidates = top_df["DI"].tolist()

    n = 100
    top_candidates = all_candidates[0:n]
    df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
    indices = list(df[df["DI"].isin(top_candidates)].index)

    for idx in indices:
        df.iloc[idx, df.columns.get_loc("test")] = "high"

    val_ind = list(random.sample(range(len(indices)), 50))
    val_indices = [x for i,x in enumerate(indices) if i in val_ind]
    indices = [x for i,x in enumerate(indices) if not i in val_ind]

    val_df = df.iloc[val_indices].copy()
    X_val = val_df[feature_cols]
    y_val = val_df["test"]
    col = X_val.columns
    X_val = pd.DataFrame(power_transform(X_val), columns=col)

    copy_df = df.iloc[indices].copy()
    copy_df = pd.concat([copy_df]*16, ignore_index=True)    

    df = df.loc[((df["class"].isin(["loss", "de_novo_loss"])) | (df["test"] == "high"))].reset_index(drop=True)

    df = pd.concat([df, copy_df], ignore_index=True)

    # Selecting train/test and validation data sets
    X = df[feature_cols]
    y = df["test"]

    col = X.columns
    X = pd.DataFrame(power_transform(X), columns=col)

  #  X, X_val, y, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

    # Testing different classifiers
    clf_names = ["logistic_regression", "svc", "random_forest", "mlp", "ada_boost", "naive_bayes"]
    clf_names = ["random_forest", "mlp", "ada_boost", "naive_bayes"]
    data_dict = dict()
    data_dict["param"] = ["validation accuracy"]
    for clf_name in clf_names:
        print(clf_name)
        data_dict[clf_name] = list()
        # setting up classifier and k-fold validation
        clf, param_grid = select_classifier(clf_name, grid_search=perform_grid_search)
        skf = StratifiedKFold(n_splits=5)
        scorers = {"accuracy_score": make_scorer(accuracy_score)}

        # perform grid search for best parameters
        grid_search = GridSearchCV(clf, param_grid, scoring=scorers, cv=skf, return_train_score=True, refit="accuracy_score")
        grid_search.fit(X, y)

        if perform_grid_search:
            print(grid_search.best_params_)

        # fit on overall model and create confusion matrix for validation set
        predicted_val = grid_search.predict(X_val)
        acc_score = f1_score(predicted_val, y_val, average="weighted")
        confusion_m = confusion_matrix(predicted_val, y_val)

        print(acc_score)
        print(confusion_m)
        data_dict[clf_name].append(acc_score)

#    o_df = pd.DataFrame(data_dict)
#    path = os.path.join(RESULTSPATH, "ML", f"{dataset_name}_means.tex")
 #   o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)

    # Loading the dataset
    df = load_pelz_dataset()["PR8"]
    df["Strain"] = "PR8"
    perform_grid_search = True

    test_classifiers(df, "Pelz", perform_grid_search)
