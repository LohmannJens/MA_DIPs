'''
    Tests different classifers on different data sets. Also tests which
    combination of features is the best.
'''
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import KFold, GridSearchCV, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix, make_scorer, precision_score, recall_score

from classifier import select_classifier

from ml_utils import load_all_sets, select_datasets, segment_ohe, junction_site_ohe, get_dirna_length, get_direct_repeat_length, get_3_to_5_ratio, get_length_proportion, full_sequence_ohe

sys.path.insert(0, "..")
from utils import RESULTSPATH
from utils import load_pelz_dataset


def test_classifiers(df: pd.DataFrame,
                     dataset_name: str,
                     label_style: str,
                     perform_grid_search: bool
                     )-> None:
    '''
        Tests three different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param dataset_name: string indicating which datasets to use as train/
                             test and validation data set
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from
        :param perform_grid_search: True if a grid search should be performed

        :return: None
    '''
    # add features
    print(df)

    feature_cols = ["Start", "End"]
    df, segment_cols = segment_ohe(df)
    df["DI_Length"] = df.apply(get_dirna_length, axis=1)
    feature_cols.append("DI_Length")
    df["Direct_repeat"] = df.apply(get_direct_repeat_length, axis=1)
    feature_cols.append("Direct_repeat")
    df, junction_start_cols = junction_site_ohe(df, "Start")
    df, junction_end_cols = junction_site_ohe(df, "End")
    feature_cols = feature_cols + segment_cols + junction_start_cols + junction_end_cols
    df["3_5_ratio"] = df.apply(get_3_to_5_ratio, axis=1)
    feature_cols.append("3_5_ratio")
    df["length_proportion"] = df.apply(get_length_proportion, axis=1)
    feature_cols.append("length_proportion")
#    df, sequence_cols = full_sequence_ohe(df)
 #   feature_cols = feature_cols + sequence_cols

    # Selecting train/test and validation data sets
    X = df[feature_cols]
    y = df["class"]

    # Testing different classifiers
    clf_names = ["logistic_regression", "svc", "random_forest", "mlp", "ada_boost", "naive_bayes"]
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
        predicted_val = grid_search.predict(X)
        acc_score = accuracy_score(predicted_val, y)
        confusion_m = confusion_matrix(predicted_val, y)

        print(acc_score)
        print(confusion_m)
        data_dict[clf_name].append(acc_score)

    o_df = pd.DataFrame(data_dict)
    path = os.path.join(RESULTSPATH, "ML", f"{dataset_name}_means.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


def test_model(df: pd.DataFrame,
               clf: object,
               f_list: list,
               d_name: str,
               label_style: str,
               )-> float:
    '''
        Fits given data to a given classifier class and returns the accuracy.
        Is used to test different feature combinations.
        :param df: pandas data frame
        :param clf: classifier class (from scikit-learn)
        :param f_list: list of features to use for testing
        :param d_name: string indicating which data sets to use
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from

        :return: Accuracy of the prediciton
    '''
    X = df[f_list]
    y = df["class"]

    clf.fit(X, y)
    y_pred = clf.predict(X)
    acc = accuracy_score(y_pred, y)
    confusion_m = confusion_matrix(y_pred, y)
    return acc


def feature_comparision(df: pd.DataFrame,
                        d_name: str,
                        label_style: str,
                        )-> None:
    '''
        Test different combinations of the given features.
        :param df: data frame containing all data sets
        :param d_name: string indicating which datasets to use as train/
                             test and validation data set
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from

        :return: None
    '''  
    # add features
    df["DI_Length"] = df.apply(get_dirna_length, axis=1)
    df["Direct_repeat"] = df.apply(get_direct_repeat_length, axis=1)
    df, segment_cols = segment_ohe(df)
    df, junction_start_cols = junction_site_ohe(df, "Start")
    df, junction_end_cols = junction_site_ohe(df, "End")
    junction_cols = junction_start_cols + junction_end_cols
    df["3_5_ratio"] = df.apply(get_3_to_5_ratio, axis=1)
    df["length_proportion"] = df.apply(get_length_proportion, axis=1)
#    df, sequence_cols = full_sequence_ohe(df)

    clf_names = ["logistic_regression", "svc", "random_forest", "mlp", "ada_boost", "naive_bayes"]
    data_dict = dict()
    comb = ["base", "DI_length", "Direct_repeat", "Segment", "Junction", "3_5_ratio", "length_proportion" ,"all"]
    data_dict["param"] = comb

    for clf_name in clf_names:
        print(clf_name)
        data_dict[clf_name] = list()
        clf, _ = select_classifier(clf_name)
        base_features = ["Start", "End"]
        single_cols = ["DI_Length", "Direct_repeat", "3_5_ratio", "length_proportion"]

        for f in comb:
            if f == "base":
                features = base_features
                if clf_name == "random_forest":
                    clf.set_params(max_features=2)
            elif f in single_cols:
                features = base_features + [f]
                if clf_name == "random_forest":
                    clf.set_params(max_features=3)
            elif f == "Segment":
                features = base_features + segment_cols
            elif f == "Junction":
                features = base_features + junction_cols
            elif f == "Sequence":
                features = base_features + sequence_cols
            elif f == "all":
                features = base_features + single_cols + segment_cols + junction_cols
            acc = test_model(df, clf, features, d_name, label_style)
            data_dict[clf_name].append(acc)

    o_df = pd.DataFrame(data_dict)
    print(o_df)
    path = os.path.join(RESULTSPATH, "ML", f"feature_testing_{d_name}.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)

    # Loading the dataset
    df = load_pelz_dataset()["PR8"]
    df["Strain"] = "PR8"
    drop_duplicates = True
    drop_duplicates = False
    label_style = "pd.cut"
    label_style = "median"
    perform_grid_search = False

    test_classifiers(df, "Pelz", label_style, perform_grid_search)
    feature_comparision(df, "Pelz", label_style)

