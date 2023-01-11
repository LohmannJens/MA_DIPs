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

from ml_utils import load_all_sets, select_datasets, segment_ohe, junction_site_ohe, get_dirna_length, get_direct_repeat_length, get_3_to_5_ratio, get_length_proportion, full_sequence_ohe

sys.path.insert(0, "..")
from utils import RESULTSPATH


def select_classifier(clf_name: str,
                      grid_search: bool=False
                      )-> (object, dict):
    '''
        Selects a scikit-learn classifier by a given name. Is implemented in an
        extra function to use the same parameters in each usage of one of the
        classifiers.
        :param clf_name: name of the classifier
        :param grid_search: Bool indicating if a grid search will be performed

        :return: 1. Selected classifier as class implemented in scikit-learn
                 2. parameter grid, if grid search is True, else empty dict
    '''
    if clf_name == "logistic_regression":
        if grid_search:
            clf = LogisticRegression(max_iter=10000, solver="saga", l1_ratio=0.5)
            param_grid = {
                "penalty": ["l1", "l2", "elasticnet"], 
                "C" : [0.01, 0.1, 1.0],
            }
        else:
            clf = LogisticRegression(penalty="l1", C=0.1, solver="saga", max_iter=10000)
            param_grid = dict()
    elif clf_name == "svc":
        if grid_search:
            clf = SVC(gamma=2, C=1)
            param_grid = {
                "C" : [0.01, 0.1, 1.0],
                "kernel" : ["linear", "rbf", "sigmoid"],
                "gamma" : ["scale", "auto", 2],
            }
        else:
            clf = SVC(gamma="auto", C=0.01, kernel="rbf") # alnaji, acc=0.36
            clf = SVC(gamma="scale", C=1.0, kernel="sigmoid") # PR8, acc=0.39
            param_grid = dict()
    elif clf_name == "random_forest":
        if grid_search:
            clf = RandomForestClassifier()
            param_grid = {
                "min_samples_split": [3, 5, 10], 
                "n_estimators" : [100, 300],
                "max_depth": [3, 5, 15, 25],
                "max_features": [3, 5, 10, 20]
            }
        else:
            clf = RandomForestClassifier(n_estimators=300, max_depth=3, min_samples_split=10, max_features=5)
            param_grid = dict()
    elif clf_name == "mlp":
        if grid_search:
            clf = MLPClassifier()
            param_grid = {
                "hidden_layer_sizes": [(50,), (100,), (250,)], 
                "alpha" : [0.001, 0.0001, 0.00001]
            }
        else:
            clf = MLPClassifier(alpha=0.001, hidden_layer_sizes=(100,), max_iter=1000)
            param_grid = dict()
    elif clf_name == "ada_boost":
        clf = AdaBoostClassifier(n_estimators=25, learning_rate=0.1)
        if grid_search:
            param_grid = {
                "n_estimators": [25, 50, 100,], 
                "learning_rate" : [0.1, 0.5, 1.0],
            }
        else:
            param_grid = dict()
    elif clf_name == "naive_bayes":
        clf = GaussianNB(var_smoothing=0.0000000001)
        if grid_search:
            param_grid = {
                "var_smoothing": [0.000000001, 0.0000000001, 0.0000000001]
            }
        else:
            param_grid = dict()
    else:
        print(f"classifier {clf_name} unknown!")
        exit()
    return clf, param_grid


def test_classifiers(df: pd.DataFrame,
                     dataset_name: str,
                     n_bins: int,
                     label_style: str,
                     y_column: str,
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
    X, y, X_val, y_val = select_datasets(df, dataset_name, feature_cols, n_bins, label_style, y_column)

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
        predicted_val = grid_search.predict(X_val)
        acc_score = accuracy_score(predicted_val, y_val)
        confusion_m = confusion_matrix(predicted_val, y_val)

        print(acc_score)
        print(confusion_m)
        data_dict[clf_name].append(acc_score)

        # if two classes given create a ROC
        if len(y.unique()) == 2:
            RocCurveDisplay.from_estimator(grid_search, X, y)
            plt.plot([0,1], [0,1])

            path = os.path.join(RESULTSPATH, "ML", f"{clf_name}_{dataset_name}_roc_curve.png")
            plt.savefig(path)
            plt.close()

    o_df = pd.DataFrame(data_dict)
    path = os.path.join(RESULTSPATH, "ML", f"{dataset_name}_means.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


def test_model(df: pd.DataFrame,
               clf: object,
               f_list: list,
               d_name: str,
               n_bins: int,
               label_style: str,
               y_column: str
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
    X, y, X_val, y_val, = select_datasets(df, d_name, f_list, n_bins, label_style, y_column)
    clf.fit(X, y)
    y_pred = clf.predict(X_val)
    acc = accuracy_score(y_pred, y_val)
    confusion_m = confusion_matrix(y_pred, y_val)
    return acc


def feature_comparision(df: pd.DataFrame,
                        d_name: str,
                        n_bins: int,
                        label_style: str,
                        y_column: str
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
            acc = test_model(df, clf, features, d_name, n_bins, label_style, y_column)
            data_dict[clf_name].append(acc)

    o_df = pd.DataFrame(data_dict)
    print(o_df)
    path = os.path.join(RESULTSPATH, "ML", f"feature_testing_{d_name}_{n_bins}.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)

# TODO: write argparse
    # drop_duplicates
    # n_bins
    # label_style
    # y_column
    # dataset
    # grid_search
    # test_classifiers <-> feature_comparision

    # Loading the dataset
    df = load_all_sets()
    drop_duplicates = True
    drop_duplicates = False
    n_bins = 2
 #   n_bins = 3
    label_style = "pd.cut"
    label_style = "median"
    datasets = ["Alnaji2019", "PR8"]
    datasets = ["Alnaji2019"]
    y_column = "comb_dup"
    y_column = "int_dup"
    y_column = "Duplicate"
    y_column = "NGS_log_norm"
    perform_grid_search = False

    if drop_duplicates:
        df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
        df.drop_duplicates("DI", keep=False, inplace=True, ignore_index=True)

    for d in datasets:
        print(f"#### {d} ####")
        test_classifiers(df, d, n_bins, label_style, y_column, perform_grid_search)
        feature_comparision(df, d, n_bins, label_style, y_column)

