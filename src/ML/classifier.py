'''
    Tests three different classifers on different data sets. Also tests which
    combination of features is the best.
'''
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import KFold 
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix

from ml_utils import load_all_sets, select_datasets, segment_ohe, junction_site_ohe, get_dirna_length, get_direct_repeat_length, get_3_to_5_ratio, get_length_proportion

sys.path.insert(0, "..")
from utils import RESULTSPATH


def select_classifier(clf_name: str)-> object:
    '''
        Selects a scikit-learn classifier by a given name. Is implemented in an
        extra function to use the same parameters in each usage of one of the
        classifiers.
        :param clf_name: name of the classifier
        :return: Selected classifier as class implemented in scikit-learn
    '''
    if clf_name == "logistic_regression":
        clf = LogisticRegression(max_iter=4000)
    elif clf_name == "svc":
        clf = SVC(gamma=2, C=1)
    elif clf_name == "random_forest":
        clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    else:
        print(f"classifier {clf_name} unknown!")
        exit()
    return clf


def test_classifiers(df: object, dataset_name: str, n_bins: int, label_style: str)-> None:
    '''
        Tests three different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param dataset_name: string indicating which datasets to use as train/
                             test and validation data set
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes

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

    # Selecting train/test and validation data sets and slice data set to get X and y
    X, y, X_val, y_val = select_datasets(df, dataset_name, feature_cols, n_bins, label_style)

    # Testing different classifiers
    clf_names = ["logistic_regression", "svc", "random_forest"]
    data_dict = dict()
    data_dict["param"] = ["k-fold avg. accurracy", "validation mean"]
    for clf_name in clf_names:
        print(clf_name)
        data_dict[clf_name] = list()
        # setting up classifier and k-fold validation
        kf = KFold(n_splits=5, random_state=None)
        clf = select_classifier(clf_name)

        # test classifier using k-fold validation
        acc_score = list()
        for train_index, test_index in kf.split(X):
            X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
            y_train, y_test = y[train_index], y[test_index]

            clf.fit(X_train,y_train)
            pred_values = clf.predict(X_test)
            acc = accuracy_score(pred_values, y_test)
            acc_score.append(acc)

        # print results of k-fold
        avg_acc_score = sum(acc_score)/len(acc_score)
        print("Avg accuracy k-fold: {}".format(avg_acc_score))
        data_dict[clf_name].append(avg_acc_score)

        # fit on overall model and create confusion matrix for validation set
        clf.fit(X, y)
        if len(X_val.index) != 0:
            predicted_val = clf.predict(X_val)
            acc_score = accuracy_score(predicted_val, y_val)
            confusion_m = confusion_matrix(predicted_val, y_val)
        else:
            predicted = clf.predict(X)
            acc_score = accuracy_score(predicted, y)
            confusion_m = confusion_matrix(predicted, y)
        print(confusion_m)
        data_dict[clf_name].append(acc_score)

        # if two classes given create a ROC
        if len(y.unique()) == 2:
            RocCurveDisplay.from_estimator(clf, X, y)
            plt.plot([0,1], [0,1])

            path = os.path.join(RESULTSPATH, "ML", f"{clf_name}_{dataset_name}_roc_curve.png")
            plt.savefig(path)
            plt.close()

    o_df = pd.DataFrame(data_dict)
    path = os.path.join(RESULTSPATH, "ML", f"{dataset_name}_means.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


def test_model(df: object, clf: object, f_list: list, d_name: str, n_bins: int, label_style: str)-> float:
    '''
        Fits given data to a given classifier class and returns the accuracy.
        Is used to test different feature combinations.
        :param df: pandas data frame
        :param clf: classifier class (from scikit-learn)
        :param f_list: list of features to use for testing
        :param d_name: string indicating which data sets to use
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes

        :return: Accuracy of the prediciton
    '''
    X, y, X_val, y_val, = select_datasets(df, d_name, f_list, n_bins, label_style)
    clf.fit(X, y)
    y_pred = clf.predict(X_val)
    acc = accuracy_score(y_pred, y_val)
    confusion_m = confusion_matrix(y_pred, y_val)
    print(confusion_m)
    return acc


def feature_comparision(df: object, d_name: str, n_bins: int, label_style: str)-> None:
    '''
        Test different combinations of the given features.
        :param df: data frame containing all data sets
        :param d_name: string indicating which datasets to use as train/
                             test and validation data set
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes

        :return: None
    '''  
    # add features
    df["DI_Length"] = df.apply(get_dirna_length, axis=1)
    df["Direct_repeat"] = df.apply(get_direct_repeat_length, axis=1)
    df, segment_cols = segment_ohe(df)
    df, junction_start_cols = junction_site_ohe(df, "Start")
    df, junction_end_cols = junction_site_ohe(df, "End")
    junction_cols = junction_start_cols + junction_end_cols

    clf_names = ["logistic_regression", "svc", "random_forest"]
    data_dict = dict()
    comb = ["base", "DI_length", "Direct_repeat", "Segment", "Junction", "all"]
    data_dict["param"] = comb

    for clf_name in clf_names:
        print(clf_name)
        data_dict[clf_name] = list()
        clf = select_classifier(clf_name)
        base_features = ["Start", "End"]

        for f in comb:
            if f == "base":
                features = base_features
            elif f == "DI_Length":
                features = base_features + [f]
            elif f == "Direct_repeat":
                features = base_features + [f]
            elif f == "Segment":
                features = base_features + segment_cols
            elif f == "Junction":
                features = base_features + junction_cols
            elif f == "all":
                features = base_features + ["DI_Length", "Direct_repeat"] + segment_cols + junction_cols
            acc = test_model(df, clf, features, d_name, n_bins, label_style)
            data_dict[clf_name].append(acc)

    o_df = pd.DataFrame(data_dict)
    print(o_df)
    path = os.path.join(RESULTSPATH, "ML", f"feature_testing_{d_name}_{n_bins}.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)
    # Loading the dataset
    df = load_all_sets()
    drop_duplicates = True
    n_bins = 2
    n_bins = 3
    label_style = "pd.cut"
    label_style = "median"
    datasets = ["Alnaji2019", "PR8"]
    if drop_duplicates:
        df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
        df.drop_duplicates("DI", keep=False, inplace=True, ignore_index=True)

    for d in datasets:
        test_classifiers(df, d, n_bins, label_style)
        feature_comparision(df, d, n_bins, label_style)

