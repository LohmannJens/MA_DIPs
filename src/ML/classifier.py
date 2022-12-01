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
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix

from ml_utils import load_all_sets, select_classifier, select_datasets

sys.path.insert(0, "..")
from utils import RESULTSPATH
from utils import get_sequence, get_seq_len

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import calculate_direct_repeat


def segment_ohe(df: object)-> (object, list):
    '''
        Converts the column with segment names into an one hot encoding.
        :param df: data frame including a row called 'Segment'

        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    ohe = OneHotEncoder()
    segment_df = pd.DataFrame(ohe.fit_transform(df[["Segment"]]).toarray())
    ohe_cols = ohe.get_feature_names_out().tolist()
    segment_df.columns = ohe_cols
    df = df.join(segment_df)
    return df, ohe_cols


def junction_site_ohe(df: object, position: str)-> (object, list):
    '''
        Gets the sequence around the start or end of a given junction site and
        converts the sequence into an one hot encoding.
        :param df: data frame including Start, End, Strain, and Segment
        :param position: is either 'Start' or 'End' to indicate which site

        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # defining static parameters
    CHARS = 'ACGU'
    CHARS_COUNT = len(CHARS)
    n = df.shape[0]
    res = np.zeros((n, CHARS_COUNT * 10), dtype=np.uint8)

    # getting sequence window for each row and do OHE
    for row in df.iterrows():
        r = row[1]
        s = r[position]
        seq = get_sequence(r["Strain"], r["Segment"])
        seq = seq[s-5:s+5]
        i = row[0]
        # Write down OHE in numpy array
        for j, char in enumerate(seq):
            pos = CHARS.rfind(char)
            res[i][j*CHARS_COUNT+pos] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{position}_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
    encoded_df.columns = col_names
    df = df.join(encoded_df)

    return df, col_names


def get_dirna_length(row: list)-> int:
    '''
        Calculates the length of the DI RNA sequence given a row of a data
        frame with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of DI RNA sequence
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    return row["Start"] + (seq_len - row["End"] + 1)


def get_direct_repeat_length(row)-> int:
    '''
        Calculates the length of the direct repeat given a row of a data frame
        with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of direct repeat
    '''
    seq = get_sequence(row["Strain"], row["Segment"])
    s = row["Start"]
    e = row["End"]
    n, _ = calculate_direct_repeat(seq, s, e, 15, 1)
    return n


def test_classifiers(df: object, dataset_name: str, n_bins: int)-> None:
    '''
        Tests three different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param dataset_name: string indicating which datasets to use as train/
                             test and validation data set
        :param n_bins:

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

    # Selecting train/test and validation data sets and slice data set to get X and y
    X, y, X_val, y_val = select_datasets(df, dataset_name, feature_cols, n_bins)

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


def test_model(df: object, clf: object, f_list: list, d_name: str, n_bins: int)-> float:
    '''
        Fits given data to a given classifier class and returns the accuracy.
        Is used to test different feature combinations.
        :param df: pandas data frame
        :param clf: classifier class (from scikit-learn)
        :param f_list: list of features to use for testing
        :param d_name: string indicating which data sets to use

        :return: Accuracy of the prediciton
    '''
    X, y, X_val, y_val, = select_datasets(df, d_name, f_list, n_bins)
    clf.fit(X, y)
    y_pred = clf.predict(X_val)
    acc = accuracy_score(y_pred, y_val)
    confusion_m = confusion_matrix(y_pred, y_val)
    print(confusion_m)
    return acc


def feature_comparision(df: object, d_name: str, n_bins: int)-> None:
    '''
        Test different combinations of the given features.
        :param df: data frame containing all data sets
        :param d_name: string indicating which datasets to use as train/
                             test and validation data set

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
    data_dict["param"] = ["base", "DI_length", "Direct_repeat", "Segment", "junction", "all"]
    for clf_name in clf_names:
        print(clf_name)
        data_dict[clf_name] = list()
        clf = select_classifier(clf_name)

        base_features = ["Start", "End"]
        data_dict[clf_name].append(test_model(df, clf, base_features, d_name, n_bins))

        for f in ["DI_Length", "Direct_repeat", "Segment", "Junction"]:
            if f == "DI_Length":
                data_dict[clf_name].append(test_model(df, clf, base_features + [f], d_name, n_bins))
            elif f == "Direct_repeat":
                data_dict[clf_name].append(test_model(df, clf, base_features + [f], d_name, n_bins))
            elif f == "Segment":
                data_dict[clf_name].append(test_model(df, clf, base_features + segment_cols, d_name, n_bins))
            elif f == "Junction":
                data_dict[clf_name].append(test_model(df, clf, base_features + junction_cols, d_name, n_bins))

        all_features = base_features + ["DI_Length", "Direct_repeat"] + segment_cols + junction_cols
        data_dict[clf_name].append(test_model(df, clf, all_features, d_name, n_bins))

    o_df = pd.DataFrame(data_dict)
    print(o_df)
    path = os.path.join(RESULTSPATH, "ML", f"feature_testing_{d_name}_{n_bins}.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)
    # Loading the dataset
    df = load_all_sets()
    n_bins = 3
    datasets = ["Alnaji2019", "PR8"]
    for d in datasets:
        test_classifiers(df, d, n_bins)
        feature_comparision(df, d, n_bins)

