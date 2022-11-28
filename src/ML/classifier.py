'''
    Tests three different classifers on different data sets.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import KFold 
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix

from ml_utils import load_all_sets

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


def get_sequence_length(row: list)-> int:
    '''
        Calculates the length of the DI RNA sequence given a row with the
        necessary data.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of DI RNA sequence
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    return row["Start"] + (seq_len - row["End"] + 1)


def get_direct_repeat_length(row)-> int:
    '''
        Calculates the length of the direct repeat given a row with the
        necessary data.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of direct repeat
    '''
    seq = get_sequence(row["Strain"], row["Segment"])
    s = row["Start"]
    e = row["End"]
    n, _ = calculate_direct_repeat(seq, s, e, 15, 1)
    return n


def test_classifiers(df: object, dataset_name: str)-> None:
    '''
        Tests three different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param dataset_name: string indicating which datasets to use as train/
                             test and validation data set

        :return: None
    '''
    # add segment and sequence information with one hot encoding
    feature_cols = ["Start", "End"]
    df, segment_cols = segment_ohe(df)
    df["DI_Length"] = df.apply(get_sequence_length, axis=1)
    feature_cols.append("DI_Length")
    df["Direct_repeat"] = df.apply(get_direct_repeat_length, axis=1)
    feature_cols.append("Direct_repeat")
    df, junction_start_cols = junction_site_ohe(df, "Start")
    df, junction_end_cols = junction_site_ohe(df, "End")
    feature_cols = feature_cols + segment_cols + junction_start_cols + junction_end_cols

    # Selecting train/test and validation data sets
    if dataset_name == "Alnaji2019":
        train_datasets = ["Alnaji2019_Cal07", "Alnaji2019_NC", "Alnaji2019_Perth"]
        val_datasets = ["Alnaji2019_BLEE"]
    elif dataset_name == "PR8":
        train_datasets = ["Pelz", "Alnaji2021"]
        val_datasets = ["Kupke"]
    else:
        train_datasets = ["Alnaji2019_Cal07", "Alnaji2019_NC", "Alnaji2019_Perth",
                          "Alnaji2019_BLEE", "Pelz", "Alnaji2021", "Kupke"]

    t_df = df.loc[df["dataset_name"].isin(train_datasets)].copy().reset_index()
    X = t_df[feature_cols]
    y = pd.cut(t_df["NGS_log_norm"], bins=3, labels=["low", "mid", "high"], ordered=False)

    if dataset_name in ["Alnaji2019", "PR8"]:
        v_df = df.loc[df["dataset_name"].isin(val_datasets)].copy().reset_index()
        X_val = v_df[feature_cols]
        y_val = pd.cut(v_df["NGS_log_norm"], bins=3, labels=["low", "mid", "high"])

    # Testing different classifiers
    clf_names = ["logistic_regression", "svc", "random_forest"]
    for clf_name in clf_names:
        print(clf_name)
        # setting up classifier and k-fold validation
        kf = KFold(n_splits=5, random_state=None)
        if clf_name == "logistic_regression":
            clf = LogisticRegression(max_iter=1500)
        elif clf_name == "svc":
            clf = SVC(gamma=2, C=1) # probably overfitting
        elif clf_name == "random_forest":
            clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
        else:
            print(f"classifier {clf_name} unknown!")
            return

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

        # fit on overall model and create confusion matrix for validation set
        clf.fit(X, y)
        if dataset_name in ["Alnaji2019", "PR8"]:   
            predicted_val = clf.predict(X_val)
            confusion_m = confusion_matrix(predicted_val, y_val)
            print(accuracy_score(predicted_val, y_val))
            print(confusion_m)

        # if two classes given create a ROC
        if len(y.unique()) == 2:
            RocCurveDisplay.from_estimator(clf, X, y)
            plt.plot([0,1], [0,1])

            path = os.path.join(RESULTSPATH, "ML", f"{clf_name}_{dataset_name}_roc_curve.png")
            plt.savefig(path)
            plt.close()
            

if __name__ == "__main__":
    # Loading the dataset
    df = load_all_sets()
    datasets = ["Alnaji2019", "PR8", "all"]
    for d in datasets:
        test_classifiers(df, d)

