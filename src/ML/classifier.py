'''
    Tests three different classifers on different data sets.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import KFold 
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix

from ml_utils import load_all_sets

sys.path.insert(0, "..")
from utils import RESULTSPATH


def segment_ohe(df: object)-> (object, list):
    '''

    '''
    ohe = OneHotEncoder()
    segment_df = pd.DataFrame(ohe.fit_transform(df[["Segment"]]).toarray())
    ohe_cols = ohe.get_feature_names_out().tolist()
    segment_df.columns = ohe_cols
    df = df.join(segment_df)
    return df, ohe_cols


def sequence_ohe(df: object)-> (object, list):
    '''

    '''   
    #print(df)
    '''
    seq_array = array(list(sequence))
    
    #integer encode the sequence
    label_encoder = LabelEncoder()
    integer_encoded_seq = label_encoder.fit_transform(seq_array)

    #one hot the sequence
    onehot_encoder = OneHotEncoder(sparse=False)
    #reshape because that's what OneHotEncoder likes
    integer_encoded_seq = integer_encoded_seq.reshape(len(integer_encoded_seq), 1)
    onehot_encoded_seq = onehot_encoder.fit_transform(integer_encoded_seq)
    '''
    return df, []


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
    df, sequence_cols = sequence_ohe(df)
    feature_cols = feature_cols + segment_cols + sequence_cols

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
            clf = LogisticRegression(max_iter=1000)
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
        print("Avg accuracy : {}".format(avg_acc_score))

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

