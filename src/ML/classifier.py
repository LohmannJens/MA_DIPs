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
from sklearn.metrics import accuracy_score, RocCurveDisplay

from ml_utils import load_all_sets

sys.path.insert(0, "..")
from utils import RESULTSPATH


def test_classifiers(df: object, dataset_name: str)-> None:
    '''
        Tests three different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param dataset_name: string indicating which datasets to use as train/
                             test and validation data set

        :return: None
    '''
    # add segment information with one hot encoding
    ohe = OneHotEncoder()
    segment_df = pd.DataFrame(ohe.fit_transform(df[["Segment"]]).toarray())
    ohe_cols = ohe.get_feature_names().tolist()
    segment_df.columns = ohe_cols
    df = df.join(segment_df)

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
    X = t_df[["Start", "End"] + ohe_cols]
    y = pd.cut(t_df["NGS_log_norm"], bins=2, labels=["low", "high"])
    if dataset_name in ["Alnaji2019", "PR8"]:
        v_df = df.loc[df["dataset_name"].isin(val_datasets)].copy().reset_index()
        X_val = v_df[["Start", "End"] + ohe_cols]
        y_val = pd.cut(v_df["NGS_log_norm"], bins=2, labels=["low", "high"])

    # Testing different classifiers
    clf_names = ["logistic_regression", "svc", "random_forest"]
    for clf_name in clf_names:
        # setting up classifier and k-fold validation
        kf = KFold(n_splits=5, random_state=None)
        if clf_name == "logistic_regression":
            clf = LogisticRegression()
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

            print(f"accuracy - {acc}")

            if dataset_name in ["Alnaji2019", "PR8"]:
                acc_val = accuracy_score(clf.predict(X_val), y_val)
                print(f"accuracy validation - {acc_val}")

        # print results
        avg_acc_score = sum(acc_score)/len(acc_score)
        print("accuracy of each fold - {}".format(acc_score))
        print("Avg accuracy : {}".format(avg_acc_score))

        # fit on overall model and plot ROC curve
        clf.fit(X, y)
        RocCurveDisplay.from_estimator(clf, X, y)
        plt.plot([0,1], [0,1])

        path = os.path.join(RESULTSPATH, "ML", f"{clf_name}_{datasets}_roc_curve.png")
        plt.savefig(path)
        plt.close()

if __name__ == "__main__":
    # Loading the dataset
    df = load_all_sets()
    datasets = ["Alnaji2019", "PR8", "all"]
    for d in datasets:
        test_classifiers(df, d)

