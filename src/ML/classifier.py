'''
    Tests different classifers on different data sets. Also tests which
    combination of features is the best.
'''
import os
import sys
import shap
import logging
import argparse
import datetime
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
from sklearn.metrics import accuracy_score, RocCurveDisplay, confusion_matrix, make_scorer, precision_score, recall_score, roc_curve, roc_auc_score

from ml_utils import load_all_sets, select_datasets, generate_features
from ml_utils import CHARS, CHARS_COUNT, MAX_LEN

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS


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
            clf = LogisticRegression(penalty="l1", C=1.0, solver="saga", max_iter=10000)
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
            clf = SVC(gamma="scale", C=1.0, kernel="rbf")
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
            clf = RandomForestClassifier(n_estimators=300, max_depth=15, min_samples_split=10, max_features=20)
            param_grid = dict()
    elif clf_name == "mlp":
        if grid_search:
            clf = MLPClassifier(max_iter=10000)
            param_grid = {
                "hidden_layer_sizes": [(50,), (100,), (250,)], 
                "alpha" : [0.001, 0.0001, 0.00001]
            }
        else:
            clf = MLPClassifier(alpha=0.0001, hidden_layer_sizes=(100,), max_iter=10000)
            param_grid = dict()
    elif clf_name == "ada_boost":
        if grid_search:
            clf = AdaBoostClassifier(n_estimators=25, learning_rate=0.1)
            param_grid = {
                "n_estimators": [25, 50, 100,], 
                "learning_rate" : [0.1, 0.5, 1.0],
            }
        else:
            clf = AdaBoostClassifier(n_estimators=100, learning_rate=1.0)
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
                     t_datasets: list,
                     v_datasets: list,
                     folder: str,
                     n_bins: int,
                     label_style: str,
                     y_column: str,
                     perform_grid_search: bool
                     )-> None:
    '''
        Tests three different classifiers on a given dataset.
        :param df: data frame containing all data sets
        :param t_datasets: list of datasets to use for training
        :param v_datasets: list of datasets to use for valdation
        :param folder: folder where to save the plots
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from
        :param perform_grid_search: True if a grid search should be performed

        :return: None
    '''
    # add features
    features = ["Segment", "DI_Length", "Direct_repeat","Junction", "3_5_ratio", "length_proportion", "delta_G", "Peptide_Length"]
    df, feature_cols = generate_features(df, features, load_precalc=True)

    features_out = "\n\t".join(features)
    logging.info(f"\nUsed features:\n\t{features_out}\n#####\n")

    # Selecting train/test and validation data sets
    X, y, X_val, y_val = select_datasets(df, t_datasets, v_datasets, feature_cols, n_bins, label_style, y_column)

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


def test_model(df: pd.DataFrame,
               clf: object,
               f_list: list,
               t_datasets: list,
               v_datasets: list,
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
        :param t_datasets: list of datasets to use for training
        :param v_datasets: list of datasets to use for valdation
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from

        :return: Accuracy of the prediciton
    '''
    X, y, X_val, y_val, = select_datasets(df, t_datasets, v_datasets, f_list, n_bins, label_style, y_column)
    clf.fit(X, y)
    y_pred = clf.predict(X_val)
    acc = accuracy_score(y_pred, y_val)
    return acc


def feature_comparision(df: pd.DataFrame,
                        t_datasets: list,
                        v_datasets: list,
                        folder: str,
                        n_bins: int,
                        label_style: str,
                        y_column: str
                        )-> None:
    '''
        Test different combinations of the given features.
        :param df: data frame containing all data sets
        :param t_datasets: list of datasets to use for training
        :param v_datasets: list of datasets to use for valdation
        :param folder: folder where to save the plots
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from

        :return: None
    '''  
    data_dict = dict()
    comb = ["base", "Segment", "DI_Length", "Direct_repeat", "Junction", "3_5_ratio", "length_proportion" ,"Peptide_Length", "all"]
    data_dict["param"] = comb
    segment_cols = [f"Segment_{s}" for s in SEGMENTS]
    start_cols = [f"Start_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
    end_cols = [f"End_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
    junction_cols = start_cols + end_cols
    sequence_cols = [f"{i}_{ch}" for i in range(1, MAX_LEN+1) for ch in CHARS]

    # add features
    df, _ = generate_features(df, comb, load_precalc=True)

    clf_names = ["logistic_regression", "svc", "random_forest", "mlp", "ada_boost", "naive_bayes"]
    for clf_name in clf_names:
        print(clf_name)
        logging.info(f"\n### {clf_name} ###")
        data_dict[clf_name] = list()
        clf, _ = select_classifier(clf_name)
        base_features = ["Start", "End"]
        single_cols = ["DI_Length", "Direct_repeat", "3_5_ratio", "length_proportion", "Peptide_Length"]

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
            acc = test_model(df, clf, features, t_datasets, v_datasets, n_bins, label_style, y_column)
            data_dict[clf_name].append(acc)

            logging.info(f"\t{f}\t{acc}")

    o_df = pd.DataFrame(data_dict)
    print(o_df)
    path = os.path.join(folder, f"feature_testing_{n_bins}.tex")
    o_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


def run_shap(df: pd.DataFrame,
             t_datasets: list,
             v_datasets: list,
             folder: str,
             n_bins: int,
             label_style: str,
             y_column: str
             )-> None:
    '''
        Calculates shap values for lin. reg. classifier.
        :param df: data frame containing all data sets
        :param t_datasets: list of datasets to use for training
        :param v_datasets: list of datasets to use for valdation
        :param folder: folder where to save the plots
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from

        :return: None
    '''  
    # add features
    features = ["Segment", "DI_Length", "Direct_repeat","Junction", "3_5_ratio", "length_proportion", "delta_G", "Peptide_Length"]
    df, feature_cols = generate_features(df, features, load_precalc=True)

    clf_name = "logistic_regression"
    print(clf_name)
    logging.info(f"\n### {clf_name} ###")

    X, y, X_val, y_val, = select_datasets(df, t_datasets, v_datasets, feature_cols, n_bins, label_style, y_column)
    clf, param_grid = select_classifier(clf_name, grid_search=perform_grid_search)
    clf.fit(X, y)

    X = pd.concat([X, X_val])

    X100 = shap.utils.sample(X, 100)
    explainer = shap.Explainer(clf, X100)
    shap_values = explainer(X)

    fig, axs = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    fig = shap.plots.beeswarm(shap_values, show=False)

    path = os.path.join(folder, f"{clf_name}_shap_beeswarm.png")
    plt.savefig(path)
    plt.close()


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)

    # argument parsing
    p = argparse.ArgumentParser(description="Run classifier to predict competitiveness of DI RNA candidates")
    p.add_argument("--n_bins", "-n", type=int, choices=[2, 3],
                   help="Number of classes to predict")
    p.add_argument("--label_style", "-l", type=str, choices=["pd.cut", "median"],
                   default="median", help="Define how the labels are assigned")
    p.add_argument("--y_column", "-y", type=str, choices=["comb_dup", "int_dup", "Duplicate", "NGS_log_norm"],
                   default="NGS_log_norm", help="Set y column to make prediction for")
    p.add_argument("--train_dataset", "-d", type=str,
                   help="Write which dataset(s) to use for training (split by ','")
    p.add_argument("--validation_dataset", "-v", type=str, default="",
                   help="Write which dataset(s) to use for validation (split by ','). If left empty train data will be split 80-20.")
    p.add_argument("--grid_search", "-g", action="store_true",
                   help="If set a grid search for best model parameters is performed")
    p.add_argument("--mode", "-m", type=str, default="test_classifier",
                   help="Define what analysis to perform: test the classifiers, analyse the features, or shap analysis")
    p.add_argument("--drop_duplicates", "-r", action="store_true",
                   help="If set the duplicates between datasets are removed from the data")
    args = p.parse_args()

    n_bins = args.n_bins
    d = args.train_dataset.split(",")
    if (len(args.validation_dataset) != 0):
        v_d = args.validation_dataset.split(",")
    else:
        v_d = list()
    label_style = args.label_style
    y_column = args.y_column
    perform_grid_search = args.grid_search

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
    df = load_all_sets()
    
    if args.drop_duplicates:
        df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
        df.drop_duplicates("DI", keep=False, inplace=True, ignore_index=True)

    if args.mode == "test_classifier":
        test_classifiers(df, d, v_d, folder, n_bins, label_style, y_column, perform_grid_search)
    elif args.mode == "feature_comparision":
        feature_comparision(df, d, v_d, folder, n_bins, label_style, y_column)
    elif args.mode == "shap":
        run_shap(df, d, v_d, folder, n_bins, label_style, y_column)
        
    logging.info(f"Script ended on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
