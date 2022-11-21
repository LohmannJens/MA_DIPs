'''

'''

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import KFold 
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, RocCurveDisplay
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder

from ml_utils import load_all_sets

# Loading the dataset
df = load_all_sets()

ohe = OneHotEncoder()
segment_df = pd.DataFrame(ohe.fit_transform(df[["Segment"]]).toarray())
ohe_cols = ohe.get_feature_names().tolist()
segment_df.columns = ohe_cols
df = df.join(segment_df)

# Selecting train/test and validation data sets
train_datasets = ["Alnaji2019_Cal07", "Alnaji2019_NC", "Alnaji2019_Perth"]
val_datasets = ["Alnaji2019_BLEE"]
#train_datasets = ["Pelz", "Alnaji2021"]
#val_datasets = ["Kupke"]

t_df = df.loc[df["dataset_name"].isin(train_datasets)].copy().reset_index()
v_df = df.loc[df["dataset_name"].isin(val_datasets)].copy().reset_index()
X = t_df[["Start", "End"] + ohe_cols]
y = pd.cut(t_df["NGS_log_norm"], bins=2, labels=["low", "high"])
X_val = v_df[["Start", "End"] + ohe_cols]
y_val = pd.cut(v_df["NGS_log_norm"], bins=2, labels=["low", "high"])

# Implementing cross validation
kf = KFold(n_splits=5, random_state=None)
clf = LogisticRegression()
clf = SVC(gamma=2, C=1) # probably overfitting
clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)

acc_score = list()

for train_index, test_index in kf.split(X):
    X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
    y_train, y_test = y[train_index], y[test_index]

    clf.fit(X_train,y_train)
    pred_values = clf.predict(X_test)

    acc = accuracy_score(pred_values, y_test)
    acc_score.append(acc)

    acc_val = accuracy_score(clf.predict(X_val), y_val)

    print(f"accuracy - {acc}")
    print(f"accuracy validation - {acc_val}")

avg_acc_score = sum(acc_score)/len(acc_score)

print("accuracy of each fold - {}".format(acc_score))
print("Avg accuracy : {}".format(avg_acc_score))

clf.fit(X, y)
RocCurveDisplay.from_estimator(clf, X, y)
plt.plot([0,1], [0,1])
plt.show()

