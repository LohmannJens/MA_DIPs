'''

'''
import os
import sys

from ydata_profiling import ProfileReport

sys.path.insert(0, "..")
from utils import RESULTSPATH
from ml_utils import load_all_sets, generate_features


if __name__ == "__main__":
    all_df = load_all_sets()

    features = ["DI_Length", "Direct_repeat", "3_5_ratio", "length_proportion", "Inframe_Deletion"]
    df, feature_cols = generate_features(all_df, features, load_precalc=False)

    print(df.describe())
    profile = ProfileReport(df, title="test")
    profile.to_file("your_report.html")
