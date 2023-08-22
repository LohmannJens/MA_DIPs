'''

'''
import os
import sys
import sweetviz

from ydata_profiling import ProfileReport

sys.path.insert(0, "..")
from utils import RESULTSPATH
from ml_utils import load_all_sets, generate_features


if __name__ == "__main__":
    all_df = load_all_sets()

    features = ["DI_Length", "Direct_repeat", "3_5_ratio", "length_proportion", "Inframe_Deletion"]
    df, feature_cols = generate_features(all_df, features, load_precalc=False)

    ### testing ydata_profiling and .describe() ###
    print(df.describe())
    profile = ProfileReport(df, title="test")
    profile.to_file("your_report.html")

    ### testing sweetviz ###
    #report = sweetviz.analyze(all_df)
    #report.show_html()
    df = df[df["dataset_name"] == "Pelz_long"]
    report = sweetviz.compare_intra(df, df["class"] == "gain", ["gain", "not gain"])
    report.show_html()