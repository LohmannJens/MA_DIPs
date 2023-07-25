
import sys

import pandas as pd
import numpy as np
#two-way anova
import statsmodels.api as sm
from statsmodels.formula.api import ols
from decimal import Decimal, ROUND_HALF_UP

#binomial test
import scipy.stats as stats

#plotting
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "..")
from utils import get_sequence, load_alnaji_excel, load_short_reads, load_full_alnaji2021, load_pelz_dataset, load_kupke, generate_sampling_data
from utils import SEGMENTS, QUANT, N_SAMPLES

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import count_direct_repeats_overall, include_correction


def create_nucleotide_ratio_matrix(df, col):
    probability_matrix = pd.DataFrame(columns=['A','C','G','U'])
    seq_matrix = df.filter([col], axis=1)
    seq_matrix = seq_matrix[col].str.split('', expand=True)
    # drop first and last column
    seq_matrix = seq_matrix.drop([0, len(seq_matrix.columns)-1], axis=1)
    
    probability_matrix['A'] = seq_matrix.apply(lambda x: dict(x.value_counts()).get('A',0)/len(x), axis=0)
    probability_matrix['C'] = seq_matrix.apply(lambda x: dict(x.value_counts()).get('C',0)/len(x), axis=0)
    probability_matrix['G'] = seq_matrix.apply(lambda x: dict(x.value_counts()).get('G',0)/len(x), axis=0)
    probability_matrix['U'] = seq_matrix.apply(lambda x: dict(x.value_counts()).get('U',0)/len(x), axis=0)
    return probability_matrix

def get_p_value_symbol(p):
    if p < 0.00001:
        return "*****"
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def two_way_anova(dependent_variable_list, factor1_list, factor2_list, dependent_variable_name, factor1_name, factor2_name):
  data = {
    dependent_variable_name: dependent_variable_list,
    factor1_name: factor1_list,
    factor2_name: factor2_list
    }
  df = pd.DataFrame(data)
  # Perform the two-way ANOVA
  model = ols(f'{dependent_variable_name} ~ {factor1_name} + {factor2_name} + {factor1_name}:{factor2_name}', data=df).fit()
  anova_table = sm.stats.anova_lm(model, typ=2)
  return anova_table

def one_way_anova(dependent_variable_list,
                  factor_list,
                  dependent_variable_name,
                  factor_name):
  data = {
    dependent_variable_name: dependent_variable_list,
    factor_name: factor_list
  }
  df = pd.DataFrame(data)
  # perform one-way ANOVA
  f_value, p_value = stats.f_oneway(*df.groupby(factor_name)[dependent_variable_name].apply(list))
  return p_value

def binomial_test_pvalue(ratio1, n_samples, ratio2):
    result = binomial_test_nucleotide_ratios(ratio1, n_samples, ratio2)
    return result.pvalue

def binomial_test_nucleotide_ratios(ratio1, n_samples, ratio2):
    """
    Perform a binomial test to compare the observed number of successes (count) 
    with the expected number of successes based on a specified success probability.

    Parameters:
        ratio1 (float): The observed success ratio.
        n_samples (int): The total number of samples.
        ratio2 (float): The hypothesized success ratio.

    Returns:
        float: The p-value resulting from the binomial test.

    Example:
        ratio1 = 0.7
        n_samples = 100
        ratio2 = 0.5
        p_value = binomial_test_nucleotide_ratios(ratio1, n_samples, ratio2)

    Null Hypothesis:
        The observed success ratio is not significantly different from the hypothesized success ratio.
    """
    count = n_samples * ratio1
    result = stats.binomtest(int(count), n_samples, ratio2)
    return result

def get_deleted_sequence(dip_id, strain):
    """return the sequence of a dip_id

    Args:
        dip_id (str): the id of the dip

    Returns:
        str: the sequence of the dip
    """
    seg, start, end = dip_id.split('_')
    seq = get_sequence(strain, seg)
    return seq[int(start):int(end)-1]

def get_dip_sequence(dip_id, strain):
    seg, start, end = dip_id.split('_')
    fl_seq = get_sequence(strain, seg)
    seq_head = fl_seq[:int(start)]
    seq_foot = fl_seq[int(end)-1:]
    del_length = int(end)-int(start)
    return seq_head + '*'*del_length + seq_foot

def sequence_df(df, strain, isize=5):
    """Generate a DataFrame with sequence information.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the DIP candidates in the 'key' column. Nomenclature: {seg}_{start}_{end}
        isize (int, optional): The size of the sequence before and after the start and end positions. Default is 5.

    Returns:
        pandas.DataFrame: A DataFrame with the following columns:
            - 'key': The original key from the input DataFrame.
            - 'seg': The segment extracted from the key.
            - 'start': The start position extracted from the key.
            - 'end': The end position extracted from the key.
            - 'seq': The dip sequence obtained from the 'key'.
            - 'deleted_sequence': The deleted sequence obtained from the 'key'.
            - 'isize': The specified size for the before and after sequences.
            - 'seq_before_start': The sequence before the start position of length 'isize'.
            - 'seq_after_start': The sequence after the start position of length 'isize'.
            - 'seq_before_end': The sequence before the end position of length 'isize'.
            - 'seq_after_end': The sequence after the end position of length 'isize'.
    """
    res_df = pd.DataFrame(columns=['key','Segment', 'Start','End','seq', 'deleted_sequence', 'isize', 'full_seq', 'Strain', 'seq_before_start', 'seq_after_start',
                                   'seq_before_end', 'seq_after_end', 'seq_around_deletion_junction'])
    for k in df.key:
        seq = get_dip_sequence(k, strain)
        start = int(k.split('_')[1].split('_')[0])
        end = int(k.split('_')[2])
        seg = k.split('_')[0]
        full_seq = get_sequence(strain, seg)
        deleted_seq = get_deleted_sequence(k, strain)
        seq_before_start = seq[start-isize:start]
        seq_after_start = deleted_seq[:isize]
        seq_before_end = deleted_seq[-isize:]
        seq_after_end = seq[end+1:end+isize+1]

        seq_around_deletion_junction = seq_before_start + seq_after_start + seq_before_end + seq_after_end
        res_df = pd.concat([res_df, pd.DataFrame({'key':k, 'Segment':seg, 'Start':start, 'End':end, 'seq':seq, 'isize':isize, 'full_seq': full_seq, 'Strain': strain,
                                'deleted_sequence':deleted_seq, 'seq_before_start':seq_before_start, 'seq_after_start':seq_after_start,
                                'seq_before_end':seq_before_end, 'seq_after_end':seq_after_end, 'seq_around_deletion_junction': seq_around_deletion_junction}, index=[0])], ignore_index=True)
    return res_df

def preprocess(strain, df):
    '''
    
    '''
    df["key"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
    return sequence_df(df, strain)

def generate_expected_data(k, v):
    '''
    
    '''
    for seg in SEGMENTS:
        df = v.loc[v["Segment"] == seg]
        if len(df) == 0:
            continue
        seq = get_sequence(k, seg)        
        s = (int(df.Start.quantile(QUANT)), int(df.Start.quantile(1-QUANT)))
        e = (int(df.End.quantile(QUANT)), int(df.End.quantile(1-QUANT)))
        if "samp_df" in locals():
            temp_df = generate_sampling_data(seq, s, e, N_SAMPLES)
            temp_df["Segment"] = seg
            samp_df = pd.concat([samp_df, temp_df], ignore_index=True)
        else:
            samp_df = generate_sampling_data(seq, s, e, N_SAMPLES)
            samp_df["Segment"] = seg
    
    return samp_df

### create sequence DF
def plot_heatmap(y,x,vals,ax, format='.2f', cmap='coolwarm', vmin=0, vmax=1, cbar=False,cbar_ax=None, cbar_kws=None):
    '''Plot heatmap for values in vals, with x (columns) and y (rows) as labels.'''
    df = pd.DataFrame({'x':x,'y':y,'vals':vals})
    df = pd.pivot_table(df, index='x', columns='y', values='vals', sort=False)
    ax = sns.heatmap(df, fmt=format, annot=True, vmin=vmin, vmax=vmax, ax=ax, cbar=cbar, cmap=cmap, cbar_ax=cbar_ax, cbar_kws=cbar_kws)
    return ax

def plot_nucleotide_ratio_around_deletion_junction_heatmaps(dfs, dfnames, col='seq_around_deletion_junction', 
                                                            height = 14, width = 10, nucleotides=['A','C','G','T']):
    """Plot heatmaps of nucleotide ratios around deletion junctions.

    Args:
        dfs (list of pandas.DataFrame): The list of DataFrames containing the data. 
                                        Each dataframe should be preprocessed with sequence_df(df)
        dfnames (list of str): The names associated with each DataFrame in `dfs`.
        col (str, optional): The column name in the DataFrames that contains the sequence segments of interest. 
                             Default is 'seq_around_deletion_junction'.
        height (float, optional): The height of the figure in inches. Default is 20.
        width (float, optional): The width of the figure in inches. Default is 16.
        nucleotides (list of str, optional): The nucleotides to be plotted. Default is ['A', 'C', 'G', 'T'].

    Returns:
        tuple: A tuple containing the figure and the axes of the subplots.
            - fig (matplotlib.figure.Figure): The generated figure.
            - axs (numpy.ndarray of matplotlib.axes.Axes): The axes of the subplots.

    """
    fig, axs = plt.subplots(figsize=(width,height), nrows=2, ncols=2)
    axs = axs.flatten()
    for i, nuc in enumerate(nucleotides):
        x = []
        y = []
        vals = []
        for dfname, df in zip(dfnames, dfs):
            #df = df.reset_index()
            probability_matrix = create_nucleotide_ratio_matrix(df, col)
            for j in probability_matrix.index:
                x.append(j)
                y.append(dfname)
                vals.append(probability_matrix.loc[j,nuc] * 100)
                
        axs[i] = plot_heatmap(x,y,vals, axs[i], vmin=0, vmax=100, cbar=False, format='.0f')
        axs[i].set_title(f'{nuc} nucleotide ratios around deletion junction [%]')
        axs[i].set_ylabel('')
        axs[i].set_yticks([ytick + 0.5 for ytick in range(len(dfnames))])
        
        axs[i].set_xlabel('position')  
        axs[i].set_xticks([xtick - 0.5 for xtick in probability_matrix.index])
        
        number_of_ticks = len(probability_matrix.index)
        quarter = number_of_ticks // 4
        negative_indexes = [pos for pos in range(-quarter,0)]
        indexes = [pos for pos in range(1,quarter+1)]
        if i % 2 == 0:
            axs[i].set_yticklabels([f'{dfname} ({len(df)})' for dfname,df in zip(dfnames,dfs)])
        else:
            axs[i].set_yticklabels([])
        if i < 2:
            # put x labels at the top
            axs[i].xaxis.set_ticks_position('top')
            axs[i].xaxis.set_label_position('top')
        axs[i].set_xticklabels(negative_indexes + indexes + negative_indexes + indexes)
        xlabels = axs[i].get_xticklabels()
        for x_idx, xlabel in enumerate(xlabels):
            if x_idx < quarter or x_idx >= quarter * 3:
                xlabel.set_color('black')
                xlabel.set_fontweight('bold')
            else:
                xlabel.set_color('grey') 
          
    fig.tight_layout()

    plt.show()

    return fig, axs

def plot_expected_vs_observed_nucleotide_enrichment_heatmaps(dfs, dfnames, expected_dfs , col='seq_around_deletion_junction', height = 14, width = 10, 
                                                 nucleotides=['A','C','G','T']):
    """Plot difference of expected vs observed nucleotide enrichment around deletion junctions as heatmap.

    Args:
        dfs (list of pandas.DataFrame): The list of DataFrames containing the data.
                                        data should be preprocessed with sequence_df(df)
        dfnames (list of str): The names associated with each DataFrame in `dfs`.
        expected_df (pandas.DataFrame): The DataFrame containing the expected data.
                                        data should be preprocessed with sequence_df(df)
        col (str, optional): The column name in the DataFrames that contains the sequences of interest. 
                             Default is 'seq_around_deletion_junction'.
        height (float, optional): The height of the figure in inches. Default is 20.
        width (float, optional): The width of the figure in inches. Default is 16.
        nucleotides (list of str, optional): The nucleotides to be plotted. Default is ['A', 'C', 'G', 'T'].

    Returns:
        tuple: A tuple containing the figure and the axes of the subplots.
            - fig (matplotlib.figure.Figure): The generated figure.
            - axs (numpy.ndarray of matplotlib.axes.Axes): The axes of the subplots.

    """
    fig, axs = plt.subplots(figsize=(width,height), nrows=2, ncols=2)
    axs = axs.flatten()

    for i, nuc in enumerate(nucleotides):
        x = []
        y = []
        vals = []
        val_labels = []
        for dfname, df, expected_df in zip(dfnames, dfs, expected_dfs):
            df = df.reset_index()
            probability_matrix = create_nucleotide_ratio_matrix(df, col)
            n_samples = len(df)
            expected_probability_matrix = create_nucleotide_ratio_matrix(expected_df, col)
            n_samples2 = len(expected_df)
            for j in probability_matrix.index:
                x.append(j)
                y.append(dfname)
                
                p1 = probability_matrix.loc[j,nuc]
                p2 = expected_probability_matrix.loc[j,nuc]
                
                if n_samples < n_samples2:
                    n_s = n_samples
                else:
                    p1, p2 = p2, p1
                    n_s = n_samples2
                pval = binomial_test_pvalue(p1, n_s, p2)
                diff = p1 - p2
                vals.append(diff)
                pval_symbol = get_p_value_symbol(pval)
                val_labels.append(pval_symbol)
        abs_vals = [abs(val) for val in vals]
        axs[i] = plot_heatmap(x,y,vals, axs[i], format='.1e', cbar=True, 
                              vmin = -0.6, vmax=0.6, cbar_kws={"pad": 0.01})
        for v_idx, val_label in enumerate(axs[i].texts):
            val_label.set_text(val_labels[v_idx])
        axs[i].set_title(f'{nuc} ratio difference (observed-expected) around deletion junction')
        axs[i].set_ylabel('')
        axs[i].set_yticks([ytick + 0.5 for ytick in range(len(dfnames))])
        
        axs[i].set_xlabel('position')  
        axs[i].set_xticks([xtick - 0.5 for xtick in probability_matrix.index])
        
        number_of_ticks = len(probability_matrix.index)
        quarter = number_of_ticks // 4
        negative_indexes = [pos for pos in range(-quarter,0)]
        indexes = [pos for pos in range(1,quarter+1)]
        if i % 2 == 0:
            axs[i].set_yticklabels([f'{dfname} ({len(df)})' for dfname,df in zip(dfnames,dfs)])
        else:
            axs[i].set_yticklabels([])
        if i < 2:
            # put x labels at the top
            axs[i].xaxis.set_ticks_position('top')
            axs[i].xaxis.set_label_position('top')
            
        axs[i].set_xticklabels(negative_indexes + indexes + negative_indexes + indexes)
        xlabels = axs[i].get_xticklabels()
        for x_idx, xlabel in enumerate(xlabels):
            if x_idx < quarter or x_idx >= quarter * 3:
                xlabel.set_color('black')
                xlabel.set_fontweight('bold')
            else:
                xlabel.set_color('grey')   
    fig.tight_layout()

    plt.show()

    return fig, axs

def plot_direct_repeat_ratio_heatmaps(dfs, dfnames, height = 14, width = 10, nucleotides=['A','C','G','T']):
    """Plot heatmaps of nucleotide ratios around deletion junctions.

    Args:
        dfs (list of pandas.DataFrame): The list of DataFrames containing the data. 
                                        Each dataframe should be preprocessed with sequence_df(df)
        dfnames (list of str): The names associated with each DataFrame in `dfs`.
        col (str, optional): The column name in the DataFrames that contains the sequence segments of interest. 
                             Default is 'seq_around_deletion_junction'.
        height (float, optional): The height of the figure in inches. Default is 20.
        width (float, optional): The width of the figure in inches. Default is 16.
        nucleotides (list of str, optional): The nucleotides to be plotted. Default is ['A', 'C', 'G', 'T'].

    Returns:
        tuple: A tuple containing the figure and the axes of the subplots.
            - fig (matplotlib.figure.Figure): The generated figure.
            - axs (numpy.ndarray of matplotlib.axes.Axes): The axes of the subplots.

    """
    fig, axs = plt.subplots(figsize=(width,height))

    x = list()
    y = list()
    vals = list()
        
    # calculate direct repeats
    for dfname, df in zip(dfnames, dfs):
        final_d = dict()

        for s in SEGMENTS:
            df_s = df[df["Segment"] == s]
            if len(df_s) == 0:
                continue
            seq = get_sequence(df_s["Strain"].unique()[0], s)
            
            counts, _ = count_direct_repeats_overall(df, seq, mode=1)
            counts = include_correction(counts)

            for k, v in counts.items():
                if k in final_d:
                    final_d[k] += v
                else:
                    final_d[k] = v

        x.extend(final_d.keys())
        y.extend([dfname for _ in range(16)])
        final = np.array(list(final_d.values()))
        vals.extend(final/final.sum())
                        
    axs = plot_heatmap(x,y,vals, axs, vmin=0, vmax=1, cbar=False, format='.5f')
    axs.set_title('direct repeat ratios around deletion junction')
    axs.set_ylabel('')
    axs.set_yticks([ytick + 0.5 for ytick in range(len(dfnames))])
    axs.set_xlabel('direct repeat length')  
          
    fig.tight_layout()
    plt.show()

    return fig, axs

def plot_expected_vs_observed_direct_repeat_heatmaps(dfs, dfnames, expected_dfs , height = 14, width = 10, 
                                                 nucleotides=['A','C','G','T']):
    """Plot difference of expected vs observed nucleotide enrichment around deletion junctions as heatmap.

    Args:
        dfs (list of pandas.DataFrame): The list of DataFrames containing the data.
                                        data should be preprocessed with sequence_df(df)
        dfnames (list of str): The names associated with each DataFrame in `dfs`.
        expected_df (pandas.DataFrame): The DataFrame containing the expected data.
                                        data should be preprocessed with sequence_df(df)
        col (str, optional): The column name in the DataFrames that contains the sequences of interest. 
                             Default is 'seq_around_deletion_junction'.
        height (float, optional): The height of the figure in inches. Default is 20.
        width (float, optional): The width of the figure in inches. Default is 16.
        nucleotides (list of str, optional): The nucleotides to be plotted. Default is ['A', 'C', 'G', 'T'].

    Returns:
        tuple: A tuple containing the figure and the axes of the subplots.
            - fig (matplotlib.figure.Figure): The generated figure.
            - axs (numpy.ndarray of matplotlib.axes.Axes): The axes of the subplots.

    """
    fig, axs = plt.subplots(figsize=(width,height))

    x = list()
    y = list()
    vals = list()
        
    # calculate direct repeats
    for dfname, df, expected_df in zip(dfnames, dfs, expected_dfs):
        final_d = dict()
        expected_final_d = dict()
        pval_labels = list()

        for s in SEGMENTS:
            df_s = df[df["Segment"] == s]
            expected_df_s = expected_df[expected_df["Segment"] == s]
            n_samples = len(df_s)
            if n_samples == 0:
                continue

            seq = get_sequence(df_s["Strain"].unique()[0], s)
            
            counts, _ = count_direct_repeats_overall(df_s, seq, mode=1)
            counts = include_correction(counts)
            for k, v in counts.items():
                if k in final_d:
                    final_d[k] += v
                else:
                    final_d[k] = v

            expected_counts, _ = count_direct_repeats_overall(expected_df_s, seq, mode=1)
            expected_counts = include_correction(expected_counts)
            for k, v in expected_counts.items():
                if k in expected_final_d:
                    expected_final_d[k] += v
                else:
                    expected_final_d[k] = v

        # test statistical significance
      #  f_obs = list()
       # f_exp = list()
        #for a in final_d.keys():
         #   f_obs.extend([a]*int(Decimal(final_d[a]).to_integral_value(rounding=ROUND_HALF_UP)))
          #  f_exp.extend([a]*int(Decimal(expected_final_d[a]).to_integral_value(rounding=ROUND_HALF_UP)))
        #stat, pvalue = stats.ks_2samp(f_obs, f_exp)

        final = np.array(list(final_d.values()))
        expected_final = np.array(list(expected_final_d.values()))
        f_obs = final/final.sum()
        f_exp = expected_final/expected_final.sum()
        stat, pvalue = stats.chisquare(f_obs, f_exp)

        symbol=get_p_value_symbol(pvalue)

        x.extend(final_d.keys())
        y.extend([f"{dfname} {symbol}" for _ in range(16)])
        final = np.array(list(final_d.values()))
        expected_final = np.array(list(expected_final_d.values()))
        vals.extend(final/final.sum() - expected_final/expected_final.sum())

    m = abs(min(vals)) if abs(min(vals)) > max(vals) else max(vals)
    axs = plot_heatmap(x,y,vals, axs, vmin=-m, vmax=m, cbar=True, format='.5f')
    axs.set_title('direct repeat ratio difference (observed-expected)')
    axs.set_ylabel('')
    axs.set_xlabel('direct repeat length')
          
    fig.tight_layout()
    plt.show()

    return fig, axs


if __name__ == "__main__":
    # load all data frames and preprocess with sequence_df(df)
    dfs = list()
    dfnames = list()
    expected_dfs = list()

    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    for k, v in all_reads_dict.items():
        dfs.append(preprocess(k, v))
        dfnames.append(k)
        expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

    df_alnaji = load_full_alnaji2021()
    dfs.append(preprocess("PR8", df_alnaji))
    dfnames.append("Alnaji2021")
    expected_dfs.append(preprocess("PR8", generate_expected_data("PR8", df_alnaji)))

    df_kupke = load_kupke(True)
    for k, v in df_kupke.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Kupke")
        expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

    df_alnaji = load_full_alnaji2021()
    for t in ["3hpi", "6hpi", "14hpi_internal", "14hpi_external", "24hpi"]:
        df = df_alnaji[df_alnaji["Timepoint"] == t].copy()
        dfs.append(preprocess("PR8", df))
        dfnames.append(f"Alnaji2021_{t}")
        expected_dfs.append(preprocess("PR8", generate_expected_data("PR8", df)))      
    
    df_pelz = load_pelz_dataset()
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz")
        expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

    df_pelz = load_pelz_dataset(long_dirna=True)
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz_long")
        expected_dfs.append(preprocess(k, generate_expected_data(k, v)))
    
    df_pelz = load_pelz_dataset(long_dirna=True, de_novo=True)
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz_denovo")
        expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

    plot_nucleotide_ratio_around_deletion_junction_heatmaps(dfs, dfnames, col='seq_around_deletion_junction', nucleotides=['A','C','G','U'])
    plot_expected_vs_observed_nucleotide_enrichment_heatmaps(dfs, dfnames, expected_dfs , col='seq_around_deletion_junction', nucleotides=['A','C','G','U'])

    plot_direct_repeat_ratio_heatmaps(dfs, dfnames, nucleotides=['A','C','G','U'])
    plot_expected_vs_observed_direct_repeat_heatmaps(dfs, dfnames, expected_dfs, nucleotides=['A','C','G','U'])