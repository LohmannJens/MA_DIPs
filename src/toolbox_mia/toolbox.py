
import os
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
from utils import get_sequence, load_alnaji_excel, load_short_reads, load_full_alnaji2021, load_pelz_dataset, load_kupke, generate_sampling_data, load_WSN_data
from utils import SEGMENTS, QUANT, N_SAMPLES, RESULTSPATH

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import count_direct_repeats_overall, include_correction


NUC = dict({"A": "Adenine", "C": "Cytosine", "G": "Guanine", "U": "Uracil"})

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
    return seq_head + '*'*del_length + seq_foot, seq_head, seq_foot

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
        seq, seq_head, seq_foot = get_dip_sequence(k, strain)
        start = int(k.split('_')[1].split('_')[0])
        end = int(k.split('_')[2])
        seg = k.split('_')[0]
        full_seq = get_sequence(strain, seg)
        deleted_seq = get_deleted_sequence(k, strain)
        seq_before_start = seq_head[-isize:]
        seq_after_start = deleted_seq[:isize]
        seq_before_end = deleted_seq[-isize:]
        seq_after_end = seq_foot[:isize]

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
        # skip if there is no range given
        # this would lead to oversampling of a single position
        if s[0] == s[1] or e[0] == e[1]:
            continue
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

def plot_nucleotide_ratio_around_deletion_junction_heatmaps(dfs, dfnames, mode):
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
    height = len(dfs)
    width = 10
    col='seq_around_deletion_junction'
    fig, axs = plt.subplots(figsize=(width,height), nrows=2, ncols=2)
    axs = axs.flatten()
    for i, nuc in enumerate(NUC.keys()):
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
                
        axs[i] = plot_heatmap(x,y,vals, axs[i], vmin=min(vals), vmax=max(vals), cbar=True, format='.0f')
        
        for val_label in axs[i].texts:
            val_label.set_size(8)
        axs[i].set_title(f'{NUC[nuc]}')
        axs[i].set_ylabel('')
        axs[i].set_yticks([ytick + 0.5 for ytick in range(len(dfnames))])
        axs[i].set_xlabel('position')  
        axs[i].set_xticks([xtick - 0.5 for xtick in probability_matrix.index])
        
        number_of_ticks = len(probability_matrix.index)
        quarter = number_of_ticks // 4
        indexes = [pos for pos in range(1,quarter*2+1)]
        if i % 2 == 0:
            axs[i].set_yticklabels([f'{dfname} ({len(df)})' for dfname,df in zip(dfnames,dfs)])
        else:
            axs[i].set_yticklabels([])

        if i < 2:
            axs[i].xaxis.set_ticks_position('top')
            axs[i].xaxis.set_label_position('top')

        axs[i].set_xticklabels(indexes + indexes)
        xlabels = axs[i].get_xticklabels()
        for x_idx, xlabel in enumerate(xlabels):
            if x_idx < quarter or x_idx >= quarter * 3:
                xlabel.set_color('black')
                xlabel.set_fontweight('bold')
            else:
                xlabel.set_color('grey') 
          
    #fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    fig.suptitle("nucleotide ratios around deletion junction [%]")

    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"nuc_occ_{mode}.png")
    plt.savefig(save_path)
    plt.close()

    return fig, axs

def plot_expected_vs_observed_nucleotide_enrichment_heatmaps(dfs, dfnames, expected_dfs, mode):
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
    height = len(dfs)
    width = 10
    col='seq_around_deletion_junction'
    fig, axs = plt.subplots(figsize=(width,height), nrows=2, ncols=2)
    axs = axs.flatten()

    for i, nuc in enumerate(NUC.keys()):
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

                n_s = min(n_s, 1000)
                n_samples2 = min(n_samples2, 1000)

                test_array = np.concatenate((np.ones(int(n_s * p1)), np.zeros(int(n_s - n_s * p1))))
                test_array2 = np.concatenate((np.ones(int(n_samples2 * p2)), np.zeros(int(n_samples2 - n_samples2 * p2))))
                # perform an ANOVA as done in Alaji2021
                pval =  stats.f_oneway(test_array, test_array2).pvalue

                diff = p1 - p2
                vals.append(diff)
                if pval < 0.00001:
                    pval_symbol = "**"
                elif pval < 0.0001:
                    pval_symbol = "*"
                else:
                    pval_symbol = ""
                val_labels.append(pval_symbol)
                
        m = abs(min(vals)) if abs(min(vals)) > max(vals) else max(vals)
        axs[i] = plot_heatmap(x,y,vals, axs[i], format='.1e', cbar=True, 
                              vmin = -m, vmax=m, cbar_kws={"pad": 0.01})
        for v_idx, val_label in enumerate(axs[i].texts):
            val_label.set_text(val_labels[v_idx])
            val_label.set_size(6)
        axs[i].set_title(f'{NUC[nuc]}')
        axs[i].set_ylabel('')
        axs[i].set_yticks([ytick + 0.5 for ytick in range(len(dfnames))])
        axs[i].set_xlabel('position')  
        axs[i].set_xticks([xtick - 0.5 for xtick in probability_matrix.index])
        
        number_of_ticks = len(probability_matrix.index)
        quarter = number_of_ticks // 4
        indexes = [pos for pos in range(1,quarter*2+1)]
        if i % 2 == 0:
            axs[i].set_yticklabels([f'{dfname} ({len(df)})' for dfname,df in zip(dfnames,dfs)])
        else:
            axs[i].set_yticklabels([])

        if i < 2:
            axs[i].xaxis.set_ticks_position('top')
            axs[i].xaxis.set_label_position('top')
            
        axs[i].set_xticklabels(indexes + indexes)
        xlabels = axs[i].get_xticklabels()
        for x_idx, xlabel in enumerate(xlabels):
            if x_idx < quarter or x_idx >= quarter * 3:
                xlabel.set_color('black')
                xlabel.set_fontweight('bold')
            else:
                xlabel.set_color('grey')   

    fig.subplots_adjust(top=0.9)
    fig.suptitle("ratio difference (observed-expected) around deletion junction")

    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"nuc_occ_diff_{mode}.png")
    plt.savefig(save_path)
    plt.close()

    return fig, axs

def plot_direct_repeat_ratio_heatmaps(dfs, dfnames, mode):
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
    height = len(dfs)/2
    width = 10
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
            if dfname == "Boussier":
                counts, _ = count_direct_repeats_overall(df_s, seq, mode=3)
            else:
                counts, _ = count_direct_repeats_overall(df_s, seq, mode=1)
            
            if dfname not in ["Mendes", "Boussier"]:
                counts = include_correction(counts)
            
            for k, v in counts.items():
                if k in final_d:
                    final_d[k] += v
                else:
                    final_d[k] = v

        x.extend(final_d.keys())
        y.extend([f"{dfname} ({len(df)})" for _ in range(6)])
        final = np.array(list(final_d.values()))
        vals.extend(final/final.sum())

    axs = plot_heatmap(x,y,vals, axs, vmin=0, vmax=1, cbar=True, format='.5f')
    axs.set_title('direct repeat ratios around deletion junction')
    axs.set_ylabel('')
    axs.set_xlabel('direct repeat length')

    x_ticks = axs.get_xticklabels()
    label = x_ticks[-2].get_text()
    x_ticks[-1].set_text(f"> {label}")
    axs.set_xticklabels(x_ticks)
    fig.tight_layout()

    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"direct_repeats_{mode}.png")
    plt.savefig(save_path)
    plt.close()

    return fig, axs

def plot_expected_vs_observed_direct_repeat_heatmaps(dfs, dfnames, expected_dfs , mode):
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
    height = len(dfs)/2
    width = 10
    fig, axs = plt.subplots(figsize=(width,height))

    x = list()
    y = list()
    vals = list()
    val_labels = list()
        
    # calculate direct repeats
    for dfname, df, expected_df in zip(dfnames, dfs, expected_dfs):
        final_d = dict()
        expected_final_d = dict()

        for s in SEGMENTS:
            df_s = df[df["Segment"] == s]
            expected_df_s = expected_df[expected_df["Segment"] == s]
            n_samples = len(df_s)
            if n_samples == 0:
                continue

            seq = get_sequence(df_s["Strain"].unique()[0], s)
            if dfname == "Boussier":
                counts, _ = count_direct_repeats_overall(df_s, seq, mode=3)
            else:
                counts, _ = count_direct_repeats_overall(df_s, seq, mode=1)
            
            if dfname not in ["Mendes", "Boussier"]:
                counts = include_correction(counts)
            for k, v in counts.items():
                if k in final_d:
                    final_d[k] += v
                else:
                    final_d[k] = v

            expected_counts, _ = count_direct_repeats_overall(expected_df_s, seq, mode=1)
            for k, v in expected_counts.items():
                if k in expected_final_d:
                    expected_final_d[k] += v
                else:
                    expected_final_d[k] = v

        final = np.array(list(final_d.values()))
        expected_final = np.array(list(expected_final_d.values()))
        f_obs = final/final.sum()
        f_exp = expected_final/expected_final.sum()
        stat, pvalue = stats.chisquare(f_obs, f_exp)

        symbol = get_p_value_symbol(pvalue)
        x.extend(final_d.keys())
        y.extend([f"{dfname} ({len(df)}) {symbol}" for _ in range(6)])
        vals.extend(final/final.sum() - expected_final/expected_final.sum())

        for f_ob, f_ex, n_samples in zip(f_obs, f_exp, final_d.values()):
            if n_samples == 0:
                pval_symbol = ""
            else:
                pvalue = binomial_test_pvalue(f_ob, int(n_samples), f_ex)
                if pvalue < 0.00001:
                    pval_symbol = "**"
                elif pvalue < 0.0001:
                    pval_symbol = "*"
                else:
                    pval_symbol = ""
            val_labels.append(pval_symbol)

    m = abs(min(vals)) if abs(min(vals)) > max(vals) else max(vals)
    axs = plot_heatmap(x,y,vals, axs, vmin=-m, vmax=m, cbar=True, format='.5f')
    axs.set_title('direct repeat ratio difference (observed-expected)')
    axs.set_ylabel('')
    axs.set_xlabel('direct repeat length')

    for v_idx, val_label in enumerate(axs.texts):
        val_label.set_text(f"{val_label.get_text()}\n{val_labels[v_idx]}")
        val_label.set_size(6)

    x_ticks = axs.get_xticklabels()
    label = x_ticks[-2].get_text()
    x_ticks[-1].set_text(f"> {label}")
    axs.set_xticklabels(x_ticks)
    fig.tight_layout()

    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"direct_repeats_diff_{mode}.png")
    plt.savefig(save_path)
    plt.close()

    return fig, axs


def plot_nucleotide_enrichment_over_time(dfs, dfnames):
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
    height = 10
    width = len(dfs) / 2
    col='seq_around_deletion_junction'
    fig, axs = plt.subplots(figsize=(width,height), nrows=2, ncols=2)
    axs = axs.flatten()
    for i, nuc in enumerate(NUC.keys()):

        relevant_indices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ,16 ,17 ,18, 19, 20]
        relevant_indices = [4, 5, 6, 7, 14, 15, 16, 17]

        x = []
        y = dict()
        for j in relevant_indices:
            y[j] = list()

        for df in dfs:
            probability_matrix = create_nucleotide_ratio_matrix(df, col)
            for j in relevant_indices:
                y[j].append(probability_matrix.loc[j,nuc])

        #x = range(len(dfnames))
        x = [0.00, 0.50, 0.99, 1.40, 3.46, 4.00, 4.47, 5.00, 5.48, 7.95, 8.96, 9.42, 12.43, 12.97, 13.50, 16.01, 16.97, 17.45, 18.00, 19.48, 19.99, 20.44, 21.00, 22.00, 26.44, 29.95, 36.42, 42.42]
        for k, v in y.items():
            # map relevant index to start/end and exact position
            if k <= 10:
                pos = "Start"
            else:
                pos = "End"
                k = k - 10
            label = f"{pos} pos {k}"
            axs[i].plot(x, v, label=label)

        axs[i].set_title(f'{nuc} nucleotide ratios around deletion junction')
        axs[i].set_ylabel('relative occurrence')
        axs[i].set_xlabel('time')

        axs[i].axvline(21.5)
        #axs[i].set_xticklabels(dfnames)
          
    fig.tight_layout()
    axs[i].legend()

    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"nuc_occ_ot.png")
    plt.savefig(save_path)

    return fig, axs


def plot_direct_repeats_over_time(dfs, dfnames):
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
    height = 6
    width = len(dfs) / 2
    fig, axs = plt.subplots(figsize=(width,height), nrows=1, ncols=1)
   # cm = plt.get_cmap('tab10')
  #  axs.set_prop_cycle('color', [cm(1.*i/8) for i in range(8)])

    y = dict({i: list() for i in range(6)})
    for df in dfs:
        all = 0
        co = dict()
        for s in SEGMENTS:
            df_s = df[df["Segment"] == s]
            if len(df_s) == 0:
                continue
                
            seq = get_sequence(df_s["Strain"].unique()[0], s)
            counts, _ = count_direct_repeats_overall(df_s, seq, mode=1)
            counts = include_correction(counts)
            for k, v in counts.items():
                if k in co:
                    co[k] += v
                else:
                    co[k] = v
                all += v

        for i in range(6):
            y[i].append(co[i] / all)
   
    
    bar_width = 0.5
    x = [0.00, 0.50, 0.99, 1.40, 3.46, 4.00, 4.47, 5.00, 5.48, 7.95, 8.96, 9.42, 12.43, 12.97, 13.50, 16.01, 16.97, 17.45, 18.00, 19.48, 19.99, 20.44, 21.00, 22.00, 26.44, 29.95, 36.42, 42.42]
    bottom = np.zeros(len(dfnames))

    for i in range(6):
        axs.bar(x, y[i], bar_width, label=i, bottom=bottom)
        bottom += y[i]

    axs.set_title(f'direct repeat ratios over time')
    axs.set_ylabel('relative occurrence')
    axs.set_xlabel('time')
    axs.axvline(21.5)
   
    fig.tight_layout()
    axs.legend()

    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"direct_repeats_ot.png")
    plt.savefig(save_path)

    return fig, axs


def plot_distribution_over_segments(dfs, dfnames, mode):
    """

    Args:
        dfs (list of pandas.DataFrame): The list of DataFrames containing the data. 
                                        Each dataframe should be preprocessed with sequence_df(df)
        dfnames (list of str): The names associated with each DataFrame in `dfs`.
        col (str, optional): The column name in the DataFrames that contains the sequence segments of interest. 
                             Default is 'seq_around_deletion_junction'.
    Returns:
        tuple: A tuple containing the figure and the axes of the subplots.
            - fig (matplotlib.figure.Figure): The generated figure.
            - axs (numpy.ndarray of matplotlib.axes.Axes): The axes of the subplots.

    """
    height = 6
    width = len(dfs) * 1.5
    fig, ax = plt.subplots(figsize=(width,height), nrows=1, ncols=1)
    cm = plt.get_cmap('tab10')
    ax.set_prop_cycle('color', [cm(1.*i/8) for i in range(8)])

    fractions_dict = dict({s: list() for s in SEGMENTS})

    for df in dfs:
        counts = df['Segment'].value_counts()
        fractions = counts / len(df)

        for s in SEGMENTS:
            if s in fractions:
                frac = fractions[s]
            else:
                frac = 0.0
            fractions_dict[s].append(frac)

    bar_width = 0.5
    x = np.arange(len(dfnames))
    bottom = np.zeros(len(dfnames))

    for s in SEGMENTS:
        ax.bar(x, fractions_dict[s], bar_width, label=s, bottom=bottom)
        bottom += fractions_dict[s]

    # Add labels, title, and legend
    ax.set_xlabel('Dataset')
    ax.set_ylabel('Relative occurrence of DIs')
    ax.set_title('Fractions of DIs from different datasets')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{dfname} ({len(df)})' for dfname,df in zip(dfnames,dfs)])
    ax.legend()

    plt.tight_layout()
    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"fraction_segments_{mode}.png")
    plt.savefig(save_path)

    return fig, ax

def calculate_deletion_shifts(dfs, dfnames, mode):
    """

    Args:
        dfs (list of pandas.DataFrame): The list of DataFrames containing the data. 
                                        Each dataframe should be preprocessed with sequence_df(df)
        dfnames (list of str): The names associated with each DataFrame in `dfs`.
        col (str, optional): The column name in the DataFrames that contains the sequence segments of interest. 
                             Default is 'seq_around_deletion_junction'.
    Returns:
        tuple: A tuple containing the figure and the axes of the subplots.
            - fig (matplotlib.figure.Figure): The generated figure.
            - axs (numpy.ndarray of matplotlib.axes.Axes): The axes of the subplots.

    """
    height = 6
    width = len(dfs) * 1.5
    fig, ax = plt.subplots(figsize=(width,height), nrows=1, ncols=1)
    cm = plt.get_cmap('tab10')
    ax.set_prop_cycle('color', [cm(1.*i/8) for i in range(8)])

    fractions_dict = dict({n: list() for n in [0, 1, 2]})

    for df in dfs:
        df["length"] = df["deleted_sequence"].apply(len)
        df["shift"] = df["length"] % 3
        
        shifts = df["shift"].value_counts()
        shifts = shifts / len(df)

        for n in [0, 1, 2]:
            fractions_dict[n].append(shifts[n])

    bar_width = 0.5
    x = np.arange(len(dfnames))
    bottom = np.zeros(len(dfnames))

    for n in [0, 1, 2]:
        ax.bar(x, fractions_dict[n], bar_width, label=n, bottom=bottom)
        bottom += fractions_dict[n]

    # Add labels, title, and legend
    ax.set_xlabel('Dataset')
    ax.set_ylabel('Relative occurrence')
    ax.set_title('Fractions of DIs with different offsets to in frame deletion')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{dfname} ({len(df)})' for dfname,df in zip(dfnames,dfs)])
    ax.legend()

    plt.tight_layout()
    save_path = os.path.join(RESULTSPATH, "toolbox_mia", f"deletion_shifts_{mode}.png")
    plt.savefig(save_path)
    
    return fig, ax


if __name__ == "__main__":
    plt.style.use('seaborn')

    # load all data frames and preprocess with sequence_df(df)
    dfs = list()
    dfnames = list()
    expected_dfs = list()

    mode = "pelz_ot"
    #mode = "all"
    #mode = "alnaji_split"
    #mode = "no_pelz"

    #mode = "test"
#    cleaned_data_dict = load_alnaji_excel()
 #   all_reads_dict = load_short_reads(cleaned_data_dict)
  #  for k, v in all_reads_dict.items():
   #     dfs.append(preprocess(k, v))
    #    dfnames.append(k)
     #   expected_dfs.append(preprocess(k, generate_expected_data(k, v)))



    if mode == "pelz_ot":
        df_pelz = load_pelz_dataset(by_time=True)
        for k, v in df_pelz.items():
            for t in ["VB3-Saat","VB3-7","VB3-8","VB3-9","VB3-13","VB3-14","VB3-15","VB3-16","VB3-17","VB3-22","VB3-24","VB3-25","VB3-31","VB3-32","VB3-33","VB3-38","VB3-40","VB3-41","VB3-42","VB3-45","VB3-46","VB3-47","VB3-48"]:
                df = v[v[t] != 0].copy()
                dfs.append(preprocess(k, df))
                dfnames.append(f"Pelz_{t}")
        
        df_pelz = load_pelz_dataset(follow_up=True)
        for k, v in df_pelz.items():
            for t in ["A1", "A2", "A3", "A4", "A5"]:
                df = v[v[t] != 0].copy()
                dfs.append(preprocess(k, df))
                dfnames.append(f"Pelz_{t}")
        plot_nucleotide_enrichment_over_time(dfs, dfnames)
        plot_direct_repeats_over_time(dfs, dfnames)
        exit()

    elif mode == "no_pelz":
        cleaned_data_dict = load_alnaji_excel()
        all_reads_dict = load_short_reads(cleaned_data_dict)
        for k, v in all_reads_dict.items():
            dfs.append(preprocess(k, v))
            dfnames.append(k)
            expected_dfs.append(preprocess(k, generate_expected_data(k, v)))
    
        alnaji_dict = load_full_alnaji2021()
        for k, v in alnaji_dict.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Alnaji2021")
            expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

        df_mendes = load_WSN_data("Mendes")
        for k, v in df_mendes.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Mendes")
            expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

        df_boussier = load_WSN_data("Boussier")
        for k, v in df_boussier.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Boussier")
            expected_dfs.append(preprocess(k, generate_expected_data(k, v)))
    
    elif mode == "alnaji_split":
        df_alnaji = load_full_alnaji2021()["PR8"]
        for t in ["3hpi", "6hpi", "14hpi_internal", "14hpi_external", "24hpi"]:
            df = df_alnaji[df_alnaji["Timepoint"] == t].copy()
            dfs.append(preprocess("PR8", df))
            dfnames.append(f"Alnaji2021_{t}")
            expected_dfs.append(preprocess("PR8", generate_expected_data("PR8", df)))

    elif mode == "all":
        cleaned_data_dict = load_alnaji_excel()
        all_reads_dict = load_short_reads(cleaned_data_dict)
        for k, v in all_reads_dict.items():
            dfs.append(preprocess(k, v))
            dfnames.append(k)
            #expected_dfs.append(preprocess(k, generate_expected_data(k, v)))
    
        alnaji_dict = load_full_alnaji2021()
        for k, v in alnaji_dict.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Alnaji2021")
            #expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

        df_pelz = load_pelz_dataset(long_dirna=True)
        for k, v in df_pelz.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Pelz_long")
            #expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

        df_mendes = load_WSN_data("Mendes")
        for k, v in df_mendes.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Mendes")
            #expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

        df_boussier = load_WSN_data("Boussier")
        for k, v in df_boussier.items():
            dfs.append(preprocess(k, v))
            dfnames.append("Boussier")
            #expected_dfs.append(preprocess(k, generate_expected_data(k, v)))

#    plot_nucleotide_ratio_around_deletion_junction_heatmaps(dfs, dfnames, mode)
 #   plot_expected_vs_observed_nucleotide_enrichment_heatmaps(dfs, dfnames, expected_dfs, mode)

  #  plot_direct_repeat_ratio_heatmaps(dfs, dfnames, mode)
   # plot_expected_vs_observed_direct_repeat_heatmaps(dfs, dfnames, expected_dfs, mode)

    plot_distribution_over_segments(dfs, dfnames, mode)
    calculate_deletion_shifts(dfs, dfnames, mode)