''' 
Module: stat.py
Author: Marc Zepeda
Created: 2024-08-26
Description: Statistics

Usage:
[Statistics]
- describe(): returns descriptive statistics for numerical columns in a DataFrame
- difference(): computes the appropriate statistical test(s) and returns the p-value(s)
- correlation(): returns a correlation matrix
- weighted_correlation(): computes the weighted correlation

[Statistics & Plotting]
- corr_line(): compute linear regression line and add it to the scatter plot 
- weighted_corr_line(): compute weighted linear regression line and add it to the scatter plot

[Comparison]
- compare(): computes FC, pval, and log transformations relative to a specified condition
- odds_ratio(): computes odds ratio relative to a specified condition (OR = (A/B)/(C/D))
- zscore(): Z-score `FC` within each `cond_col` using stats computed from rows where `var_col == var`.
'''

# Import packages
import itertools
import pandas as pd
import numpy as np
from scipy.stats import skew, kurtosis, ttest_ind, ttest_rel, f_oneway, ttest_ind, ttest_rel, mannwhitneyu, wilcoxon, rankdata, fisher_exact
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests
from . import io, tidy as t, plot as p

# Statistics
def describe(df: pd.DataFrame | str, cols:list=[], group:str='', dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    describe(): returns descriptive statistics for numerical columns in a DataFrame
    
    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    cols (list, optional): list of numerical columns to compute statistics
    group (str, optional): column name to split tidy dataframes
    dir (str, optional): save directory
    file (str, optional): save file
    
    Dependencies: pandas, numpy, scipy.stats, & io
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    if group!='': df = df.pivot(columns=group) # Splits tidy dataframe
    if len(cols)>0: df = df[cols] # Isolate specified columns

    # Compute descriptive statistics
    descriptive = pd.DataFrame()
    if group !='': descriptive[group] = [multindex for multindex in df.mean().keys()] # Group
    descriptive['mean'] = df.mean().reset_index(drop=True) # Mean
    descriptive['median'] = df.median().reset_index(drop=True) # Median
    descriptive['variance'] = df.var().reset_index(drop=True) # Variance
    descriptive['std_dev'] = df.std().reset_index(drop=True) # Standard Deviation
    descriptive['mad'] = df.apply(lambda x: np.median(np.abs(x - x.median()))).reset_index(drop=True) # Median Absolute Deviation
    descriptive['min'] = df.min().reset_index(drop=True) # Minimum
    descriptive['max'] = df.max().reset_index(drop=True) # Maximum
    descriptive['range'] = [ma - mi for ma,mi in t.zip_cols(descriptive,['max','min'])] # Range
    descriptive['skewness'] = df.apply(lambda x: skew(x, nan_policy='omit')).reset_index(drop=True) # Skewness
    descriptive['kurtosis'] = df.apply(lambda x: kurtosis(x, nan_policy='omit')).reset_index(drop=True) # Kurtosis
    descriptive['count'] = df.count().reset_index(drop=True) # Count (non-missing observations)
    descriptive['sum'] = df.sum().reset_index(drop=True) # Sum
    descriptive['25%'] = df.quantile(0.25).reset_index(drop=True)  # Quantiles (25%, 50%, 75%)
    descriptive['50%'] = df.quantile(0.50).reset_index(drop=True)
    descriptive['75%'] = df.quantile(0.75).reset_index(drop=True)

    # Save & return descriptive statistics
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=descriptive)  
    return descriptive

def difference(df: pd.DataFrame | str, data_col: str, compare_col: str, compare: list,
               same: bool=False, para: bool=True, alpha: float=0.05, within_cols:list=[], method:str='holm',
               dir: str=None, file: str=None) -> pd.DataFrame:
    ''' 
    difference(): computes the appropriate statistical test(s) and returns the p-value(s)
    
    Parameters:
    df (dataframe | str): tidy dataframe (or file path)
    data_col (str): data column name
    compare_col (str): comparisons column name
    compare (list): list of comparisons
    same (bool, optional): same (True) or different (False) subjects
    para (bool, optional): parameteric (True) or nonparametric (False)
    alpha (float, optional): significance level for the test (Default: 0.05)
    within_cols (list, optional): list of column names corresponding to different conditions or time points with the same subject (optional; para=True, same=True)
    method (str, optional): multiple hypothesis testing correction method (Default: holm)
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: pandas, scipy.stats, statsmodels.stats & io
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    if not same: # Different samples

        if para==True: # Data that follows a normal distribution

            if len(compare)==2: # Comparing 2 samples

                # Student's T Test Calculation
                print('Statistical Test: Student\'s T-test \n - Only use if you have two groups or if you are comparing just two of the n groups and are not concerned about inflating the Type I error rate (e.g., False Positives).')
                t_stat, p_value = ttest_ind(*(df[df[compare_col]==comp][data_col] for comp in compare))
                if p_value<alpha: null='reject'
                else: null='fail to reject'
                inference = pd.DataFrame({'test':['Student\'s T-test'],
                                          'comparison': [','.join(compare)], 
                                          'p_value': [p_value],
                                          'null_hypothesis':[null]})
            
            elif len(compare)>2: # Comparing more than 2 samples
                
                # 1-way Anova Test
                print('Statistical Test: 1-way Anova\n - If you want to compare all three groups simultaneously to determine if there is a significant difference in their means.\n - If the ANOVA finds a significant difference, then perform post-hoc tests to determine which specific groups differ from each other.')
                f_stat, p_value = f_oneway(*(df[df[compare_col]==comp][data_col] for comp in compare))
                if p_value<alpha: null='reject'
                else: null='fail to reject'
                inference = pd.DataFrame({'test':['1-way Anova Test'],
                                            'comparison': [','.join(compare)],
                                            'p_value':[p_value],
                                            'null_hypothesis': [null]})
                
                # Tukey's Honestly Significant Difference Test
                print('Statistical Test: Tukey\'s Honestly Significant Difference (HSD) Test\n - To control for Type-I error, adjust the significance threshold (α) to account for the number of tests.\n - Use with equal group sizes.')
                df2 = pd.concat([df[df[compare_col]==comp] for comp in compare]).reset_index(drop=True) # Only include specified data
                tukey_result = pairwise_tukeyhsd(df2[data_col], df2[compare_col], alpha=alpha)
                tukey_df = pd.DataFrame(data=tukey_result._results_table.data[1:], columns=tukey_result._results_table.data[0])
                inference = pd.concat([inference,
                                       pd.DataFrame({'test':['Tukey\'s HSD Test']*len(tukey_df['group1']),
                                                     'comparison': [','.join([tukey_df.iloc[i]['group1'],tukey_df.iloc[i]['group2']]) for i in range(len(tukey_df))], 
                                                     'p_value': tukey_df['p-adj'],
                                                     'null_hypothesis': ['reject' if p_adj<alpha else 'fail to reject' for p_adj in tukey_df['p-adj']]})]).reset_index(drop=True)
                
                # Multiple Hypothesis Correction for Student's T Test
                print(f'Statistical Test: Student\'s Paired T-test ({method} correction) \n - To control for Type-I error, adjust the significance threshold (α) to account for the number of tests.\n - Use with unequal group sizes.')
                tests = list(itertools.combinations(df2[compare_col].unique(), 2)) # Generate all comparisons from unique conditions
                p_values = [] # Perform pairwise t-tests
                comparisons = []
                for (cond1, cond2) in tests:
                    t_stat, p_value = ttest_ind(df2[df2[compare_col] == cond1][data_col], df2[df2[compare_col] == cond2][data_col])
                    p_values.append(p_value)
                    comparisons.append(','.join([cond1,cond2]))
                corrected_p_values = list(multipletests(p_values, method=method)[0]) # Apply Bonferroni correction
                inference = pd.concat([inference,
                                        pd.DataFrame({'test':[f'Student\'s T-test ({method} correction)']*len(tests),
                                                        'comparison': comparisons, 
                                                        'p_value': p_values,
                                                        'null_hypothesis':['reject' if corrected_p_value else 'fail to reject' for corrected_p_value in corrected_p_values]})]).reset_index(drop=True)
            
            else: print('Error: Invalid compare. List needs to contain 2 or more stings')
        
        else: # Data that does not follow a normal distribution
            
            # Mann Whitney U Test
            print(f'Statistical Test: Mann Whitney U Test ({method} correction)\n - The Mann-Whitney U test is a non-parametric test used to compare differences between two independent groups when the dependent variable is either ordinal or continuous, but not normally distributed.')
            df2 = pd.concat([df[df[compare_col]==comp] for comp in compare]).reset_index(drop=True) # Only include specified data
            tests = list(itertools.combinations(df2[compare_col].unique(), 2)) # Generate all test comparisons from unique conditions
            p_values = []
            comparisons = []
            for (cond1, cond2) in tests:
                u_stat, p_value = mannwhitneyu(df2[df2[compare_col] == cond1][data_col], df2[df2[compare_col] == cond2][data_col], alternative='two-sided')
                p_values.append(p_value)
                comparisons.append(','.join([cond1,cond2]))
            corrected_p_values = list(multipletests(p_values, method=method)[0]) # Apply Bonferroni correction
            inference = pd.DataFrame({'test':['Mann Whitney U Test' if len(tests)==1 else f'Mann Whitney U Test ({method} correction)']*len(tests),
                                        'comparison': comparisons, 
                                        'p_value': p_values,
                                        'null_hypothesis':['reject' if corrected_p_value else 'fail to reject' for corrected_p_value in corrected_p_values]})

    else: # Same sample

        if para==True:  # Data that follows a normal distribution

            if len(compare)==2: # Comparing 2 related (paired) samples
                
                # Student's Paired T-test
                print('Statistical Test: Paired Student\'s T-test \n - Only use if you have two groups or if you are comparing just two of the n groups and are not concerned about inflating the Type I error rate (e.g., False Positives).\n - Assumes the two groups are related or paired, meaning that each data point in one group has a corresponding data point in the other group.')
                t_stat, p_value = ttest_rel(*(df[df[compare_col]==comp][data_col] for comp in compare))
                if p_value<alpha: null='reject'
                else: null='fail to reject'
                inference = pd.DataFrame({'test':['Paired Student\'s T-test'],
                                            'comparison': [','.join(compare)], 
                                            'p_value': [p_value],
                                            'null_hypothesis':[null]})

            elif len(compare)>2: # Comparing 2 or more related (paired) samples
                
                # Repeated Anova
                print('Statistical Test: Repeated 1-way Anova\n - Use repeated measures ANOVA when you have multiple measurements from the same subjects or units under different conditions or time points.\n - It is particularly useful when the goal is to reduce error variability and account for within-subject correlations, thereby increasing the power of the test.')
                df2 = pd.concat([df[df[compare_col]==comp] for comp in compare]).reset_index(drop=True) # Only include specified data
                anova = AnovaRM(data=df2, depvar=data_col, subject=compare_col,within=within_cols)
                anova_result = anova.fit()
                anova_df = anova_result.summary().tables[0]
                inference = pd.DataFrame({'test':['Repeated Anova']*len(within_cols),
                                            'comparison': [','.join(compare) + f'; Within = {w}' for w in within_cols], 
                                            'p_value': [p_val for p_val in anova_df['Pr > F']],
                                            'null_hypothesis':['reject' if p_val<alpha else 'fail to reject' for p_val in anova_df['Pr > F']]})


                # Multiple Hypotheis Correction for Student's Paired T-test
                print(f'Statistical Test: Student\'s Paired T-test ({method} correction) \n - To control for Type-I error, adjust the significance threshold (α) to account for the number of tests.')
                for w in within_cols: # Iterate through within subject columns
                    tests = list(itertools.combinations(df2[w].unique(), 2)) # Generate all pairwise comparisons from unique conditions
                    p_values = [] # Perform pairwise t-tests
                    comparisons = []
                    for (cond1, cond2) in tests:
                        t_stat, p_value = ttest_rel(df2[df2[w] == cond1][data_col], df2[df2[w] == cond2][data_col])
                        p_values.append(p_value)
                        comparisons.append(','.join(compare) + f'; Within = {w}; ' + ','.join([cond1,cond2]))
                    corrected_p_values = list(multipletests(p_values, method=method)[0]) # Apply Bonferroni correction
                    inference = pd.concat([inference,
                                            pd.DataFrame({'test':[f'Student\'s Paired T-test ({method} correction)']*len(tests),
                                                            'comparison': comparisons, 
                                                            'p_value': p_values,
                                                            'null_hypothesis':['reject' if corrected_p_value else 'fail to reject' for corrected_p_value in corrected_p_values]})]).reset_index(drop=True)

            else: print('Error: Invalid compare. List needs to contain 2 or more stings')

        else: # Data that does not follow a normal distribution
            
            # Wilcoxon Signed-Rank Test
            print(f'Statistical Test: Wilcoxon Signed-Rank Test ({method} correction) \n - To control for Type-I error, adjust the significance threshold (α) to account for the number of tests.')
            df2 = pd.concat([df[df[compare_col]==comp] for comp in compare]).reset_index(drop=True) # Only include specified data
            tests = list(itertools.combinations(df2[compare_col].unique(), 2)) # Generate all pairwise comparisons from unique conditions
            p_values = [] # Perform pairwise t-tests
            comparisons = []
            for (cond1, cond2) in tests:
                w_stat, p_value = wilcoxon(df2[df2[compare_col] == cond1][data_col], df2[df2[compare_col] == cond2][data_col])
                p_values.append(p_value)
                comparisons.append(','.join([cond1,cond2]))
            corrected_p_values = list(multipletests(p_values, method=method)[0]) # Apply Bonferroni correction
            inference = pd.DataFrame({'test':['Wilcoxon Signed-Rank Test' if len(tests)==1 else f'Wilcoxon Signed-Rank Test ({method} correction)']*len(tests),
                                        'comparison': comparisons, 
                                        'p_value': p_values,
                                        'null_hypothesis':['reject' if corrected_p_value else 'fail to reject' for corrected_p_value in corrected_p_values]})
    
    # Save & return descriptive statistics
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=inference)  
    return inference

def correlation(df: pd.DataFrame | str, var_cols: list=[], value_cols: list=[], method: str='pearson', numeric_only: bool=True,
                plot: bool=True, dir:str=None, file_data:str=None, file_plot:str=None, **kwargs_plot) -> pd.DataFrame:
    ''' 
    correlation(): returns a correlation matrix & plot
    
    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    var_cols (list, optional): list of 2 variable column names for tidy dataframe (optional; pivot table index & column)
    value_cols (list, optional): list of numerical column names to compute statistics; single column name for tidy dataframe (optional)
    method (str, optional): pearson, spearman, or kendall (Default: pearson)
    numeric_only (bool, optional): only calculates correlations for numeric columns (Default: True)
    plot (bool, optional): generate correlation matrix plot (Default: True)
    dir (str, optional): save directory
    file_data (str, optional): save data file (e.g., csv)
    file_plot (str, optional): save plot file (e.g., pdf)
    kwargs_plot (dict, optional): plotting keyword arguments
    
    Depedencies: pandas, io
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    if (len(var_cols)==2)&(len(value_cols)==1): df = df.pivot(index=var_cols[0],columns=var_cols[1],values=value_cols[0]) # Splits tidy dataframe
    elif len(value_cols)>=1: df = df[value_cols] # Isolate specified columns for non-tidy dataframe
    df_corr = df.corr(method=method,numeric_only=numeric_only) # Correlation matrix with specified method

    # Plot, save & return correlation matrix
    if plot == True:
        p.heat(df=df_corr, cbar_label=method, dir=dir, file=file_plot, **kwargs_plot)
    if dir is not None and file_data is not None:
        io.save(dir=dir,file=file_data,obj=df_corr,id=True)
    return df_corr

def weighted_correlation(df: pd.DataFrame | str, x: str, y: str, weight: str=None,
                        method: str='pearson') -> float:
    ''' 
    weighted_correlation(): computes the weighted correlation
    
    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    x (str): first variable column name
    y (str): second variable column name
    weight (str, optional): weights column name (Default: None)
    method (str, optional): pearson, spearman, or kendall (Default: pearson)
    
    Dependencies: numpy
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    if weight is None: # Unweighted correlation
        weight ='__weight__'
        df[weight] = [1]*len(df)

    def weighted_pearson(x, y, w):
        ''' 
        weighted_pearson() Compute weighted Pearson correlation coefficient
        '''
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        w = np.asarray(w, dtype=float)

        # normalize weights (optional but numerically nicer)
        w = w / w.sum()

        mx = np.sum(w * x)
        my = np.sum(w * y)

        cov_xy = np.sum(w * (x - mx) * (y - my))
        sx = np.sqrt(np.sum(w * (x - mx)**2))
        sy = np.sqrt(np.sum(w * (y - my)**2))

        return cov_xy / (sx * sy)

    if method == 'pearson':
        return weighted_pearson(df[x], df[y], df[weight])

    elif method == 'spearman':
        x = np.asarray(df[x], dtype=float)
        y = np.asarray(df[y], dtype=float)
        w = np.asarray(df[weight], dtype=float)

        # rank the data (average ranks for ties, like standard Spearman)
        rx = rankdata(x, method="average")
        ry = rankdata(y, method="average")

        return weighted_pearson(rx, ry, w)
    
    elif method == 'kendall':
        x = np.asarray(df[x], dtype=float)
        y = np.asarray(df[y], dtype=float)
        n = len(x)
        w = np.asarray(df[weight], dtype=float)

        num = 0.0
        den_x = 0.0
        den_y = 0.0

        for i in range(n - 1):
            for j in range(i + 1, n):
                wxw = w[i] * w[j]

                sx = np.sign(x[i] - x[j])
                sy = np.sign(y[i] - y[j])

                # skip pairs tied in x or y (no ordering information)
                if sx == 0 or sy == 0:
                    continue

                num += wxw * sx * sy      # +wxw for concordant, -wxw for discordant
                den_x += wxw * (sx**2)    # == wxw
                den_y += wxw * (sy**2)    # == wxw

        if den_x == 0 or den_y == 0:
            return np.nan

        return num / np.sqrt(den_x * den_y)

    else:
        raise ValueError(f"Invalid method: {method}. Choose 'pearson', 'spearman', or 'kendall'.")

# Statistics & Plotting
def corr_line(df: pd.DataFrame, x: str, y: str, ax=None, 
                color='black', linestyle='--', linewidth=1, alpha=0.5):
    ''' 
    corr_line(): compute linear regression line and add it to the scatter plot 

    Parameters:
    df (pandas.DataFrame): dataframe with x and y data
    x (str): x data column
    y (str): y data column
    ax (matplotlib.axes.Axes, optional): Axes object to plot on
    color (str, optional): line color (Default: 'black')
    linestyle (str, optional): line style (Default: '--')
    linewidth (int, optional): line width (Default: 1)
    alpha (float, optional): transparency level (Default: 0.5)
    '''
    x = np.asarray(df[x], dtype=float)
    y = np.asarray(df[y], dtype=float)

    # Fit slope (b) and intercept (a)
    b, a = np.polyfit(x, y, 1)

    # Make line
    xx = np.linspace(x.min(), x.max(), 200)
    yy = a + b * xx

    ax.plot(xx, yy, color=color, linestyle=linestyle, linewidth=linewidth, alpha=alpha)
    return a, b

def weighted_corr_line(df: pd.DataFrame, x: str, y: str, weight: str, ax=None, 
                        color='black', linestyle='--', linewidth=1, alpha=0.5):
    '''
    weighted_corr_line(): compute weighted linear regression line and add it to the scatter plot

    Parameters:
    df (pandas.DataFrame): dataframe with x, y, and weight data
    x (str): x data column
    y (str): y data column
    weight (str): weight data column
    ax (matplotlib.axes.Axes, optional): Axes object to plot on
    color (str, optional): line color (Default: 'black')
    linestyle (str, optional): line style (Default: '--')
    linewidth (int, optional): line width (Default: 1)
    alpha (float, optional): transparency level (Default: 0.5)
    '''
    x = np.asarray(df[x], float)
    y = np.asarray(df[y], float)
    w = np.asarray(df[weight], float)

    # Normalize weights
    w = w / w.sum()

    # Weighted means
    mx = np.sum(w * x)
    my = np.sum(w * y)

    # Weighted slope
    b = np.sum(w * (x - mx) * (y - my)) / np.sum(w * (x - mx)**2)

    # Weighted intercept
    a = my - b * mx

    # Plot the line
    xx = np.linspace(x.min(), x.max(), 200)
    yy = a + b * xx

    ax.plot(xx, yy, color=color, linestyle=linestyle, linewidth=linewidth, alpha=alpha)
    return a, b

# Comparison
def compare(df: pd.DataFrame | str, sample: str, cond: str, cond_comp: str, 
            var: str, count: str, pseudocount: int=1, alternative: str='two-sided', replicate: str = None,
            dir: str = None, file: str = None, verbose: bool = False) -> pd.DataFrame:
    ''' 
    compare(): computes FC, pval, and log transformations relative to a specified condition

    Parameters:
    df (dataframe | str): tidy dataframe (or file path)
    sample (str): sample column name
    cond (str): condition column name
    cond_comp (str): condition for comparison group
    var (str): variable column name
    count (str): count column name
    pseudocount (int, optional): pseudocount to avoid log(0) & /0 (Default: 1)
    alternative (str, optional): alternative hypothesis for statistical test ('two-sided', 'less', or 'greater'; Default: 'two-sided')
    dir (str, optional): save directory
    file (str, optional): save file
    verbose (bool, optional): print progress to console (Default: False)
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Validate replicate column
    if replicate is not None:
        if replicate not in df.columns: # Replicate column exists
            raise ValueError(f"replicate column '{replicate}' not found in df")
        
        # Even number of replicates per condition
        df_replicate_vc = df[replicate].value_counts().reset_index()
        df_replicate_vc_1 = df_replicate_vc.iloc[0]['count']
        if any(df_replicate_vc['count'] != df_replicate_vc_1):
            raise ValueError(f"Unequal number of replicates per condition in replicate column '{replicate}'")    

    # Get metadata
    meta_cols = list(df.columns)
    if count not in meta_cols: # Validate count column
        raise ValueError(f"count column '{count}' not found in df")
    meta_cols.remove(count)

    if verbose: 
        print(f'Add pseudocount ({pseudocount}) to avoid log(0) & compute fraction per million (FPM)')
    df[f'{count}+{pseudocount}'] = df[count] + pseudocount

    # Compute per-sample library-normalized values
    dc = t.split(df=df, key=sample)
    for s, df_sample in dc.items():
        count_total = float(df_sample[count].sum())
        n_rows = int(df_sample.shape[0])
        count_total_pc = count_total + pseudocount * n_rows

        # Classic FPM; avoid divide-by-zero
        if count_total > 0:
            dc[s]['FPM'] = [c / count_total * 1_000_000 for c in df_sample[count]]
        else:
            dc[s]['FPM'] = [0.0] * n_rows

        # Pseudocount-adjusted FPM for robust FC (adds pseudocount to each observation and to the library size)
        if count_total_pc > 0:
            dc[s]['FPM_pc'] = [(c + pseudocount) / count_total_pc * 1_000_000 for c in df_sample[count]]
        else:
            dc[s]['FPM_pc'] = [0.0] * n_rows

    if replicate is None: # Aggregate version

        if verbose: 
            print('Aggregate version:\nGroup samples by condition & compute averages')

        df_cond_stat = pd.DataFrame()
        # Define aggregation dictionary dynamically
        agg_dict = {}
        for col in meta_cols: # Include metadata columns as sets
            agg_dict[col] = lambda x: set(x)
        # Include both mean and list of original values
        agg_dict[count] = 'mean'
        agg_dict[f'{count}+{pseudocount}'] = 'mean'
        agg_dict['FPM'] = ['mean', list]
        agg_dict['FPM_pc'] = ['mean', list]

        # Join samples back into 1 dataframe and split by condition
        for c, df_cond in t.split(df=t.join(dc=dc, col=sample), key=cond).items():
            if verbose: 
                print(c)

            # Group by variable and aggregate; keep grouping key as a column
            df_cond_agg = df_cond.groupby(by=var, as_index=False).agg(agg_dict)
            df_cond_agg.columns = ['_'.join(col).strip('_') for col in df_cond_agg.columns]
            for col in df_cond_agg.columns:
                if '_<lambda>' in col:
                    if any(isinstance(item, str) for item in df_cond_agg.iloc[0][col]):
                        ls = []
                        for s in df_cond_agg[col]: 
                            ls.append([s_ if isinstance(s_, str) else str(s_) for s_ in s])
                        df_cond_agg[col] = [','.join(sorted(l)) for l in ls]
                    else:
                        df_cond_agg[col] = [sorted(s) for s in df_cond_agg[col]]
            df_cond_agg.columns = df_cond_agg.columns.str.replace('_<lambda>', '', regex=True)
            df_cond_agg[cond] = c
            df_cond_stat = pd.concat([df_cond_stat, df_cond_agg]).reset_index(drop=True)
        
        # Fold change & p-value relative comparison group
        if verbose: 
            print(f'Compute FC & pval relative to {cond_comp}:')
        df_stat = pd.DataFrame()
        df_comp = df_cond_stat[df_cond_stat[cond] == cond_comp]
        df_other = df_cond_stat[df_cond_stat[cond] != cond_comp]
        
        # Only evaluate variables present in both the comparison group and the other groups
        vars_in_both = sorted(set(df_other[var].tolist()).intersection(set(df_comp[var].tolist())))
        for v in vars_in_both:
            if verbose: 
                print(f'{v}')
            df_other_edit = df_other[df_other[var] == v].copy()
            df_comp_edit = df_comp[df_comp[var] == v]

            # If the comparison group has no rows for this variable, skip
            if df_comp_edit.empty:
                continue

            # Carry through comparison means for transparency
            df_other_edit[f'{count}_mean_compare'] = [df_comp_edit.iloc[0][f'{count}_mean']] * df_other_edit.shape[0]
            df_other_edit[f'{count}+{pseudocount}_mean_compare'] = [df_comp_edit.iloc[0][f'{count}+{pseudocount}_mean']] * df_other_edit.shape[0]
            df_other_edit['FPM_pc_mean_compare'] = [df_comp_edit.iloc[0]['FPM_pc_mean']] * df_other_edit.shape[0]
            df_other_edit['FPM_mean_compare'] = [df_comp_edit.iloc[0]['FPM_mean']] * df_other_edit.shape[0]

            # Robust FC using pseudocount-adjusted FPM; guard against zero denominator
            _den = float(df_comp_edit.iloc[0]['FPM_pc_mean'])
            if _den <= 0:
                _den = np.finfo(float).eps
            df_other_edit['FC'] = df_other_edit['FPM_pc_mean'] / _den

            # Two-sample t-tests using the adjusted per-sample values
            ttests = [
                ttest_ind(list(other_vals), list(df_comp_edit.iloc[0]['FPM_pc_list']), alternative=alternative)
                for other_vals in df_other_edit['FPM_pc_list']
            ]
            df_other_edit['pval'] = [ttest[1] for ttest in ttests]
            df_other_edit['tstat'] = [ttest[0] for ttest in ttests]
            df_stat = pd.concat([df_stat, df_other_edit])
        df_stat['compare'] = [cond_comp] * df_stat.shape[0]
        df_stat = df_stat.sort_values(by=[cond, var]).reset_index(drop=True)

    else: # Replicate-based version
        
        if verbose: 
            print(f'Replicate-based version:\nCompute FC within replicates relative to {cond_comp}:')

        df_stat = pd.DataFrame()
        df_comp = t.join(dc=dc, col=sample)
        df_comp = df_comp[df_comp[cond] == cond_comp]
        df_other = t.join(dc=dc, col=sample)
        df_other = df_other[df_other[cond] != cond_comp]

        # Pair within replicate: match on (var, replicate)
        paired = df_other.merge(
            df_comp[[var, replicate, 'FPM_pc', 'FPM']],
            on=[var, replicate],
            how='inner',
            suffixes=('', '_compare')
        )

        # Compute FC within replicate
        eps = np.finfo(float).eps
        paired['compare'] = cond_comp
        paired['FPM_pc_compare'] = paired['FPM_pc_compare'].astype(float)
        paired['FC'] = paired['FPM_pc'].astype(float) / np.maximum(paired['FPM_pc_compare'], eps)
        paired['log2FC'] = np.log2(np.maximum(paired['FC'].astype(float), eps))
        df_stat = paired.sort_values(by=[cond, var, replicate]).reset_index(drop=True)
        
    # Save & return statistics dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df_stat) 
    return df_stat

def odds_ratio(df: pd.DataFrame | str, cond: str, cond_comp: str, 
            var: str, var_comp: str, count: str, pseudocount: int=1,
            alternative: str='two-sided', dir:str=None, file:str=None, verbose:bool=False) -> pd.DataFrame:
    ''' 
    odds_ratio(): computes odds ratios relative to a specified condition & variable (e.g., unedited & WT)

    Parameters:
    df (dataframe | str): tidy dataframe (or file path)
    cond (str): condition column name
    cond_comp (str): condition for comparison group
    var (str): variable column name for odds ratio computation (e.g., mut1, mut2, ..., WT)
    var_comp (str): variable name for comparison (e.g., WT)
    count (str): count column name
    pseudocount (int, optional): pseudocount to avoid /0 (Default: 1)
    alternative (str, optional): 'two-sided', 'less', or 'greater' (Default: 'two-sided')
    dir (str, optional): save directory
    file (str, optional): save file
    verbose (bool, optional): print progress to console (Default: False)
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Get metadata
    meta_cols = list(df.columns)
    meta_cols.remove(count)

    if verbose: print(f'Combine counts across samples in the same condition & add pseudocount ({pseudocount}) to avoid /0')
    # Define aggregation dictionary dynamically
    agg_dict = {}
    for col in meta_cols: # Include metadata columns as sets
        agg_dict[col] = lambda x: set(x)
    # Include count sum
    agg_dict[count] = 'sum'

    # Group by condition and variable; aggregate; keep grouping key as a column
    df = df.groupby(by=[cond, var], as_index=False).agg(agg_dict)
    for col in meta_cols:
        if any(isinstance(item, str) for item in df.iloc[0][col]):
            ls = []
            for s in df[col]: ls.append([s_ if isinstance(s_, str) else str(s_) for s_ in s])
            df[col] = [','.join(sorted(l)) for l in ls]
        else:
            df[col] = [sorted(s) for s in df[col]]
    
    # Add pseudocount to counts
    df[f'{count}+{pseudocount}'] = df[count] + pseudocount

    if verbose: print(f'Split comparison groups: condition from ({var_comp}) from other variables')
    # Get comparison variable count and count+pseudocount for comparison condition
    df_cond_comp = df[df[cond] == cond_comp].reset_index(drop=True)
    df_cond_comp_var_comp = df_cond_comp[df_cond_comp[var] == var_comp]
    if len(df_cond_comp_var_comp) == 0:
        raise ValueError(f'No data found for condition "{cond_comp}" and variable "{var_comp}". Cannot compute odds ratio.')
    elif len(df_cond_comp_var_comp) > 1:
        raise ValueError(f'Multiple entries found for condition "{cond_comp}" and variable "{var_comp}". Ensure unique entries for odds ratio computation.')
    else:
        cond_comp_var_comp_count = df_cond_comp_var_comp[count].values[0]
        cond_comp_var_comp_count_pc = df_cond_comp_var_comp[f'{count}+{pseudocount}'].values[0]

    # Remove comparison condition from main dataframe
    df = df[df[cond] != cond_comp].reset_index(drop=True)

    if verbose: print(f'Split dataframe by condition & compute odds ratio:')
    # Iterate through conditions 
    df_stat = pd.DataFrame()
    for c in df[cond].unique():
        
        # Isolate condition comparison dataframe
        df_cond = df[df[cond] == c].reset_index(drop=True)
        df_cond_var_comp = df_cond[df_cond[var] == var_comp]
        if len(df_cond_var_comp) == 0:
            raise ValueError(f'No data found for condition "{c}" and variable "{var_comp}". Cannot compute odds ratio.')
        elif len(df_cond_var_comp) > 1:
            raise ValueError(f'Multiple entries found for condition "{c}" and variable "{var_comp}". Ensure unique entries for odds ratio computation.')
        else:
            cond_var_comp_count = df_cond_var_comp[count].values[0]
            cond_var_comp_count_pc = df_cond_var_comp[f'{count}+{pseudocount}'].values[0]
        
        # Remove comparison variable from condition dataframe
        df_cond = df_cond[df_cond[var] != var_comp].reset_index(drop=True)

        # Include comparison counts for transparency
        df_cond[f'{count}_compare_{var_comp}'] = [cond_var_comp_count] * df_cond.shape[0]
        df_cond[f'{count}+{pseudocount}_compare_{var_comp}'] = [cond_var_comp_count_pc] * df_cond.shape[0]
        df_cond[f'{count}_compare_{cond_comp}_{var_comp}'] = [cond_comp_var_comp_count] * df_cond.shape[0]
        df_cond[f'{count}+{pseudocount}_compare_{cond_comp}_{var_comp}'] = [cond_comp_var_comp_count_pc] * df_cond.shape[0]
        df_cond[f'compare_{cond}'] = [cond_comp] * df_cond.shape[0]
        df_cond[f'compare_{var}'] = [var_comp] * df_cond.shape[0]

        # Iterate through variables within condition; compute odds ratio & pval from fisher exact test
        fe_or_ls = []
        pval_ls = []

        for v in df_cond[var].unique(): # Iterate through variables
            # Isolate variable & variable comparison counts
            cond_var_count_pc = df_cond[df_cond[var] == v][f'{count}+{pseudocount}'].values[0]
            try:
                cond_comp_var_count_pc = df_cond_comp[df_cond_comp[var] == v][f'{count}+{pseudocount}'].values[0]
            except IndexError:
                raise ValueError(f'No data found for condition "{cond_comp}" and variable "{v}". Cannot compute odds ratio.')
            
            # Compute odds ratio & fisher exact test
            fe_or, pval = fisher_exact(table=[[cond_var_count_pc, cond_var_comp_count_pc], [cond_comp_var_count_pc, cond_comp_var_comp_count_pc]],
                                        alternative=alternative)
            fe_or_ls.append(fe_or)
            pval_ls.append(pval)
        
        # Add results to condition dataframe
        df_cond['fisher_exact_odds_ratio'] = fe_or_ls
        df_cond['fisher_exact_pval'] = pval_ls
        
        # Append to statistics dataframe
        df_stat = pd.concat([df_stat, df_cond]).reset_index(drop=True)
    
    # Save & return statistics dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df_stat) 
    return df_stat

def zscore(
    df: pd.DataFrame | str,
    val: str,
    cond_col: str,
    var_col: str,
    var: str,
    out_col: str | None = None,
    ddof: int = 1,
    dir:str=None,
    file:str=None,
) -> pd.DataFrame:
    """
    zscore(): Z-score `val` within each `cond_col` using stats computed from rows where `var_col == var`.

    For each condition:
      - baseline = rows with df[var_col] == var
      - mu, sigma computed from baseline[val]
      - z = (val - mu) / sigma applied to *all* rows in that condition

    Parameters:
    df (pd.DataFrame): Input dataframe.
    cond_col (str): Column name defining condition groups.
    val (str): Column name of values to z-score.
    var_col (str): Column name defining variable groups (used to pick baseline rows).
    var (str): The value in `var_col` used as the baseline set for mean/std.
    out_col (str, optional): Output column name (Default: f"{val}_z").
    ddof (int): Delta degrees of freedom for std (0 = population, 1 = sample).
    dir (str, optional): save directory
    file (str, optional): save file
    
    Returns: Copy of df (pd.DataFrame) with added z-score column.

    Notes:
    - If a condition has no baseline rows, output z-scores are NaN for that condition.
    - If baseline std is 0 or NaN, output z-scores are NaN for that condition.
    """
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    if out_col is None:
        out_col = f"{val}_z"

    required = {val, cond_col, var_col}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    out = df.copy()

    # Compute baseline mean/std per condition using only rows where var_col == var
    baseline = out.loc[out[var_col].eq(var), [cond_col, val]].copy()

    stats = (
        baseline.groupby(cond_col, dropna=False)[val]
        .agg(mu="mean", sigma=lambda s: s.std(ddof=ddof))
    )

    # Map per-condition mu/sigma back onto all rows
    out["_mu__"] = out[cond_col].map(stats["mu"])
    out["_sigma__"] = out[cond_col].map(stats["sigma"])

    denom = out["_sigma__"].replace(0, np.nan)
    out[out_col] = (out[val] - out["_mu__"]) / denom
    out.drop(columns=["_mu__", "_sigma__"], inplace=True)
    
    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=out)  
    return out