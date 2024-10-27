### stat.py ###
# Author: Marc Zepeda
# Date: 2024-08-26

import itertools
import pandas as pd
import numpy as np
from scipy.stats import skew, kurtosis, ttest_ind, ttest_rel, f_oneway, ttest_ind, ttest_rel, mannwhitneyu, wilcoxon
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

# Statistics methods
def describe(df: pd.DataFrame, cols=[], group=''):
    ''' 
    describe(): returns descriptive statistics for numerical columns in a DataFrame
    
    Parameters:
    df (dataframe): pandas dataframe
    cols (list, optional): list of numerical columns to compute statistics
    group (str, optional): column name to split tidy dataframes
    
    Dependencies: pandas, numpy, & scipy.stats
    '''
    if group!='': df = df.pivot(columns=group) # Splits tidy dataframe
    if len(cols)>0: df = df[cols] # Isolate specified columns

    descriptive = pd.DataFrame()    
    descriptive['mean'] = df.mean() # Mean
    descriptive['median'] = df.median() # Median
    descriptive['variance'] = df.var() # Variance
    descriptive['std_dev'] = df.std() # Standard Deviation
    descriptive['mad'] = df.apply(lambda x: np.median(np.abs(x - x.median()))) # Median Absolute Deviation
    descriptive['min'] = df.min() # Minimum
    descriptive['max'] = df.max() # Maximum
    descriptive['range'] = df.max() - df.min() # Range
    descriptive['skewness'] = df.apply(lambda x: skew(x, nan_policy='omit')) # Skewness
    descriptive['kurtosis'] = df.apply(lambda x: kurtosis(x, nan_policy='omit')) # Kurtosis
    descriptive['count'] = df.count() # Count (non-missing observations)
    descriptive['sum'] = df.sum() # Sum
    descriptive['25%'] = df.quantile(0.25)  # Quantiles (25%, 50%, 75%)
    descriptive['50%'] = df.quantile(0.50)
    descriptive['75%'] = df.quantile(0.75)

    return descriptive

def difference(df: pd.DataFrame,data_col: str,compare_col: str,compare: list,same=False,para=True,alpha=0.05,within_cols=[],method='holm'):
    ''' 
    difference(): computes the appropriate statistical test(s) and returns the p-value(s).
    
    Parameters:
    df: tidy DataFrame
    data_col: data column name
    compare_col: comparisons column name
    compare: list of comparisons
    same: same (True) or different (False) subjects
    compare: list of comparisons
    para: parameteric (True) or nonparametric (False)
    alpha: significance level for the test (Default: 0.05)
    within_cols: list of column names corresponding to different conditions or time points with the same subject (optional; para=True, same=True)
    method: multiple hypothesis testing correction method (Default: holm)

    Dependencies: pandas, scipy.stats, & statsmodels.stats
    '''
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
                
                return inference
            
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
                                        pd.DataFrame({'test':['Tukey\'s HSD Test']*len(compare),
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

                return inference
            
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

            return inference

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
                
                return inference

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
                
                return inference

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
            
            return inference

def correlation(df: pd.DataFrame, var_cols=[], value_cols=[], method='pearson',numeric_only=True):
    ''' 
    correlation(): returns a correlation matrix
    
    Parameters:
    df (dataframe): pandas dataframe
    var_cols (list, optional): list of 2 variable column names for tidy dataframe (optional; pivot table index & column)
    value_cols (list, optional): list of numerical column names to compute statistics; single column name for tidy dataframe (optional)
    method (str, optional): pearson, spearman, or kendall (Default: pearson)
    numeric_only (bool, optional): only calculates correlations for numeric columns (Default: True)
    
    Depedencies: pandas
    '''
    if (len(var_cols)==2)&(len(value_cols)==1): df = df.pivot(index=var_cols[0],columns=var_cols[1],values=value_cols[0]) # Splits tidy dataframe
    elif len(value_cols)>=1: df = df[value_cols] # Isolate specified columns for non-tidy dataframe
    return df.corr(method=method,numeric_only=numeric_only) # Correlation matrix with specified method