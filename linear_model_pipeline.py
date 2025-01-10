import numpy as np 
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
from random import choices
from sklearn.linear_model import LinearRegression
import pickle
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import random
import seaborn as sns
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu
from scipy.stats import bootstrap
import statsmodels.api as sm
from scipy.stats import normaltest
from scipy.stats import ttest_1samp
from scipy.stats import ttest_ind
import statsmodels.stats.api as sms

def compute_betas_and_rsquared(adult_data,infant_data,channels):
    adult_data_mean = adult_data.groupby(['channel']).mean()
    x = np.array(adult_data_mean['tau'])

    infant_data_mean = infant_data.groupby(['subject','channel']).mean()
    infant_data_mean['subject'] = [item[0] for item in infant_data_mean.index]
    all_y = infant_data_mean[['subject','tau']]

    all_sub_list = list(all_y['subject'].unique())
    betas_out = np.zeros(len(all_sub_list))
    rsquared_out = np.zeros(len(all_sub_list))
    y_matrix = np.zeros((len(all_sub_list),channels))
    for i,sub_code in enumerate(all_sub_list):
        y_matrix[i,:] = all_y[all_y['subject']==sub_code]['tau'].values
    
    if np.isnan(y_matrix).any():
        imputer = IterativeImputer(random_state=0)
        imputer.fit(y_matrix)
        y_matrix = imputer.transform(y_matrix)

    for count,sub in enumerate(all_sub_list):
        y = y_matrix[count,:]
        model = LinearRegression().fit(x.reshape(-1,1),y)
        betas_out[count] = np.mean(model.coef_)
        rsquared_out[count] = model.score(x.reshape(-1,1),y)

    return betas_out,rsquared_out


def compute_stats(betas_dict, rsquared_dict, alpha, group_flag):
    combinations = ['video-eyes_open','video-eyes_closed','eyes_open-eyes_closed']
    timepoints= [6,9,16]
    conditions= ['video','eyes_open','eyes_closed']
    with open(f'/eeg/EEG_timescales/Results/nobootstrapping_linear_model_adult_comparison_allstats_pipeline{group_flag}.txt','a') as f:
        for t in timepoints:
            f.write(f'\n ##### {t} \n')  
            f.write(f'Corrected alpha = {alpha/3} \n')
            betas=betas_dict[t]
            rsquared=rsquared_dict[t]
            f.write(f'\n ##### Age: {t} \n')
            for cond in conditions:
                print(f'normality {t} {cond}:{shapiro(betas[cond])}')
                t_1sample,p_t_1sample= ttest_1samp(betas[cond],0)
                f.write(f'{cond} (M={np.mean(betas[cond])}, SD={np.std(betas[cond])}, CI=[{sms.DescrStatsW(betas[cond]).tconfint_mean()}]) different from zero: t={round(t_1sample,3)}, p={p_t_1sample} \n')
                rsquared_mean = np.mean(rsquared[cond])
                f.write(f'Average R2 for {t} and {cond}: {np.round(rsquared_mean,3)} \n')
            f.write('## Comparisons \n')
            for combination in combinations:
                key1 = combination.split('-')[0]
                key2 = combination.split('-')[1]
                u,p=ttest_ind(betas[key1],betas[key2])
                f.write(f'\n {key1} vs {key2}: t={round(u,3)}, p={p} \n')



