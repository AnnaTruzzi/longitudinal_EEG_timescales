import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mne
import os
from scipy import optimize
import scipy 
import glob
import random
import pickle
import tau_estimate
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
import linear_model_pipeline
from scipy.stats import shapiro
import seaborn as sns
from scipy.stats import ttest_1samp

# function to calculate Cohen's d for independent samples
# from: https://machinelearningmastery.com/effect-size-measures-in-python/
def cohend(d1, d2):
	# calculate the size of samples
	n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
	s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
	# calculate the pooled standard deviation
	s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
	u1, u2 = np.mean(d1), np.mean(d2)
	# calculate the effect size
	return (u1 - u2) / s


if __name__ == '__main__':

    #group='base' # or 'validation'
    group='validation'
    comparison_across_groups=False
    comparison_within_adults=False
    linear_model_analysis=True

    timepoints = [6,9,16,'adults']
    adult_conditions = ['video','eyes_open','eyes_closed']
    rootpth = '/eeg/ITS_Project_NewPreproc'
    

    if group=='base':
        group_flag = ''
    else:
        group_flag = '_validation'


    if comparison_across_groups:
        tau_dict = {}
        infant_data = pd.read_csv(f'/eeg/EEG_timescales/Data/develop_tau{group_flag}_v3.csv')
        infant_data = infant_data[(infant_data['pernan']<0.25) & (infant_data['event']!= 9999)]
        infant_data_mean = infant_data.groupby(['ses_age','channel']).mean()[['tauimp']]
        ses_age = [i[0] for i in infant_data_mean.index]
        infant_data_mean['ses_age'] = ses_age
        channel = [i[1] for i in infant_data_mean.index]
        infant_data_mean['channel'] = channel
        adult_data = pd.read_csv(f'/eeg/EEG_timescales/Data/adult_tau_v3.csv')
        adult_data = adult_data[adult_data['pernan']<0.25]
        adult_data.loc[adult_data.event=='2222','event'] = 'video'
        for cond in adult_conditions:
            adult_data_cond=adult_data[adult_data['event']==cond]
            adult_data_mean = adult_data_cond.groupby('channel').mean()[['tauimp']]
            adult_data_mean['ses_age']=np.repeat('adults',adult_data_mean.shape[0])
            channel_adult = [i[0] for i in adult_data_mean.index]
            adult_data_mean['channel'] = channel_adult
            all_tau = pd.concat([infant_data_mean, adult_data_mean], ignore_index=True)
            all_tau = all_tau.rename(columns={"tauimp": "tau"})

            with open(f'//eeg/EEG_timescales/Results/tau_across_groups_analysis_adult{cond}{group_flag}.txt', 'w') as f:
                statistic,p = kruskal(all_tau[all_tau['ses_age']==6]['tau'],all_tau[all_tau['ses_age']==9]['tau'],all_tau[all_tau['ses_age']==16]['tau'],all_tau[all_tau['ses_age']=='adults']['tau'])
                print(statistic,p)
                f.write(f'\n Kruskal test: chi = {statistic}, p = {round(p,3)} \n')
            
                combinations = [(6,9),(6,16),(6,'adults'),(9,16),(9,'adults'),(16,'adults')]
                f.write(f'ACROSS GROUP COMPARISONS - alpha={0.05/6} \n')

                for combination in combinations:
                    u,p = mannwhitneyu(all_tau[all_tau['ses_age']==combination[0]]['tau'],all_tau[all_tau['ses_age']==combination[1]]['tau'])
                    d = cohend(all_tau[all_tau['ses_age']==combination[0]]['tau'],all_tau[all_tau['ses_age']==combination[1]]['tau'])
                    print(f'{combination[0]} vs {combination[1]} - U = {u}, p = {round(p,4)}, d = {d}')
                    mean1 = np.nanmean(np.array(all_tau[all_tau['ses_age']==combination[0]]['tau']))
                    mean2= np.nanmean(np.array(all_tau[all_tau['ses_age']==combination[1]]['tau']))
                    sd1 = np.nanstd(np.array(all_tau[all_tau['ses_age']==combination[0]]['tau']))
                    sd2= np.nanstd(np.array(all_tau[all_tau['ses_age']==combination[1]]['tau']))
                    f.write(f'###### {combination[0]} (M={mean1}, SD={sd1}) vs {combination[1]} (M={mean2}, SD={sd2})\n')
                    f.write(f'Mann Whitney U - U = {u}, p = {p}, d = {d} \n')


    if comparison_within_adults:
        tau_list=[]
        event_list=[]
        data = pd.read_csv(f'/eeg/EEG_timescales/Data/adult_tau_v3.csv')
        data = data[data['pernan']<0.25]
        data.loc[data.event=='2222','event'] = 'video'
        data = data.rename(columns={"tauimp": "tau"})
        for cond in adult_conditions:
            event_data=data[data['event']==cond]
            data_mean = event_data.groupby('channel').mean()
            tau_list.extend(list(data_mean['tau']))
            event_list.extend(np.repeat(cond,len(data_mean['tau'])))
        
        all_tau_dict = {'cond_name':event_list,'tau':tau_list}
        all_tau=pd.DataFrame(all_tau_dict)

        with open(f'//eeg/EEG_timescales/Results/tau_across_adult conditions.txt', 'w') as f:
            statistic,p = kruskal(all_tau[all_tau['cond_name']=='video']['tau'],all_tau[all_tau['cond_name']=='eyes_open']['tau'],all_tau[all_tau['cond_name']=='eyes_closed']['tau'])
            print(statistic,p)
            f.write(f'\n Kruskal test: chi = {statistic}, p = {round(p,4)} \n')
            f.write(f'ACROSS ADULT CONDITIONS COMPARISONS - alpha={0.05/3} \n')
            combinations = [('video','eyes_open'),('video','eyes_closed'),('eyes_open','eyes_closed')]

            for combination in combinations:
                u,p = mannwhitneyu(all_tau[all_tau['cond_name']==combination[0]]['tau'],all_tau[all_tau['cond_name']==combination[1]]['tau'])
                d = cohend(all_tau[all_tau['cond_name']==combination[0]]['tau'],all_tau[all_tau['cond_name']==combination[1]]['tau'])
                print(f'{combination[0]} vs {combination[1]} - U = {u}, p = {p}, d = {d}')
                mean1 = np.mean(np.array(all_tau[all_tau['cond_name']==combination[0]]['tau']))
                mean2 = np.mean(np.array(all_tau[all_tau['cond_name']==combination[1]]['tau']))
                sd1 = np.std(np.array(all_tau[all_tau['cond_name']==combination[0]]['tau']))
                sd2 = np.std(np.array(all_tau[all_tau['cond_name']==combination[1]]['tau']))
                f.write(f'###### {combination[0]} (M={mean1}, SD={sd1}) vs {combination[1]} (M={mean2}, SD={sd2})\n')
                f.write(f'Mann Whitney U- U = {u}, p = {p}, d = {d} \n')
        
        
    if linear_model_analysis:
        alpha=0.05
        betas_dict = {}
        rsquared_dict = {}
        adult_data = pd.read_csv(f'/eeg/EEG_timescales/Data/adult_tau_v3.csv')
        adult_data = adult_data[adult_data['pernan']<0.25]
        adult_data.loc[adult_data.event=='2222','event'] = 'video'
        channels=len(adult_data['channel'].unique())
        adult_data = adult_data.rename(columns={"tauimp": "tau"})
        infant_data_all = pd.read_csv(f'/eeg/EEG_timescales/Data/develop_tau{group_flag}_v3.csv')
        infant_data_all = infant_data_all[(infant_data_all['pernan']<0.25) & (infant_data_all['event']!= 9999)]
        infant_data_all = infant_data_all.rename(columns={"tauimp": "tau"})       
        for t in timepoints:
            if t!='adults':
                betas_timepoint = {}
                rsquared_timepoint = {}
                infant_data = infant_data_all[infant_data_all['ses_age']==t]
                for cond in adult_conditions:
                    adult_data_cond=adult_data[adult_data['event']==cond]

                    betas_cond,rsquared = linear_model_pipeline.compute_betas_and_rsquared(adult_data_cond,infant_data,channels)
                    betas_timepoint[cond] = betas_cond
                    rsquared_timepoint[cond] = rsquared

                    s,p_norm=shapiro(np.array(betas_cond))

                    if p_norm<0.05:
                        ttest,p=mannwhitneyu(np.array(betas_cond), np.zeros(np.array(betas_cond).shape[0]))
                    else:
                        ttest,p=ttest_1samp(np.array(betas_cond),0)

                    sns.displot(np.array(betas_cond))
                    plt.suptitle(f'betas for adult {cond} vs infant {t} \n normality: s={round(s,3)}, p={round(p_norm,4)} - ttest={round(ttest,2)}, p={round(p,4)}')
                    plt.savefig(f'/eeg/EEG_timescales/Figure/nobootstrapping_betas_distribution_{cond}_{t}.png')
                    plt.close()

                betas_dict[t] = betas_timepoint
                rsquared_dict[t] = rsquared_timepoint

        with open (f'/eeg/EEG_timescales/Results/betas_linear_model_adult_simpleregression{group_flag}.pickle','wb') as handle:
            pickle.dump(betas_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
        month_column = np.concatenate((np.repeat(list(betas_dict.keys())[0],betas_dict[list(betas_dict.keys())[0]]['video'].shape[0]*3),np.repeat(list(betas_dict.keys())[1],betas_dict[list(betas_dict.keys())[1]]['video'].shape[0]*3),np.repeat(list(betas_dict.keys())[2],betas_dict[list(betas_dict.keys())[2]]['video'].shape[0]*3)))
        cond_column = np.concatenate((np.repeat(['video','eyes_open','eyes_closed'],betas_dict[list(betas_dict.keys())[0]]['video'].shape[0]),np.repeat(['video','eyes_open','eyes_closed'],betas_dict[list(betas_dict.keys())[1]]['video'].shape[0]),np.repeat(['video','eyes_open','eyes_closed'],betas_dict[list(betas_dict.keys())[2]]['video'].shape[0])))
        betas_column = np.concatenate((betas_dict[list(betas_dict.keys())[0]]['video'],betas_dict[list(betas_dict.keys())[0]]['eyes_open'],betas_dict[list(betas_dict.keys())[0]]['eyes_closed'],betas_dict[list(betas_dict.keys())[1]]['video'],betas_dict[list(betas_dict.keys())[1]]['eyes_open'],betas_dict[list(betas_dict.keys())[1]]['eyes_closed'],betas_dict[list(betas_dict.keys())[2]]['video'],betas_dict[list(betas_dict.keys())[2]]['eyes_open'],betas_dict[list(betas_dict.keys())[2]]['eyes_closed']))
        rsquared_column = np.concatenate((rsquared_dict[list(rsquared_dict.keys())[0]]['video'],rsquared_dict[list(rsquared_dict.keys())[0]]['eyes_open'],rsquared_dict[list(rsquared_dict.keys())[0]]['eyes_closed'],rsquared_dict[list(rsquared_dict.keys())[1]]['video'],rsquared_dict[list(rsquared_dict.keys())[1]]['eyes_open'],rsquared_dict[list(rsquared_dict.keys())[1]]['eyes_closed'],rsquared_dict[list(rsquared_dict.keys())[2]]['video'],rsquared_dict[list(rsquared_dict.keys())[2]]['eyes_open'],rsquared_dict[list(rsquared_dict.keys())[2]]['eyes_closed']))
        betas_df_dict = {'month': month_column,
                        'adult_cond': cond_column,
                        'beta': betas_column,
                        'rsquared':rsquared_column}
        
        betas_df = pd.DataFrame(betas_df_dict)
        betas_df.to_csv(f'/eeg/EEG_timescales/Results/betas_linear_model_adult_simpleregression{group_flag}.csv')

        linear_model_pipeline.compute_stats(betas_dict,rsquared_dict,alpha,group_flag)                   

