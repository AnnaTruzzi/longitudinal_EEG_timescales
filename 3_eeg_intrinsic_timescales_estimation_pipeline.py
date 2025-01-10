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


if __name__ == '__main__':

    #group='base' # or 'validation'
    group=['base','validation']

    timepoints = ['6mo','9mo','16mo','adults']
    adult_conditions = ['video','eyes_open','eyes_closed']
    rootpth = '/eeg/ITS_Project_NewPreproc'
    nlags = 1000
    removelag0=1
    
    infant_sample_df = pd.read_csv('/eeg/EEG_timescales/Data/Datasets_per_group.csv')

    for g in group:
        if g=='base':
            half_sample = list(infant_sample_df['half'])
            group_flag = ''
        else:
            half_sample = list(infant_sample_df['rep'])
            group_flag = '_validation'

        half_sample = [sub.split('_')[0] for sub in half_sample]
        channels = 109

        estimate_timescales=False

        if estimate_timescales:
            for t in timepoints:
                print(f'########### RUNNING MODEL FOR {t}')
                datapth = f'/eeg/ITS_Project_NewPreproc/{t}/10s/'
                resultspth = f'/eeg/EEG_timescales/Results/{t}/10s'
                if not os.path.isdir(resultspth):
                    os.mkdir(resultspth)
                sub_list = tau_estimate.get_subj_list(datapth,t)

                tau_estimate.estimate_eeg_tau(sub_list,t,half_sample,datapth,channels,nlags,removelag0,resultspth,group_flag)
            
            

        for t in timepoints:
            print(f'{g} - {t}')
            resultspth = f'/eeg/EEG_timescales/Results/{t}/10s'
            tau = pd.read_csv(os.path.join(resultspth,f'EEGdata_tau_{t}_10s_numpy{nlags}_nooutliers{group_flag}.csv'))
        
            epochs_all=tau.shape[0]
            print(f'all epochs={epochs_all}')
            sample_all=len(tau['sub'].unique())
            print(sample_all)

            # Exclude 9999 epochs (too much movement)
            tau = tau[tau['cond_name']!='movement']
            low_movement_epochs = tau.shape[0]
            print(f'low movement epochs={low_movement_epochs}')
            sample_low_movement=len(tau['sub'].unique())

            # Exclude epochs with too many nans
            highnan_idx = []
            for i,row in tau.iterrows():
                tau_row = row[4:]
                nan_percentage = (np.count_nonzero(pd.isnull(tau_row))/tau_row.shape[0])*100
                if nan_percentage>25:
                    highnan_idx.append(i)
            tau_nohighnan = tau.drop(highnan_idx)
            tau_nohighnan.to_csv(os.path.join(resultspth,f'EEGdata_tau_{t}_10s_numpy{nlags}_nooutliers_nohighnan{group_flag}.csv'))
            
            epochs_nohighnan=tau_nohighnan.shape[0]
            print(f'no high nan epochs={epochs_nohighnan}')
            sample_nohighnan=len(tau_nohighnan['sub'].unique())
            print(sample_nohighnan)

            ch_list=[str(ch) for ch in np.arange(0,109)]
            tau_nohighnan_long_data=pd.melt(tau_nohighnan, id_vars=['sub','cond_name'], value_vars=ch_list)

            tau_nohighnan_long_data.columns=['sub','cond_name','channel','tau']

            tau_nohighnan_long_data.to_csv(f'/eeg/EEG_timescales/Results/{t}/10s/EEGdata_tau_{t}_10s_numpy1000_longformat_seconds_nooutliers_lowmovement_nohighnans{group_flag}.csv')

            # summarize epochs with too many nans
            samplesize_filename = f'/eeg/EEG_timescales/Results/samplesize_10s_numpy{nlags}_nooutliers.txt'

            with open(samplesize_filename, 'a') as f:
                if group_flag=='':
                    group='control'
                else:
                    group='validation'
                f.write(f'### Age: {t} - Group: {group} \n')
                f.write(f'all sample = {sample_all} - all epochs = {epochs_all} \n')
                f.write(f'low movement sample = {sample_low_movement} - low movement epochs = {low_movement_epochs} \n')
                f.write(f'nohighnans sample = {sample_nohighnan} - nohighnans epochs = {epochs_nohighnan} \n')
                f.close()
            
            # plot Nans
            tau_array = np.array(tau.iloc[:,4:])
            nan_total = np.count_nonzero(np.isnan(tau_array))
            nan_total_perc = (nan_total*100)/(tau_array.shape[0]*tau_array.shape[1])
            nan_per_electrode = np.zeros((109))
            for electrode in range(0,tau_array.shape[1]):
                nan_count = np.count_nonzero(np.isnan(tau_array[:,electrode]))
                nan_percentage = (nan_count*100)/(tau_array[:,electrode].shape[0])
                nan_per_electrode[electrode] = nan_percentage

            electrodes_area = pd.read_csv('/eeg/EEG_timescales/Data/electrodes_area.csv',header=None,index_col=None,)
            #areas = electrodes_area.iloc[:,0].unique()
            areas = ['Central','Frontal','Frontalpole','Occipital','Parietal','Temporal']

            areas_idx = {}
            nan_reordered = []
            colours = ["#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352"]
            sd_bars = True
            my_plot_colours = []
            for n,area in enumerate(areas):
                idx=np.where(electrodes_area.iloc[:,0]==area)[0]
                nan_reordered.extend(nan_per_electrode[idx])
                areas_idx[area]=idx
                print(f'{area} has {len(idx)} electrodes')
                my_plot_colours.extend(np.repeat(colours[n],len(idx)))

            if sd_bars:
                plt.figure(figsize=(12, 6))
                plt.bar(range(0,tau_array.shape[1]),np.array(nan_reordered),color=my_plot_colours,label=[areas], yerr=np.nanstd(nan_reordered), capsize=0.5, linewidth=0.2)
                flag = '_withsdbars'
            else:
                plt.bar(range(0,tau_array.shape[1]),np.array(nan_reordered),color=my_plot_colours,label=[areas])
                flag = ''
            plt.xlim((0,110))
            plt.ylim((0,13))
            #plt.legend(areas)
            print(areas)
            plt.suptitle(f'Nan percentage per electrode in {t} - {g} \n Total nan % = {nan_total_perc}, Max % = {np.max(nan_per_electrode)}')
            plt.savefig(f'/eeg/EEG_timescales/Figure/nan_count_{t}_{group_flag}{flag}.png')
            plt.savefig(f'/eeg/EEG_timescales/Figure/nan_count_{t}_{group_flag}{flag}.pdf')
            plt.close()

        ## THE TAU DATASETS ARE THEN PASSED ON TO MATLAB TO ESTIMATE NANS INTERPOLATING FROM NEARBY ELECTRODES AND CREATING THE DATAFRAMES
        ## SHARED IN THE DATA FOLDER. THESE DATASETS CREATED IN MATLAB ARE THE ONES USED TO RUN THE ANALYSIS.

