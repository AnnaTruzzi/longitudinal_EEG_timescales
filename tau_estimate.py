import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mne
import os
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.graphics.tsaplots import plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.diagnostic import het_breuschpagan
from scipy import optimize
import scipy 
import glob
import random
import pickle


def get_subj_list(datapth,t):
    sub_list = []
    for filename in glob.glob(os.path.join(datapth,'*.set')):
        if t!='adults':
            subname = filename.split('/')[-1].split('_')[0]
            sub_list.append(subname)
        else:
            subname = filename.split('/')[-1]
            sub_list.append(subname)

    return sub_list


def remap_epochs_labels(events,event_id):
    for i,epoch in enumerate(events):
        label = epoch[2]
        new_label = list(event_id.keys())[list(event_id.values()).index(label)]
        if '/' in new_label:
            new_label = new_label.split('/')[0]
        if '.' in new_label:
            new_label = new_label.split('.')[0]
        events[i][2] = new_label
        events_out = [epoch[2] for epoch in events]
    return events_out


def get_events_list(datapth,sub_list):
    allevents_list = []
    sub_epochs_list = []
    for sub in sub_list:
        eegfile = mne.io.read_epochs_eeglab(os.path.join(datapth,sub+'.set'))
        events = [epoch[2] for epoch in eegfile.events]
        remap_events = remap_epochs_labels(eegfile.events, eegfile.event_id)
        allevents_list.extend(remap_events)
        sub_epochs_list.extend([sub+'_'+event for event in remap_events])
    
    return allevents_list,sub_epochs_list


def load_eeg_data(datapth,sub, remap_labels=False):
    eegfile = mne.io.read_epochs_eeglab(os.path.join(datapth,sub))
    data = eegfile.get_data()
    events = [epoch[2] for epoch in eegfile.events]
    if remap_labels:
        events = remap_epochs_labels(eegfile.events, eegfile.event_id)

    return data,events


def autocorr_decay(dk,A,tau,B):
    return A*(np.exp(-(dk/tau))+B)


def estimate_eeg_tau(sub_list, timepoint, half_sample,datapth,channels,nlags,removelag0,resultspth,group_flag):
    events_tot = 0
    sub_epoch_list = [[],[]]
    alltau_list = []
    for sub in sub_list:
        if timepoint=='adults' or sub in half_sample:
            if (sub=='B186' or sub=='B602') and timepoint=='6mo':
                sub_name = f'{sub}_6m2.set'
            elif (sub=='B079' or sub=='B100' or sub=='B194') and timepoint=='9mo':
                sub_name = f'{sub}_9m2.set'
            elif timepoint=='adults':
                sub_name = sub
            else:
                sub_name = f'{sub}_{timepoint}.set'
            
            print(sub)
            data,events = load_eeg_data(datapth, sub_name, remap_labels=True)
            sub_name_tot = np.repeat(sub,len(events))
            sub_epoch_list[0].extend(sub_name_tot)
            sub_epoch_list[1].extend(events)
            subtau = np.zeros((len(events),channels))
            for e,epoch in enumerate(events):
                events_tot+=1
                xdata=np.arange(nlags)
                autocorr_values = np.zeros((channels, nlags))
                for c in range(0,channels):
                    xc=data[e,c,:]-np.mean(data[e,c,:])
                    fullcorr=np.correlate(xc, xc, mode='full')
                    fullcorr=fullcorr / np.max(fullcorr)
                    start=len(fullcorr) // 2
                    stop=start+nlags
                    autocorr_values[c,:]=fullcorr[start:stop] 


                fig, ax = plt.subplots(1, 1, figsize=(10,20))
                for c in range(autocorr_values.shape[0]):
                    ax.plot(range(nlags),autocorr_values[c,:])
                ax.set_xlabel("Lag (TRs)")
                ax.set_ylabel("Autocorrelation")
                plt.suptitle(f"Sub {sub} - Epoch {e}")
                plt.show()
                fig.savefig(f'/eeg/EEG_timescales/Figure/{timepoint}_autocorr/{sub}_{e}_lag{nlags}{group_flag}.png')
                plt.close()

                for c in range(0,channels):
                    if np.all((data[e,c,:] == 0)):
                        subtau[e,c] = np.nan
                    else: 
                        try:
                            A, tau, B = optimize.curve_fit(autocorr_decay,xdata[removelag0:],autocorr_values[c,removelag0:],p0=[0,np.random.rand(1)[0]+0.01,0],bounds=(([0,0,-np.inf],[np.inf,np.inf,np.inf])),method='trf',maxfev=2000)[0]
                            subtau[e,c] = tau
                        except:
                            subtau[e,c] = np.nan


            alltau_list.append(subtau)

    alltau = np.concatenate(alltau_list,axis=0)

    p95 = np.nanpercentile(alltau[0:,0:],95)
    alltau[np.where(alltau>p95)] = np.nan
    alltau[np.where(alltau<0)] = np.nan
    alltau = alltau * 0.001 
    np.savetxt(os.path.join(resultspth,f'EEGdata_tau_{timepoint}_10s_numpy{nlags}{group_flag}_nooutliers.txt'),alltau)


    tau_df=pd.DataFrame(data=alltau[0:,0:],
            index=[i for i in range(alltau.shape[0])],
            columns=[i for i in range(alltau.shape[1])])
    tau_df.insert(loc=0, column='sub', value=sub_epoch_list[0])
    tau_df.insert(loc=1, column='event', value=sub_epoch_list[1])
    tau_df.insert(loc=2,column='cond_name',value=np.repeat('None',len(sub_epoch_list[0])))
    if timepoint=='adults':
        tau_df.loc[tau_df.event=='2222','cond_name'] = 'video'
        tau_df.loc[tau_df.event==3333,'cond_name'] = 'eyes_open'
        tau_df.loc[tau_df.event==4444,'cond_name'] = 'eyes_closed'
    else:
        tau_df.loc[tau_df.event==1111,'cond_name'] = 'bubbles'
        tau_df.loc[tau_df.event==2222,'cond_name'] = 'video'
        tau_df.loc[tau_df.event==9999,'cond_name'] = 'movement'


    tau_df.to_csv(os.path.join(resultspth,f'EEGdata_tau_{timepoint}_10s_numpy{nlags}_nooutliers{group_flag}.csv'))

    ch_list=[int(ch) for ch in np.arange(0,109)]
    long_data=pd.melt(tau_df, id_vars=['sub','cond_name'], value_vars=ch_list)

    long_data.columns=['sub','cond_name','channel','tau']
    long_data.to_csv(f'/eeg/EEG_timescales/Results/{timepoint}/10s/EEGdata_tau_{timepoint}_10s_numpy1000_longformat_seconds_nooutliers{group_flag}.csv')

