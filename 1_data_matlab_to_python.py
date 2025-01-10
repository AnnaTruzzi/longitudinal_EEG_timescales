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
import scipy 

def remap_epochs_labels(events,event_id):
    for i,epoch in enumerate(events):
        label = epoch[2]
        new_label = list(event_id.keys())[list(event_id.values()).index(label)]
        if '/' in new_label:
            new_label = new_label.split('/')[0]
        events[i][2] = new_label
    return events


rootpth = 'C:\\Users\\Anna\\Documents\\Research\\Projects\\ONGOING\\Project dHCP_Autocorrelation\\EEG_autocorrelations'
eegfile = mne.io.read_epochs_eeglab(os.path.join(rootpth,'Data','B003_6mo_processed_data.set'))
#eegfile.resample(200)
#x = data[10,:]
data = eegfile.get_data()
events_list = [epoch[2] for epoch in eegfile.events]
remapped_events = remap_epochs_labels(eegfile.events, eegfile.event_id)
#events1list_idx = np.where(np.array(events_list) == 1)

x = data[10,10,:]
fig,ax = plt.subplots(nrows=2, ncols=1, figsize=(12,8))
plot_acf(x,ax=ax[0],lags=100,title=None,zero=False)
ax[0].set_ylabel('ACF')
plot_pacf(x,ax=ax[1],lags=100,title=None,zero=False)
ax[1].set_ylabel('PACF')
plt.suptitle('ACF & PACF EEG data')
plt.savefig(os.path.join(rootpth,'Results','ACF_PACF_EEGsignal.png'))
#plt.show()
plt.close()


mod = ARIMA(endog=x, order=(1,0,1) ,enforce_stationarity=False)
m = mod.fit()
res = m.resid
print(scipy.stats.shapiro(res))
res_lbtest = acorr_ljungbox(res,lags=[100],model_df=1+1, return_df=True)
print(res_lbtest)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
ax.plot(x)
ax.set_ylabel('Signal')
ax.plot(m.fittedvalues)
#plt.show()
plt.close()
