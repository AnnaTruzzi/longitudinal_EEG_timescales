import glob
import random
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def get_subj_list(datapth):
    sub_list = []
    for filename in glob.glob(os.path.join(datapth,'*.set')):
        subname = filename.split('/')[-1].split('_')[0]
        sub_list.append(subname)

    return sub_list


datapth = '/eeg/ITS_Project_NewPreproc/6mo/MADE_APICE_Combined_Pipelines/Narrow_Filter/10s_Epochs_NoOver'
sub_list = get_subj_list(datapth)
random.seed(23)

half_sample = random.sample(sub_list, len(sub_list)//2)
replication_sample = [item for item in sub_list if item not in half_sample]
subj_dict = {'half':half_sample,
            'rep':replication_sample}
subj_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in subj_dict.items()]))
subj_df.to_csv('/eeg/EEG_timescales/Data/half_and_replication_subj_EEGdata.csv')
