import numpy as np
import pandas as pd
import numbers
import scipy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from phi import *
import random
import os
import math
from helpers import *
import sys

from skewness import good_dict
from skewness import args0
os.system('python skewness.py')


from phi import *

import math
import numpy as np
def gsl_stats_mean(data,stride,size):
    data_0 = data;stride_0 = stride;size_0 = size;
    mean_0=None;mean_2=None;mean_1=None;mean_3=None;

    gen_bad = random() < probability
    mean_0=0 
    phi0 = Phi()
    for i_0 in range(0,size_0):
        phi0.set()
        mean_2 = phi0.phiEntry(mean_0,mean_1)

        mean_1 = fuzzy(mean_2+(data_0[i_0*stride_0]-mean_2)/(i_0+1), gen_bad)
    mean_3 = phi0.phiExit(mean_0,mean_1)
    lo = locals()
    record_locals(lo, test_counter)
    return mean_3

def compute_variance(data,stride,n,mean):
    data_1 = data;stride_1 = stride;n_0 = n;mean_4 = mean;
    variance_0=None;variance_2=None;variance_1=None;variance_3=None;delta_1=None;delta_0=None;delta_2=None;

    variance_0=0 
    phi0 = Phi()
    for i_1 in range(0,n_0):
        phi0.set()
        variance_2 = phi0.phiEntry(variance_0,variance_1)
        delta_1 = phi0.phiEntry(None,delta_0)

        delta_0=(data_1[i_1*stride_1]-mean_4) 
        variance_1 = variance_2+(delta_0*delta_0-variance_2)/(i_1+1)
    variance_3 = phi0.phiExit(variance_0,variance_1)
    delta_2 = phi0.phiExit(None,delta_0)
    lo = locals()
    record_locals(lo, test_counter)
    return variance_3

def gsl_stats_sd_m(data,stride,n,mean):
    data_2 = data;stride_2 = stride;n_1 = n;mean_5 = mean;
    sd_0=None;variance_4=None;

    variance_4=compute_variance(data_2,stride_2,n_1,mean_5) 
    sd_0=math.sqrt(variance_4*(n_1/(n_1-1)))
    lo = locals()
    record_locals(lo, test_counter)
    return sd_0

def gsl_stats_skew_m_sd(data,stride,n,mean,sd):
    data_3 = data;stride_3 = stride;n_2 = n;mean_6 = mean;sd_1 = sd;
    x_1=None;x_0=None;x_2=None;skew_0=None;skew_2=None;skew_1=None;skew_3=None;

    skew_0=0 
    phi0 = Phi()
    for i_2 in range(0,n_2):
        phi0.set()
        x_1 = phi0.phiEntry(None,x_0)
        skew_2 = phi0.phiEntry(skew_0,skew_1)

        x_0=(data_3[i_2*stride_3]-mean_6)/sd_1
        skew_1 = skew_2+(x_0*x_0*x_0-skew_2)/(i_2+1) 
    x_2 = phi0.phiExit(None,x_0)
    skew_3 = phi0.phiExit(skew_0,skew_1)
    lo = locals()
    record_locals(lo, test_counter)
    return skew_3

def gsl_stats_skew(data,stride,n):
    data_4 = data;stride_4 = stride;n_3 = n;
    sd_2=None;mean_7=None;skewness_0=None;

    mean_7=gsl_stats_mean(data_4,stride_4,n_3) 
    sd_2=gsl_stats_sd_m(data_4,stride_4,n_3,mean_7) 
    skewness_0 = gsl_stats_skew_m_sd(data_4,stride_4,n_3,mean_7,sd_2)
    lo = locals()
    record_locals(lo, test_counter)
    return skewness_0



#generate python causal map
causal_map = dict(mean_7=['data_4','stride_4','n_3'],mean_3=['mean_0','mean_1'],mean_2=['mean_0','mean_1'],x_0=['data_3','i_2','stride_3','mean_6','sd_1'],mean_1=['mean_2','data_0','i_0','stride_0','mean_2','i_0'],x_2=['x_0'],mean_0=[],skew_0=[],x_1=['x_0'],skew_1=['skew_2','x_0','x_0','x_0','skew_2','i_2'],skew_2=['skew_0','skew_1'],skew_3=['skew_0','skew_1'],sd_2=['data_4','stride_4','n_3','mean_7'],sd_0=['variance_4','n_1','n_1'],delta_0=['data_1','i_1','stride_1','mean_4'],variance_4=['data_2','stride_2','n_1','mean_5'],delta_1=['delta_0'],variance_3=['variance_0','variance_1'],variance_2=['variance_0','variance_1'],delta_2=['delta_0'],variance_1=['variance_2','delta_0','delta_0','variance_2','i_1'],skewness_0=['data_4','stride_4','n_3','mean_7','sd_2'],variance_0=[],)

#added phi names
phi_names_set = {'mean_2','mean_3','variance_2','delta_1','variance_3','delta_2','x_1','skew_2','x_2','skew_3',}




#---------end of the faulty program------------
def record_locals(lo, i):
    for name in lo:
        if '_IV' in name:
            continue
        if isinstance(lo[name], numbers.Number) and name in causal_map:
            if name not in global_value_dict:
                columns = [name]
                global_value_dict[name] = pd.DataFrame(columns=columns)
            new_row = [np.float64(lo[name])]
            global_value_dict[name].loc[i] = new_row
            


global_value_dict = {}

bad_dict = {}
global_value_dict = {}
test_counter = 0
args1 = args0
insertion_count = 0
probability = float(sys.argv[1])/100.0
for arg1 in args1:
    sk = gsl_stats_skew(arg1, 1, len(arg1))
    bad_dict[test_counter] = sk
    test_counter += 1

diff_dict = {index : 0.0 if bad_dict[index] == good_dict[index] else 1.0 for index in bad_dict }
total_failed = sum(1 for index in diff_dict if diff_dict[index] == 1.0)


def label_predicate(df):
    if df[key] == mean:
        label = 0
    if (df[key] < mean) and (df[key] >= mean - sd):
        label = -1
    if (df[key] < mean - sd) and (df[key] >= mean - 2 * sd):
        label = -2
    if (df[key] < mean - 2 * sd) and (df[key] >= mean - 3 * sd):
        label = -3
    if (df[key] < mean - 3 * sd):
        label = -4
    if (df[key] > mean) and (df[key] <= mean + sd):
        label = 1
    if (df[key] > mean + sd) and (df[key] <= mean + 2 * sd):
        label = 2
    if (df[key] > mean + 2 * sd) and (df[key] <= mean + 3 * sd):
        label = 3
    if (df[key] > mean + 3 * sd):
        label = 4
    return label

for key in global_value_dict:
    df = global_value_dict[key]
    rows = df.index
    outcome_list = [diff_dict[i] for i in rows]
    df['outcome'] = outcome_list
    sd = df[key].std()
    mean = df[key].mean()

    df['label'] = df.apply(label_predicate, axis = 1)



suspicious_df = pd.DataFrame(columns=['variable_name', 'importance_score', 'p_label'])

for key in global_value_dict:
    df = global_value_dict[key]
    grouped = df.groupby('label')
    group_dict = grouped.groups
    max_importance = 0
    #initalization
    p_label_max = -5
    
    F_p_obs = df['outcome'].value_counts()[1.0] if 1.0 in df['outcome'].value_counts() else 0
    S_p_obs = df['outcome'].value_counts()[0.0] if 0.0 in df['outcome'].value_counts() else 0
    for p_label in group_dict:
        #key-->label  value-->indexed matches label
        matched_rows = group_dict[p_label]
        df_p = df.loc[matched_rows]
        F_p = df_p['outcome'].value_counts()[1.0] if 1.0 in df_p['outcome'].value_counts() else 0
        S_p = df_p['outcome'].value_counts()[0.0] if 0.0 in df_p['outcome'].value_counts() else 0
        increase_p = F_p/(S_p + F_p) + F_p_obs/(S_p_obs + F_p_obs)
        importance_p = 2 / ((1 / increase_p) + (1/(math.log(F_p + 0.0001) / math.log(total_failed + 0.0001))))
        if key == 't_1':
            print(key, p_label, increase_p, importance_p, F_p, S_p, F_p_obs, S_p_obs)
        if importance_p > max_importance:
            max_importance = importance_p
            p_label_max = p_label
    row = [key, max_importance, p_label_max]
    suspicious_df.loc[len(suspicious_df)] = row

def filter_phi_rows(suspicious_df, phi_names_set):
    return suspicious_df[~suspicious_df['variable_name'].isin(phi_names_set)]

suspicious_df = suspicious_df.sort_values(by='importance_score', ascending=False)

suspicious_final_rank = filter_phi_rows(suspicious_df, phi_names_set)
print('*************Target variables in total: ', len(suspicious_final_rank),'*************')
print(suspicious_final_rank)
        

    
with open(os.path.basename(__file__)[:-3] + "-" + sys.argv[1] + "-Trial" + sys.argv[2] + ".txt", "w") as f:
    f.write('*************Target variables in total: ' + str(len(suspicious_final_rank)) + '*************\n')
    f.write(str(suspicious_final_rank.to_csv()))