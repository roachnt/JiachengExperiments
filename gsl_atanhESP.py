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

from gsl_atanh import good_dict
os.system('python gsl_atanh.py')
from phi import *
a_0=None;s_0=None;to_return_0=None;to_return_1=None;to_return_2=None;to_return_3=None;to_return_4=None;x_0=None;x_1=None;y_0=None;z_0=None


GSL_DBL_EPSILON=2.2204460492503131e-16 
GSL_POSINF=math.inf 
GSL_NEGINF=-math.inf 
GSL_NAN=math.nan 
def gsl_log1p(x):
    x_0 = x;
    to_return_0=None;y_0=None;z_0=None;

    y_0=1+x_0 
    z_0=y_0-1 
    to_return_0=math.log(y_0)-(z_0-x_0)/y_0 + bug
    lo = locals()
    record_locals(lo, test_counter)
    return to_return_0

def gsl_atanh(x):
    x_1 = x;
    a_0=None;s_0=None;to_return_1=None;to_return_2=None;to_return_3=None;to_return_4=None;

    a_0=abs(x_1) 
    s_0=-1 if x_1<0 else 1 
    if a_0>1:
        lo = locals()
        record_locals(lo, test_counter)
        return GSL_NAN
    elif a_0==1:
        to_return_1=GSL_NEGINF if x_1<0 else GSL_POSINF 
        lo = locals()
        record_locals(lo, test_counter)
        return to_return_1
    elif a_0>=0.5:
        to_return_2=s_0*0.5*gsl_log1p(2*a_0/(1-a_0)) 
        lo = locals()
        record_locals(lo, test_counter)
        return to_return_2
    elif a_0>GSL_DBL_EPSILON:
        to_return_3=s_0*0.5*gsl_log1p(2*a_0+2*a_0*a_0/(1-a_0)) 
        lo = locals()
        record_locals(lo, test_counter)
        return to_return_3
    else:
        lo = locals()
        record_locals(lo, test_counter)
        return x_1
    phiPreds = [a_0>1,a_0==1,a_0>=0.5,a_0>GSL_DBL_EPSILON]
    phiNames = [None,to_return_1,to_return_2,to_return_3,None]
    to_return_4= phiIf(phiPreds, phiNames)



#generate python causal map
causal_map = {'a_0':['x_1'],'to_return_0':['y_0','z_0','x_0','y_0'],'s_0':['x_1'],'to_return_2':['s_0','a_0','a_0'],'to_return_1':['x_1'],'to_return_4':['to_return_1','to_return_2','to_return_3'],'to_return_3':['s_0','a_0','a_0','a_0','a_0'],'z_0':['y_0'],'y_0':['x_0'],}

#added phi names
phi_names_set = {'to_return_4',}


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
            

def fluky(good_val, bad_val, p):
        r = random.random()
        if r <= p:
            return bad_val
        else:
            return good_val


bad_dict = {}
global_value_dict = {}
test_counter = 0
arg1s = np.arange(0, 1, 0.001)
bug = 0
for arg1 in arg1s:
    bug = fluky(0, 0.138 , 0.95)
    bad_outcome = gsl_atanh(arg1)
    bad_dict[test_counter] = bad_outcome
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
        importance_p = 2 / ((1 / increase_p) + (1/(math.log(F_p + 0.00001) / math.log(total_failed+0.00001))))
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