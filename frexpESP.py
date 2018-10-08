import numpy as np
import pandas as pd
import numbers
import scipy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import Lasso
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from phi import *
import random
import os
import math
from helpers import *
import sys


from frexp import good_dict
os.system('python frexp.py')

from phi import *
ex_0=None;ex_1=None;ei_0=None;ei_1=None;ei_2=None;ei_3=None;ei_4=None;ei_6=None;ei_5=None;ei_7=None;ei_9=None;ei_8=None;ei_10=None;ei_11=None;e_0=None;e_1=None;e_2=None;e_3=None;e_4=None;e_5=None;e_6=None;e_7=None;f_0=None;f_2=None;f_1=None;f_3=None;f_5=None;f_4=None;f_6=None;f_7=None;x_0=None;x_1=None


M_LN2=0.69314718055994530941723212146 
DBL_MIN_EXP=-1021 
def gsl_finite(x):
    x_0 = x;
    lo = locals()
    record_locals(lo, test_counter)
    return np.isfinite(x_0)

def gsl_frexp(x,e):
    x_1 = x;e_0 = e;
    ex_0=None;ex_1=None;ei_0=None;ei_1=None;ei_2=None;ei_3=None;ei_4=None;ei_6=None;ei_5=None;ei_7=None;ei_9=None;ei_8=None;ei_10=None;ei_11=None;e_1=None;e_2=None;e_3=None;e_4=None;e_5=None;e_6=None;e_7=None;f_0=None;f_2=None;f_1=None;f_3=None;f_5=None;f_4=None;f_6=None;f_7=None;

    gen_bad = random() < probability
    if x_1==0.0:
        e_1=0 
        lo = locals()
        record_locals(lo, test_counter)
        return 0.0,e_1
    elif  not gsl_finite(x_1):
        e_2=0
        lo = locals()
        record_locals(lo, test_counter) 
        return x_1,e_2
    elif abs(x_1)>=0.5 and abs(x_1)<1:
        e_3=0 
        lo = locals()
        record_locals(lo, test_counter)
        return x_1,e_3
    else:
        ex_0=math.ceil(math.log(abs(x_1))/M_LN2) 
        ei_0=ex_0 
        if ei_0<DBL_MIN_EXP:
            ei_1=DBL_MIN_EXP 
        phiPreds = [ei_0<DBL_MIN_EXP]
        phiNames = [ei_1,ei_0]
        ei_2= phiIf(phiPreds, phiNames)
        if ei_2>-DBL_MIN_EXP:
            ei_3=-DBL_MIN_EXP 
        phiPreds = [ei_2>-DBL_MIN_EXP]
        phiNames = [ei_3,ei_2]
        ei_4= phiIf(phiPreds, phiNames)
        f_0=fuzzy(x_1*pow(2.0,-ei_4), gen_bad)
        if  not gsl_finite(f_0):
            e_4=0 
            lo = locals()
            record_locals(lo, test_counter)
            return f_0,e_4
        phiPreds = [ not gsl_finite(f_0)]
        phiNames = [e_4,e_0]
        e_5= phiIf(phiPreds, phiNames)
        phi0 = Phi()
        while abs(phi0.phiLoopTest(f_0,f_1))>=1.0:
            phi0.set()
            ei_6 = phi0.phiEntry(ei_4,ei_5)
            f_2 = phi0.phiEntry(f_0,f_1)

            ei_5 = ei_6+1
            f_1 = f_2/2.0
        ei_7 = phi0.phiExit(ei_4,ei_5)
        f_3 = phi0.phiExit(f_0,f_1)
        phi0 = Phi()
        while abs(phi0.phiLoopTest(f_3,f_4))>0 and abs(phi0.phiLoopTest(f_3,f_4))<0.5:
            phi0.set()
            ei_9 = phi0.phiEntry(ei_7,ei_8)
            f_5 = phi0.phiEntry(f_3,f_4)

            ei_8 = ei_9-1
            f_4 = f_5*2.0
        ei_10 = phi0.phiExit(ei_7,ei_8)
        f_6 = phi0.phiExit(f_3,f_4)
        e_6=ei_10 
        lo = locals()
        record_locals(lo, test_counter)
        return f_6,e_6
    phiPreds = [x_1==0.0, not gsl_finite(x_1),abs(x_1)>=0.5 and abs(x_1)<1]
    phiNames = [None,None,None,ex_0]
    ex_1= phiIf(phiPreds, phiNames)
    phiPreds = [x_1==0.0, not gsl_finite(x_1),abs(x_1)>=0.5 and abs(x_1)<1]
    phiNames = [None,None,None,ei_10]
    ei_11= phiIf(phiPreds, phiNames)
    phiPreds = [x_1==0.0, not gsl_finite(x_1),abs(x_1)>=0.5 and abs(x_1)<1]
    phiNames = [e_1,e_2,e_3,e_6]
    e_7= phiIf(phiPreds, phiNames)
    phiPreds = [x_1==0.0, not gsl_finite(x_1),abs(x_1)>=0.5 and abs(x_1)<1]
    phiNames = [None,None,None,f_6]
    f_7= phiIf(phiPreds, phiNames)



#generate python causal map
causal_map = {'ex_1':['ex_0'],'ex_0':['x_1'],'e_1':[],'f_0':['x_1','ei_4'],'e_3':[],'f_2':['f_0','f_1'],'e_2':[],'f_1':['f_2'],'e_5':['e_4','e_0'],'f_4':['f_5'],'ei_10':['ei_7','ei_8'],'e_4':[],'f_3':['f_0','f_1'],'ei_11':['ei_10'],'f_6':['f_3','f_4'],'e_7':['e_1','e_2','e_3','e_6'],'f_5':['f_3','f_4'],'e_6':['ei_10'],'f_7':['f_6'],'ei_9':['ei_7','ei_8'],'ei_8':['ei_9'],'ei_7':['ei_4','ei_5'],'ei_6':['ei_4','ei_5'],'ei_5':['ei_6'],'ei_4':['ei_3','ei_2'],'ei_3':[],'ei_2':['ei_1','ei_0'],'ei_1':[],'ei_0':['ex_0'],}

#added phi names
phi_names_set = {'ei_2','ei_4','e_5','ei_6','f_2','ei_7','f_3','ei_9','f_5','ei_10','f_6','ex_1','ei_11','e_7','f_7',}

#------end of program---------------------------
def record_locals(lo, i):
    for name in lo:
        if '_IV' in name:
            continue
        if isinstance(lo[name], numbers.Number) and name in causal_map:
            if name not in global_value_dict:
                columns = causal_map[name].copy()
                columns.insert(0, name)
                global_value_dict[name] = pd.DataFrame(columns=columns)
            new_row = [np.float64(lo[name])]

            for pa in causal_map[name]:
                if isinstance(lo[pa], numbers.Number):
                    new_row.append(np.float64(lo[pa]))
                else:
                    new_row.append(lo[pa])
            global_value_dict[name].loc[i] = new_row


bad_dict = {}
global_value_dict = {}
arg1s = np.arange(0, 1000)
test_counter = 0
insertion_count = 0
probability = float(sys.argv[1])/100.0
for arg1 in arg1s:
    e = 0.0
    bad_outcome = gsl_frexp(arg1, e)
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
        importance_p = 2 / ((1 / increase_p) + (1/(math.log(F_p + 0.00001) / math.log(total_failed + 0.00001))))
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