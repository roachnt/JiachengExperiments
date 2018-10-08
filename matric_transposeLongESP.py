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

insertion_count = 0

from matric_transposeLong import good_dict
from matric_transposeLong import args0
os.system('python matric_transposeLong.py')


from phi import *
m_block_0=None;n1_0=None;n2_0=None;view_0=None;array_0=None;m_data_0=None;m_size2_0=None;m_size2_1=None;m_owner_0=None;size1_0=None;tmp_5=None;tmp_3=None;tmp_1=None;tmp_0=None;tmp_2=None;tmp_4=None;tmp_6=None;m_size1_0=None;m_size1_1=None;size2_0=None;view_matrix_block_0=None;view_matrix_owner_0=None;view_matrix_data_0=None;view_matrix_0=None;view_matrix_1=None;m_data_e2_6=None;m_data_e2_4=None;m_data_e2_2=None;m_data_e2_0=None;m_data_e2_1=None;m_data_e2_3=None;m_data_e2_5=None;m_data_e2_7=None;i_0=None;j_0=None;k_0=None;e1_5=None;e1_3=None;e1_1=None;e1_0=None;e1_2=None;e1_4=None;e1_6=None;m_0=None;m_1=None;e2_5=None;e2_3=None;e2_1=None;e2_0=None;e2_2=None;e2_4=None;e2_6=None;view_matrix_size1_0=None;view_matrix_size2_0=None;m_tda_0=None;m_tda_7=None;m_tda_5=None;m_tda_3=None;m_tda_1=None;m_tda_2=None;m_tda_4=None;m_tda_6=None;m_tda_8=None;m_data_e1_6=None;m_data_e1_4=None;m_data_e1_2=None;m_data_e1_0=None;m_data_e1_1=None;m_data_e1_3=None;m_data_e1_5=None;m_data_e1_7=None;view_matrix_tda_0=None

import numpy as np
class gsl_block_struct(object):
    __slots__=['size','data'] 
gsl_block=gsl_block_struct 
NULL_MATRIX_VIEW=[[0,0,0,0,0,0]] 
NULL_MATRIX=[0,0,0,0,0,0] 
GSL_SUCCESS=1 
MULTIPLICITY=1 
class gsl_matrix(object):
    __slots__=['size1','size2','tda','data','block','owner'] 
class _gsl_matrix_view(object):
    __slots__=['matrix'] 
gsl_matrix_view=_gsl_matrix_view 
def gsl_matrix_view_array(array,n1,n2):
    array_0 = array;n1_0 = n1;n2_0 = n2;
    m_block_0=None;view_matrix_owner_0=None;view_matrix_data_0=None;view_matrix_0=None;view_matrix_1=None;m_0=None;view_0=None;view_matrix_size1_0=None;view_matrix_size2_0=None;m_tda_0=None;m_data_0=None;m_size2_0=None;m_owner_0=None;m_size1_0=None;view_matrix_tda_0=None;view_matrix_block_0=None;

    view_0=_gsl_matrix_view() 
    view_matrix_0=gsl_matrix() 
    view_0.matrix=view_matrix_0 
    view_matrix_size1_0=0 
    view_0.matrix.size1=view_matrix_size1_0 
    view_matrix_size2_0=0 
    view_0.matrix.size2=view_matrix_size2_0 
    view_matrix_tda_0=0 
    view_0.matrix.tda=view_matrix_tda_0 
    view_matrix_data_0=0 
    view_0.matrix.data=view_matrix_data_0 
    view_matrix_block_0=0 
    view_0.matrix.block=view_matrix_block_0 
    view_matrix_owner_0=0 
    view_0.matrix.owner=view_matrix_owner_0 
    m_0=gsl_matrix() 
    m_0.size1,m_0.size2,m_0.tda,m_0.data,m_0.block,m_0.owner=NULL_MATRIX 
    m_data_0=array_0 
    m_0.data=m_data_0 
    m_size1_0=n1_0 
    m_0.size1=m_size1_0 
    m_size2_0=n2_0 
    m_0.size2=m_size2_0 
    m_tda_0=n2_0 
    m_0.tda=m_tda_0 
    m_block_0=0 
    m_0.block=m_block_0 
    m_owner_0=0 
    m_0.owner=m_owner_0 
    view_matrix_1=m_0 
    view_0.matrix=view_matrix_1 
    lo = locals()
    record_locals(lo, test_counter)
    return view_0

def gsl_matrix_transpose(m):
    m_1 = m;
    m_tda_7=None;m_tda_5=None;m_tda_3=None;m_tda_1=None;m_tda_2=None;m_tda_4=None;m_tda_6=None;m_tda_8=None;size1_0=None;m_size2_1=None;tmp_5=None;tmp_3=None;tmp_1=None;tmp_0=None;tmp_2=None;tmp_4=None;tmp_6=None;m_size1_1=None;size2_0=None;m_data_e1_6=None;m_data_e1_4=None;m_data_e1_2=None;m_data_e1_0=None;m_data_e1_1=None;m_data_e1_3=None;m_data_e1_5=None;m_data_e1_7=None;m_data_e2_6=None;m_data_e2_4=None;m_data_e2_2=None;m_data_e2_0=None;m_data_e2_1=None;m_data_e2_3=None;m_data_e2_5=None;m_data_e2_7=None;e1_5=None;e1_3=None;e1_1=None;e1_0=None;e1_2=None;e1_4=None;e1_6=None;e2_5=None;e2_3=None;e2_1=None;e2_0=None;e2_2=None;e2_4=None;e2_6=None;

    gen_bad = random() < probability
    global insertion_count
    if gen_bad:
        insertion_count += 1
        
    m_size1_1=m_1.size1 
    size1_0=m_size1_1 
    m_size2_1=m_1.size2 
    size2_0=m_size2_1 
    if size1_0!=size2_0:
        print("matrix must be square to take transpose") 
    phi0 = Phi()
    for i_0 in range(0,size1_0):
        phi0.set()
        m_tda_7 = phi0.phiEntry(None,m_tda_6)
        tmp_5 = phi0.phiEntry(None,tmp_4)
        m_data_e1_6 = phi0.phiEntry(None,m_data_e1_5)
        m_data_e2_6 = phi0.phiEntry(None,m_data_e2_5)
        e1_5 = phi0.phiEntry(None,e1_4)
        e2_5 = phi0.phiEntry(None,e2_4)

        phi1 = Phi()
        for j_0 in range(i_0+1,size2_0):
            phi1.set()
            m_tda_5 = phi1.phiEntry(None,m_tda_4)
            tmp_3 = phi1.phiEntry(None,tmp_2)
            m_data_e1_4 = phi1.phiEntry(None,m_data_e1_3)
            m_data_e2_4 = phi1.phiEntry(None,m_data_e2_3)
            e1_3 = phi1.phiEntry(None,e1_2)
            e2_3 = phi1.phiEntry(None,e2_2)

            phi2 = Phi()
            for k_0 in range(0,MULTIPLICITY):
                phi2.set()
                m_tda_3 = phi2.phiEntry(None,m_tda_2)
                tmp_1 = phi2.phiEntry(None,tmp_0)
                m_data_e1_2 = phi2.phiEntry(None,m_data_e1_1)
                m_data_e2_2 = phi2.phiEntry(None,m_data_e2_1)
                e1_1 = phi2.phiEntry(None,e1_0)
                e2_1 = phi2.phiEntry(None,e2_0)

                m_tda_1=m_1.tda 
                e1_0=(i_0*m_tda_1+j_0)*MULTIPLICITY+k_0 
                m_tda_2=m_1.tda
                e2_0=(j_0*m_tda_2+i_0)*MULTIPLICITY+k_0 
                m_data_e1_0=m_1.data[e1_0] 
                tmp_0=m_data_e1_0
                m_data_e2_0=m_1.data[e2_0] 
                m_data_e1_1=m_data_e2_0 
                m_1.data[e1_0]=m_data_e1_1 
                m_data_e2_1=fuzzy(tmp_0, gen_bad)
                m_1.data[e2_0]=m_data_e2_1 
            m_tda_4 = phi2.phiExit(None,m_tda_2)
            tmp_2 = phi2.phiExit(None,tmp_0)
            m_data_e1_3 = phi2.phiExit(None,m_data_e1_1)
            m_data_e2_3 = phi2.phiExit(None,m_data_e2_1)
            e1_2 = phi2.phiExit(None,e1_0)
            e2_2 = phi2.phiExit(None,e2_0)
        m_tda_6 = phi1.phiExit(None,m_tda_4)
        tmp_4 = phi1.phiExit(None,tmp_2)
        m_data_e1_5 = phi1.phiExit(None,m_data_e1_3)
        m_data_e2_5 = phi1.phiExit(None,m_data_e2_3)
        e1_4 = phi1.phiExit(None,e1_2)
        e2_4 = phi1.phiExit(None,e2_2)
    m_tda_8 = phi0.phiExit(None,m_tda_6)
    tmp_6 = phi0.phiExit(None,tmp_4)
    m_data_e1_7 = phi0.phiExit(None,m_data_e1_5)
    m_data_e2_7 = phi0.phiExit(None,m_data_e2_5)
    e1_6 = phi0.phiExit(None,e1_4)
    e2_6 = phi0.phiExit(None,e2_4)
    lo = locals()
    record_locals(lo, test_counter)
    return GSL_SUCCESS



#generate python causal map
causal_map = dict(m_size1_1=['m_1'],view_matrix_1=['m_0'],view_matrix_0=[],m_owner_0=[],e1_5=['e1_4'],view_matrix_data_0=[],e1_6=['e1_4'],e1_3=['e1_2'],e1_4=['e1_2'],view_matrix_size1_0=[],tmp_5=['tmp_4'],tmp_4=['tmp_2'],tmp_6=['tmp_4'],tmp_1=['tmp_0'],m_size1_0=['n1_0'],tmp_0=['m_data_e1_0'],tmp_3=['tmp_2'],tmp_2=['tmp_0'],m_tda_5=['m_tda_4'],m_tda_4=['m_tda_2'],m_tda_7=['m_tda_6'],m_tda_6=['m_tda_4'],m_tda_1=['m_1'],m_tda_0=['n2_0'],m_tda_3=['m_tda_2'],m_tda_2=['m_1'],m_tda_8=['m_tda_6'],e1_1=['e1_0'],m_data_e1_0=['m_1','e1_0'],m_data_e1_1=['m_data_e2_0'],e1_2=['e1_0'],m_data_e1_2=['m_data_e1_1'],e1_0=['i_0','m_tda_1','j_0','k_0'],m_data_e1_3=['m_data_e1_1'],m_data_e1_4=['m_data_e1_3'],m_data_e1_5=['m_data_e1_3'],m_data_e1_6=['m_data_e1_5'],m_data_e1_7=['m_data_e1_5'],size2_0=['m_size2_1'],m_size2_0=['n2_0'],m_size2_1=['m_1'],view_0=[],m_block_0=[],e2_6=['e2_4'],e2_4=['e2_2'],e2_5=['e2_4'],e2_2=['e2_0'],e2_3=['e2_2'],view_matrix_tda_0=[],view_matrix_size2_0=[],view_matrix_block_0=[],m_0=[],m_data_0=['array_0'],view_matrix_owner_0=[],e2_0=['j_0','m_tda_2','i_0','k_0'],e2_1=['e2_0'],m_data_e2_0=['m_1','e2_0'],m_data_e2_1=['tmp_0'],m_data_e2_2=['m_data_e2_1'],m_data_e2_3=['m_data_e2_1'],m_data_e2_4=['m_data_e2_3'],m_data_e2_5=['m_data_e2_3'],m_data_e2_6=['m_data_e2_5'],m_data_e2_7=['m_data_e2_5'],size1_0=['m_size1_1'],)

#added phi names
phi_names_set = {'m_tda_7','tmp_5','m_data_e1_6','m_data_e2_6','e1_5','e2_5','m_tda_5','tmp_3','m_data_e1_4','m_data_e2_4','e1_3','e2_3','m_tda_3','tmp_1','m_data_e1_2','m_data_e2_2','e1_1','e2_1','m_tda_4','tmp_2','m_data_e1_3','m_data_e2_3','e1_2','e2_2','m_tda_6','tmp_4','m_data_e1_5','m_data_e2_5','e1_4','e2_4','m_tda_8','tmp_6','m_data_e1_7','m_data_e2_7','e1_6','e2_6',}

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
probability = float(sys.argv[1])/100.0
for arg1 in args1:
    m = gsl_matrix_view_array(arg1.copy(), 8, 8)
    gsl_matrix_transpose(m.matrix)
    bad_dict[test_counter] = (m.matrix.data[0], m.matrix.data[1],m.matrix.data[8], m.matrix.data[63])
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
        importance_p = 2 / ((1 / increase_p) + (1/(math.log(F_p + 0.001) / math.log(total_failed + 0.001))))
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
    bad_runs, good_runs = get_run_ratio(bad_dict, good_dict)
    f.write("Number of Fault Insertions: " + str(insertion_count) + "\n")
    f.write("Number of Faulty Executions: " + str(bad_runs) + "\n")
    f.write(str(suspicious_final_rank.to_csv()))