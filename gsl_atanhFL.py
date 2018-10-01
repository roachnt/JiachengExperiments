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
    bug = fluky(0, 1.38 , 0.99)
    bad_outcome = gsl_atanh(arg1)
    bad_dict[test_counter] = bad_outcome
    test_counter += 1


diff_dict = {index : 0.0 if bad_dict[index] == good_dict[index] else 1.0 for index in bad_dict }


for key in global_value_dict:
    rows = global_value_dict[key].index
    outcome_list = [diff_dict[i] for i in rows]
    global_value_dict[key]['outcome'] = outcome_list


def get_quantiled_tr(W):
    # 10 quantiles from 0.05 to 0.95
    quantile_list = []
    for i in np.arange(0.05, 1.05, 0.1):
        quantile_list.append(W.quantile(i))
    return quantile_list


def predict_causal_risk_list(train_set_X, quantiles, model):

    risk_list = []
    print(train_set_X.columns[0] + " being treatment...")
    X_with_quantile = train_set_X.drop(train_set_X.columns[0], axis=1)

    for quantile in quantiles:
        X_with_quantile.insert(loc=0, column=train_set_X.columns[0],
                               value=np.full((len(X_with_quantile), 1), quantile))
        # X_with_quantile[train_set_X.columns[col_index_todrop]] = np.full((len(X_with_quantile), 1), quantile)
        # print(X_with_quantile.describe())
        risk_list.append(model.predict(X_with_quantile).mean())
        X_with_quantile = X_with_quantile.drop(train_set_X.columns[0], axis=1)
    return risk_list



def suspicious_ranking(global_value_dict, model_to_use):

    suspicious_df = pd.DataFrame(columns=['variable_name', 'max_risk_diff', 'quantile1', 'quantile2'])
    for name in global_value_dict:

        #df cleaning
        #df = global_value_dict[name].select_dtypes(include=[np.number]).dropna(axis=1, how='all')
        df = global_value_dict[name].select_dtypes(include=[np.number]).dropna(axis=1, how='any')
        train_set = df
        #train_set, test_set = train_test_split(df, test_size=0.2, random_state=42)
        train_set_X = train_set.drop(['outcome'], axis=1)
        train_set_Y = train_set['outcome']
        if model_to_use == 0:
            model = RandomForestClassifier(n_estimators=500, max_leaf_nodes=16, n_jobs=-1)
        if model_to_use == 1:
            model = Lasso(alpha=0.1)

        
        model.fit(train_set_X, train_set_Y)

        W = df.iloc[:, 0].to_frame()
        quantiles = get_quantiled_tr(W)
        risk_list = predict_causal_risk_list(train_set_X, quantiles, model)
        max_risk = max(risk_list)
        min_risk = min(risk_list)
        row = [df.columns[0], max_risk - min_risk, risk_list.index(max_risk),
               risk_list.index(min_risk)]
        suspicious_df.loc[len(suspicious_df)] = row
    suspicious_df = suspicious_df.sort_values(by='max_risk_diff', ascending=False)
    return filter_phi_rows(suspicious_df, phi_names_set)

def filter_phi_rows(suspicious_df, phi_names_set):
    return suspicious_df[~suspicious_df['variable_name'].isin(phi_names_set)]


# 0-> random forest  1 -> lasso
result = suspicious_ranking(global_value_dict, 0)
pd.set_option("display.precision", 8)
print('*************Target variables in total: ', len(result),'*************')
print(result)