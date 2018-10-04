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
from math import *
from collections import namedtuple
import sys

from J0Long import good_dict
os.system('python J0Long.py') 

from phi import *
dd_0=None;dd_2=None;dd_1=None;dd_3=None;reason_0=None;ampl_0=None;ampl_1=None;abs_sum_0=None;cs_c_cs_order_IV_0=None;result_val_0=None;result_val_1=None;result_val_2=None;result_val_3=None;ca_val_IV_0=None;ca_val_IV_1=None;cs_a_IV_0=None;result_val_IV_0=None;result_val_IV_1=None;result_val_IV_2=None;order_0=None;val_0=None;eps_0=None;result_err_0=None;result_err_1=None;result_err_2=None;result_err_3=None;result_err_4=None;result_err_5=None;result_err_6=None;result_err_7=None;result_err_8=None;result_err_9=None;result_err_10=None;result_err_11=None;cs_b_IV_0=None;status_0=None;stat_ca_0=None;stat_ca_1=None;GSL_NAN_0=None;ca_err_IV_0=None;ca_err_IV_1=None;cs_c_j_IV_1=None;cs_c_j_IV_0=None;cs_c_j_IV_2=None;fn_0=None;fn_1=None;return_val_IV_0=None;ceps_0=None;ceps_1=None;ceps_2=None;ct_val_IV_0=None;ct_val_IV_1=None;result_0=None;result_1=None;result_2=None;result_3=None;result_4=None;sqrty_0=None;sqrty_1=None;cp_val_IV_0=None;cp_val_IV_1=None;cp_val_IV_2=None;stat_cp_0=None;stat_cp_1=None;gsl_errno_0=None;stat_ct_0=None;stat_ct_1=None;y2_0=None;ca_0=None;ca_1=None;seps_0=None;seps_1=None;seps_2=None;a_0=None;a_1=None;a_2=None;a_3=None;a_4=None;b_0=None;b_1=None;b_2=None;b_3=None;b_4=None;temp_1=None;temp_0=None;temp_2=None;temp_3=None;c_0=None;c_1=None;c_2=None;c_3=None;err_0=None;d_0=None;d_2=None;d_1=None;d_3=None;d_4=None;d_5=None;d_6=None;d_7=None;e_0=None;e_2=None;e_1=None;e_3=None;e_4=None;e_5=None;sy_0=None;cp_err_IV_0=None;cp_err_IV_1=None;cs_c_0_IV_0=None;j_0=None;e2_0=None;e2_1=None;cp_0=None;cp_1=None;cs_0=None;ct_0=None;ct_1=None;s_0=None;cy_0=None;x_0=None;x_1=None;x_2=None;y_0=None;y_1=None;y_2=None;z_0=None;z_1=None;order_sp_0=None


GSL_SUCCESS=0 
GSL_FAILURE=-1 
GSL_CONTINUE=-2 
GSL_EDOM=1 
GSL_ERANGE=2 
GSL_EFAULT=3 
GSL_EINVAL=4 
GSL_EFAILED=5 
GSL_EFACTOR=6 
GSL_ESANITY=7 
GSL_ENOMEM=8 
GSL_EBADFUNC=9 
GSL_ERUNAWAY=1 
GSL_EMAXITER=11 
GSL_EZERODIV=12 
GSL_EBADTOL=13 
GSL_ETOL=14 
GSL_EUNDRFLW=15 
GSL_EOVRFLW=16 
GSL_ELOSS=17 
GSL_EROUND=18 
GSL_EBADLEN=19 
GSL_ENOTSQR=20 
GSL_ESING=21 
GSL_EDIVERGE=22 
GSL_EUNSUP=23 
GSL_EUNIMPL=24 
GSL_ECACHE=25 
GSL_ETABLE=26 
GSL_ENOPROG=27 
GSL_ENOPROGJ=28 
GSL_ETOLF=29 
GSL_ETOLX=30 
GSL_ETOLG=31 
GSL_EOF=32 
GSL_DBL_EPSILON=2.2204460492503131e-16 
class gsl_sf_result:
    def __init__(self,val,err):
        val_0 = val;err_0 = err;
        

        self.val=val_0 
        self.err=err_0 

bj0_data=[0.100254161968939137,-0.665223007764405132,0.248983703498281314,-0.0332527231700357697,0.0023114179304694015,-0.0000991127741995080,0.0000028916708643998,-0.0000000612108586630,0.0000000009838650793,-0.0000000000124235515,0.0000000000001265433,-0.0000000000000010619,0.0000000000000000074] 
class cheb_series:
    def __init__(self,c,order,a,b,order_sp):
        c_0 = c;order_0 = order;a_0 = a;b_0 = b;order_sp_0 = order_sp;
        

        self.c=c_0 
        self.order=order_0 
        self.a=a_0 
        self.b=b_0 
        self.order_sp=order_sp_0 

bj0_cs=cheb_series(bj0_data,12,-1,1,9) 
def cheb_eval_e(cs,x,result):
    cs_0 = cs;x_0 = x;result_0 = result;
    dd_0=None;dd_2=None;dd_1=None;dd_3=None;temp_1=None;temp_0=None;temp_2=None;temp_3=None;d_0=None;d_2=None;d_1=None;d_3=None;d_4=None;e_0=None;e_2=None;e_1=None;e_3=None;e_4=None;cs_c_j_IV_1=None;cs_c_j_IV_0=None;cs_c_j_IV_2=None;cs_c_cs_order_IV_0=None;cs_c_0_IV_0=None;result_err_0=None;cs_a_IV_0=None;result_val_IV_0=None;cs_b_IV_0=None;y_0=None;y2_0=None;

    d_0=0.0 
    dd_0=0.0 
    cs_a_IV_0=cs_0.a 
    cs_b_IV_0=cs_0.b 
    y_0=(2.0*x_0-cs_a_IV_0-cs_b_IV_0)/(cs_b_IV_0-cs_a_IV_0) 
    y2_0=2.0*y_0 
    e_0=0.0 
    phi0 = Phi()
    for j_0 in range(cs_0.order,0,-1):
        phi0.set()
        dd_2 = phi0.phiEntry(dd_0,dd_1)
        temp_1 = phi0.phiEntry(None,temp_0)
        d_2 = phi0.phiEntry(d_0,d_1)
        e_2 = phi0.phiEntry(e_0,e_1)
        cs_c_j_IV_1 = phi0.phiEntry(None,cs_c_j_IV_0)

        temp_0=d_2 
        cs_c_j_IV_0=cs_0.c[j_0] 
        d_1=y2_0*d_2-dd_2+cs_c_j_IV_0 
        e_1 = e_2+fabs(y2_0*temp_0)+fabs(dd_2)+fabs(cs_c_j_IV_0)
        dd_1=temp_0 
    dd_3 = phi0.phiExit(dd_0,dd_1)
    temp_2 = phi0.phiExit(None,temp_0)
    d_3 = phi0.phiExit(d_0,d_1)
    e_3 = phi0.phiExit(e_0,e_1)
    cs_c_j_IV_2 = phi0.phiExit(None,cs_c_j_IV_0)
    temp_3=d_3 
    cs_c_0_IV_0=cs_0.c[0] 
    d_4=y_0*d_3-dd_3+0.5*cs_c_0_IV_0 
    e_4 = e_3+fabs(y_0*temp_3)+fabs(dd_3)+0.5*fabs(cs_c_0_IV_0)
    result_val_IV_0=d_4 
    result_0.val=result_val_IV_0 
    cs_c_cs_order_IV_0=cs_0.c[cs_0.order] 
    result_err_0=GSL_DBL_EPSILON*e_4+fabs(cs_c_cs_order_IV_0) 
    result_0.err=result_err_0 
    lo = locals()
    record_locals(lo, test_counter)
    return GSL_SUCCESS

GSL_SQRT_DBL_EPSILON=1.4901161193847656e-08 
GSL_ROOT5_DBL_EPSILON=7.4009597974140505e-04 
bm0_data=[0.09284961637381644,-0.00142987707403484,0.00002830579271257,-0.00000143300611424,0.00000012028628046,-0.00000001397113013,0.00000000204076188,-0.00000000035399669,0.00000000007024759,-0.00000000001554107,0.00000000000376226,-0.00000000000098282,0.00000000000027408,-0.00000000000008091,0.00000000000002511,-0.00000000000000814,0.00000000000000275,-0.00000000000000096,0.00000000000000034,-0.00000000000000012,0.00000000000000004] 
_gsl_sf_bessel_amp_phase_bm0_cs=cheb_series(bm0_data,20,-1,1,10) 
bth0_data=[-0.24639163774300119,0.001737098307508963,-0.000062183633402968,0.000004368050165742,-0.000000456093019869,0.000000062197400101,-0.000000010300442889,0.000000001979526776,-0.000000000428198396,0.000000000102035840,-0.000000000026363898,0.000000000007297935,-0.000000000002144188,0.000000000000663693,-0.000000000000215126,0.000000000000072659,-0.000000000000025465,0.000000000000009229,-0.000000000000003448,0.000000000000001325,-0.000000000000000522,0.000000000000000210,-0.000000000000000087,0.000000000000000036] 
_gsl_sf_bessel_amp_phase_bth0_cs=cheb_series(bth0_data,23,-1,1,12) 
def gsl_sf_bessel_cos_pi4_e(y,eps,result):
    y_1 = y;eps_0 = eps;result_1 = result;
    seps_0=None;seps_1=None;seps_2=None;d_5=None;sy_0=None;abs_sum_0=None;ceps_0=None;ceps_1=None;ceps_2=None;result_err_1=None;result_err_2=None;result_err_3=None;result_err_4=None;result_err_5=None;result_err_6=None;e2_0=None;e2_1=None;result_val_0=None;s_0=None;cy_0=None;

    sy_0=sin(y_1) 
    cy_0=cos(y_1) 
    s_0=sy_0+cy_0 
    d_5=sy_0-cy_0 
    abs_sum_0=fabs(cy_0)+fabs(sy_0) 
    if fabs(eps_0)<GSL_ROOT5_DBL_EPSILON:
        e2_0=eps_0*eps_0 
        seps_0=eps_0*(1.0-e2_0/6.0*(1.0-e2_0/20.0)) 
        ceps_0=1.0-e2_0/2.0*(1.0-e2_0/12.0) + bug
    else:
        seps_1=sin(eps_0) 
        ceps_1=cos(eps_0) 
    phiPreds = [fabs(eps_0)<GSL_ROOT5_DBL_EPSILON]
    phiNames = [seps_0,seps_1]
    seps_2= phiIf(phiPreds, phiNames)
    phiPreds = [fabs(eps_0)<GSL_ROOT5_DBL_EPSILON]
    phiNames = [ceps_0,ceps_1]
    ceps_2= phiIf(phiPreds, phiNames)
    phiPreds = [fabs(eps_0)<GSL_ROOT5_DBL_EPSILON]
    phiNames = [e2_0,None]
    e2_1= phiIf(phiPreds, phiNames)
    result_val_0=(ceps_2*s_0-seps_2*d_5)/sqrt(2) 
    result_1.val=result_val_0 
    result_err_1=2.0*GSL_DBL_EPSILON*(fabs(ceps_2)+fabs(seps_2))*abs_sum_0/sqrt(2) 
    result_1.err=result_err_1 
    if y_1>1.0/GSL_DBL_EPSILON:
        result_err_2=result_1.err 
        result_err_3 = result_err_2*0.5*y_1
        result_1.err=result_err_3 
    elif y_1>1.0/GSL_SQRT_DBL_EPSILON:
        result_err_4=result_1.err 
        result_err_5 = result_err_4*256.0*y_1*GSL_SQRT_DBL_EPSILON
        result_1.err=result_err_5 
    phiPreds = [y_1>1.0/GSL_DBL_EPSILON,y_1>1.0/GSL_SQRT_DBL_EPSILON]
    phiNames = [result_err_3,result_err_5,result_err_1]
    result_err_6= phiIf(phiPreds, phiNames)
    lo = locals()
    record_locals(lo, test_counter)
    return GSL_SUCCESS

def GSL_ERROR_SELECT_2(a,b):
    a_1 = a;b_1 = b;
    

    return a_1 if a_1!=GSL_SUCCESS else (b_1 if b_1!=GSL_SUCCESS else GSL_SUCCESS)

def GSL_ERROR_SELECT_3(a,b,c):
    a_2 = a;b_2 = b;c_1 = c;
    

    return a_2 if a_2!=GSL_SUCCESS else GSL_ERROR_SELECT_2(b_2,c_1)

def GSL_ERROR_SELECT_4(a,b,c,d):
    a_3 = a;b_3 = b;c_2 = c;d_6 = d;
    

    return a_3 if a_3!=GSL_SUCCESS else GSL_ERROR_SELECT_3(b_3,c_2,d_6)

def GSL_ERROR_SELECT_5(a,b,c,d,e):
    a_4 = a;b_4 = b;c_3 = c;d_7 = d;e_5 = e;
    

    return a_4 if a_4!=GSL_SUCCESS else GSL_ERROR_SELECT_4(b_4,c_3,d_7,e_5)

def gsl_sf_bessel_J0_e(x,result):
    x_1 = x;result_2 = result;
    stat_ca_0=None;stat_ca_1=None;ampl_0=None;ampl_1=None;ca_err_IV_0=None;ca_err_IV_1=None;cp_err_IV_0=None;cp_err_IV_1=None;result_err_7=None;result_err_8=None;result_err_9=None;result_err_10=None;result_err_11=None;cp_0=None;cp_1=None;ct_val_IV_0=None;ct_val_IV_1=None;result_val_1=None;result_val_2=None;result_val_3=None;sqrty_0=None;sqrty_1=None;ct_0=None;ct_1=None;ca_val_IV_0=None;ca_val_IV_1=None;cp_val_IV_0=None;cp_val_IV_1=None;cp_val_IV_2=None;stat_cp_0=None;stat_cp_1=None;result_val_IV_1=None;result_val_IV_2=None;stat_ct_0=None;stat_ct_1=None;y_2=None;z_0=None;z_1=None;ca_0=None;ca_1=None;

    y_2=fabs(x_1) 
    if y_2<2.0*GSL_SQRT_DBL_EPSILON:
        result_val_1=1.0 
        result_2.val=result_val_1 
        result_err_7=y_2*y_2 
        result_2.err=result_err_7 
        lo = locals()
        record_locals(lo, test_counter)
        return GSL_SUCCESS
    elif y_2<=4.0:
        lo = locals()
        record_locals(lo, test_counter)
        return cheb_eval_e(bj0_cs,0.125*y_2*y_2-1.0,result_2)
    else:
        ca_0=gsl_sf_result(0.0,0.0) 
        ct_0=gsl_sf_result(0.0,0.0) 
        cp_0=gsl_sf_result(0.0,0.0) 
        z_0=32.0/(y_2*y_2)-1.0 
        stat_ca_0=cheb_eval_e(_gsl_sf_bessel_amp_phase_bm0_cs,z_0,ca_0) 
        stat_ct_0=cheb_eval_e(_gsl_sf_bessel_amp_phase_bth0_cs,z_0,ct_0) 
        ct_val_IV_0=ct_0.val 
        stat_cp_0=gsl_sf_bessel_cos_pi4_e(y_2,ct_val_IV_0/y_2,cp_0) 
        sqrty_0=sqrt(y_2) 
        ca_val_IV_0=ca_0.val 
        ampl_0=(0.75+ca_val_IV_0)/sqrty_0
        cp_val_IV_0=cp_0.val 
        result_val_2=ampl_0*cp_val_IV_0 
        result_2.val=result_val_2 
        cp_val_IV_1=cp_0.val 
        ca_err_IV_0=ca_0.err 
        cp_err_IV_0=cp_0.err 
        result_err_8=fabs(cp_val_IV_1)*ca_err_IV_0/sqrty_0+fabs(ampl_0)*cp_err_IV_0 
        result_2.err=result_err_8 
        result_val_IV_1=result_2.val 
        result_err_9=result_2.err 
        result_err_10 = result_err_9+GSL_DBL_EPSILON*fabs(result_val_2)
        result_2.err=result_err_10 
        lo = locals()
        record_locals(lo, test_counter)
        return GSL_ERROR_SELECT_3(stat_ca_0,stat_ct_0,stat_cp_0)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,stat_ca_0]
    stat_ca_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,ampl_0]
    ampl_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,ca_err_IV_0]
    ca_err_IV_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,cp_err_IV_0]
    cp_err_IV_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [result_err_7,None,result_err_10]
    result_err_11= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,cp_0]
    cp_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,ct_val_IV_0]
    ct_val_IV_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [result_val_1,None,result_val_2]
    result_val_3= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,sqrty_0]
    sqrty_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,ct_0]
    ct_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,ca_val_IV_0]
    ca_val_IV_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,cp_val_IV_1]
    cp_val_IV_2= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,stat_cp_0]
    stat_cp_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,result_val_IV_1]
    result_val_IV_2= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,stat_ct_0]
    stat_ct_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,z_0]
    z_1= phiIf(phiPreds, phiNames)
    phiPreds = [y_2<2.0*GSL_SQRT_DBL_EPSILON,y_2<=4.0]
    phiNames = [None,None,ca_0]
    ca_1= phiIf(phiPreds, phiNames)

def GSL_ERROR_VAL(reason,gsl_errno,GSL_NAN):
    reason_0 = reason;gsl_errno_0 = gsl_errno;GSL_NAN_0 = GSL_NAN;
    

    return 

def EVAL_RESULT(fn,result):
    fn_0 = fn;result_3 = result;
    return_val_IV_0=None;status_0=None;

    status_0=fn_0 
    if status_0!=GSL_SUCCESS:
        GSL_ERROR_VAL(fn_0,status_0,result_3.val) 
    return_val_IV_0=result_3.val 
    lo = locals()
    record_locals(lo, test_counter)
    return return_val_IV_0

def gsl_sf_bessel_J0(x):
    x_2 = x;
    result_4=None;fn_1=None;

    result_4=gsl_sf_result(0.0,0.0) 
    fn_1=gsl_sf_bessel_J0_e(x_2,result_4) 
    lo = locals()
    record_locals(lo, test_counter)
    return EVAL_RESULT(fn_1,result_4)



#generate python causal map
causal_map = {'cs_c_0_IV_0':['cs_0'],'cs_c_j_IV_1':['cs_c_j_IV_0'],'cs_c_j_IV_0':['cs_0','j_0'],'cy_0':['y_1'],'d_0':[],'ca_err_IV_0':['ca_0'],'d_2':['d_0','d_1'],'ca_err_IV_1':['ca_err_IV_0'],'d_1':['y2_0','d_2','dd_2','cs_c_j_IV_0'],'d_4':['y_0','d_3','dd_3','cs_c_0_IV_0'],'d_3':['d_0','d_1'],'cs_a_IV_0':['cs_0'],'dd_3':['dd_0','dd_1'],'abs_sum_0':['cy_0','sy_0'],'d_5':['sy_0','cy_0'],'dd_1':['temp_0'],'dd_2':['dd_0','dd_1'],'dd_0':[],'sy_0':['y_1'],'ca_0':[],'cs_c_cs_order_IV_0':['cs_0','cs_0'],'cp_val_IV_1':['cp_0'],'cp_val_IV_0':['cp_0'],'z_0':['y_2','y_2'],'cs_b_IV_0':['cs_0'],'z_1':['z_0'],'stat_ct_1':['stat_ct_0'],'temp_0':['d_2'],'temp_1':['temp_0'],'temp_2':['temp_0'],'stat_ct_0':['z_0','ct_0'],'stat_cp_1':['stat_cp_0'],'cs_c_j_IV_2':['cs_c_j_IV_0'],'stat_cp_0':['y_2','ct_val_IV_0','y_2','cp_0'],'cp_val_IV_2':['cp_val_IV_1'],'status_0':['fn_0'],'temp_3':['d_3'],'ca_1':['ca_0'],'result_val_1':[],'result_val_0':['ceps_2','s_0','seps_2','d_5'],'seps_1':['eps_0'],'result_val_IV_1':['result_2'],'stat_ca_1':['stat_ca_0'],'ct_val_IV_1':['ct_val_IV_0'],'result_val_3':['result_val_1','result_val_2'],'result_val_IV_0':['d_4'],'seps_0':['eps_0','e2_0','e2_0'],'stat_ca_0':['z_0','ca_0'],'ct_val_IV_0':['ct_0'],'result_val_2':['ampl_0','cp_val_IV_0'],'seps_2':['seps_0','seps_1'],'fn_1':['x_2','result_4'],'ct_0':[],'ct_1':['ct_0'],'e_1':['e_2','y2_0','temp_0','dd_2','cs_c_j_IV_0'],'e_0':[],'e_3':['e_0','e_1'],'e_2':['e_0','e_1'],'e_4':['e_3','y_0','temp_3','dd_3','cs_c_0_IV_0'],'ampl_1':['ampl_0'],'ampl_0':['ca_val_IV_0','sqrty_0'],'ceps_1':['eps_0'],'ca_val_IV_0':['ca_0'],'s_0':['sy_0','cy_0'],'ceps_0':['e2_0','e2_0'],'ca_val_IV_1':['ca_val_IV_0'],'y_0':['x_0','cs_a_IV_0','cs_b_IV_0','cs_b_IV_0','cs_a_IV_0'],'ceps_2':['ceps_0','ceps_1'],'return_val_IV_0':['result_3'],'y_2':['x_1'],'result_4':[],'e2_0':['eps_0','eps_0'],'cp_0':[],'cp_err_IV_0':['cp_0'],'result_err_9':['result_2'],'e2_1':['e2_0'],'result_err_8':['cp_val_IV_1','ca_err_IV_0','sqrty_0','ampl_0','cp_err_IV_0'],'cp_err_IV_1':['cp_err_IV_0'],'cp_1':['cp_0'],'result_err_10':['result_err_9','result_val_2'],'result_err_5':['result_err_4','y_1'],'result_err_11':['result_err_7','result_err_10'],'result_err_4':['result_1'],'result_err_7':['y_2','y_2'],'sqrty_1':['sqrty_0'],'result_err_6':['result_err_3','result_err_5','result_err_1'],'sqrty_0':['y_2'],'result_err_1':['ceps_2','seps_2','abs_sum_0'],'result_err_0':['e_4','cs_c_cs_order_IV_0'],'result_val_IV_2':['result_val_IV_1'],'result_err_3':['result_err_2','y_1'],'result_err_2':['result_1'],'y2_0':['y_0'],}

#added phi names
phi_names_set = {'dd_2','temp_1','d_2','e_2','cs_c_j_IV_1','dd_3','temp_2','d_3','e_3','cs_c_j_IV_2','seps_2','ceps_2','e2_1','result_err_6','stat_ca_1','ampl_1','ca_err_IV_1','cp_err_IV_1','result_err_11','cp_1','ct_val_IV_1','result_val_3','sqrty_1','ct_1','ca_val_IV_1','cp_val_IV_2','stat_cp_1','result_val_IV_2','stat_ct_1','z_1','ca_1',}


#--------------------end of progarm-----------

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
arg1s = np.arange(0, 1000)
test_counter = 0
bug = 0
probability = float(sys.argv[1])/100.0
for arg1 in arg1s:
    bug = fluky(0, -0.00027, probability)
    bad_outcome = gsl_sf_bessel_J0(arg1)
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
        importance_p = 2 / ((1 / increase_p) + (1/(log(F_p + 0.00001) / log(total_failed + 0.00001))))
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
    

with open(os.path.basename(__file__)[:-3] + str(probability) + ".txt", "w") as f:
    f.write('*************Target variables in total: ' + str(len(suspicious_final_rank)) + '*************\n')
    f.write(str(suspicious_final_rank))