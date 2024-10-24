###########################################################################################
# To predict T2DM-PaC comorbidity
# Input: pre-processed and normalized gene expression data
# OUTPUT: Score file with 3 columns, for T2DM prediction, PaC prediction, and comorbidity prediction.  
# OUTPT: '1' indicates the presence and '0' indicates the absence of respective disease
###########################################################################################
# load the expression values of 67 genes features (upload in google colab drive) 
# upload the model files gnb_pac.sav, xgb_pac.sav, svm_t2d.sav, lr_t2d.sav, (upload in google colab drive)


#install required packages 
pip install pickle # usually not required in google colab
pip install pandas # usually not required in google colab
pip install numpy # usually not required in google colab
pip install repeat # usually not required in google colab
pip install os # usually not required in google colab
pip install dill
pip install scikit-optimize

# Data processing
import pickle
import pandas as pd
import numpy as np
from itertools import repeat
import os
import dill
import skopt

#data loading
os.chdir("./T2DM-PaC_comorbidity_predictor-main")#change working directory
input_file_name = "./pre_processing_example/input_data.csv" #change file name with your input file
input_data = pd.read_csv(input_file_name)
input_data = input_data.transpose()
input_data.columns = input_data.iloc[0,]
input_data.drop(["genes"], inplace = True)
sample_index = input_data.index
input_data_2 =np.array(input_data)

####MODEL LOADING####
#load PaC model
pac_gnb_model = pickle.load(open('./Models/gnb_pac.sav', 'rb'))
pac_xgb_model = pickle.load(open('./Models/xgb_pac.sav', 'rb'))
cutoff_threshold_pac = 0.48
#load t2D model
t2d_svm_model = pickle.load(open('./Models/svm_t2d.sav', 'rb'))
t2d_lr_model = pickle.load(open('./Models/lr_t2d.sav', 'rb'))
cutoff_threshold_t2d = 0.52

####MAKING PREDICTION####
#t2d prediction
t2d_lr_como = pd.DataFrame(t2d_lr_model.predict_proba(input_data),columns=["0","1"])
t2d_svm_como = pd.DataFrame(t2d_svm_model.predict_proba(input_data),columns=["0","1"])
#pac prediction
pac_xgb_como = pd.DataFrame(pac_xgb_model.predict_proba(input_data_2),columns=["0","1"])
pac_gnb_como = pd.DataFrame(pac_gnb_model.predict_proba(input_data),columns=["0","1"])

####WEIGHTED MODEL COMBINATION####
#pac combination
w_pac_gnb_co = 0.5161821449336343*pac_gnb_como.iloc[:,1]
w_pac_xgb_co = 0.9708451791592324*pac_xgb_como.iloc[:,1]
average_predictions_pac = (w_pac_gnb_co + w_pac_xgb_co)/(0.5161821449336343+0.9708451791592324)
results_pac_test = (average_predictions_pac >= cutoff_threshold_pac).astype(int)
#t2d combination
w_t2d_svm_co = 0.8262620868032353*t2d_svm_como.iloc[:,1]
w_t2d_lr_co = 0.5484805963516988*t2d_lr_como.iloc[:,1]
average_predictions_t2d = (w_t2d_svm_co + w_t2d_lr_co)/(0.8262620868032353+0.5484805963516988)
results_t2d_test = (average_predictions_t2d >= cutoff_threshold_t2d).astype(int)

#comorbidity voting
average_prediction = (results_t2d_test + results_pac_test)/2
results_como = (average_prediction == 1).astype(int)

df_c = pd.concat([results_pac_test,results_t2d_test,results_como], axis=1)
df_c.columns = ["PaC_prediction","T2D_prediction","comorbidity_prediction"]
df_c.index = sample_index
df_c.to_csv("output_data.csv")#,index = False)
