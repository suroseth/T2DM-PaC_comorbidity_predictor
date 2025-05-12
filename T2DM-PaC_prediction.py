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
from scipy import stats

#function to fit into distribution of training set
def match_quantiles(counts_sub, old_mu, old_phi, new_mu, new_phi):
 
    new_counts_sub = np.zeros_like(counts_sub, dtype=float)
    for a in range(counts_sub.shape[0]):
        for b in range(counts_sub.shape[1]):
            if counts_sub[a, b] <= 1:
                new_counts_sub[a, b] = counts_sub[a, b]
            else:
                r_old = 1/old_phi[a]  # size parameter (r in negative binomial)
                p_old = r_old/(r_old + old_mu[a, b])  # convert mu to probability parameter
                
                # Calculate CDF (equivalent to pnbinom in R)
                tmp_p = stats.nbinom.cdf(counts_sub[a, b]-1, r_old, p_old)
                
                if abs(tmp_p - 1) < 1e-4:
                    new_counts_sub[a, b] = counts_sub[a, b]  # for outlier count
                else:
                    r_new = 1/new_phi[a]  # new size parameter
                    p_new = r_new/(r_new + new_mu[a, b])  # new probability parameter
                    
                    # Calculate quantile (equivalent to qnbinom in R)
                    new_counts_sub[a, b] = 1 + stats.nbinom.ppf(tmp_p, r_new, p_new)
    
    return new_counts_sub


def adding_parameters(new_counts, parameters, batch_number, gene_list):

    if not isinstance(parameters, pd.DataFrame):
        parameters = pd.DataFrame(parameters)
    if not isinstance(new_counts, pd.DataFrame):
        new_counts = pd.DataFrame(new_counts)
    
    # Sort parameters by row names
    parameters = parameters.sort_index()
    
    # Filter new_counts to include only genes in gene_list
    new_counts = new_counts.loc[new_counts.index.isin(gene_list)]
    
    # Sort new_counts by row names
    new_counts = new_counts.sort_index()
    
    # Extract parameters for the specific batch
    col1 = f"old_mu_avg{batch_number}"
    old_mu_avg = parameters.loc[parameters.index.isin(gene_list), col1].values.reshape(-1, 1)
    
    col2 = f"old_phi_avg{batch_number}"
    old_phi_avg = parameters.loc[parameters.index.isin(gene_list), col2].values
    
    col3 = f"new_mu_avg{batch_number}"
    new_mu_avg = parameters.loc[parameters.index.isin(gene_list), col3].values.reshape(-1, 1)
    
    col4 = f"new_phi_avg{batch_number}"
    new_phi_avg = parameters.loc[parameters.index.isin(gene_list), col4].values
    
    # Initialize result matrix
    result_x = np.zeros((new_counts.shape[0], 0))
    
    # Process each column
    for i in range(new_counts.shape[1]):
        x = match_quantiles(
            counts_sub=new_counts.iloc[:, i].values.reshape(-1, 1),
            old_mu=old_mu_avg,
            old_phi=old_phi_avg,
            new_mu=new_mu_avg,
            new_phi=new_phi_avg
        )
        result_x = np.column_stack((result_x, x))
    
    # Create DataFrame with proper row and column names
    result_df = pd.DataFrame(result_x, index=new_counts.index, columns=new_counts.columns)
    
    return result_df



#data loading
#os.chdir("./T2DM-PaC_comorbidity_predictor-main")#change working directory
input_file_name = "./pre_processing_example/input_data.csv" #change file name with your input file
input_data = pd.read_csv(input_file_name)
#input pac model
input_data_pac = adding_parameters(input_data, pac_combat_param, batch_number=8, gene_list=genes_67_features['Genes'])
input_data_pac = input_data_pac.transpose()
input_data_pac.columns = input_data_pac.iloc[0,]
input_data_pac.drop(["genes"], inplace = True)
sample_index = input_data_pac.index
input_data_pac_2 =np.array(input_data_pac)

#input t2d model
input_data_t2d = adding_parameters(input_data, t2d_combat_param, batch_number=7, gene_list=genes_67_features['Genes'])
input_data_t2d = input_data_t2d.transpose()
input_data_t2d.columns = input_data_t2d.iloc[0,]
input_data_t2d.drop(["genes"], inplace = True)
sample_index = input_data_t2d.index
input_data_t2d_2 =np.array(input_data_t2d)

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
t2d_lr_como = pd.DataFrame(t2d_lr_model.predict_proba(input_data_t2d),columns=["0","1"])
t2d_svm_como = pd.DataFrame(t2d_svm_model.predict_proba(input_data_t2d),columns=["0","1"])
#pac prediction
pac_xgb_como = pd.DataFrame(pac_xgb_model.predict_proba(input_data_pac_2),columns=["0","1"])
pac_gnb_como = pd.DataFrame(pac_gnb_model.predict_proba(input_data_pac),columns=["0","1"])

####WEIGHTED MODEL COMBINATION####
#pac combination
w_pac_gnb_co = 0.5161821449336343*pac_gnb_como.iloc[:,1]
w_pac_xgb_co = 0.9856469185969381*pac_xgb_como.iloc[:,1]
average_predictions_pac = (w_pac_gnb_co + w_pac_xgb_co)/(0.5161821449336343+0.9856469185969381)
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
