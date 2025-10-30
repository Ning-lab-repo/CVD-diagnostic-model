import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.utils import resample
from joblib import Parallel, delayed

np.random.seed(42)  

cvds = ['I10','I25','I20','I48','I21','I50','I47','I73','I26','G45','I51','I63','I49','I44','I34','I35','I45','I67','I77','I12','I64','I42','I69','I60','I70','I22',
    'I65','I74','I24','I31','I08','I71','I61','I78','I46','I27','I05','I25.1','I20.9','I25.9','I25.8','I20.0','I25.2','I21.9','I21.1','I50.1','I26.9','I47.1',
    'I73.9','I21.0','G45.9','I63.9','I34.0','I51.7','I20.8','I21.4','I12.0','I50.0','I35.0','I44.7','I45.1','I47.2','I49.9','I65.2','I73.0','I77.1','I69.4',
    'I60.9','I70.2','I42.0','I35.1','I50.9','I67.9','I34.1','I22.9','I31.9','I74.3','I24.8','I49.8','I78.1','I44.2','I21.2','I44.0','I51.8','I42.9','I71.4','I61.9']

n_bootstraps = 1000
confidence_level = 0.95
n_jobs = 30

def bootstrap_auc(y_true, y_pred_proba, random_state=None):
  
    indices = resample(np.arange(len(y_true)), replace=True, random_state=random_state)
    if len(np.unique(y_true[indices])) < 2:
        return None
    fpr_b, tpr_b, _ = roc_curve(y_true[indices], y_pred_proba[indices])
    return auc(fpr_b, tpr_b)


for cvd in cvds:
    
    cvd_roc_curves = {}
    cvd_conf_intervals = {}
    
    if '.' in cvd:
        
        model_files_cvd = {
            'Logistic Regression': f'./LR_Z_DJ_patient_predictions_{cvd}.csv',
            'Random Forest': f'./RF_Z_DJ_patient_predictions_{cvd}.csv',
            'XGBoost': f'./X_Z_DJ_patient_predictions_{cvd}.csv'
        }
    else:
        
        model_files_cvd = {
            'Logistic Regression': f'./LR_Z_DJ_patient_predictions_{cvd}.csv',
            'Random Forest': f'./RF_Z_DJ_patient_predictions_{cvd}.csv',
            'XGBoost': f'./X_Z_DJ_patient_predictions_{cvd}.csv'
        }
    
    

    for model_name, file_path in model_files_cvd.items():
        try:
            df = pd.read_csv(file_path)
            y_true = df['Actual_Label']
            y_pred_proba = df['Predicted_Score']
            
   
            fpr, tpr, thresholds = roc_curve(y_true, y_pred_proba)
            roc_auc = auc(fpr, tpr)
            
       
            bootstrapped_scores = Parallel(n_jobs=n_jobs)(
                delayed(bootstrap_auc)(y_true, y_pred_proba, i) for i in range(n_bootstraps))
            bootstrapped_scores = [score for score in bootstrapped_scores if score is not None]
            
 
            sorted_scores = np.array(bootstrapped_scores)
            sorted_scores.sort()
            lower_bound = np.percentile(sorted_scores, (1-confidence_level)/2 * 100)
            upper_bound = np.percentile(sorted_scores, (1+confidence_level)/2 * 100)
            
       
            cvd_roc_curves[model_name] = (fpr, tpr, roc_auc)
            cvd_conf_intervals[model_name] = (lower_bound, upper_bound)
            
        except FileNotFoundError:
            print(f"文件未找到: {file_path}")
            continue
    

    plt.rcParams['font.family'] = 'Arial'
    plt.figure(figsize=(5, 5))
    
    for model_name, (fpr, tpr, auc_value) in cvd_roc_curves.items():
        lower_bound, upper_bound = cvd_conf_intervals[model_name]
        plt.plot(fpr, tpr, label=f'{model_name} (AUC = {auc_value:.3f} [{lower_bound:.3f}-{upper_bound:.3f}])')
    
    plt.plot([0, 1], [0, 1], linestyle='--', color='lightgray', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.01])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.legend(loc='lower right')
    
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.grid(False)
    
    plt.savefig(f'./z_dj_roc_325/roc_dj_{cvd}.pdf', format='pdf')
    plt.close()  
