import joblib
import shap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed

if 'DISPLAY' in os.environ:
    del os.environ['DISPLAY']

def process_cvd(cvd):
    print(f"Processing cvd: {cvd}")
    

    scaler_path = os.path.join(base_path, f'X_Z_D_final_scaler_{cvd}.joblib')
    imputer_path = os.path.join(base_path, f'X_Z_D_final_simple_imputer_{cvd}.joblib')
    model_path = os.path.join(base_path, f'X_Z_D_final_xgb_model_{cvd}.joblib')
    
 
    if not (os.path.exists(scaler_path) and os.path.exists(imputer_path) and os.path.exists(model_path)):
        print(f"One or more files not found for cvd {cvd}. Skipping this cvd.")
        return
    

    final_scaler = joblib.load(scaler_path)
    final_imputer = joblib.load(imputer_path)
    model = joblib.load(model_path)
    

    X = final_imputer.transform(features)
    X = final_scaler.transform(X)
    X_sample = pd.DataFrame(X, columns=features.columns)
    

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_sample)
    if isinstance(shap_values, list):
        shap_values = shap_values[1] 
    print(shap_values)
    

    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_values, X_sample, max_display=10, show=False)
    output_pdf = os.path.join(out_path, f'xgb_{cvd}_Dshap.pdf')
    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()
    

    importance = pd.DataFrame({
        'feature': features.columns,
        'mean_abs_shap': np.abs(shap_values).mean(axis=0)
    }).nlargest(50, 'mean_abs_shap')
    

    importance_csv = os.path.join(out_path, f'xgb_Dshap_importance_{cvd}.csv')
    importance.to_csv(importance_csv, index=False)
    

    top_features = importance['feature'].tolist()
    top_feature_indices = [features.columns.get_loc(feature) for feature in top_features]
    top_shap_values = shap_values[:, top_feature_indices]
    top_shap_df = pd.DataFrame(top_shap_values, columns=top_features)
    

    top_shap_csv = os.path.join(out_path, f'xgb_{cvd}_Dtop50_shap_values.csv')
    top_shap_df.to_csv(top_shap_csv, index=False)
    
    print(f"Finished processing cvd: {cvd}\n")

if __name__ == '__main__':

    base_path = './xgb_zukbnmr_d'
    out_path = './xgbSHAP/xgb_d_re1'
 
    plt.rcParams['font.family'] = 'Arial'
    
   
    cvd_ids = [
        'I10','I25','I20','I48','I21','I50','I47','I73','I26','G45','I51','I63','I49','I44','I34','I35',
        'I45','I67','I77','I12','I64','I42','I69','I60','I70','I22','I65','I74','I24','I31','I08','I71',
        'I61','I78','I46','I27','I05','I25.1','I20.9','I25.9','I25.8','I20.0','I25.2','I21.9','I21.1',
        'I50.1','I26.9','I47.1','I73.9','I21.0','G45.9','I63.9','I34.0','I51.7','I20.8','I21.4','I12.0',
        'I50.0','I35.0','I44.7','I45.1','I47.2','I49.9','I65.2','I73.0','I77.1','I69.4','I60.9','I70.2',
        'I42.0','I35.1','I50.9','I67.9','I34.1','I22.9','I31.9','I74.3','I24.8','I49.8','I78.1','I44.2',
        'I21.2','I44.0','I51.8','I42.9','I71.4','I61.9'
    ]
    
    
    dataraw = pd.read_csv('./ukbnmr_325_data.csv')
    
   
    assessment_centers = [11012, 11021, 11011, 11008, 11009, 11024, 11020, 11018, 11010, 
                          11016, 11001, 11017, 11013, 11002, 11007, 11014, 10003, 11006, 
                          11025, 11026, 11027, 11028]
    data = dataraw[dataraw['UK Biobank assessment centre | Instance 0'].isin(assessment_centers)].copy()
    
   
    features = data.iloc[:, 67:392]
    
  
    Parallel(n_jobs=os.cpu_count() // 2)(
        delayed(process_cvd)(cvd) for cvd in cvd_ids
    )
    
    print("Processing completed for all cvd IDs.")
