import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, f1_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import joblib
import os
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    roc_curve, auc, accuracy_score, f1_score,
    matthews_corrcoef, confusion_matrix, precision_score, recall_score
)
from joblib import Parallel, delayed

if 'DISPLAY' in os.environ:
    del os.environ['DISPLAY']

base_path = 'D:/BaiduNetdiskDownload'
if not os.path.exists(base_path):
    os.makedirs(base_path)
dataraw = pd.read_csv('./ukbnmr_325_data.csv.csv')


cvd_id_list = ['I10','I25','I20','I48','I21','I50','I47','I73','I26','G45','I51','I63','I49','I44','I34','I35','I45','I67','I77','I12','I64','I42','I69','I60','I70','I22',
    'I65','I74','I24','I31','I08','I71','I61','I78','I46','I27','I05','I25.1','I20.9','I25.9','I25.8','I20.0','I25.2','I21.9','I21.1','I50.1','I26.9','I47.1',
    'I73.9','I21.0','G45.9','I63.9','I34.0','I51.7','I20.8','I21.4','I12.0','I50.0','I35.0','I44.7','I45.1','I47.2','I49.9','I65.2','I73.0','I77.1','I69.4',
    'I60.9','I70.2','I42.0','I35.1','I50.9','I67.9','I34.1','I22.9','I31.9','I74.3','I24.8','I49.8','I78.1','I44.2','I21.2','I44.0','I51.8','I42.9','I71.4','I61.9'] 

assessment_centers = [11012, 11021, 11011, 11008, 11009,11024, 11020, 11018, 11010, 11016, 11001, 11017, 11013, 11002, 11007, 11014, 10003, 11006, 11025, 11026, 11027, 11028]
data = dataraw[dataraw['UK Biobank assessment centre | Instance 0'].isin(assessment_centers)].copy()
X_all = data.iloc[:, 67:392]

data['_id_list'] = data['diagnosis_cvd_id'].astype(str).str.split('|')

def process_cvd(cvd, n_threads):
    if '.' not in cvd:
        y = data['diagnosis_cvd_id_list'].apply(lambda codes: int(any(code.startswith(cvd) for code in codes)))
    else:
        y = data['diagnosis_cvd_id_list'].apply(lambda codes: int(any(code == cvd for code in codes)))
  
    
    model = LogisticRegression(C=0.1, penalty='l2', solver='lbfgs', class_weight='balanced', max_iter=10000, random_state=42)
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_pred_prob = np.zeros(len(y), dtype=float)
    
    for train_idx, val_idx in cv.split(X_all, y):
        X_train, X_val = X_all.iloc[train_idx], X_all.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]
        imputer = SimpleImputer(strategy='median')
        X_train_imputed = imputer.fit_transform(X_train)
        X_val_imputed = imputer.transform(X_val)
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train_imputed)
        X_val_scaled = scaler.transform(X_val_imputed)
        model.fit(X_train_scaled, y_train)
        y_pred_prob[val_idx] = model.predict_proba(X_val_scaled)[:, 1]
    
    fpr, tpr, _ = roc_curve(y, y_pred_prob)
    roc_auc = auc(fpr, tpr)
    f1 = f1_score(y, (y_pred_prob > 0.5).astype(int))
    
    try:
        plt.rcParams['font.sans-serif'] = 'Arial'
        plt.figure()
        plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'LR (AUC = {roc_auc:.4f})')
        plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate', fontsize=16)
        plt.ylabel('True Positive Rate', fontsize=16)
        plt.title('ROC Curve', fontsize=16, pad=10)
        plt.legend(loc='lower right', fontsize=14)
        plt.savefig(os.path.join(base_path, f'LR_Z_D_roc_curve_{cvd}.png'))
        plt.close()
    except Exception as e:
        print(f"Error in plotting the ROC curve ({cvd}): {e}")
    

    
    results_df = pd.DataFrame({
        'Patient_ID': data['Participant ID'],
        'Predicted_Score': y_pred_prob,
        'Actual_Label': y
    })
    results_df.to_csv(os.path.join(base_path, f'LR_Z_D_patient_predictions_{cvd}.csv'), index=False)
    X_all_values = X_all.values
    y_all = y.values
    
 
    final_imputer = SimpleImputer(strategy='median')
    X_all_imputed = final_imputer.fit_transform(X_all_values)
    

    final_scaler = StandardScaler()
    X_all_scaled = final_scaler.fit_transform(X_all_imputed)
    

    final_model = LogisticRegression(C=0.1, penalty='l2', solver='lbfgs', class_weight='balanced', max_iter=10000, random_state=42)
    final_model.fit(X_all_scaled, y_all)
    
    # Save final model, imputer, and scaler
    joblib.dump(final_model, os.path.join(base_path, f'LR_Z_D_final_LR_model_{cvd}.joblib'))
    joblib.dump(final_imputer, os.path.join(base_path, f'LR_Z_D_final_simple_imputer_{cvd}.joblib'))
    joblib.dump(final_scaler, os.path.join(base_path, f'LR_Z_D_final_scaler_{cvd}.joblib'))
    return {'cvd_id': cvd, 'AUC': roc_auc, 'F1 Score': f1}


n_cores = os.cpu_count()  # 192
n_tasks = len(cvd_id_list)  # 87
n_threads_per_task = max(1, n_cores // n_tasks)  # 192 // 87 = 2

results_list = Parallel(n_jobs=n_tasks)(
    delayed(process_cvd)(cvd, n_threads_per_task) for cvd in cvd_id_list
)

results = pd.DataFrame(results_list)
results.to_csv(os.path.join(base_path, 'LR_Z_D_cvd_model_results.csv'), index=False)
print(results.head())




np.random.seed(42)


BOOTSTRAP_ITERATIONS = 1000
CONFIDENCE_LEVEL = 0.95 



def format_value(value):
    if np.isnan(value):
        return 'N/A'
    else:
        return f"{value:.4f}"


def compute_metrics(y_true, y_pred, sample_weight=None):
    acc = accuracy_score(y_true, y_pred, sample_weight=sample_weight)
    sn = recall_score(y_true, y_pred, sample_weight=sample_weight)  
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, sample_weight=sample_weight).ravel()
    sp = tn / (tn + fp) if (tn + fp) > 0 else np.nan  # Specificity
    ppv = precision_score(y_true, y_pred, sample_weight=sample_weight, zero_division=0)  # Positive Predictive Value
    npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan  # Negative Predictive Value
    mcc = matthews_corrcoef(y_true, y_pred)
    f1 = f1_score(y_true, y_pred, sample_weight=sample_weight)
    return acc, sn, sp, ppv, npv, mcc, f1


def process_single_cvd(cvd):

    file_path = os.path.join(base_path, f'LR_Z_D_patient_predictions_{cvd}.csv')
    df = pd.read_csv(file_path)
    label = df['Actual_Label']
    score = df['Predicted_Score']
    numy = (label == 1).sum()


    class_weights = label.value_counts(normalize=True)
    weights = label.map({0: 1 / class_weights[0], 1: 1 / class_weights[1]})


    fpr, tpr, thresholds = roc_curve(label, score)
    roc_auc = auc(fpr, tpr)


    youden_index = tpr - fpr
    optimal_idx = np.argmax(youden_index)
    optimal_threshold = thresholds[optimal_idx]


    predictions = (score >= optimal_threshold).astype(int)


    acc, sn, sp, ppv, npv, mcc, f1 = compute_metrics(label, predictions, sample_weight=weights)


    bootstrap_metrics = {
        'AUC': [],
        'Accuracy': [],
        'Sensitivity': [],
        'Specificity': [],
        'PPV': [],
        'NPV': [],
        'MCC': [],
        'F1 Score': []
    }


    for i in range(BOOTSTRAP_ITERATIONS):
 
        indices = np.random.randint(0, len(label), len(label))
        label_sample = label.iloc[indices]
        score_sample = score.iloc[indices]
        weights_sample = weights.iloc[indices]


        if len(np.unique(label_sample)) < 2:
            continue


        fpr_sample, tpr_sample, thresholds_sample = roc_curve(label_sample, score_sample)
        roc_auc_sample = auc(fpr_sample, tpr_sample)
        bootstrap_metrics['AUC'].append(roc_auc_sample)

        youden_index_sample = tpr_sample - fpr_sample
        optimal_idx_sample = np.argmax(youden_index_sample)
        if optimal_idx_sample >= len(thresholds_sample):
            continue
        optimal_threshold_sample = thresholds_sample[optimal_idx_sample]


        predictions_sample = (score_sample >= optimal_threshold_sample).astype(int)


        try:
            acc_sample, sn_sample, sp_sample, ppv_sample, npv_sample, mcc_sample, f1_sample = compute_metrics(
                label_sample, predictions_sample, sample_weight=weights_sample)
            bootstrap_metrics['Accuracy'].append(acc_sample)
            bootstrap_metrics['Sensitivity'].append(sn_sample)
            bootstrap_metrics['Specificity'].append(sp_sample)
            bootstrap_metrics['PPV'].append(ppv_sample)
            bootstrap_metrics['NPV'].append(npv_sample)
            bootstrap_metrics['MCC'].append(mcc_sample)
            bootstrap_metrics['F1 Score'].append(f1_sample)
        except:
            continue


    def get_confidence_interval(data, confidence=CONFIDENCE_LEVEL):
        lower = np.percentile(data, (1 - confidence) / 2 * 100)
        upper = np.percentile(data, (1 + confidence) / 2 * 100)
        return lower, upper


    metrics_ci = {}
    for metric, values in bootstrap_metrics.items():
        if len(values) == 0:
            metrics_ci[metric] = (np.nan, np.nan)
        else:
            lower, upper = get_confidence_interval(values)
            metrics_ci[metric] = (lower, upper)


    results = {
        'cvd': cvd,
        'numy': numy,
        'AUC': (roc_auc, metrics_ci['AUC'][0], metrics_ci['AUC'][1]),
        'Accuracy': (acc, metrics_ci['Accuracy'][0], metrics_ci['Accuracy'][1]),
        'Sensitivity (Sn)': (sn, metrics_ci['Sensitivity'][0], metrics_ci['Sensitivity'][1]),
        'Specificity (Sp)': (sp, metrics_ci['Specificity'][0], metrics_ci['Specificity'][1]),
        'PPV': (ppv, metrics_ci['PPV'][0], metrics_ci['PPV'][1]),
        'NPV': (npv, metrics_ci['NPV'][0], metrics_ci['NPV'][1]),
        'MCC': (mcc, metrics_ci['MCC'][0], metrics_ci['MCC'][1]),
        'F1 Score': (f1, metrics_ci['F1 Score'][0], metrics_ci['F1 Score'][1])
    }


    for metric in ['AUC', 'Accuracy', 'Sensitivity (Sn)', 'Specificity (Sp)', 'PPV', 'NPV', 'MCC', 'F1 Score']:
        point, ci_lower, ci_upper = results[metric]
        print(f"{metric}: {format_value(point)} (95% CI: {format_value(ci_lower)} - {format_value(ci_upper)})")


    auc_str = f"{format_value(roc_auc)} ({format_value(metrics_ci['AUC'][0])} - {format_value(metrics_ci['AUC'][1])})"
    acc_str = f"{format_value(acc)} ({format_value(metrics_ci['Accuracy'][0])} - {format_value(metrics_ci['Accuracy'][1])})"
    sn_str = f"{format_value(sn)} ({format_value(metrics_ci['Sensitivity'][0])} - {format_value(metrics_ci['Sensitivity'][1])})"
    sp_str = f"{format_value(sp)} ({format_value(metrics_ci['Specificity'][0])} - {format_value(metrics_ci['Specificity'][1])})"
    ppv_str = f"{format_value(ppv)} ({format_value(metrics_ci['PPV'][0])} - {format_value(metrics_ci['PPV'][1])})"
    npv_str = f"{format_value(npv)} ({format_value(metrics_ci['NPV'][0])} - {format_value(metrics_ci['NPV'][1])})"
    mcc_str = f"{format_value(mcc)} ({format_value(metrics_ci['MCC'][0])} - {format_value(metrics_ci['MCC'][1])})"
    f1_str = f"{format_value(f1)} ({format_value(metrics_ci['F1 Score'][0])} - {format_value(metrics_ci['F1 Score'][1])})"

    return {
        'cvd': cvd,
        'numy': numy,
        'AUC': auc_str,
        'Accuracy': acc_str,
        'Sensitivity (Sn)': sn_str,
        'Specificity (Sp)': sp_str,
        'PPV': ppv_str,
        'NPV': npv_str,
        'MCC': mcc_str,
        'F1 Score': f1_str
    }

final_results = Parallel(n_jobs=os.cpu_count() // 2)(
    delayed(process_single_cvd)(cvd) for cvd in cvd_id_list
)


final_df = pd.DataFrame(final_results)


print(final_df)


final_df.to_csv(os.path.join(base_path, 'LR_auc_Statistical_Results.csv'), index=False)
print(f"Results saved to {base_path}")
