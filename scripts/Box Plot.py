import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


data = pd.read_csv('D:/BaiduNetdiskDownload/ukbnmr_325_data.csv')



id_list = ['I12.0','I50.0','I21.2','I21.0','I25.8','I25.1','I21.1','I21.9','I50.1','I22.9','I71.4','I25.2','I73.9','I24.8','I77.1','I21.4','I25.9','I20.9','I65.2',
'I20.8','I20.0','I42.0','I50.9','I74.3','I70.2','I35.0','I63.9','I51.8','I69.4','I51.7','I42.9','I67.9']


excluded_codes = ['I79', 'I80', 'I81', 'I82', 'I83', 'I84', 
                  'I85', 'I86', 'I87', 'I88', 'I89', 'I95', 'I97', 'I98', 'I99']


data['y'] = data.apply(
    lambda row: 1 if (
        any(cvd.startswith('I') and not any(cvd.startswith(code) for code in excluded_codes) for cvd in str(row['zhen_cvd_id']).split('|')) or
        any(cvd.startswith('G45') for cvd in str(row['zhen_cvd_id']).split('|'))
    ) else 0, axis=1
)

metabolite_columns = ['IDL_CE','IDL_C','LA_pct','LDL_C_pct','LA']
#metabolite_names = ['LA/Total FA%','TG/Total Lipids-S-VLDL%','Chol/Total Lipids-M-VLDL%','TG/Total Lipids-M-VLDL%','CE/Total Lipids-M-VLDL%']
metabolite_names =['IDL_CE','IDL_C','LA_pct','LDL_C_pct','LA']
fig, axes = plt.subplots(5, 1, figsize=(10, 15))


colors = colors = [
    'lightblue', 'lightgreen', 'lightcoral', 'lightsalmon', 'lightpink', 
    'lightgray', 'skyblue', 'mediumseagreen', 'gold', 'purple', 'steelblue', 
    'dodgerblue', 'seagreen', 'darkorange', 'tomato', 'violet', 'blue', 
    'green', 'red', 'darkblue', 'darkgreen', 'darkgray', 'maroon', 
    'darkviolet', 'saddlebrown', 'mediumpurple', 'mediumslateblue', 'indianred',
    'peachpuff', 'lavender', 'palevioletred', 'beige', 'mintcream'
]



for idx, metabolite_column in enumerate(metabolite_columns):
    boxplot_data = {}  
    p_values = []  

  
    Non_CVD_data = data[data['y'] == 0].loc[:, metabolite_column].dropna()
    Non_CVD_median = Non_CVD_data.median() if len(Non_CVD_data) > 0 else None

    for cvd in id_list:

        matches = data['zhen_cvd_id'].apply(lambda x: any(part == cvd for part in str(x).split('|')))

     
        cvd_data = data[matches].loc[:, metabolite_column].dropna()
        
        if len(cvd_data) > 0:
            boxplot_data[cvd] = cvd_data.tolist()
        
        if len(cvd_data) > 0 and len(Non_CVD_data) > 0:
            _, p_value = mannwhitneyu(Non_CVD_data, cvd_data, alternative='two-sided')
        else:
            p_value = 1  
        p_values.append(p_value)


    boxplot_data['Non_CVD'] = Non_CVD_data.tolist()


    _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')
    print(p_values_corrected)

    ax = axes[idx]
    bplot = ax.boxplot(
        list(boxplot_data.values()),
        vert=True,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(color='gray'),
        whiskerprops=dict(color='gray'),
        capprops=dict(color='gray'),
        medianprops=dict(color='black')
    )


    for patch, color in zip(bplot['boxes'], colors[:len(boxplot_data)]):
        patch.set_facecolor(color)


    if Non_CVD_median is not None:
        ax.axhline(y=Non_CVD_median, color='lightgray', linestyle='--', label='Non-CVD Median')


    ax.set_ylabel(metabolite_names[idx], fontsize=10, color='black', fontname='Arial')


    if idx == len(metabolite_columns) - 1:
        ax.set_xticks(range(1, len(boxplot_data) + 1))
        ax.set_xticklabels(boxplot_data.keys(), rotation=45, fontname='Arial')
    else:
        ax.set_xticks([])


    for i, p_value in enumerate(p_values_corrected):
        if i >= len(boxplot_data) - 1: 
            break
      
        whiskers = bplot['whiskers'][2 * i + 1].get_ydata()  
        cap = bplot['caps'][2 * i + 1].get_ydata()[1] 
      
        marker_y = cap + 0.000001 * (ax.get_ylim()[1] - ax.get_ylim()[0])  
        if p_value < 0.001:
            ax.text(i + 1, marker_y, '***', horizontalalignment='center', fontsize=8, color='black')
        elif p_value < 0.01:
            ax.text(i + 1, marker_y, '**', horizontalalignment='center', fontsize=8, color='black')
        elif p_value < 0.05:
            ax.text(i + 1, marker_y, '*', horizontalalignment='center', fontsize=8, color='black')


plt.tight_layout()
plt.subplots_adjust(top=0.95)


output_path = './Box Plot by Subcategory.pdf'
plt.savefig(output_path, format='pdf')
plt.close()

