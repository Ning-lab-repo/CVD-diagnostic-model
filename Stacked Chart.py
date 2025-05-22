import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

diseases = [
    "Chronic rheumatic heart diseases",
    "Hypertensive diseases",
    "Ischaemic heart diseases",
    "Other forms of heart disease",
    "Cerebrovascular diseases",
    "Diseases of arteries, arterioles and capillaries",
  #  "Transient Cerebral lschemic Attacks and Related Syndromes",
    "Pulmonary heart disease and diseases of pulmonary circulation",
]

colors = [
    '#EA4335',  
    '#FFD700',  
    '#34A853',  
    '#4169E1',  
    '#BA55D3',  
    '#FF4081',  
  #  '#03A9F4',  
    '#FF9800',  
]


file_paths = [
    './nmr_Category.csv',
    './nmr_Subcategory.csv',
    './baseline_Category.csv',
    './baseline_Subcategory.csv',

]

fig, axs = plt.subplots(4, 1, figsize=(16, 20)) 
plt.subplots_adjust(hspace=0.4)  

for idx, file_path in enumerate(file_paths):
   
    data = pd.read_csv(file_path)

   
    all_columns = data.columns.tolist()
    valid_diseases = [disease for disease in diseases if disease in all_columns]
    valid_colors = [colors[diseases.index(disease)] for disease in valid_diseases]

  
    disease_data = data[valid_diseases]
    total_heights = disease_data.sum(axis=1)

   
    data['Total'] = total_heights
    data = data.sort_values(by='Total', ascending=False).reset_index(drop=True)

   
    elements = data['Element']
    disease_data = data[valid_diseases]

  
    bar_width = 0.8  
    bottom = None  
    for col_idx, disease in enumerate(valid_diseases):
        axs[idx].bar(
            range(len(elements)),
            disease_data[disease],
            label=disease if idx == 0 else None,  
            bottom=bottom,
            color=valid_colors[col_idx],
            width=bar_width,  
        )
        bottom = disease_data[disease] if bottom is None else bottom + disease_data[disease]

 
    for i, total in enumerate(data['Total']):
        axs[idx].text(
            i, total * 1.01, f'{int(total)}', ha='center', va='bottom', fontsize=10, color='black'
        )

    axs[idx].set_ylim(0, data['Total'].max() * 1.1)

   
    axs[idx].set_xticks(range(len(elements)))
    axs[idx].set_xticklabels(elements, rotation=45, ha='right')

  
    axs[idx].spines['top'].set_visible(False)
    axs[idx].spines['right'].set_visible(False)
    axs[idx].spines['left'].set_linewidth(1.5)
    axs[idx].spines['bottom'].set_linewidth(1.5)

   
    axs[idx].set_xlim(-0.5, len(elements) - 0.5)  

   
    axs[idx].set_ylabel('Roles in diagnosing CVD', fontsize=12)


handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors]
fig.legend(
    handles,
    diseases,
   # title='Diseases',
    loc='lower center',
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    frameon=False,
    fontsize=12,
    title_fontsize=12,
)


output_file = './Stacked_Chart.pdf'  
plt.tight_layout()
plt.savefig(output_file, format='pdf', bbox_inches='tight')
plt.close()
