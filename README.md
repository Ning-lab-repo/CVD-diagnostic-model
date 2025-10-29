  This study maps the metabolomic landscape of cardiovascular disease (CVD) in the UK Biobank to understand shared versus disease‑specific biology and to prioritise robust cross‑CVD biomarkers. Using 244,567 participants with 325 Nightingale NMR biomarkers and rich clinical data, we defined CVD classes/subclasses by ICD‑10, performed leakage‑free preprocessing (fold‑wise imputation/standardisation), addressed imbalance with SMOTE and class weighting, and trained nested‑CV models (logistic regression, random forest, XGBoost) to quantify metabolomic distinctiveness (AUC) rather than build diagnostic tools. We interpreted models with SHAP to stabilise key features and built a diagnosis‑independent disease similarity space that combines overlap of significantly shifted metabolites and rank concordance of mean Z‑profiles. Sensitivity analyses (statin adjustment) and external validation in Scotland/Wales confirmed the robustness and generalisability of patterns. Results highlight pronounced, reproducible signatures for ischaemic and hypertensive–renal entities and a lipid/fatty‑acid–centred signal (IDL/LDL fractions, LA metrics) with complementary glycaemic, renal, hepatic, and inflammatory markers; a compact cross‑CVD metabolite panel recapitulates much of the separation across diseases. Overall, the work delivers an empirical atlas of CVD metabolomic heterogeneity, prioritises stable cross‑CVD biomarkers, and provides a reproducible framework and codebase to guide mechanistic follow‑up and future translational studies.

Data structure:
Required columns:
Participant identifier: Participant ID
CVD labels
Feature matrix: 325 NMR features + 66 clinical features
cvd_id_list contains ICD-10 codes. 
Train cohort filter:
Only England assessment centers are kept:  [11012, 11021, 11011, 11008, 11009,11024, 11020, 11018, 11010, 11016, 11001, 11017, 11013, 11002, 11007, 11014, 10003, 11006, 11025, 11026, 11027, 11028].

How the script runs (pipeline):
5-fold Stratified CV:
Within each fold it fits a SimpleImputer(median) and StandardScaler only on the training fold (no leakage) and applies to the validation fold.
Trains an XGBClassifier with class imbalance handling via scale_pos_weight = Nneg/Npos in that fold.
Collects out-of-fold predicted probabilities for all subjects to compute CV ROC/AUC and F1 (threshold 0.5 for this quick metric).



# box plot.py  
Compare the differences in metabolite distribution between different cardiovascular disease (CVD) subtypes and the non-CVD group, draw box plots and conduct statistical tests (Mann-Whitney U test), and perform FDR correction for multiple test results. Finally, save the box plots as PDF files.

# Heatmap(cvd vs noncvd).R 
This R code implements the metabolite data based on the UK Biobank samples, and creates a heatmap for the metabolite expression of both cardiovascular disease (CVD) and non-CVD populations.

# Heatmap(Sampling 100).R 
This R code's function is to conduct sample screening and sampling (100) based on the NMR metabolite data from UK Biobank, and to draw a metabolite heatmap with group annotations related to cardiovascular disease (CVD) diagnostic information.

# Horizontal tree diagram.R 
Construct and draw a hierarchical tree diagram (tree structure diagram), which can visually display the relationships and quantities of various diseases and their subcategories.

# Lollipop Chart.R 
Based on the AUC indicators in the cardiovascular disease (CVD) data, a sorted lollipop chart is drawn to visually present the magnitudes and differences of the AUC values for different diseases.

# model.py 
Model training, evaluation and calculation of statistical confidence intervals based on the list of diagnostic codes for cardiovascular diseases (CVD), as well as parallel processing of multiple disease codes.

# pie.R 
Segment the data by AUC values, count the quantities, and draw a pie chart.

# roc.py 
Estimate the confidence interval of AUC using the bootstrap method, and plot the ROC curves for each disease and save the figures as PDF files.

# Sankey Diagram.R 
Based on three groups of nodes (the left feature column, the middle disease ID, and the right feature column in the middle), a sorted Sankey Diagram is constructed to visualize the flow relationships among the three layers of variables. Finally, the diagram is saved as HTML and PNG files.

# shap.py 
Generate SHAP summary charts and numerical values for key features, process multiple CVD disease codes in batches, and save the results as PDF and CSV files.

# Stacked Chart.py 
Draw a stacked bar chart to show the diagnostic role of "different elements" in each file for these diseases, and finally generate a 4-row 1-column multi-subplot PDF chart.
