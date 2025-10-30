UK Biobank NMR Metabolic Biomarkers (325 variables) — Data Source

Source：
Ritchie SC, Surendran P, Karthikeyan S, et al. Quality control and removal of technical variation of NMR metabolic biomarker data in ~120,000 UK Biobank participants. Scientific Data. 2023;10:64. doi:10.1038/s41597-023-01949-y.

The original Nightingale Health NMR release contains 249 variables:

168 absolute concentrations = 107 “non-derived” (directly measured) + 61 “composite” (sums of non-derived), and

81 ratios provided in the release. 

Ritchie et al. removed technical variation and re-computed the 61 composites and 81 ratios from the 107 base markers, then added 76 new biomarker ratios (e.g., lipid/cholesterol/fatty-acid percentages for VLDL/LDL/HDL classes and total serum), yielding 249 + 76 = 325 biomarkers that are widely used in downstream studies. Formulae and variable lists are in their Supplementary Tables; code is in the ukbnmr R package. 
