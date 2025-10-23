# longitudinal_EEG_timescales
This repository contains the data and code used to produce the results published in XXX

## 📂 Repository structure

| Folder / File | Purpose |
|---|---|
| `Data/` | Contains processed datasets and intermediate files used by the analysis scripts |
| `1_data_matlab_to_python.py` | Converts MATLAB outputs into Python-friendly formats (e.g. `.csv` or `.pkl`) |
| `2_select_subj_random.py` | Splits participants into exploratory and validation cohorts (random assignment) |
| `3_eeg_intrinsic_timescales_estimation_pipeline.py` | Main pipeline to estimate intrinsic timescales per subject, handle outliers, and aggregate results |
| `4_eeg_ITS_GroupComparisons&LinearModels_pipeline.py` | Runs group comparisons, mixed models, and associations with rhythm features |
| `linear_model_pipeline.py` | Helper script for more general linear modeling |
| `tau_estimate.py` | Functions relating to estimating tau / intrinsic timescales metrics |
| `LICENSE` | MIT license for reuse and redistribution |

---

## 🧪 How to run analyses / reproduce results

Below is a suggested order. Make sure all dependencies (e.g. `numpy`, `pandas`, `statsmodels`, `scipy`) are installed in your environment.

1. **Convert MATLAB outputs**  
   python 1_data_matlab_to_python.py
   This creates data files in a Python format that downstream scripts will use.

2. **Split participants into random exploratory and validation groups**
    python 2_select_subj_random.py
    Random seed is fixed to ensure reproducibility

2. **Estimate timescales per subject**
    python 3_eeg_intrinsic_timescales_estimation_pipeline.py
    This is the heart of the pipeline: it computes intrinsic timescales, filters out outliers, and creates summary data tables.
    Bespoke functions used in this pipeline are defined in the script tau_estimate.py

    The long format data in output are passed on to Matlab for inputation (see script XXX). Processed timescales data used in further anlaysis are shared in the compressed file /Data/Longitudinal Data.zip.

3. **Modeling developmental trajectory**
    LLM XXXXX

4. **Correlation with alpha properties**
    XXXX

5. **Comparison with adults**
    python 4_eeg_ITS_GroupComparisons&LinearModels_pipeline.py
    This takes the precomputed subject-level metrics and runs the comparisons and linear models between infants and adults.


📋 Requirements & environment

Python 3.x (recommended ≥ 3.8)

Common packages: numpy, pandas, scipy, statsmodels, matplotlib / seaborn (for plotting)

(If any MATLAB-to-python conversion uses scipy.io or similar, ensure that is installed)


🧠 Notes

Random seeds are fixed to ensure reproducibility of subject group assignments.

File paths can be updated directly in the scripts if you wish to adapt the workflow for your own data.

All scripts are written to be modular, so you can import specific functions for your own analyses (e.g., from tau_estimate.py).

📄 Licensing & citation

This code is released under the MIT license (see LICENSE), meaning you’re free to use, modify, and share it as long as you include attribution.

If you use this repository or adapt the code in your work, please cite the associated manuscript (or preprint) and include a reference to this repository (e.g. via DOI or its GitHub URL).