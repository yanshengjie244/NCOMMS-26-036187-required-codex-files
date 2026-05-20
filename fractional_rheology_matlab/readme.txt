Fractional Fit Results - MATLAB Code

This folder contains MATLAB scripts and data files for fitting and visualizing
fractional-order and variable-order rheological models.

Folder structure
----------------
Main_codex
  Main analysis scripts and figure-generation data for the fractional-order
  models used in the main results.

Supp_codex
  Supplementary fitting scripts and data for grouped rheology datasets.

Main scripts
------------
Main_codex/Fit_fractional_order.m
  Fits experimental rheology data using IFDM, VFDM, FDM, and SKM models, and
  plots the fitted curves together with the order trend alpha(t).

Main_codex/Fit_data_IFDM_VFDM_FDM.m
  Fits and compares IFDM, VFDM, and FDM models for the selected dataset.

Main_codex/Fit_data_ML_model.m
  Fits and compares IFDM, VMLM, VFDM, FDM, and SKM model variants.

Main_codex/Model_evolution_fractional_order.m
  Generates model-evolution curves for different fractional orders.

Main_codex/Prediction_Fig.4a-b_data_IFDM_HMLM.m
Main_codex/Prediction_Fig.4c-h_data_IFDM_HMLM.m
  Fits the early-time region and predicts the later-time region for Fig. 4
  datasets, reporting prediction accuracy by MAPE.

Supplementary scripts
---------------------
Supp_codex/fit_Supp_group*_HBM_VFDM_HBJM_HBCM.m
Supp_codex/fit_Supp_data2_group*_HBM_VFDM_HBJM_HBCM.m
  Group-specific launch scripts for supplementary datasets.

Supp_codex/run_group_HBM_VFDM_HBJM_HBCM.m
  Shared fitting routine for HBM, VFDM, HBJM, and HBCM models.

Data files
----------
Main_codex/Fig.4a-b_data.xlsx
Main_codex/Fig.4c-d_data.csv
Main_codex/Fig.4e-f_data.csv
Main_codex/Fig.4g-h_data.csv
Supp_codex/Supp_data_1.xlsx
Supp_codex/Supp_data_2.xlsx

Notes
-----
1. Run each script from the folder that contains its required data file, or add
   that folder to the MATLAB path before running.
2. The scripts use MATLAB Optimization Toolbox functions such as lsqcurvefit.
3. Figures are generated in MATLAB figure windows. Most scripts do not save
   figures automatically.
4. Code comments have been converted to English. Model equations, variables,
   data ranges, and plotting logic were kept unchanged.
