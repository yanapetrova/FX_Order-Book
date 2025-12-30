# Output Files and Intermediate Results

The replication uses the following intermediate estimation result files (available in the Release section):

Case 1:
- Case_1_12-Nov-2025troll_7200freq_minutely_samp_size=201037_lasso_pca_spca_ada_lassorf.mat
- Case_1_19-Nov-2025troll_1440freq_hourly_samp_size=201037_lasso_pca_spca_ada_lassorf.mat

Case 2:
- Case_2_13-Nov-2025troll_7200freq_minutely_samp_size=201037_lasso_pca_spca_ada_lassorf.mat
- Case_2_19-Nov-2025troll_1440freq_hourly_samp_size=201037_lasso_pca_spca_ada_lassorf.mat

Case 3:
- Case_3_14-Nov-2025troll_7200freq_minutely_samp_size=201037_lasso_pca_spca_ada_lassorf.mat
- Case_3_19-Nov-2025troll_1440freq_hourly_samp_size=201037_lasso_pca_spca_ada_lassorf.mat

Case 4:
- Case_4_15-Nov-2025troll_7200freq_minutely_samp_size=201037_lasso_pca_spca_ada_lassorf.mat
- Case_4_19-Nov-2025troll_1440freq_hourly_samp_size=201037_lasso_pca_spca_ada_lassorf.mat

Case 5:
- Case_5_16-Nov-2025troll_7200freq_minutely_samp_size=201037_lasso_pca_spca_ada_lassorf.mat
- Case_5_19-Nov-2025troll_1440freq_hourly_samp_size=201037_lasso_pca_spca_ada_lassorf.mat

These files are available through the Release page of this repository:

https://github.com/yanapetrova/FX_Order-Book/releases

In the same reciprocity you can find the PC_variation.explained, which provides the inputs and the final graphs shown in Figures 1 and 2 of the paper.

### How to use these files
1. Download all Case_*.mat files from the Release page.
2. Place them into the folder you will use as `savepath`.
3. Set `savepath` in `use_main_functionReplication.m` accordingly.
4. Set `reEstimateResults = 0` to use the intermediate files.
