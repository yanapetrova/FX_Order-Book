# Data Documentation

## 1. datasetFX.mat (Available via GitHub Release)

The file `datasetFX.mat` contains one-minute predictors constructed from LMAX order book tick data.  
It is required to reproduce the empirical results in the paper.

Because the file exceeds GitHub’s upload size limits, it is stored in the Release section of this repository:

https://github.com/yanapetrova/FX_Order-Book/releases

### How to use this file
1. Download `datasetFX.mat` from the Release page.
2. Place it in a folder on your computer.
3. Set the variable `datalocation` in `use_main_functionReplication.m` to the path of this folder.

### About the data
- Data is aggregated from proprietary LMAX tick-level order book data (first five bid/ask levels).
- Aggregated at a one-minute frequency following the methodology in Table 3 of the paper.
- The aggregated dataset may be used for non-commercial research purposes.

## 2. Non-sharable data used in the Appendix

The EBS historical order book data used for certain appendix results cannot be included in this replication package.

For details on access and restrictions, see the main repository README file.
