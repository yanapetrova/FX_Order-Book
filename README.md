# Replication Package for  
## "Assessing Cross-Currency Predictability in Forex Markets: Insights from Limit Order Book Data"

Authors:  
Anders Vilhelmsson (Lund University)  
Yana Petrova (Aarhus University)

Replication package assembled: 17 November 2025

---

## 1. Repository Structure

```
FX_Order-Book/
├── code/                 # MATLAB scripts for running the replication
├── data/                 # Documentation for datasetFX.mat
└── output/               # Documentation for Case_*.mat intermediate files
```

### /code
Contains all MATLAB code required to replicate the results:
- `use_main_functionReplication.m`
- `main_function_replication.m`
- `SPCAdr.m`

### /data
Contains `README_DATA.md`, which describes how to obtain and use `datasetFX.mat`.

### /output
Contains `README_OUTPUT.md`, which explains how to obtain and use all `Case_*.mat` files.

---

## 2. Computing Environment

The replication was performed using the following software and hardware:

**Software**
- MATLAB 2022b

**MATLAB Toolboxes**
- System Identification Toolbox 10.0  
- Statistics and Machine Learning Toolbox 12.4  

**Hardware**
- Windows 10  
- AMD Ryzen 9 5900X (12 cores)  
- NVIDIA RTX 2080  
- 64 GB RAM  

**Expected runtimes**
- Hourly-frequency estimation: approximately 45 minutes  
- One-minute-frequency estimation: approximately 85 hours  
- Using intermediate `Case_*.mat` files: a few minutes  

---

## 3. Data

### 3.1 datasetFX.mat  
Contains one-minute aggregated predictors derived from LMAX order book tick data.  
Because the file exceeds GitHub’s size limits, it is available through the Release page:

https://github.com/yanapetrova/FX_Order-Book/releases

### 3.2 Case_*.mat files  
Ten intermediate model estimation files used to reproduce the results without recomputing lengthy estimations.  
Available through the Release page.

### 3.3 Non-sharable data  
Some appendix results use EBS historical order book data, which cannot be distributed due to licensing restrictions.

---

## 4. Running the Replication

1. Download this repository.
2. Download `datasetFX.mat` from the Release page and specify its location by setting `datalocation` in:
   ```
   use_main_functionReplication.m
   ```
3. (Optional) To skip long estimation runs, download all `Case_*.mat` files and set the folder path using `savepath`.
4. In MATLAB, run:
   ```
   use_main_functionReplication
   ```
5. Output tables and figures will be saved automatically in the folder specified by `savepath`.

---

## 5. Notes on Reproducibility

- All scripts required to reproduce results are included in `/code`.
- Large data files are stored in the Release section.
- The repository is organized according to IJF reproducibility standards.
- All tables and figures generated correspond to the manuscript results.

---

## 6. Contact

For questions regarding the replication package, data, or methodology, please contact:

- Yana Petrova, Aarhus University — ypetrova@econ.au.dk  
- Anders Vilhelmsson, Lund University

