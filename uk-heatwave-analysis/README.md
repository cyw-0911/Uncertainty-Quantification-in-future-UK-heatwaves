# UNCERTAINTY QUANTIFICATION IN Future UK HEATWAVES
Python preprocessing pipeline for UK heatwave detection and Bayesian modelling (OpenBUGS).
This repository contains the code developed for my dissertation **“Uncertainty Quantification in Future UK Heatwaves.”**  
The Python scripts extract and process gridded climate datasets (HadUK-Grid & UKCORDEX), compute region-level heatwave metrics, and generate inputs for a hierarchical Bayesian model implemented in **OpenBUGS**.

> ## Raw data availability

The raw climate datasets used in this project (HadUK-Grid & UKCORDEX daily tasmax/tasmin)
are extremely large (tens of GB) and are distributed under the CEDA licence.
Therefore, they are **not included in this repository**.

## Repository Structure

```
src/
  demo_hadukgrid.py        # Example: reading HadUK-Grid & producing maps
  regions_and_means.py     # Region masks, spatial means, example figures
  heatwave_detection.py    # Core heatwave detection (grid → region → year)
  run_ukcordex_batch.py    # Batch pipeline for all UKCORDEX model folders

data/
  raw/                     # (Not included) CEDA-licensed NetCDF files
  processed/               # Small demonstration outputs

model/
  BUGS model files (for OpenBUGS)
```

## Pipeline Overview

### **1. Read large NetCDF climate files**
Supports:
- Daily tasmax / tasmin  
- Monthly tas  
- OSGB projection coordinates  
- Precomputed 16-region OSGB mask  

### **2. Map grid cells to UK regions**
Region masks directly provide grid indices, avoiding shapefiles.

### **3. Detect heatwaves for each grid cell**
Since the heatwave has no universal definition, I use the definition introduced by UK Met Office: 
- tasmax ≥ regional threshold  
- tasmin ≥ regional threshold  
- ≥ 3 consecutive days  
- intensity = (mean(max) + peak(max)) / 2  

### **4. Aggregate into region × year metrics**
Outputs include: (3 metrics important to indicate the heatwave performance)
- Heatwave_Count  
- Mean_Duration  
- Mean_Intensity  

### **5. Batch-process entire UKCORDEX ensemble**
Automated handling of:
- All GCM–RCM–RunID combinations  
- Multiple period files per model  
- Matching tasmax/tasmin pairs  
- Producing one summary CSV per model  

## Example Outputs Included
```
data/processed/
```
These illustrate the expected structure of temperature means and heatwave summaries.


## Bayesian Modelling in OpenBUGS

All statistical modelling for this project — including the **Bayesian hierarchical model**,  
posterior inference, uncertainty quantification, and model comparison — is performed in **OpenBUGS**.

Python is used **only for data extraction, cleaning, regional aggregation and early EDA**.  
The full Bayesian analysis begins after the Python pipeline outputs the processed regional  
heatwave metrics into `data/processed/`.

## Software Requirements

### **1. OpenBUGS (required)**
Download from:  
https://www.mrc-bsu.cam.ac.uk/software/bugs-project

OpenBUGS must be installed separately.

### **2. R + Required Packages**

Run the following in R:

```r
install.packages("R2OpenBUGS")
```
All R code used for running the Bayesian models (including the scripts that call OpenBUGS) is stored inside the `model/` directory.

## How to Run

Place raw NetCDF data under:

```
data/raw/
```

Then run:

```bash
pip install -r requirements.txt

cd src/
python demo_hadukgrid.py
python regions_and_means.py
python run_ukcordex_batch.py
```

Processed results will be saved to:

```
data/processed/
figs/
```

---

## Data Licensing

HadUK-Grid and UKCORDEX datasets are protected under CEDA licence terms.  
They are **not redistributed** here.  
Only derived outputs (plots, small CSV samples) are included.

---

## Contact
Please open an issue for questions, suggestions, or discussion.


