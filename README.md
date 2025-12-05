# coralline_10years
Title: Local Resilience Amid Climate-Driven Stress: Divergent Coralline Barren Trajectories in the Context of Korea’s Coastal Restoration Initiatives

This repository contains all Python scripts, processed datasets, and analytical outputs used in the study:

**“Standardized Panel Indicators for Assessing Coastal Barren Grounds Around the Korean Peninsula.”**

The repository is organized to allow full reproducibility of the analyses, from raw data processing to final statistical outputs and figures.
Data and scripts are structured into **three stages of upload**:

1. **1st upload:** Raw and base-level structured survey datasets
2. **2nd upload:** Depth-adjusted environmental datasets & merged panel data
3. **3rd upload:** Admin-level trend analyses, correlation matrices, decoupling frameworks, and integrated visualization outputs

---

# **1. Repository Structure**

```
├── code/                     # Python scripts and notebooks for analysis
├── data/                     # All CSV data files (raw → processed → final)
├── figures/                  # All generated PNG figures from step 3
├── docs/ (optional)          # Manuscript drafts, supplementary notes
└── README.md                 # This file
```

### **Main scripts (examples)**

* `01_load_and_clean.ipynb` – Cleaning raw survey data
* `02_merge_env_and_depth.ipynb` – Merging depth-adjusted environmental variables
* `03_panel_and_slope_analysis.ipynb` – Constructing panel datasets, trend estimation, correlations
* `04_visualization_step3.ipynb` – Generating decoupling matrices, integrated trend maps, and heatmaps

---

# **2. Data Uploads and Contents**

## **2.1. First Upload — Raw & Base Survey Data**

### **`sc_coralline_covers_col.csv`**

Standardized field-survey dataset containing quadrat-level or station-level observations.

Includes:

* Station ID
* Year/quarter
* Depth
* Cover of coralline algae, macroalgae, barren ground classification (CB stage)
* Notes on ecological conditions

Used to build the foundation of the **CB3 panel dataset**.

---

### **`sc_filtered_stid_under30m_b.csv`**

Filtered set of stations meeting analysis criteria (e.g., ≤ 30 m depth).

Includes:

* Station metadata (lat/lon, depth, admin unit)
* Temporal coverage

Serves as the **panel skeleton** for subsequent merging with environmental variables.

---

## **2.2. Second Upload — Depth-Adjusted Environmental Dataset (Step 1)**

### **`step1_stid_depth_adjusted_env.csv`**

Depth-corrected environmental time series for each station.

Variables include:

* Temperature
* Salinity
* DO
* Chl-a
* DIN / DIP
* NH3-N, NO2-, NO3-
* SS, SAS
* TN, TP
* pH

This dataset is used to construct station-level panels for CB3–environment linkage analysis.

---

## **2.3. Third Upload — Trend Analyses, Correlations & Integrated Frameworks (Step 2/3 Outputs)**

### **Panel-level merged data**

**`step2_stid_quarter_adj_panel.csv`**
Integrated station × quarter panel dataset for CB3 and all environmental variables.

---

### **Trend & correlation summaries**

| File                                            | Description                                                          |
| ----------------------------------------------- | -------------------------------------------------------------------- |
| `step2_admin_CB_slope_summary.csv`              | Admin-level CB3 trend (slope) summary                                |
| `step2_admin_env_slope_summary.csv`             | Admin-level environmental trend summary                              |
| `step2_admin_env_trend_summary.csv`             | Categorized (increase/decrease/stable) environmental trends          |
| `step2_admin_CB_env_slope_long.csv`             | Long-format slope dataset for CB3 vs each environmental variable     |
| `step2_admin_CB_env_timeseries_correlation.csv` | Within-admin time series correlations (CB3 vs environment)           |
| `step2_CB3_env_correlation_table.csv`           | Correlation matrix (CB3 vs environmental variables)                  |
| `step2_stid_quarter_slope_summary.csv`          | Station-level slope summaries                                        |
| `step3_env_trend_slope_matrix.csv`              | Raw & normalized slope matrices for environmental trends             |
| `step3_env_CB3_combined_table.csv`              | Integrated trend combination matrix (env↑/↓ × CB3↑/↓ classification) |
| `step3_region_env_slope_summary.csv`            | Regional (East, South, West, Jeju) mean environmental trends         |

---

# **3. Visualization Outputs (Step 3 Figures)**

All figures are generated directly from repository scripts and stored under `figures/`.
Below is a catalog of the major outputs.

---

## **3.1. Admin-level CB3 trends (barren-ground intensification / alleviation)**

**`step3_CB3_admin_trend_barplot.png`**

Shows 2014–2024 CB3 slope per admin unit.
Blue = East, Orange = Jeju, Green = South, Red = West.

* Negative slope → reduction of severe barren ground (improvement)
* Positive slope → expansion of barren ground (intensification)

---

## **3.2. CB3–environment time-series correlations**

**`step3_CB3_env_correlation_heatmap.png`**

Heatmap of correlations between CB3 and key environmental variables.
Blue = negative correlation, Red = positive correlation.

Reveals whether barren ground dynamics are synchronous with environmental forcing.

---

## **3.3. Slope–slope scatter matrix (CB3 trend vs. environmental trend)**

**`step3_CB3_vs_env_slope_scatter_matrix.png`**

Each subplot represents the relationship between CB3 slope and the slope of a single environmental variable.
Used to identify consistent pressure–response patterns (or lack thereof).

---

## **3.4. Environment–Ecology Decoupling Matrix**

**`step3_Decoupling_Matrix_Temperature.png`**

Plots CB3 slope vs. temperature trend to reveal four zones:

* **Intensification:** env↑ & CB3↑
* **Decoupling (mitigation resilience):** env↑ & CB3↓
* **Joint improvement:** env↓ & CB3↓
* **Counter-intuitive:** env↓ & CB3↑

This framework quantifies whether ecological conditions track or diverge from environmental stress.

---

## **3.5. Integrated environmental–CB3 trend combination matrix**

**`step3_env_CB3_integrated_matrix.png`**

Categorizes each admin unit × variable pair into:

* env↑ & CB3↑ (intensification)
* env↑ & CB3↓ (decoupling)
* env↓ & CB3↓ (joint improvement)
* env↓ & CB3↑ (counter-intuitive)
* no trend / no data

This matrix provides a comprehensive overview of ecological–environmental coherence.

---

## **3.6. Environmental trend heatmaps (raw & normalized slopes)**

### **`step3_env_trend_slope_heatmap_norm.png`**

Normalized (0–1) slopes for all environmental variables → compares relative stress intensities across regions.

### **`step3_env_trend_slope_heatmap_raw.png`**

Raw slopes retaining physical units → compares absolute magnitudes of change.

---

## **3.7. Regional mean environmental trend summary**

**`step3_region_env_slope_barplot.png`**

Shows average slope of each environmental variable across four regions:

* East
* Jeju
* South
* West

Highlights broad-scale regional environmental patterns (e.g., nutrient declines, salinity trends, warming rates).

---

# **4. Replication Guide**

### **Clone the repository**

```bash
git clone https://github.com/<YOUR_ID>/<REPO_NAME>.git
cd <REPO_NAME>
```

### **Set up environment**

```bash
conda create -n gall_sea python=3.11
conda activate gall_sea
pip install -r requirements.txt
```

### **Run analyses (recommended order)**

1. **Raw data pre-processing**
   `01_load_and_clean.ipynb`

2. **Merge depth-adjusted environmental data**
   `02_merge_env_and_depth.ipynb`

3. **Construct panel dataset & compute slopes**
   `03_panel_and_slope_analysis.ipynb`

4. **Generate all figures and integrated matrices**
   `04_visualization_step3.ipynb`

---

# **5. Contact & Support**

* For questions, please open an **Issue** in this repository.
* For academic inquiries, contact the corresponding author listed in the manuscript.
* Additional documentation will be added progressively under `docs/`.


