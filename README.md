# **Revised README Content**

## **Repository: coralline_10years**

### **Title:**

**Local Resilience Amid Climate-Driven Stress: Divergent Coralline Barren Trajectories in the Context of Korea’s Coastal Restoration Initiatives**

This repository contains all Python scripts, processed datasets, and analytical outputs used in the study:

**“Standardized Panel Indicators for Assessing Coastal Barren Grounds Around the Korean Peninsula.”**

The repository is structured for **full reproducibility**, covering all steps from raw data processing to final statistical outputs and visualizations.
All datasets are uploaded in three chronological stages:

1. **First upload:** Raw and base-level field/survey datasets
2. **Second upload:** Depth-adjusted environmental datasets (Step 1 outputs)
3. **Third upload:** Admin-level trend analyses, correlations, decoupling matrices, and visualization outputs (Step 2–3)

---

# **1. Repository Structure**

```
📁 repository_root/
│
├── run_all_steps.py
│
├── step1_stid_depth_adjusted_env.csv
├── step2_tab2_admin_CB3_slope_summary.csv
├── step2_CB3_env_correlation_table.csv
├── env_trend_slope_matrix.csv
├── env_CB3_combined_table.csv
├── region_env_slope_summary.csv
│
├── fig3_Normalized_slopes_of_env._variables.png
├── fig4_Trends_of_CB3_by_municipality_2014-2024.png
├── fig5_Relationships_between_environmental_variables_and_CB3_change.png
├── fig6_Environment-ecology_decoupling_matrix_for_temperature.png
├── fig7_Heatmap_for_municipality-env._variable_in_CB3.png
├── fig8_Integrated_trend_matrix_of_environmental_change_and_CB3_intensification.png
│
└── 📁 src_data/
      ├── src_filtered_stid_under30m_b.csv
      └── src_coralline_covers_col.csv
```

### **Main Scripts (Examples)**

run_all_steps.py

This is the complete analysis pipeline integrating all three stages:

Step 1 – Depth Adjustment
Performs per-station regression to remove depth effects from environmental variables.
Output: step1_stid_depth_adjusted_env.csv

Step 2 – Trend & Correlation Analysis
Builds quarterly panels, computes environmental/CB3 slopes, and generates CB3–environment correlation matrices.
Output examples:

step2_tab2_admin_CB3_slope_summary.csv

step2_CB3_env_correlation_table.csv

Step 3 – Visualization & Integrated Matrices
Generates manuscript-ready figures and integrated ecological–environmental frameworks.
Output examples:

fig3_Normalized_slopes_of_env._variables.png

fig4_Trends_of_CB3_by_municipality_2014-2024.png

fig6_Environment-ecology_decoupling_matrix_for_temperature.png

env_CB3_combined_table.csv

---

# **2. Data Uploads and Contents**

## **2.1. First Upload — Raw & Structured Survey Data**

### **`sc_coralline_covers_col.csv`**

Field-survey dataset containing quadrat- or station-level observations.

**Includes:**

* Station ID
* Year / quarter
* Depth
* Coralline algae cover, macroalgae cover
* CB stage classification (normal / progressive / critical)
* Notes on ecological conditions

Used to construct the long-term **CB1/CB2/CB3 panel dataset**.

---

### **`sc_filtered_stid_under30m_b.csv`**

Filtered stations satisfying analysis criteria (≤ 30 m depth).

**Includes:**

* lat/lon, depth, admin unit
* survey period
  Forms the **panel backbone for merging with environmental variables**.

---

## **2.2. Second Upload — Depth-Adjusted Environmental Dataset (Step 1 Output)**

### **`step1_stid_depth_adjusted_env.csv`**

Station-level environmental time series with depth effects removed via per-station linear regression.

**Variables include (all depth-adjusted):**

* Temperature
* Salinity
* DO
* Chl-a
* pH
* Nutrients: DIN, DIP, NO₂⁻, NO₃⁻, NH₃-N
* TN, TP
* SS, SAS

This dataset is used to build quarterly environmental panels used in **CB3–environment linkage analysis**.

---

## **2.3. Third Upload — Trend Analyses, Correlations & Integrated Frameworks (Step 2–3 Outputs)**

### **Panel-level merged data**

> (*Note: In the current script, this file is **not written to disk**; panel merging occurs entirely in memory.*)

* `step2_stid_quarter_adj_panel.csv`
  *(Disabled in script — not generated.)*

---

### **Trend & Correlation Summaries**

| File                                                             | Description                                                                      |
| ---------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `step2_tab2_admin_CB3_slope_summary.csv`                         | **CB3 trend slopes** (admin-level), used in main analyses                        |
| `step2_CB3_env_correlation_table.csv`                            | **CB3 vs environmental variables correlation matrix** (admin-level, yearly mean) |
| *(Not exported)* `step2_admin_env_slope_summary.csv`             | Admin-level environmental trend slopes (in memory only)                          |
| *(Not exported)* `step2_admin_env_trend_summary.csv`             | Categorical environmental trends                                                 |
| *(Not exported)* `step2_admin_CB_env_slope_long.csv`             | Env slope × CB slope long-format table                                           |
| *(Not exported)* `step2_admin_CB_env_timeseries_correlation.csv` | Yearly within-admin correlations                                                 |
| *(Not exported)* `step2_stid_quarter_slope_summary.csv`          | Station-level environmental slopes                                               |

---

### **Step 3 Outputs (Saved Files)**

| File                           | Description                                         |
| ------------------------------ | --------------------------------------------------- |
| `env_trend_slope_matrix.csv`   | Admin-level environmental slope matrix (raw slopes) |
| `env_CB3_combined_table.csv`   | Integrated env–CB3 trend + correlation summary      |
| `region_env_slope_summary.csv` | Regional mean environmental slopes                  |

---

# **3. Visualization Outputs (Step 3 Figures)**

All figures are generated by the scripts and stored under `figures/`.

## **3.1. Admin-level CB3 Trends**

**`fig4_Trends_of_CB3_by_municipality_2014-2024.png`**
Barplot of CB3 slopes (2014–2024).

Interpretation:

* **Negative slope → improvement** (reduction in critical barren ground)
* **Positive slope → intensification**

---

## **3.2. CB3–Environment Correlation Heatmap**

**`fig7_Heatmap_for_municipality-env._variable_in_CB3.png`**
Correlation of CB3 vs each environmental variable by admin unit.

*Zero correlations are masked in the visualization.*

---

## **3.3. Slope–Slope Scatter Matrix**

**`fig5_Relationships_between_environmental_variables_and_CB3_change.png`**
Subplots showing CB3 slope vs environmental slope (one variable per panel).

Reveals pressure–response relationships.

---

## **3.4. Environment–Ecology Decoupling Matrix**

**`fig6_Environment-ecology_decoupling_matrix_for_temperature.png`**

Quadrant interpretation:

* **Intensification:** env↑ & CB3↑
* **Decoupling (resilience):** env↑ & CB3↓
* **Joint improvement:** env↓ & CB3↓
* **Counter-intuitive:** env↓ & CB3↑

---

## **3.5. Integrated Env–CB3 Trend Combination Matrix**

**`fig8_Integrated_trend_matrix_of_environmental_change_and_CB3_intensification.png`**

Categorizes each (admin × variable) pair into:

* env↑ & CB3↑
* env↑ & CB3↓
* env↓ & CB3↓
* env↓ & CB3↑
* no trend / no data

---

## **3.6. Environmental Trend Heatmaps**

**`fig3_Normalized_slopes_of_env._variables.png`**
Sign-preserving normalized slopes (−1~1) for all environmental variables.

*(Raw and Z-score heatmaps exist in code but are disabled for now.)*

---

## **3.7. Regional Mean Environmental Trends**

**`Regional_mean_slopes_of_environmental_variables.png`**

Shows mean environmental trend per region (East, South, West, Jeju).

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

1. **Raw data preprocessing**
   `01_load_and_clean.ipynb`

2. **Depth-adjusted environmental variables**
   `02_merge_env_and_depth.ipynb`

3. **Panel construction & slope estimations**
   `03_panel_and_slope_analysis.ipynb`

4. **Visualizations & integrated frameworks**
   `04_visualization_step3.ipynb`

---

# **5. Contact & Support**

* For repository support, please open an **Issue**.
* For academic inquiries, contact the **corresponding author** in the manuscript.
* Additional documentation will be uploaded progressively under `docs/`.

