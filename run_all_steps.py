# -*- coding: utf-8 -*-
"""
Last modified on Dec 2 12:51:25 2025

@author: konjoon and chatgpt

This script runs a complete 3-step pipeline to analyze how
environmental variables are related to the expansion of
coralline barren (CB) areas.

------------------------------------------------------------
STEP 1: Depth adjustment (station-level, monthly)
------------------------------------------------------------
Input:
    - sc_filtered_stid_under30m_b.csv
        * Station-level environmental observations (<= 30 m)
        * Contains depth column and raw env variables (temp, salinity, etc.)

Output (all prefixed with "step1_"):
    - step1_stid_depth_adjusted_env.csv
        * Same table but with depth effect removed for selected variables
        * New columns: <var>_adj
        * Original env columns (temp_b, salinity_b, ...) are removed

Note:
    In this integrated script, STEP 2 directly uses env_final
    (the depth-adjusted DataFrame) passed in memory from STEP 1.
    The CSV file is written only as an optional backup/log.

------------------------------------------------------------
STEP 2: Env–CB analysis (stid/admin-level, quarterly & yearly)
------------------------------------------------------------
Main external input:
    - sc_coralline_covers_col.csv (admin-level CB area by year)

Intermediate outputs (all prefixed with "step2_"):
    - step2_stid_quarter_adj_panel.csv
    - step2_stid_quarter_slope_summary.csv
    - step2_admin_env_quarter_panel.csv

Main outputs (reported in console as STEP 2 results):
    - step2_admin_env_slope_summary.csv
    - step2_admin_env_trend_summary.csv
    - step2_admin_CB_slope_summary.csv
    - step2_admin_CB_env_slope_long.csv
    - step2_admin_CB_env_slope_correlation.csv
    - step2_admin_CB_env_timeseries_correlation.csv
    - step2_CB3_env_correlation_table.csv

(Advanced CB3 ranking/driver tables are currently disabled.)

------------------------------------------------------------
STEP 3: Report tables & plots (manuscript-ready summaries)
------------------------------------------------------------
Inputs (all passed in memory from STEP 2):
    - admin_env_trend_summary
    - admin_cb_slope
    - admin_CB_env_cor_admin
    - admin_CB_env_long
    - cb3_matrix (CB3 vs env correlation matrix)

Outputs (all prefixed with "step3_"):
    - step3_env_trend_slope_matrix.csv
    - step3_env_trend_slope_heatmap_raw.png
    - step3_env_trend_slope_heatmap_norm.png
    - step3_CB3_env_correlation_heatmap.png
    - step3_env_CB3_combined_table.csv
    - step3_CB3_admin_trend_barplot.png
    - step3_env_CB3_integrated_matrix.png
    - step3_region_env_slope_summary.csv
    - step3_region_env_slope_barplot.png
    - step3_CB3_vs_env_slope_scatter_matrix.png
    - step3_Decoupling_Matrix_Temperature.png

Note on notation:
    - In all regression summaries, the coefficient of determination
      is stored under the key/column name "r^2" instead of "r_squared".
"""

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import math

# ======================================================================
# Global switches / settings
# ======================================================================

# Output filename prefixes per step
step1_prefix = "step1_"
step2_prefix = "step2_"
step3_prefix = "step3_"

# ----------------------------------------------------------------------
# Global CSV float precision setting
# ----------------------------------------------------------------------
CSV_FLOAT_DECIMALS = 3
CSV_FLOAT_FORMAT = f"%.{CSV_FLOAT_DECIMALS}f"

# Whether to save advanced CB3 ranking & driver tables in STEP 2.
SAVE_ADVANCED_CB3_TABLES = False


###############################################################
# Step 1. Removing depth effect (station-level, monthly)
###############################################################

print("=== RUN STEP 1: Depth adjustment ===")

# ------------------------------------------------------------
# 1) Load raw station-level environmental data
# ------------------------------------------------------------
env_path = "sc_filtered_stid_under30m_b.csv"
env = pd.read_csv(env_path)

# ------------------------------------------------------------
# 2) Specify variables to be depth-adjusted
# ------------------------------------------------------------
vars_to_adjust = [
    "temp_b", "salinity_b", "pH_b", "DO_b",
    "COD_b", "NH3-N_b", "NO2-_b", "NO3-_b",
    "DIN_b", "TN_b", "DIP_b", "TP_b",
    "SAS_b", "SS_b", "Chla_b"
]
vars_to_adjust = [v for v in vars_to_adjust if v in env.columns]


# ------------------------------------------------------------
# 3) Helper: per-station linear regression of var ~ depth
# ------------------------------------------------------------
def depth_adjust_by_lm(group: pd.DataFrame, var: str) -> pd.Series:
    y = group[var].astype(float)
    x = group["depth"].astype(float)

    mask = (~y.isna()) & (~x.isna())
    if mask.sum() < 2:
        return pd.Series(np.nan, index=group.index)

    X = sm.add_constant(x[mask])
    model = sm.OLS(y[mask], X).fit()

    resid = pd.Series(np.nan, index=group.index)
    resid.loc[mask] = model.resid
    return resid


# ------------------------------------------------------------
# 4) Apply per-station regression and create *_adj variables
# ------------------------------------------------------------
env_adj = env.copy()

for var in vars_to_adjust:
    adj_series = (
        env_adj.groupby("stid", group_keys=False)[["depth", var]]
        .apply(lambda g: depth_adjust_by_lm(g, var))
    )
    env_adj[f"{var}_adj"] = adj_series

# ------------------------------------------------------------
# 5) Drop original (non-adjusted) variables and save result
# ------------------------------------------------------------
env_final = env_adj.drop(columns=vars_to_adjust, errors="ignore")

output_path = step1_prefix + "stid_depth_adjusted_env.csv"
env_final.to_csv(
    output_path,
    index=False,
    encoding="utf-8-sig",
    float_format=CSV_FLOAT_FORMAT
)

print(f"[STEP1] depth-adjusted file saved → {output_path}")


###############################################################
# Step 2. Data processing & Env–CB analysis
###############################################################

print("\n=== RUN STEP 2: Env–CB analysis ===")

# -------------------------------------------------------------------
# 0. Set working directory to script directory (optional)
# -------------------------------------------------------------------
try:
    THIS_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(THIS_DIR)
except NameError:
    pass  # e.g., Jupyter environment


# -------------------------------------------------------------------
# Utility: simple linear regression y ~ x
# -------------------------------------------------------------------
def linear_regression_summary(x, y):
    """
    Run a simple linear regression y ~ x and return key statistics.

    Returns
    -------
    dict with keys:
        intercept, slope, p_value, r^2, n
        where r^2 is the coefficient of determination (R-squared).
    """
    x = pd.Series(x)
    y = pd.Series(y)

    mask = (~x.isna()) & (~y.isna())
    x = x[mask].to_numpy(dtype=float)
    y = y[mask].to_numpy(dtype=float)

    if x.size < 3:
        return {
            "intercept": np.nan,
            "slope": np.nan,
            "p_value": np.nan,
            "r^2": np.nan,
            "n": int(x.size),
        }

    X = sm.add_constant(x)
    model = sm.OLS(y, X).fit()

    return {
        "intercept": float(model.params[0]),
        "slope": float(model.params[1]),
        "p_value": float(model.pvalues[1]),
        "r^2": float(model.rsquared),
        "n": int(x.size),
    }


# ===================================================================
# 1. env_final → stid-level quarterly panel & slopes
# ===================================================================

# Use in-memory env_final from STEP 1
env = env_final.copy()

required_env_cols = ["stid", "region", "admin", "year", "month"]
if not all(col in env.columns for col in required_env_cols):
    raise ValueError(
        "env_final must contain columns: "
        + ", ".join(required_env_cols)
    )

env["stid"] = env["stid"].astype(str)
env["region"] = env["region"].astype(str)
env["admin"] = env["admin"].astype(str)
env["year"] = pd.to_numeric(env["year"], errors="coerce")
env["month"] = pd.to_numeric(env["month"], errors="coerce")

adj_vars = [c for c in env.columns if c.endswith("_adj")]
if len(adj_vars) == 0:
    raise ValueError("No adjusted variables ('*_adj') found in depth-adjusted env data.")


def month_to_quarter(m):
    """Map month(1–12) to quarter(1–4)."""
    if 1 <= m <= 3:
        return 1
    elif 4 <= m <= 6:
        return 2
    elif 7 <= m <= 9:
        return 3
    elif 10 <= m <= 12:
        return 4
    return np.nan


env["quarter"] = env["month"].apply(month_to_quarter)

env_q_stid = (
    env.loc[~env["quarter"].isna()]
    .groupby(["stid", "region", "admin", "year", "quarter"], as_index=False)[adj_vars]
    .mean()
)

env_q_stid["time_num"] = env_q_stid["year"] + (env_q_stid["quarter"] - 1) / 4.0
env_q_stid = env_q_stid.sort_values(["stid", "year", "quarter"])

stid_quarter_panel_path = step2_prefix + "stid_quarter_adj_panel.csv"
env_q_stid.to_csv(
    stid_quarter_panel_path,
    index=False,
    float_format=CSV_FLOAT_FORMAT
)

# ---- stid-level slopes ----
env_long_stid = env_q_stid.melt(
    id_vars=["stid", "region", "admin", "year", "quarter", "time_num"],
    value_vars=adj_vars,
    var_name="variable",
    value_name="value",
).dropna(subset=["value"])


def fit_stid_slope(group: pd.DataFrame) -> pd.Series:
    res = linear_regression_summary(group["time_num"].values, group["value"].values)
    return pd.Series(
        {
            "n": res["n"],
            "start_year": group["year"].min(),
            "end_year": group["year"].max(),
            "intercept": res["intercept"],
            "slope": res["slope"],
            "p_value": res["p_value"],
            "r^2": res["r^2"],
        }
    )


stid_slope_summary = (
    env_long_stid
    .groupby(["stid", "region", "admin", "variable"])
    .apply(fit_stid_slope, include_groups=False)
    .reset_index()
)

stid_slope_path = step2_prefix + "stid_quarter_slope_summary.csv"
stid_slope_summary.to_csv(
    stid_slope_path,
    index=False,
    float_format=CSV_FLOAT_FORMAT
)


# ===================================================================
# 2. Admin-level quarterly environmental panel (weighted)
# ===================================================================

slope_stid = stid_slope_summary.copy()
slope_stid["n"] = pd.to_numeric(slope_stid["n"], errors="coerce")
slope_stid["r^2"] = pd.to_numeric(slope_stid["r^2"], errors="coerce")

# weight_raw = n * max(r^2, 0)
slope_stid["weight_raw"] = slope_stid["n"] * slope_stid["r^2"].clip(lower=0).fillna(0)


def compute_weight(row):
    """
    Station weight:
        - weight_raw if > 0
        - else n if n > 0
        - else 1
    """
    if pd.notna(row["weight_raw"]) and row["weight_raw"] > 0:
        return row["weight_raw"]
    elif pd.notna(row["n"]) and row["n"] > 0:
        return row["n"]
    else:
        return 1.0


slope_stid["weight"] = slope_stid.apply(compute_weight, axis=1)
slope_stid["variable"] = slope_stid["variable"].astype(str)
slope_stid["stid"] = slope_stid["stid"].astype(str)

env_q_long = env_q_stid.melt(
    id_vars=["stid", "region", "admin", "year", "quarter", "time_num"],
    value_vars=adj_vars,
    var_name="variable",
    value_name="value",
).dropna(subset=["value"])

env_q_long_w = env_q_long.merge(
    slope_stid[["stid", "variable", "weight"]],
    on=["stid", "variable"],
    how="left",
)

env_q_long_w["weight"] = env_q_long_w["weight"].where(
    (env_q_long_w["weight"] > 0) & (~env_q_long_w["weight"].isna()),
    1.0,
)


def weighted_admin_env(group: pd.DataFrame) -> pd.Series:
    """
    Weighted mean across stations → admin-level value per (region, admin, year, quarter, variable).
    """
    valid = group["value"].notna()
    if valid.sum() == 0:
        return pd.Series({"value_admin": np.nan, "n_stid": 0})

    v = group.loc[valid, "value"]
    w = group.loc[valid, "weight"]
    denom = w.sum()
    if denom > 0:
        value_admin = float((v * w).sum() / denom)
    else:
        value_admin = np.nan
    n_stid = int(group.loc[valid, "stid"].nunique())

    return pd.Series({"value_admin": value_admin, "n_stid": n_stid})


admin_env_q = (
    env_q_long_w
    .groupby(["region", "admin", "year", "quarter", "variable"])
    .apply(weighted_admin_env, include_groups=False)
    .reset_index()
)

admin_env_q["time_num"] = admin_env_q["year"] + (admin_env_q["quarter"] - 1) / 4.0
admin_env_q = admin_env_q.sort_values(["region", "admin", "variable", "year", "quarter"])

# admin_env_panel_path = step2_prefix + "admin_env_quarter_panel.csv"
# admin_env_q.to_csv(
#     admin_env_panel_path,
#     index=False,
#     float_format=CSV_FLOAT_FORMAT
# )


# ===================================================================
# 3. Admin-level environmental trends (slopes)
# ===================================================================

def fit_admin_env_slope(group: pd.DataFrame) -> pd.Series:
    res = linear_regression_summary(group["time_num"].values, group["value_admin"].values)
    return pd.Series(
        {
            "n_points": res["n"],
            "start_year": group["year"].min(),
            "end_year": group["year"].max(),
            "intercept": res["intercept"],
            "slope": res["slope"],
            "p_value": res["p_value"],
            "r^2": res["r^2"],
        }
    )


admin_env_slope = (
    admin_env_q
    .groupby(["region", "admin", "variable"])
    .apply(fit_admin_env_slope, include_groups=False)
    .reset_index()
)

admin_env_slope_path = step2_prefix + "admin_env_slope_summary.csv"
admin_env_slope.to_csv(
    admin_env_slope_path,
    index=False,
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP2] admin_env_slope_summary saved → {admin_env_slope_path}")

# ---- categorical trend summary ----
records = []
for _, row in admin_env_slope.iterrows():
    region = row["region"]
    admin = row["admin"]
    var_raw = str(row["variable"])
    slope = row["slope"]
    pval = row["p_value"]

    if pd.isna(slope):
        trend = "no_data"
    elif slope > 0:
        trend = "increasing"
    elif slope < 0:
        trend = "decreasing"
    else:
        trend = "no_change"

    if pd.isna(pval):
        signif = "no_data"
    elif pval < 0.05:
        signif = "significant"
    elif pval < 0.1:
        signif = "marginal"
    else:
        signif = "not_significant"

    records.append({
        "region": region,
        "admin": admin,
        "variable": var_raw,
        "raw_variable": var_raw,
        "slope": slope,
        "p_value": pval,
        "trend": trend,
        "significance": signif,
    })

admin_env_trend_summary = pd.DataFrame(records)

admin_env_trend_path = step2_prefix + "admin_env_trend_summary.csv"
admin_env_trend_summary.to_csv(
    admin_env_trend_path,
    index=False,
    encoding="utf-8-sig",
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP2] admin_env_trend_summary saved → {admin_env_trend_path}")


# ===================================================================
# 4. Admin-level CB trends (sc_coralline_covers_col.csv)
# ===================================================================

cb_panel_path = "sc_coralline_covers_col.csv"
cb_panel = pd.read_csv(cb_panel_path)

required_cb_cols = [
    "region", "admin", "year",
    "normal_km2", "prog_km2", "critical_km2"
]
if not all(col in cb_panel.columns for col in required_cb_cols):
    raise ValueError(
        "sc_coralline_covers_col.csv must contain: "
        "region, admin, year, normal_km2, prog_km2, critical_km2."
    )

cb_panel["region"] = cb_panel["region"].astype(str)
cb_panel["admin"] = cb_panel["admin"].astype(str)
cb_panel["year"] = pd.to_numeric(cb_panel["year"], errors="coerce")

cb_admin = cb_panel.assign(
    CB_1=cb_panel["normal_km2"],
    CB_2=cb_panel["prog_km2"],
    CB_3=cb_panel["critical_km2"],
)[["region", "admin", "year", "CB_1", "CB_2", "CB_3"]].drop_duplicates()

cb_cols = ["CB_1", "CB_2", "CB_3"]

cb_long = cb_admin.melt(
    id_vars=["region", "admin", "year"],
    value_vars=cb_cols,
    var_name="CB_type",
    value_name="CB_value",
).dropna(subset=["CB_value"])


def fit_admin_cb_slope(group: pd.DataFrame) -> pd.Series:
    res = linear_regression_summary(group["year"].values, group["CB_value"].values)
    return pd.Series(
        {
            "n_points": res["n"],
            "start_year": group["year"].min(),
            "end_year": group["year"].max(),
            "cb_intercept": res["intercept"],
            "cb_slope": res["slope"],
            "cb_p_value": res["p_value"],
            "cb_r^2": res["r^2"],
        }
    )


admin_cb_slope = (
    cb_long
    .groupby(["region", "admin", "CB_type"])
    .apply(fit_admin_cb_slope, include_groups=False)
    .reset_index()
)

admin_cb_slope_path = step2_prefix + "admin_CB_slope_summary.csv"
admin_cb_slope.to_csv(
    admin_cb_slope_path,
    index=False,
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP2] admin_CB_slope_summary saved → {admin_cb_slope_path}")


# ===================================================================
# 5. Cross-admin CB–environment slope correlation
# ===================================================================

admin_env_slope_renamed = admin_env_slope.rename(
    columns={
        "n_points": "env_n_points",
        "start_year": "env_start_year",
        "end_year": "env_end_year",
        "intercept": "env_intercept",
        "slope": "env_slope",
        "p_value": "env_p_value",
        "r^2": "env_r^2",
    }
)

admin_cb_slope_renamed = admin_cb_slope.rename(
    columns={
        "n_points": "cb_n_points",
        "start_year": "cb_start_year",
        "end_year": "cb_end_year",
    }
)

admin_CB_env_long = admin_env_slope_renamed.merge(
    admin_cb_slope_renamed,
    on=["region", "admin"],
    how="left",
)

admin_CB_env_long = admin_CB_env_long.dropna(subset=["env_slope", "cb_slope"])

admin_CB_env_long_path = step2_prefix + "admin_CB_env_slope_long.csv"
admin_CB_env_long.to_csv(
    admin_CB_env_long_path,
    index=False,
    float_format=CSV_FLOAT_FORMAT
)


def slope_correlation(group: pd.DataFrame) -> pd.Series:
    df = group.dropna(subset=["env_slope", "cb_slope"])
    n_admin = df.shape[0]
    if n_admin >= 3:
        x = df["cb_slope"].values
        y = df["env_slope"].values
        cor_slope = float(np.corrcoef(x, y)[0, 1])
        X = sm.add_constant(x)
        model = sm.OLS(y, X).fit()
        lm_p_value = float(model.pvalues[1])
    else:
        cor_slope = np.nan
        lm_p_value = np.nan
    return pd.Series({"n_admin": n_admin, "cor_slope": cor_slope, "lm_p_value": lm_p_value})


admin_CB_env_cor = (
    admin_CB_env_long
    .dropna(subset=["env_slope", "cb_slope"])
    .groupby(["CB_type", "variable"])
    .apply(slope_correlation, include_groups=False)
    .reset_index()
)

# admin_CB_env_cor_path = step2_prefix + "admin_CB_env_slope_correlation.csv"
# admin_CB_env_cor.to_csv(
#     admin_CB_env_cor_path,
#     index=False,
#     float_format=CSV_FLOAT_FORMAT
# )


# ===================================================================
# 6. Within-admin CB–environment time-series correlation
# ===================================================================

admin_env_year = (
    admin_env_q
    .groupby(["region", "admin", "year", "variable"], as_index=False)["value_admin"]
    .mean()
    .rename(columns={"value_admin": "env_year_mean"})
)

cb_year_long = cb_admin.melt(
    id_vars=["region", "admin", "year"],
    value_vars=cb_cols,
    var_name="CB_type",
    value_name="CB_value",
).dropna(subset=["CB_value"])

admin_env_year["year"] = pd.to_numeric(admin_env_year["year"], errors="coerce")
cb_year_long["year"] = pd.to_numeric(cb_year_long["year"], errors="coerce")

admin_env_year["year"] = admin_env_year["year"].round().astype("Int64")
cb_year_long["year"] = cb_year_long["year"].round().astype("Int64")

admin_CB_env_year = admin_env_year.merge(
    cb_year_long,
    on=["region", "admin", "year"],
    how="left",
)


def admin_ts_correlation(group: pd.DataFrame) -> pd.Series:
    df = group.dropna(subset=["env_year_mean", "CB_value"])
    n_years = df.shape[0]
    if n_years >= 3:
        x = df["env_year_mean"].values
        y = df["CB_value"].values
        cor_ts = float(np.corrcoef(x, y)[0, 1])
        X = sm.add_constant(x)
        model = sm.OLS(y, X).fit()
        lm_slope = float(model.params[1])
        lm_p_value = float(model.pvalues[1])
    else:
        cor_ts = np.nan
        lm_slope = np.nan
        lm_p_value = np.nan

    return pd.Series(
        {
            "n_years": n_years,
            "cor_ts": cor_ts,
            "lm_slope": lm_slope,
            "lm_p_value": lm_p_value,
        }
    )


admin_CB_env_cor_admin = (
    admin_CB_env_year
    .dropna(subset=["env_year_mean", "CB_value"])
    .groupby(["region", "admin", "CB_type", "variable"])
    .apply(admin_ts_correlation, include_groups=False)
    .reset_index()
)

admin_CB_env_cor_admin_path = step2_prefix + "admin_CB_env_timeseries_correlation.csv"
admin_CB_env_cor_admin.to_csv(
    admin_CB_env_cor_admin_path,
    index=False,
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP2] admin_CB_env_timeseries_correlation saved → {admin_CB_env_cor_admin_path}")


# ===================================================================
# 7. CB3 vs environmental correlation heatmap data
# ===================================================================

df_corr = admin_CB_env_cor_admin.copy()
cb3 = df_corr[df_corr["CB_type"] == "CB_3"].copy()

if cb3.empty:
    cb3_matrix = pd.DataFrame()
    print("Warning: No CB3 rows found in admin_CB_env_timeseries_correlation.")
else:
    cb3_matrix = cb3.pivot(index=["region", "admin"], columns="variable", values="cor_ts")
    cb3_matrix = cb3_matrix.fillna(0)

    cb3_table_path = step2_prefix + "CB3_env_correlation_table.csv"
    cb3_matrix.reset_index().to_csv(
        cb3_table_path,
        index=False,
        float_format=CSV_FLOAT_FORMAT
    )
    print(f"[STEP2] CB3_env_correlation_table saved → {cb3_table_path}")

    sns.set_theme(style="white", context="talk")
    plt.figure(figsize=(14, 10))

    ax = sns.heatmap(
        cb3_matrix,
        cmap="bwr",
        vmin=-1,
        vmax=1,
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "Correlation (CB_3 vs Env)"},
    )

    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=90,
        ha="center",
        fontsize=13
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=13
    )

    ax.set_xlabel("Environmental Variables", fontsize=16)
    ax.set_ylabel("Admin", fontsize=16)
    ax.set_title("Correlation between CB and Environmental Variables",
                 fontsize=20, pad=20)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("Correlation (CB vs Env)", fontsize=16)

    plt.tight_layout()
    plt.show()


# ===================================================================
# 8. Advanced CB3 analysis (disabled by default)
# ===================================================================

def _p_to_sig_label(p):
    if pd.isna(p):
        return "no_data"
    elif p < 0.05:
        return "significant"
    elif p < 0.1:
        return "marginal"
    else:
        return "not_significant"


if SAVE_ADVANCED_CB3_TABLES:
    # Advanced analysis code could be placed here if needed.
    pass
else:
    print("[STEP2] Advanced CB3 ranking/driver tables are not generated (SAVE_ADVANCED_CB3_TABLES=False).")


# ===================================================================
# STEP 3. Report tables & plots
# ===================================================================

print("\n=== RUN STEP 3: Report tables & plots ===")

core_vars = ["Temp", "DIN", "TN", "SS", "pH", "Salinity"]

# Use in-memory DataFrames from STEP 2
env_trend_core = admin_env_trend_summary.copy()
cb_slope = admin_cb_slope.copy()
cb_ts_corr = admin_CB_env_cor_admin.copy()
admin_cb_env_long = admin_CB_env_long.copy()

# -------------------------------------------------------------------
# 2. Env slope matrix (wide) + heatmaps
# -------------------------------------------------------------------

env_slope_matrix = env_trend_core.pivot_table(
    index=["region", "admin"],
    columns="raw_variable",
    values="slope"
)

env_slope_matrix = env_slope_matrix.reindex(sorted(env_slope_matrix.columns), axis=1)

env_slope_matrix_csv_path = step3_prefix + "env_trend_slope_matrix.csv"
env_slope_matrix.reset_index().to_csv(
    env_slope_matrix_csv_path,
    index=False,
    encoding="utf-8-sig",
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP3] env_trend_slope_matrix saved → {env_slope_matrix_csv_path}")

# 2-A. Raw slope heatmap
if not env_slope_matrix.empty:
    sns.set_theme(style="white", context="talk")

    n_vars = env_slope_matrix.shape[1]
    n_admins = env_slope_matrix.shape[0]
    fig_width = max(10, 0.5 * n_vars)
    fig_height = max(6, 0.4 * n_admins)

    plt.figure(figsize=(fig_width, fig_height))

    vmax = np.nanmax(np.abs(env_slope_matrix.values))
    if np.isnan(vmax) or vmax == 0:
        vmax = 1.0

    ax = sns.heatmap(
        env_slope_matrix,
        cmap="bwr",
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "Slope of environmental variable"}
    )

    ax.set_xlabel("Environmental variables", fontsize=14)
    ax.set_ylabel("Region / Admin", fontsize=14)
    ax.set_title("Admin-level environmental trends (raw slopes)",
                 fontsize=18, pad=20)

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()

    env_slope_heatmap_raw_path = step3_prefix + "env_trend_slope_heatmap_raw.png"
    plt.savefig(env_slope_heatmap_raw_path, dpi=300)
    plt.close()
    print(f"[STEP3] env_trend_slope_heatmap_raw saved → {env_slope_heatmap_raw_path}")
else:
    print("Warning: No data for env_slope_matrix heatmap (raw).")

# 2-B. Normalized slope heatmap
if not env_slope_matrix.empty:
    env_slope_norm = env_slope_matrix.copy()

    for col in env_slope_norm.columns:
        col_min = env_slope_norm[col].min()
        col_max = env_slope_norm[col].max()
        if pd.notna(col_min) and pd.notna(col_max) and col_max > col_min:
            env_slope_norm[col] = (env_slope_norm[col] - col_min) / (col_max - col_min)
        else:
            env_slope_norm[col] = 0.5

    sns.set_theme(style="white", context="talk")

    n_vars = env_slope_norm.shape[1]
    n_admins = env_slope_norm.shape[0]
    fig_width = max(10, 0.5 * n_vars)
    fig_height = max(6, 0.4 * n_admins)

    plt.figure(figsize=(fig_width, fig_height))

    ax = sns.heatmap(
        env_slope_norm,
        cmap="bwr",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "Normalized slope (0-1, per variable)"}
    )

    ax.set_xlabel("Environmental variables", fontsize=14)
    ax.set_ylabel("Region / Admin", fontsize=14)
    ax.set_title("Admin-level environmental trends (normalized slopes)",
                 fontsize=18, pad=20)

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()

    env_slope_heatmap_norm_path = step3_prefix + "env_trend_slope_heatmap_norm.png"
    plt.savefig(env_slope_heatmap_norm_path, dpi=300)
    plt.close()
    print(f"[STEP3] env_trend_slope_heatmap_norm saved → {env_slope_heatmap_norm_path}")
else:
    print("Warning: No data for env_slope_matrix heatmap (normalized).")


# -------------------------------------------------------------------
# 3. CB3 vs environment correlation heatmap (report version)
# -------------------------------------------------------------------

if "cb3_matrix" not in locals():
    print("Warning: cb3_matrix not found in memory; skipping CB3 heatmap.")
else:
    cb3_plot = cb3_matrix.copy()

    if set(cb3_plot.index.names) == {"region", "admin"}:
        cb3_plot = cb3_plot.reset_index()
        cb3_plot["region_admin"] = (
            cb3_plot["region"].astype(str) + " / " + cb3_plot["admin"].astype(str)
        )
        cb3_plot = cb3_plot.drop(columns=["region", "admin"]).set_index("region_admin")
    elif {"region", "admin"}.issubset(cb3_plot.columns):
        cb3_plot["region_admin"] = (
            cb3_plot["region"].astype(str) + " / " + cb3_plot["admin"].astype(str)
        )
        cb3_plot = cb3_plot.drop(columns=["region", "admin"]).set_index("region_admin")

    if not cb3_plot.empty:
        sns.set_theme(style="white", context="talk")

        plt.figure(figsize=(10, max(6, 0.4 * cb3_plot.shape[0])))

        ax = sns.heatmap(
            cb3_plot,
            cmap="bwr",
            vmin=-1,
            vmax=1,
            linewidths=0.5,
            linecolor="gray",
            cbar_kws={"label": "Correlation (CB3 vs Env)"}
        )

        ax.set_xlabel("Environmental variables", fontsize=14)
        ax.set_ylabel("Region / Admin", fontsize=14)
        ax.set_title("Within-admin correlation: CB3 vs environmental variables",
                     fontsize=18, pad=20)

        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        plt.tight_layout()

        cb3_heatmap_path = step3_prefix + "CB3_env_correlation_heatmap.png"
        plt.savefig(cb3_heatmap_path, dpi=300)
        plt.close()
        print(f"[STEP3] CB3_env_correlation_heatmap saved → {cb3_heatmap_path}")
    else:
        print("Warning: CB_3 correlation matrix is empty; no report heatmap created.")


# -------------------------------------------------------------------
# 4. Combined Env + CB3 summary table
# -------------------------------------------------------------------
cb3_slope_only = cb_slope[cb_slope["CB_type"] == "CB_3"].copy()
cb3_slope_only = cb3_slope_only[[
    "region", "admin", "CB_type",
    "n_points", "start_year", "end_year",
    "cb_slope", "cb_p_value", "cb_r^2"
]]

cb3_ts = cb_ts_corr[cb_ts_corr["CB_type"] == "CB_3"].copy()

combined = env_trend_core.merge(
    cb3_ts,
    left_on=["region", "admin", "raw_variable"],
    right_on=["region", "admin", "variable"],
    how="left",
    suffixes=("_env", "_ts")
)

combined = combined.merge(
    cb3_slope_only,
    on=["region", "admin"],
    how="left"
)


def encode_cb_trend(row):
    slope = row["cb_slope"]
    pval = row["cb_p_value"]
    if pd.isna(slope) or pd.isna(pval):
        return "no_data"
    if slope > 0:
        direction = "increasing"
    elif slope < 0:
        direction = "decreasing"
    else:
        direction = "no_change"
    if pval < 0.05:
        signif = "significant"
    elif pval < 0.1:
        signif = "marginal"
    else:
        signif = "not_significant"
    return f"{direction} ({signif})"


def encode_corr_pattern(row):
    cor = row["cor_ts"]
    p = row["lm_p_value"]
    if pd.isna(cor):
        return "no_data"
    if p is not None and p < 0.05:
        signif = "significant"
    elif p is not None and p < 0.1:
        signif = "marginal"
    else:
        signif = "not_significant"
    if cor > 0:
        dir_str = "positive"
    elif cor < 0:
        dir_str = "negative"
    else:
        dir_str = "zero"
    return f"{dir_str} ({signif})"


combined["CB3_trend_summary"] = combined.apply(encode_cb_trend, axis=1)
combined["CB3_env_corr_summary"] = combined.apply(encode_corr_pattern, axis=1)

# Remove rows where both CB3 trend and correlation are "no_data"
mask_valid = ~(
    (combined["CB3_trend_summary"] == "no_data") &
    (combined["CB3_env_corr_summary"] == "no_data")
)
combined = combined[mask_valid].copy()

# 최종 출력용 컬럼 선택 (env_variable 하나만, n_points 제외)
combined_out = combined[[
    "region",
    "admin",
    "raw_variable",   # → env_variable 로 rename
    "start_year",     # env_variable 바로 뒤에 위치
    "end_year",
    "slope",
    "p_value",
    "trend",
    "significance",
    "n_years",
    "cor_ts",
    "lm_p_value",
    "CB3_env_corr_summary",
    "cb_slope",
    "cb_p_value",
    "cb_r^2",
    "CB3_trend_summary"
]].copy()

# 컬럼 이름 표준화
combined_out = combined_out.rename(columns={
    "raw_variable": "env_variable",   # env_core_variable + env_raw_variable → env_variable 하나
    "slope": "env_slope",
    "p_value": "env_p_value",
    "n_years": "n_year"               # n_years → n_year
})

combined_out_path = step3_prefix + "env_CB3_combined_table.csv"
combined_out.to_csv(
    combined_out_path,
    index=False,
    encoding="utf-8-sig",
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP3] env_CB3_combined_table saved → {combined_out_path}")


# -------------------------------------------------------------------
# Figure: CB3 admin trend barplot
# -------------------------------------------------------------------

cb3_plot_df = cb_slope[cb_slope["CB_type"] == "CB_3"].copy()
cb3_plot_df = cb3_plot_df.dropna(subset=["region", "admin", "cb_slope"])

if cb3_plot_df.empty:
    print("Warning: No valid CB_3 rows (region/admin/cb_slope). Figure 7 will not be generated.")
else:
    cb3_plot_df["region_admin"] = cb3_plot_df["region"].astype(str) + " / " + cb3_plot_df["admin"].astype(str)
    cb3_plot_df = cb3_plot_df.sort_values("cb_slope", ascending=True)

    admin_order = pd.unique(cb3_plot_df["region_admin"])
    cb3_plot_df["region_admin"] = pd.Categorical(cb3_plot_df["region_admin"], categories=admin_order, ordered=True)

    sns.set_theme(style="whitegrid", context="talk")

    plt.figure(figsize=(12, 6))
    ax = sns.barplot(
        data=cb3_plot_df,
        x="region_admin",
        y="cb_slope",
        hue="region",
        dodge=False
    )

    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Region / Admin", fontsize=12)
    ax.set_ylabel("CB3 area trend ($km^2$/yr)", fontsize=12)
    ax.set_title("Admin-level trends in CB3 (critical-stage) barren area, 2014-2024",
                 fontsize=16, pad=16)

    plt.xticks(rotation=60, ha="right")
    plt.legend(title="Region", loc="upper left", bbox_to_anchor=(1.02, 1.0))
    plt.tight_layout()

    fig7_path = step3_prefix + "CB3_admin_trend_barplot.png"
    plt.savefig(fig7_path, dpi=300)
    plt.close()
    print(f"[STEP3] CB3_admin_trend_barplot saved → {fig7_path}")


# -------------------------------------------------------------------
# Figure: Integrated Env–CB3 trend matrix
# -------------------------------------------------------------------

if 'combined_out' not in locals():
    try:
        combined_out = pd.read_csv(step3_prefix + "env_CB3_combined_table.csv")
        print("Loaded env_CB3_combined_table.csv for integrated matrix.")
    except FileNotFoundError:
        combined_out = None
        print("Warning: env_CB3_combined_table.csv not found; Figure 8 will not be generated.")

if combined_out is None or combined_out.empty:
    print("Warning: No data in combined_out; skipping integrated environmental-CB3 matrix (Figure 8).")
else:
    df_int = combined_out.copy()

    def slope_dir(val, eps=1e-6):
        if pd.isna(val):
            return 0
        if abs(val) < eps:
            return 0
        return 1 if val > 0 else -1

    categories = []
    for _, row in df_int.iterrows():
        ed = slope_dir(row["env_slope"])
        cd = slope_dir(row["cb_slope"])
        if pd.isna(row["env_slope"]) or pd.isna(row["cb_slope"]):
            categories.append("no_data")
        elif ed == 0 or cd == 0:
            categories.append("no_trend")
        elif ed > 0 and cd > 0:
            categories.append("env↑, CB3↑ (intensification)")
        elif ed > 0 and cd < 0:
            categories.append("env↑, CB3↓ (decoupling)")
        elif ed < 0 and cd > 0:
            categories.append("env↓, CB3↑ (counter-intuitive)")
        elif ed < 0 and cd < 0:
            categories.append("env↓, CB3↓ (joint improvement)")
        else:
            categories.append("no_trend")

    df_int["category"] = categories
    df_int["region_admin"] = df_int["region"].astype(str) + " / " + df_int["admin"].astype(str)

    cat_order = [
        "env↑, CB3↑ (intensification)",
        "env↑, CB3↓ (decoupling)",
        "env↓, CB3↑ (counter-intuitive)",
        "env↓, CB3↓ (joint improvement)",
        "no_trend",
        "no_data",
    ]
    cat_to_code = {c: i for i, c in enumerate(cat_order)}

    df_int["cat_code"] = df_int["category"].map(cat_to_code)

    # 여기서 env_core_variable → env_variable 사용
    plot_df = df_int.dropna(subset=["env_variable", "cat_code"]).copy()

    if plot_df.empty:
        print("Warning: No valid rows for integrated matrix; skipping Figure 8.")
    else:
        matrix = plot_df.pivot_table(
            index="region_admin",
            columns="env_variable",
            values="cat_code",
            aggfunc="first"
        )

        matrix = matrix.dropna(how="all")

        no_data_code = cat_to_code["no_data"]

        def all_nodata(row):
            non_na = row.dropna()
            if non_na.empty:
                return True
            return non_na.isin([no_data_code]).all()

        mask_all_nodata = matrix.apply(all_nodata, axis=1)
        matrix = matrix[~mask_all_nodata]

        colors = [
            "#d73027",
            "#fc8d59",
            "#91bfdb",
            "#1a9850",
            "#dddddd",
            "#aaaaaa",
        ]

        cmap = ListedColormap(colors)

        sns.set_theme(style="white", context="talk")

        fig_width = 15
        fig_height = max(8, 0.4 * matrix.shape[0])
        plt.figure(figsize=(fig_width, fig_height))

        ax = sns.heatmap(
            matrix,
            cmap=cmap,
            vmin=-0.5,
            vmax=len(cat_order) - 0.5,
            cbar=False,
            linewidths=0.5,
            linecolor="gray"
        )

        ax.set_xlabel("Environmental variables", fontsize=16)
        ax.set_ylabel("Region / Admin", fontsize=16)
        ax.set_title("Integrated environmental–CB3 trend matrix", fontsize=22, pad=20)

        plt.xticks(rotation=45, ha="right", fontsize=12)
        plt.yticks(rotation=0, fontsize=12)

        legend_patches = [
            mpatches.Patch(color=colors[i], label=cat_order[i])
            for i in range(len(cat_order))
        ]

        ax.legend(
            handles=legend_patches,
            title="Trend combination",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
            fontsize=12,
            title_fontsize=14
        )

        plt.tight_layout(rect=[0, 0, 0.95, 1])

        fig8_path = step3_prefix + "env_CB3_integrated_matrix.png"
        plt.savefig(fig8_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[STEP3] env_CB3_integrated_matrix saved → {fig8_path}")


# -------------------------------------------------------------------
# 5. Region-level mean slope summary + barplot
# -------------------------------------------------------------------

region_group = (
    env_trend_core
    .groupby(["region", "raw_variable"], as_index=False)
    .agg(
        mean_slope=("slope", "mean"),
        std_slope=("slope", "std"),
        n_admin=("slope", "count"),
    )
)

region_group_csv_path = step3_prefix + "region_env_slope_summary.csv"
region_group.to_csv(
    region_group_csv_path,
    index=False,
    encoding="utf-8-sig",
    float_format=CSV_FLOAT_FORMAT
)
print(f"[STEP3] region_env_slope_summary saved → {region_group_csv_path}")

# Barplot: mean slope per variable, grouped by region
if not region_group.empty:
    sns.set_theme(style="whitegrid", context="talk")

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=region_group,
        x="raw_variable",
        y="mean_slope",
        hue="region",
        errorbar=None,
    )

    ax.set_xlabel("Environmental variables", fontsize=14)
    ax.set_ylabel("Mean slope (by region)", fontsize=14)
    ax.set_title(
        "Regional mean trends of environmental variables",
        fontsize=18,
        pad=20,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    region_bar_path = step3_prefix + "region_env_slope_barplot.png"
    plt.savefig(region_bar_path, dpi=300)
    plt.close()
    print(f"[STEP3] region_env_slope_barplot saved → {region_bar_path}")
else:
    print("Warning: No data for region_group barplot.")


# -------------------------------------------------------------------
# 6. Multi-panel scatter: env_slope vs CB3 slope
# -------------------------------------------------------------------

cb3_slope_long = admin_cb_env_long[admin_cb_env_long["CB_type"] == "CB_3"].copy()

if not cb3_slope_long.empty:
    cb3_slope_long["region_admin"] = (
            cb3_slope_long["region"].astype(str)
            + " / "
            + cb3_slope_long["admin"].astype(str)
    )
        
    all_vars = sorted(cb3_slope_long["variable"].dropna().unique())
    panel_tags = [f"({chr(97 + i)})" for i in range(len(all_vars))]

    nvars = len(all_vars)
    ncols = 3
    nrows = math.ceil(nvars / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5 * ncols, 3.2 * nrows),
        sharey=True
    )
    axes = np.array(axes).reshape(-1)

    for i, var_name in enumerate(all_vars):
        ax = axes[i]
        sub = cb3_slope_long[
            (cb3_slope_long["variable"] == var_name) &
            (~cb3_slope_long["env_slope"].isna()) &
            (~cb3_slope_long["cb_slope"].isna())
        ]

        if sub.shape[0] < 3:
            ax.set_visible(False)
            continue

        sns.regplot(
            data=sub,
            x="env_slope",
            y="cb_slope",
            scatter_kws={"s": 40, "alpha": 0.7},
            line_kws={"linewidth": 2},
            ci=None,
            ax=ax
        )

        for _, row in sub.iterrows():
            ax.text(
                row["env_slope"],
                row["cb_slope"],
                row["region_admin"],   # 예: "Jeju / Seogwipo"
                fontsize=6,
                alpha=0.6
            )
            
        ax.set_xlabel(f"{panel_tags[i]} {var_name} slope")
        ax.set_ylabel("CB3 slope")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Slope–Slope Relationships: CB3 vs Environmental Variables",
        fontsize=16,
        y=0.99
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    scatter_path = step3_prefix + "CB3_vs_env_slope_scatter_plots.png"
    plt.savefig(scatter_path, dpi=300)
    plt.close()
    print(f"[STEP3] CB3_vs_env_slope_scatter_plots saved → {scatter_path}")


# -------------------------------------------------------------------
# 7. 2x2 Decoupling Matrix Plot (Temperature vs CB3)
# -------------------------------------------------------------------

target_var = "temp_b_adj"

df_matrix = admin_cb_env_long[
    (admin_cb_env_long["CB_type"] == "CB_3") &
    (admin_cb_env_long["variable"] == target_var)
].copy()

if df_matrix.empty:
    print(f"Warning: No data found for {target_var} in 2x2 Matrix Plot.")
else:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.scatterplot(
        data=df_matrix,
        x="env_slope",
        y="cb_slope",
        hue="region",
        style="region",
        s=100,
        alpha=0.8,
        ax=ax
    )

    ax.legend(
        title="Region",
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=12,
        title_fontsize=10,
        frameon=True
    )

    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.axvline(0, color="gray", linestyle="--", linewidth=1)

    for i, row in df_matrix.iterrows():
        ax.text(
            row["env_slope"] + 0.002,
            row["cb_slope"] + 0.005,
            row["admin"],
            fontsize=12,
            color='black',
            alpha=0.7
        )

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()

    ax.text(x_max * 0.9, y_max * 0.9,
            "Intensification\n(High Risk)",
            ha='right', va='top',
            fontsize=13, color='red', weight='bold')

    ax.text(x_max * 0.9, y_min * 0.9,
            "Decoupling Zone\n(Policy Resilience)",
            ha='right', va='bottom',
            fontsize=13, color='green', weight='bold')

    ax.text(x_min * 0.9, y_min * 0.9,
            "Joint Improvement",
            ha='left', va='bottom',
            fontsize=13, color='blue', alpha=0.6)

    ax.text(x_min * 0.9, y_max * 0.9,
            "Counter-intuitive",
            ha='left', va='top',
            fontsize=13, color='gray', alpha=0.6)

    ax.set_xlabel(f"Environmental Trend Slope ({target_var})", fontsize=12)
    ax.set_ylabel("CB3 Trend Slope (km²/year)", fontsize=12)
    ax.set_title(f"Environment-Ecology Decoupling Matrix\n(Variable: {target_var})",
                 fontsize=10, pad=20)

    ax.fill_between([0, x_max], 0, y_min, color='green', alpha=0.05)

    plt.tight_layout()

    matrix_plot_path = step3_prefix + "Decoupling_Matrix_Temperature.png"
    plt.savefig(matrix_plot_path, dpi=300)
    plt.close()
    print(f"[STEP3] Decoupling_Matrix_Temperature saved → {matrix_plot_path}")

print("[STEP3] All steps completed!")
