# %%

import joblib
import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV

import shap

# %%

# 1. Load your saved model

rf_model = joblib.load("models/best_randomforest.pkl")
xtra_model = joblib.load("models/best_extratrees.pkl")
xgb_model = joblib.load("models/best_xgboost.pkl")
lgbm_model = joblib.load("models/best_lightgbm.pkl")

# %%

#  Load data

#load original data used for training
X = pd.read_csv("data/predictors_python_saline_bulk_ASV.tsv", sep=" ")
y = pd.read_csv("data/response_python_saline_bulk_ASV.tsv", sep=" ", header= 0).values.ravel() 

le = LabelEncoder()
y_encoded = le.fit_transform(y)  # E.G. 'Acidic' -> 0, 'Alkaline' -> 1 or 'Saline' -> 1

# %%

#SHAP rank consensus
def get_shap_importance(model, X, model_name):
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)


    # shap_values is a list (XGB/LGBM multiclass) ---
    if isinstance(shap_values, list):
        # List of arrays, each (n_samples, n_features)
        shap_importance = np.mean([np.abs(sv) for sv in shap_values], axis=0) # mean over samples

    # shap_values is 3D (RandomForest multiclass) ---
    elif shap_values.ndim == 3:
        # Shape = (n_samples, n_features, n_classes)
        shap_importance = np.abs(shap_values).mean(axis=(0, 2))  # mean over samples + classes

    # shap_values is 2D (binary/regression) ---
    else:
        shap_importance = np.abs(shap_values).mean(axis=0) # mean over samples

    # Final check
    assert shap_importance.shape[0] == X.shape[1], \
        f"Shape mismatch: {shap_importance.shape[0]} SHAP features vs {X.shape[1]} input features"
    

    return pd.DataFrame({
        "feature": X.columns,
        f"shap_{model_name}": shap_importance
    }).sort_values(by=f"shap_{model_name}", ascending=False)

rf_shap  = get_shap_importance(rf_model,  X, "RF")
xtra_shap  = get_shap_importance(xtra_model,  X, "ExtraTrees")
xgb_shap = get_shap_importance(xgb_model, X, "XGB")
lgbm_shap= get_shap_importance(lgbm_model,X, "LGBM")

# %%

# Merge SHAP importances
all_shap = rf_shap.merge(xtra_shap, on="feature").merge(lgbm_shap, on="feature").merge(xgb_shap, on="feature")

# Compute SHAP ranks
for col in ["shap_RF", "shap_XGB", "shap_LGBM","shap_ExtraTrees"]:
    all_shap[col.replace("shap", "rank")] = all_shap[col].rank(ascending=False)

all_shap["mean_rank_shap"] = all_shap[["rank_RF", "rank_ExtraTrees", "rank_LGBM","rank_XGB"]].mean(axis=1)

# %%

# Permutation rank consensus

def get_perm_importance(model, X, y, model_name):
    result = permutation_importance(model, X, y, n_repeats=100, random_state=1998, n_jobs=-1, scoring='balanced_accuracy')
    importance = pd.DataFrame({
        "feature": X.columns,
        f"perm_{model_name}": result.importances_mean
    }).sort_values(by=f"perm_{model_name}", ascending=False)
    return importance

rf_perm  = get_perm_importance(rf_model,  X, y_encoded, "RF")
xtra_perm  = get_perm_importance(rf_model,  X, y_encoded, "ExtraTrees")
xgb_perm = get_perm_importance(xgb_model, X, y_encoded, "XGB")
lgbm_perm= get_perm_importance(lgbm_model,X, y_encoded, "LGBM")

# %%

# Merge permutation importances
all_perm = rf_perm.merge(xtra_perm, on="feature").merge(lgbm_perm, on="feature").merge(xgb_perm, on="feature")

# Compute permutation ranks
for col in ["perm_RF", "perm_ExtraTrees", "perm_LGBM","perm_XGB"]:
    all_perm[col.replace("perm", "rank")] = all_perm[col].rank(ascending=False)

all_perm["mean_rank_perm"] = all_perm[["rank_RF", "rank_XGB", "rank_LGBM","rank_ExtraTrees"]].mean(axis=1)

# %%
#Combine SHAP + Permutation Consensus

consensus = all_shap[["feature", "mean_rank_shap"]].merge(
    all_perm[["feature", "mean_rank_perm"]], on="feature"
)

# Add a final "consensus_score" (average of SHAP + Perm ranks)
consensus["final_consensus_rank"] = consensus[["mean_rank_shap", "mean_rank_perm"]].mean(axis=1)

# Sort by the consensus
consensus = consensus.sort_values(by="final_consensus_rank", ascending=True)

print(consensus.head(20))

# %%
import matplotlib.pyplot as plt
import seaborn as sns

# Heatmap of top 20 features by consensus
# %%

shap_ranks = all_shap[["feature","rank_RF", "rank_ExtraTrees","rank_XGB","rank_LGBM","mean_rank_shap"]].copy()
shap_ranks = shap_ranks.sort_values(by="mean_rank_shap", ascending=True)
top_features_shap = shap_ranks.head(20)
top_features_shap = top_features_shap.set_index("feature")

plt.figure(figsize=(8, 10))
ax = sns.heatmap(
    top_features_shap.iloc[:, :4],
    annot=True, fmt=".0f",
    cmap="rocket_r",  # reverse so dark = important (low rank)
    cbar_kws={'label': 'Rank (1 = most important)'}
)
plt.xticks(rotation=45, ha='right')
plt.title("Feature Importance Ranks by Model")
plt.ylabel("Feature")
plt.xlabel("Model & Method")

# Left-align y-axis labels
#ax.set_yticklabels(ax.get_yticklabels(), ha="left", rotation=0)
#ax.yaxis.set_tick_params(pad=250)

plt.tight_layout()
plt.show()