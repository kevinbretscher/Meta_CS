# %%

import numpy as np
import joblib
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV, cross_validate
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from scipy.stats import randint, uniform
from scipy.stats import ttest_rel
import xgboost as xgb
import lightgbm as lgb

import pandas as pd
from sklearn.preprocessing import LabelEncoder

# %%
#Load data frame and transform response variable

X = pd.read_csv("data/predictors_meta_SC_rhizo_genus.tsv", sep=" ")
y = pd.read_csv("data/response_meta_SC_rhizo_genus.tsv", sep=" ", header= 0).values.ravel() 

le = LabelEncoder()
y_encoded = le.fit_transform(y)  # transform data to 1/0

# %%
# CONFIG
N_SPLITS = 3
N_ITER = 1000  # Number of iterations for RandomizedSearchCV
RANDOM_SEED = 1998
metrics = ['accuracy', 'f1', 'roc_auc']

skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_SEED)

# PARAMETER SPACES
param_spaces = {
    'RandomForest': {
        'n_estimators': randint(100, 500),
        'max_depth': randint(3, 20),
        'min_samples_split': randint(2, 20),
        'min_samples_leaf': randint(1, 10),
        'max_features': ['sqrt', 'log2', None]
    },
    'ExtraTrees': {  
        'n_estimators': randint(100, 500),
        'max_depth': randint(3, 20),
        'min_samples_split': randint(2, 20),
        'min_samples_leaf': randint(1, 10),
        'max_features': ['sqrt', 'log2', None]
    },
    'XGBoost': {
        'n_estimators': randint(30, 500),
        'max_depth': randint(3, 15),
        'learning_rate': uniform(0.01, 0.3),
        'subsample': uniform(0.5, 0.5),
        'colsample_bytree': uniform(0.5, 0.5)
    },
    'LightGBM': {
        'n_estimators': randint(30, 500),
        'max_depth': randint(3, 15),
        'learning_rate': uniform(0.01, 0.3),
        'subsample': uniform(0.5, 0.5),
        'colsample_bytree': uniform(0.5, 0.5)
    },
    'LogisticRegression': {
        'C': uniform(0.01, 10)  # inverse regularization strength
    },
    'LinearSVC': {
        'C': uniform(0.01, 10)
    }
}

# %%

# MODEL DEFINITIONS
model_defs = {
    'RandomForest': RandomForestClassifier(random_state=RANDOM_SEED),
    'ExtraTrees': ExtraTreesClassifier(random_state=RANDOM_SEED),
    'XGBoost': xgb.XGBClassifier(eval_metric='error', random_state=RANDOM_SEED),
    'LightGBM': lgb.LGBMClassifier(random_state=RANDOM_SEED, feature_pre_filter=False, min_data_in_bin=1),
    'LogisticRegression': LogisticRegression(penalty='l1', solver='liblinear', random_state=RANDOM_SEED),
    'LinearSVC': LinearSVC(random_state=RANDOM_SEED, max_iter=5000)
}

# %%

# RUN TUNING + CV
results = {}
best_params = {}

for name, model in model_defs.items():
    print(f"\n Tuning {name}...")
    
    search = RandomizedSearchCV(
        model, param_spaces[name], n_iter=N_ITER, scoring='accuracy',
        cv=skf, random_state=RANDOM_SEED, n_jobs=-1
    )
    search.fit(X, y_encoded)
    best_model = search.best_estimator_
    best_params[name] = search.best_params_
    
    print(f"Best params for {name}: {search.best_params_}")
    
    # Final CV with tuned model
    cv_results = cross_validate(best_model, X, y_encoded, cv=skf, scoring=metrics)
    results[name] = {m: cv_results[f'test_{m}'] for m in metrics}
    
    # RETRAIN BEST MODEL ON FULL DATA
    best_model.fit(X, y_encoded)
    save_path = f"best_{name.replace(' ', '_').lower()}.pkl"
    joblib.dump(best_model, save_path)
    print(f"Saved {name} model to {save_path}")

# %%
#CREATE LEADERBOARD
leaderboard_data = []
for model_name, scores in results.items():
    row = {'Model': model_name}
    for metric, values in scores.items():
        row[f"{metric}_mean"] = values.mean()
        row[f"{metric}_std"] = values.std()
    leaderboard_data.append(row)

leaderboard_df = pd.DataFrame(leaderboard_data)

# Save leaderboard
leaderboard_df.to_csv("model_leaderboard.csv", index=False)
print("\n Leaderboard saved as CSV.")

print("\n Model Performance (mean ± std):")
print(leaderboard_df)
