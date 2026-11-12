import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.model_selection import cross_val_score, train_test_split

# Data loading
print("Loading reference dataset...")
df_ref = pd.read_csv(
    "data/brca1_clinvar_2026-10-05_reference_labels.tsv", sep="\t", comment="#"
)

# 1 = Pathogenic/Likely Pathogenic, 0 = Benign/Likely Benign
df_ref["target"] = df_ref["germline_classification"].apply(
    lambda x: (
        1
        if any(
            term in str(x).lower()
            for term in ["pathogenic", "likely pathogenic"]
        )
        else 0
    )
)

df_ref["aa_ref_code"] = df_ref["aa_ref"].astype("category").cat.codes
df_ref["aa_alt_code"] = df_ref["aa_alt"].astype("category").cat.codes

features = ["aa_pos", "aa_ref_code", "aa_alt_code"]
X = df_ref[features].fillna(0)
y = df_ref["target"]

# Splitting dataset into training and testing
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Model training
print("Training random forest model...")
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Model evaluation
y_pred_proba = model.predict_proba(X_test)[:, 1]
roc_auc = roc_auc_score(y_test, y_pred_proba)
print(f"\nROC-AUC score on test set: {roc_auc:.4f}")

# Application to the VUS dataset
print("\nLoading prediction dataset...")
df_pred = pd.read_csv(
    "data/brca1_clinvar_2026-10-05_predict.tsv", sep="\t", comment="#"
)
df_pred["aa_ref_code"] = (
    df_pred["aa_ref"].astype("category").cat.codes
)
df_pred["aa_alt_code"] = df_pred["aa_alt"].astype("category").cat.codes

X_pred = df_pred[features].fillna(0)

df_pred["predicted_pathogenic_probability"] = model.predict_proba(X_pred)[:, 1]

print("\nOverview of the 10 first predictions of VUS :")
print(
    df_pred[["variation_id", "hgvs_p", "predicted_pathogenic_probability"]].head(
        10
    )
)

#output file
output_filename = "brca1_vus_predictino_results_random_forest.csv"

output_df  =df_pred[
    [
        "variation_id", "hgvs_c", "hgvs_p",
                     "aa_ref", "aa_pos", "aa_alt",
                     "predicted_pathogenic_probability"
    ]
].copy()

output_df.to_csv(output_filename, index=False)
