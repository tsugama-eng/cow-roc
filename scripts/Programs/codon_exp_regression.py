#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import argparse
from sklearn.linear_model import RidgeCV, LassoCV, LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.decomposition import PCA
from scipy import odr
import scipy.linalg as la

# --- CLR変換 ---
def clr_transform(matrix: pd.DataFrame) -> pd.DataFrame:
    gm = np.exp(np.mean(np.log(matrix), axis=1))
    clr = np.log(matrix.div(gm, axis=0))
    return clr

# --- ILR変換 ---
def ilr_transform(matrix: pd.DataFrame) -> pd.DataFrame:
    log_matrix = np.log(matrix)
    gm = np.exp(np.mean(log_matrix, axis=1))
    clr = log_matrix - np.log(gm.to_numpy())[:, None]   # ← 修正
    # QR分解で直交基底を構築
    Q, _ = la.qr(np.eye(clr.shape[1]) - np.ones((clr.shape[1], clr.shape[1]))/clr.shape[1])
    ilr = clr @ Q[:, :clr.shape[1]-1]
    return pd.DataFrame(ilr, index=matrix.index)

def run_model(codon_file, tpm_file, tpm_col=None, model_type="ridge", auto_alpha=True,
              log_transform="log1p", feature_transform="clr",
              model_out="model_coeffs.tsv", pred_out="predictions.tsv",
              tpm_threshold=0.0):

    codon_df = pd.read_csv(codon_file, sep="\t", index_col=0)
    tpm_df = pd.read_csv(tpm_file, sep="\t", index_col=0)

    df = codon_df.join(tpm_df, how="inner")

    # --- Codon features ---
    codon_features_raw = df.iloc[:, 0:64].copy()
    codon_features_raw = codon_features_raw.replace(0, 1e-6)
    codon_features_raw = codon_features_raw.reindex(sorted(codon_features_raw.columns), axis=1)

    if feature_transform == "clr":
        codon_features = clr_transform(codon_features_raw)
    elif feature_transform == "ilr":
        codon_features = ilr_transform(codon_features_raw)
    else:
        codon_features = codon_features_raw

    extra_features = df.iloc[:, 64:codon_df.shape[1]].copy()

    parts = [
        pd.DataFrame(codon_features, index=df.index, columns=codon_features.columns),
        extra_features
    ]
    # 空でないものだけ残す
    parts = [p for p in parts if p.shape[1] > 0]
    features = pd.concat(parts, axis=1)

    # --- TPM処理 ---
    if tpm_col is not None:
        y_raw = df.iloc[:, codon_df.shape[1] + tpm_col]
    else:
        y_raw = df.iloc[:, codon_df.shape[1]:].mean(axis=1)

    mask_threshold = y_raw > tpm_threshold
    y_raw = y_raw.loc[mask_threshold]
    features = features.loc[mask_threshold, :]

    if log_transform == "log1p":
        y = np.log1p(y_raw.values)
        X = features.values
    elif log_transform == "log10":
        y_raw_safe = y_raw.replace(0, 1e-6)
        y = np.log10(y_raw_safe)
        X = features.values
        y_raw = y_raw_safe
    else:
        y = y_raw.values
        X = features.values

    print("Number of samples after preprocessing:", len(y))
    if len(y) == 0:
        raise ValueError("No samples left after preprocessing. Check input files or threshold.")

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # --- モデル選択 ---
    if model_type == "deming":
        def f(B, x):
            return B[0]*x + B[1]
        linear = odr.Model(f)
        data = odr.Data(X_train[:,0], y_train)
        odr_inst = odr.ODR(data, linear, beta0=[1., 0.])
        out = odr_inst.run()
        model = out
        y_pred = f(out.beta, X_test[:,0])

    elif model_type == "ridge":
        model = RidgeCV(alphas=np.logspace(-3, 3, 20)) if auto_alpha else RidgeCV(alphas=[1.0])
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

    elif model_type == "lasso":
        model = LassoCV(alphas=np.logspace(-3, 3, 20), max_iter=5000) if auto_alpha else LassoCV(alphas=[1.0], max_iter=5000)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

    elif model_type == "linear":
        model = LinearRegression()
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

    elif model_type == "weighted":
        weights = np.sqrt(np.abs(y_train))
        model = LinearRegression()
        model.fit(X_train, y_train, sample_weight=weights)
        y_pred = model.predict(X_test)

    elif model_type == "pcr":
        pca = PCA(n_components=min(X_train.shape[1], 10))
        X_train_pca = pca.fit_transform(X_train)
        X_test_pca = pca.transform(X_test)
        model = LinearRegression()
        model.fit(X_train_pca, y_train)
        y_pred = model.predict(X_test_pca)

    else:
        raise ValueError("Unknown model_type. Choose from 'ridge', 'lasso', 'linear', 'weighted', 'pcr', 'deming'.")

    print(f"Model: {model_type}")
    print("R^2:", r2_score(y_test, y_pred))
    if hasattr(model, "alpha_"):
        print("Best alpha:", model.alpha_)

    # --- 係数出力 ---
    if model_type in ["ridge", "lasso", "linear", "weighted", "pcr"]:
        coef_df = pd.DataFrame({
            "variable": features.columns[:len(model.coef_)],
            "coef": model.coef_
        })
    elif model_type == "deming":
        coef_df = pd.DataFrame({
            "variable": ["x"],
            "coef": [model.beta[0]],
            "intercept": [model.beta[1]]
        })

    coef_df.to_csv(model_out, sep="\t", index=False)

    # --- 予測出力 ---
    if model_type == "deming":
        all_pred = f(model.beta, X[:,0])
    elif model_type == "pcr":
        all_pred = model.predict(pca.transform(X))
    else:
        all_pred = model.predict(X)

    if log_transform == "log1p":
        observed_transformed = np.log1p(y_raw.values)
    elif log_transform == "log10":
        observed_transformed = np.log10(y_raw.values)
    else:
        observed_transformed = None

    pred_df = pd.DataFrame({
        "gene_id": features.index,
        "observed_TPM": y_raw.values,
        "observed_transformed": observed_transformed,
        "predicted": all_pred
    })
    pred_df.to_csv(pred_out, sep="\t", index=False)

    return model

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Codon usage vs TPM regression")
    parser.add_argument("codon_file", help="Codon feature file (TSV)")
    parser.add_argument("tpm_file", help="TPM file (TSV)")
    parser.add_argument("--tpm_col", type=int, default=None, help="Column index of TPM to use (default: mean of all)")
    parser.add_argument("--method", choices=["deming", "ridge", "lasso", "linear", "weighted", "pcr"], default="linear", help="Regression method")
    parser.add_argument("--auto_alpha", action="store_true", help="Enable automatic alpha search (default: True)")
    parser.add_argument("--no_auto_alpha", action="store_true", help="Disable automatic alpha search")
    parser.add_argument("--log_transform", choices=["log1p", "log10", "none"], default="log1p", help="Transformation for TPM (default: log1p)")
    parser.add_argument("--feature_transform", choices=["clr", "ilr", "none"], default="clr", help="Codon feature transform (default: clr)")
    parser.add_argument("--model_out", default=None, help="Output file for coefficients")
    parser.add_argument("--pred_out", default=None, help="Output file for predictions")
    parser.add_argument("--threshold", type=float, default=0.0, help="TPM threshold (default: 0.0)")
    parser.add_argument("--model_pkl", default="trained_model.pkl",
                    help="Output pickle file for trained model (default: trained_model.pkl)")
                    
    args = parser.parse_args()

    auto_alpha = True
    if args.no_auto_alpha:
        auto_alpha = False

    if args.model_out is None:
        base = os.path.splitext(os.path.basename(args.tpm_file))[0]
        col_str = str(args.tpm_col) if args.tpm_col is not None else "mean"
        args.model_out = f"{base}.{col_str}.{args.method}.coeffs.tsv"

    if args.pred_out is None:
        base = os.path.splitext(os.path.basename(args.tpm_file))[0]
        col_str = str(args.tpm_col) if args.tpm_col is not None else "mean"
        args.pred_out = f"{base}.{col_str}.{args.method}.pred.tsv"

    model = run_model(args.codon_file, args.tpm_file, args.tpm_col, args.method,
                  auto_alpha, args.log_transform, args.feature_transform,
                  args.model_out, args.pred_out, args.threshold)

    # --- Save trained model ---
    import pickle
    with open(args.model_pkl, "wb") as f:
        pickle.dump(model, f)

    print(f"Saved trained model to: {args.model_pkl}")