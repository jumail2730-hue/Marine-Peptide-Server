#!/usr/bin/env python3
"""train_extratrees.py - train and validate the Extra Trees anticancer-peptide model.

Usage: python3 train_extratrees.py training_balanced_props.tsv [--outdir data/final/model]

1. Holds out 20% of the peptides (stratified) as a test set.
2. 10-fold cross-validation on the other 80%, then one score on the held-out 20%.
3. Refits on ALL peptides and saves the model to <outdir>/model_et.joblib.
"""
import argparse
import os
import joblib
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict, train_test_split
import et_common as ec


def make(seed):
    return ExtraTreesClassifier(n_estimators=500, min_samples_leaf=1, max_features="sqrt",
                                n_jobs=-1, random_state=seed)


def show(name, m):
    print(f"{name:<22} Sn {m['Sn']:.3f}  Sp {m['Sp']:.3f}  Acc {m['Acc']:.3f}  "
          f"MCC {m['MCC']:.3f}  AUC {m['AUC']:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("props")
    ap.add_argument("--outdir", default="model")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    df, X, y = ec.load(a.props)
    print(f"Loaded {len(y)} peptides ({y.sum()} ACP, {(y == 0).sum()} non-ACP)")
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, stratify=y, random_state=a.seed)
    print(f"Train/CV: {len(ytr)} peptides; held-out test: {len(yte)} peptides\n")

    cv = StratifiedKFold(10, shuffle=True, random_state=a.seed)
    p_cv = cross_val_predict(make(a.seed), Xtr, ytr, cv=cv, method="predict_proba")[:, 1]
    show("10-fold cross-val", ec.metrics(ytr, p_cv))
    p_te = make(a.seed).fit(Xtr, ytr).predict_proba(Xte)[:, 1]
    show("Held-out test", ec.metrics(yte, p_te))

    final = make(a.seed).fit(X, y)
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "model_et.joblib")
    joblib.dump(final, out)
    print(f"\nSaved model (refitted on all {len(y)} peptides): {out}")

    top = sorted(zip(final.feature_importances_, ec.NAMES), reverse=True)[:12]
    print("\nMost important features:")
    for imp, name in top:
        print(f"  {name:<20} {imp:.3f}")


if __name__ == "__main__":
    main()
