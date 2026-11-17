"""et_common.py - shared feature code for the Extra Trees model (Python side of Student 3).

Features = the same three sets as the Perl logistic model:
  aac   amino acid composition                      20
  nct   composition of first 5 / last 5 residues    40
  prop  Student 2 properties (computed by BioPerl)   9
Length and molecular weight are left out on purpose.
"""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, matthews_corrcoef

AA = "ACDEFGHIKLMNPQRSTVWY"
PROP = ["pI", "net_charge", "gravy", "hydrophobic_pct", "cationic_pct",
        "anionic_pct", "boman", "aliphatic_index", "hydrophobic_moment"]
NAMES = ([f"aac_{a}" for a in AA] + [f"nterm_{a}" for a in AA] +
         [f"cterm_{a}" for a in AA] + PROP)


def aac(s):
    return [s.count(a) / len(s) for a in AA]


def features(seq, row):
    """seq = peptide string, row = dict-like with the Student 2 property columns."""
    return aac(seq) + aac(seq[:5]) + aac(seq[-5:]) + [float(row[p]) for p in PROP]


def load(path, labelled=True):
    """Read a properties table made by compute_props.pl. Returns (dataframe, X, y)."""
    df = pd.read_csv(path, sep="\t", dtype={"pep_id": str, "label": str, "sequence": str})
    if labelled:
        df = df[df["label"].isin(["ACP", "nonACP"])].reset_index(drop=True)
    X = np.array([features(r["sequence"], r) for _, r in df.iterrows()], dtype=float)
    y = (df["label"] == "ACP").astype(int).values if labelled else None
    return df, X, y


def metrics(y, p, thr=0.5):
    call = p >= thr
    tp = int(((y == 1) & call).sum()); fn = int(((y == 1) & ~call).sum())
    tn = int(((y == 0) & ~call).sum()); fp = int(((y == 0) & call).sum())
    return {"Sn": tp / max(tp + fn, 1), "Sp": tn / max(tn + fp, 1), "Acc": (tp + tn) / len(y),
            "MCC": matthews_corrcoef(y, call), "AUC": roc_auc_score(y, p)}
