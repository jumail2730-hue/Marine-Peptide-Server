#!/usr/bin/env python3
"""predict_et.py - score peptides with the saved Extra Trees model.

Usage: python3 predict_et.py model_et.joblib props.tsv > scores.tsv

props.tsv comes from Student 2's compute_props.pl (BioPerl), for example:
  perl student2/compute_props.pl data/final/discovery_marine.fasta > data/final/discovery_marine_props.tsv

Output (best first): rank, pep_id, label, sequence, anticancer_probability, call
"""
import sys
import joblib
import et_common as ec


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: predict_et.py model_et.joblib props.tsv > scores.tsv")
    model = joblib.load(sys.argv[1])
    df, X, _ = ec.load(sys.argv[2], labelled=False)
    df["p"] = model.predict_proba(X)[:, 1]
    df = df.sort_values("p", ascending=False).reset_index(drop=True)
    print("rank\tpep_id\tlabel\tsequence\tanticancer_probability\tcall")
    for i, r in df.iterrows():
        label = r["label"] if isinstance(r["label"], str) else ""
        call = "ACP" if r["p"] >= 0.5 else "non-ACP"
        print(f"{i + 1}\t{r['pep_id']}\t{label}\t{r['sequence']}\t{r['p']:.4f}\t{call}")
    print(f"Scored {len(df)} peptides", file=sys.stderr)


if __name__ == "__main__":
    main()
