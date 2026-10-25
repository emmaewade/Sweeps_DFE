#!/usr/bin/env python
"""Collect finished gamma-DFE fits (single and recurrent sweeps, genic-only resimulation):

  202610     single sweep        (dadi_analysis[_folded])
  recurrent  recurrent sweeps    (dadi_analysis_recurrent[_folded])

Per fit: gamma shape, mean |s|, Na and the five |s| bins.

Writes working-data/202610/comparisons/dfe_by_fit.csv
Usage (dadi_env):  python compare_dfe_single_recurrent_202610.py
"""
import ast
import csv
import glob
import os
import re

import numpy as np
from scipy.stats import gamma

W = "/u/home/e/eewade/project-awbigham/data/working-data"
OUT = f"{W}/202610/comparisons"
SOURCES = {
    "202610": f"{W}/202610/dadi_analysis",
    "recurrent": f"{W}/202610/dadi_analysis_recurrent",
}
EDGES = [0, 1e-5, 1e-4, 1e-3, 1e-2, np.inf]
BIN_LABELS = ["-1e-5<s<=0", "-1e-4<s<=-1e-5", "-1e-3<s<=-1e-4",
              "-1e-2<s<=-1e-3", "s<=-1e-2"]
RUN_RE = re.compile(r"(humansim_(sweep|neutral)|recurrentsweep)_dfehuman_mu[^_]+_bens([0-9.]+)"
                    r"(?:_benfrac([0-9.]+))?_rec([0-9e.-]+)_seq\d+_(?:rep|pool\d+_reps)(\d+)(?:-\d+)?/n(\d+)_d(\d+)mb/(\w+)/dfe")


def dfe_bins(shape, mean):
    cdf = gamma.cdf(EDGES, shape, scale=mean / shape)
    return np.diff(cdf)


dfe_rows = []
for source, root in SOURCES.items():
    for folded in (False, True):
        r = root + ("_folded" if folded else "")
        for f in glob.glob(f"{r}/*/n*_d*mb/*/dfe/dfe_results.csv"):
            m = RUN_RE.search(f)
            if not m:
                continue
            kind = "recurrent" if m.group(1) == "recurrentsweep" else m.group(2)
            ben_s, benfrac, rec, rep, n, d, model = m.group(3, 4, 5, 6, 7, 8, 9)
            with open(f) as fh:
                res = [x for x in csv.DictReader(fh) if x["dfe_model"] == "gamma"][0]
            shape, scale = ast.literal_eval(res["unscaled_params"])
            key = dict(source=source, fold="folded" if folded else "unfolded", kind=kind,
                       ben_s=ben_s, benfrac=benfrac or "", rec=rec, rep=int(rep),
                       n=int(n), d=int(d), model=model)
            row = dict(key, shape=shape, mean_s=shape * scale,
                       Na=float(res["ancestral_Na_estimate"]), ns_s_ratio=float(res["ns_s_ratio"]),
                       loglik=float(res["max_log_likelihood"]))
            row.update(zip(BIN_LABELS, dfe_bins(shape, shape * scale)))
            dfe_rows.append(row)

os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/dfe_by_fit.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(dfe_rows[0]))
    w.writeheader()
    w.writerows(dfe_rows)
print(f"{len(dfe_rows)} rows -> {OUT}/dfe_by_fit.csv")
