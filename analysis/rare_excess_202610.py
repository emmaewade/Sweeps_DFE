#!/usr/bin/env python3
"""Relative excess of rare variants, 1 - pi/theta_W, for syn and nonsyn SFS.

Reads the empirical folded SFS used by the dadi fits (n = 100):
  single sweep: dadi_analysis_folded/<run>/n100_d{1,3,5}mb/two_epoch/dfe/empirical_{syn,nonsyn}_sfs.fs
  recurrent:    dadi_analysis_recurrent_folded/<run>/n100_d5mb/...
For a folded SFS eta_i (i = 1..n/2), pi = sum i(n-i) eta_i / C(n,2) and theta_W = S / a_n.
1 - pi/theta_W is 0 in expectation under the standard neutral model and > 0 when rare
variants are overrepresented. Tajima's D is also reported for reference.

Outputs (working-data/202610/rare_excess/):
  rare_excess_by_rep.csv        one row per replicate x window x class
  rare_excess_summary.csv       per condition: pooled-SFS value and per-replicate median/IQR
"""
import os
import re

import numpy as np
import pandas as pd

BASE = "/u/home/e/eewade/project-awbigham/data/working-data/202610"
OUT = os.path.join(BASE, "rare_excess")
N = 100
RUN_RE = re.compile(r"^(humansim_(sweep|neutral)|recurrentsweep)_dfehuman_mu[^_]+_bens([0-9.]+)"
                    r"(?:_benfrac([0-9.]+))?_rec([0-9e.-]+)_seq[0-9]+_rep([0-9]+)$")


def read_fs(f):
    lines = open(f).read().split("\n")
    v = np.array(lines[1].split(), dtype=float)
    mask = np.array(lines[2].split(), dtype=int) == 1
    v[mask] = 0.0
    return v


def stats(fs, n=N):
    i = np.arange(len(fs))
    poly = (i > 0) & (i < n)
    eta, i = fs[poly], i[poly]
    S = eta.sum()
    a1 = np.sum(1 / np.arange(1, n))
    a2 = np.sum(1 / np.arange(1, n) ** 2)
    pi = np.sum(i * (n - i) * eta) / (n * (n - 1) / 2)
    thw = S / a1
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / a1 ** 2
    e1, e2 = c1 / a1, c2 / (a1 ** 2 + a2)
    D = (pi - thw) / np.sqrt(e1 * S + e2 * S * (S - 1)) if S > 1 else np.nan
    return S, pi, thw, (1 - pi / thw) if thw > 0 else np.nan, D


rows, pooled = [], {}
for root, kind_default, dists in [("dadi_analysis_folded", None, (1, 3, 5)),
                                  ("dadi_analysis_recurrent_folded", "recurrent", (5,))]:
    for run in sorted(os.listdir(os.path.join(BASE, root))):
        m = RUN_RE.match(run)
        if not m:
            continue
        kind = kind_default or m.group(2)
        for d in dists:
            for cls in ("syn", "nonsyn"):
                f = os.path.join(BASE, root, run, f"n{N}_d{d}mb", "two_epoch", "dfe", f"empirical_{cls}_sfs.fs")
                if not os.path.exists(f):
                    continue
                fs = read_fs(f)
                S, pi, thw, rx, D = stats(fs)
                key = (kind, float(m.group(3)), m.group(4) or "", m.group(5), d, cls)
                rows.append(dict(zip(["kind", "ben_s", "benfrac", "rec", "d_mb", "class"], key),
                                 rep=int(m.group(6)), S=S, pi=pi, theta_w=thw, rare_excess=rx, tajima_d=D))
                pooled[key] = pooled.get(key, 0) + fs

os.makedirs(OUT, exist_ok=True)
by_rep = pd.DataFrame(rows).sort_values(["kind", "ben_s", "benfrac", "rec", "d_mb", "class", "rep"])
by_rep.to_csv(os.path.join(OUT, "rare_excess_by_rep.csv"), index=False)

keys = ["kind", "ben_s", "benfrac", "rec", "d_mb", "class"]
summ = by_rep.groupby(keys).agg(n_reps=("rep", "size"),
                                median=("rare_excess", "median"),
                                q25=("rare_excess", lambda x: x.quantile(0.25)),
                                q75=("rare_excess", lambda x: x.quantile(0.75)),
                                tajima_d_median=("tajima_d", "median")).reset_index()
summ["pooled"] = [stats(pooled[tuple(r[k] for k in keys)])[3] for _, r in summ.iterrows()]
summ.to_csv(os.path.join(OUT, "rare_excess_summary.csv"), index=False)
print(f"{len(by_rep)} rows -> {OUT}")
