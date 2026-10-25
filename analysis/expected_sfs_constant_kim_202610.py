#!/usr/bin/env python
"""Expected nonsynonymous SFS under constant size (N = 10,000, the simulated
census size) and the true human gamma DFE of Kim et al. 2017, for the manuscript
figures (make_manuscript_figures_202610.R, Fig 2B).

Built the same way as the fitdadi gamma fits (fit_dfe_poisson_test.py): dadi
Cache1D over Demographics1D.snm, gamma_bounds (1e-5, 500), 100 gamma points.
SLiM s (fitness 1 + s for homozygotes, h = 0.5) = 2 x dadi s, so the dadi-scale
gamma scale = 2 * N * (mean_s_slim / 2) / shape.

Writes one proportion-normalised spectrum per (n, fold) to
working-data/202610/manuscript_figures/data/constant_kim_nonsyn_sfs.csv
(columns n, fold, allele_count, prop). Theta cancels in proportions.

Usage (dadi_env):  python expected_sfs_constant_kim_202610.py [n ...]
"""
import csv
import os
import sys

import dadi
import dadi.DFE as DFE
import numpy as np

N = 10000
SHAPE = 0.186
MEAN_S_SLIM = 0.01314833
scale_gamma = 2 * N * (MEAN_S_SLIM / 2) / SHAPE
OUT = "/u/home/e/eewade/project-awbigham/data/working-data/202610/manuscript_figures/data"

ns = [int(x) for x in sys.argv[1:]] or [100, 1000]
rows = []
for n in ns:
    pts = [1200, 1600, 2000]  # same grid as the fits (fit_demo_dfe_human_202610.sh PTS)
    cache = DFE.Cache1D([], [n], DFE.DemogSelModels.equil, pts=pts,
                        gamma_bounds=(1e-5, 500), gamma_pts=100, cpus=int(os.environ.get("NSLOTS", 4)))
    fs = cache.integrate([SHAPE, scale_gamma], None, DFE.PDFs.gamma, 1.0, None)  # (params, ns, sel_dist, theta, pts)
    for fold, spec in (("unfolded", fs), ("folded", fs.fold())):
        v = np.asarray(spec.data, float)
        mask = np.asarray(spec.mask, bool)
        keep = [i for i in range(1, n) if not mask[i]]
        tot = sum(v[i] for i in keep)
        rows += [dict(n=n, fold=fold, allele_count=i, prop=v[i] / tot) for i in keep]
    print(f"n = {n} done", flush=True)

os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/constant_kim_nonsyn_sfs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["n", "fold", "allele_count", "prop"])
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} rows -> {OUT}/constant_kim_nonsyn_sfs.csv")
