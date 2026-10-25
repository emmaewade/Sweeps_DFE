# sweeps-dfe

Scripts for  Wade EE, Lohmueller KE. Robustness of the inference of the distribution of fitness effects to selective sweep. *BioRxiv*. (2026)

1. `simulation/`: SLiM single-sweep (`ConstantPopSizeOneSweep_v4_genic.slim`) and recurrent-sweep
   (`RecurrentSweep_v2_genic.slim`) arrays; the task-argument files list every simulated condition
2. `fitting/`: dadi demography (`fit_demo_fixedgrid_test.py`) then fitdadi gamma DFE
   (`fit_dfe_poisson_test.py`) 
3. `analysis/`: collect DFE fits (`compare_dfe_single_recurrent_202610.py`), rare-variant
   excess, and the expected nonsynonymous SFS under constant size and the Kim et al. 2017 DFE
4. `figures/`: main and supplementary figures and tables.


