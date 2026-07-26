# Vendored scDesign2 (patched fork)

Source: /mnt/storage/jaquino/scDesign2/scDesign2 (jaquino's working tree)
Base:   JSB-UCLA/scDesign2 @ f32cfe0 (master HEAD, 2024-11-01)
Patches (uncommitted in the origin tree, vendored here as-is):
  - R/data_simulation.R  -> data_simulation.diff.R  (adds fold_change= for DE/DTU injection)
  - R/model_fitting.R    -> model_fitting.fastglm.R (fastglm-based marginal fitting)
Installed into r.sif via `R CMD INSTALL vendor/scDesign2`. Do NOT replace with the
public scDesign2 -- it lacks simulate_count_scDesign2(fold_change=).
