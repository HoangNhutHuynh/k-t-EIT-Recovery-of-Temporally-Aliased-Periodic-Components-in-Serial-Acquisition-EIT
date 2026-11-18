# k-t EIT: recovery of temporally aliased periodic components in serial-acquisition EIT

MATLAB code (with EIDORS) for the manuscript

> H. N. Huynh, "k-t EIT: Recovery of Temporally Aliased Periodic Components in Serial-Acquisition Electrical Impedance Tomography," submitted to *IEEE Transactions on Medical Imaging*.

Author: Hoang Nhut Huynh, National Chung Hsing University. License: MIT (see `LICENSE`).

## What the method does

Most EIT systems acquire the measurements of a frame one after another, so measurement *m* of frame *f* is taken at `t = f/fs + tau_m`. A periodic component above half the frame rate folds onto a lower frequency `beta`, but each measurement still carries the phase `exp(i 2 pi nu tau_m)` of the **true** frequency `nu`. k-t EIT puts the acquisition time of every measurement into the forward model:

```
Y(beta) = sum_k D(nu_k) J X_k,      D(nu) = diag(exp(i 2 pi nu tau_m))
```

- one component per frequency bin: phase derotation with the true `nu`, then the usual one-step reconstruction (same cost);
- components that fold onto the same bin: joint MAP solve, possible only because the acquisition is serial;
- the acquisition schedule `(duty, shift, direction)` is estimated from the data by MAP fitting and checked by two-fold cross-validation.

All methods (naive, phase correction, k-t, EIDORS `frame_timing_correct`, temporal prior) share the same Jacobian, prior and regularization, so differences come from the treatment of time only.

## Requirements

- MATLAB (Signal Processing Toolbox not required).
- EIDORS 3.12 (https://eidors3d.sourceforge.net). EIDORS is distributed under its own license (GPL) and is **not** included here. Either run its `startup.m` before using this code or set `cfg.eidors_startup` in `kt_config.m`.
- Python 3 with `numpy`, `scipy`, `matplotlib` (only for redrawing the manuscript figures, `figures/requirements.txt`).

## Data (not included)

The real records are public in the EIDORS contributed-data repository: https://eidors3d.sourceforge.net/data_contrib.shtml (sections *Animal Data* and *Human Data*). They are not redistributed here. Each contributor may set its own license; attribution is required, so cite the source publication of every record you use.

Download the records and place them as follows (folder and file names are those used by the scripts):

| Folder under `Data/` | Files | Record | Source publication |
|---|---|---|---|
| `chb_hfov_porcine/` | `CHB_pig06-55.get` | HFOV pig, 44 frames/s, 104 stored measurements | Wolf et al., *Pediatr. Crit. Care Med.* 11:610, 2010 |
| `pig_2003/` | `p1130107.get`, `p1130122.get` | PEEP-trial pig after lavage and after surfactant, 13 frames/s | Frerichs et al., *Intensive Care Med.* 29:2312, 2003 |
| `if-neonate-ventilated/` | `nm18.get` | ventilated infant, 13 frames/s | Heinrich et al., *Intensive Care Med.* 32:1392, 2006 |
| `Pig_torso/` | `g08-c003.get` ... `g10-c029.get` | two piglets, stepwise vital-capacity manoeuvre, 13 frames/s | Dargaville et al., *Intensive Care Med.* 36:1953, 2010 |

All records are Goe-MF II raw files (`.get`, 256 float32 values per frame) and are read by `eidors/kt_read_goe_raw.m`. The thorax contours come from the EIDORS `shape_library` (`pig_23kg`, `neonate`).

## Quick start

```matlab
% from the repository root, with EIDORS on the path
run_E0_selftest        % checks the k-t core on a synthetic problem (no EIDORS needed); all lines PASS
run_E1_sim             % main simulation, a few seconds
run_FINAL_E7           % frame-removal test on the real records (~10 min)
```

Every run writes a folder `runs/<date>_<time>_<experiment>_<tag>/` with `config.txt`, `console.txt` and `metrics.csv` (long format: `run_id, experiment, tag, param, pvalue, component, method, metric, value`), and appends to `runs/all_metrics.csv`, `runs/run_log.csv` and `runs/RUN_LOG.md`. Figures are written only when `cfg.save_figs = true`.

All settings are in `kt_config.m`; a script changes a copy of `cfg` and calls one experiment function, so a variant is a new block with a new `cfg.tag`.

## Reproducing the manuscript

| Manuscript item | Script | Output used |
|---|---|---|
| Simulation, oracle and 90° phase-lag cases | `run_FINAL.m` | `*_E1sim_FINAL_E1a_oracle`, `*_FINAL_E1d_lag90` |
| Fig. 2(a), principal angles | `run_FINAL.m` | `*_E3collision_FINAL_E3` |
| Fig. 2(b, c), nonlinear forward model, 20 noise seeds × 3 SNR | `run_R1_revision.m` (section X2, ~45 min) | `runs/R1_X2_nonlinear_sim.csv` |
| Table II, schedule estimation with cross-validation | `run_FINAL.m` (E4 to E6), `run_E10b_pig_torso_timing.m` | `*_E4real_FINAL_*`, `*_E4real_E10b_*` |
| Table III, HFOV regional measures (5-s blocks) | `run_FINAL.m` | `*_E4real_FINAL_E4_hfov` |
| Table IV, frame-removal results | `run_FINAL_E7.m`, `run_E10_pig_torso.m` | `*_E7decim_FINAL_E7_*`, `*_E7decim_E10_*` |
| Table V, regularization sweep; factor-3 removal; slow line | `run_E9_extras.m` | `*_E7decim_E9a_*`, `*_E9b_*`, `*_E9c_*` |
| Table VI, priors, FTC linear mode, f0 error; alias-order bound | `run_R1_revision.m` (X4, X5, X6, X7) | `*_R1_X*` |
| Figs. 3 to 5, images and error maps | `run_FINAL_E7.m`, `run_E10b_pig_torso_timing.m` (figure section), then `run_export_fig_data.m` | `*_data.mat` next to the `.fig` files |

Sections of `run_R1_revision.m` can be selected without editing the file:

```matlab
RUN = {'X7'}; run_R1_revision; clear RUN
```

### Redrawing the figures (Python)

After the MATLAB runs above:

```bash
pip install -r figures/requirements.txt
python figures/fig_sim.py              # Fig. 2
python figures/fig_images.py           # Figs. 3 to 5
python figures/graphical_abstract.py   # graphical abstract
```

Each script uses the newest matching run in `runs/` and writes PDF, SVG and PNG into `figures/`.

## Repository structure

```
kt_config.m, kt_addpaths.m   defaults and path set-up
run_*.m                      experiment scripts (E0 to E10, FINAL, R1); each block is one logged run
core/     acquisition-time model, alias map, spectra, noise weights, reconstruction matrix,
          k-t solve (single and joint), schedule estimation (MAP, free, cross-validation),
          separability (principal angles), temporal-prior baseline, comparison metrics
eidors/   EIDORS interface: models and thorax warping, phantom and simulation,
          Goe-MF II raw reader, row matching, reciprocity completion
exp/      experiment drivers: simulation, sweeps, collision, real data, frame removal,
          FTC baseline, nonlinearity index, image grids
log/      run logging, metrics export, figure export, figure-data export
figures/  Python scripts for the manuscript figures
```

Development scripts `run_E1` to `run_E8` document the path to the final configuration (for example, the coherence-based schedule estimator in `core/kt_estimate_timing.m` was replaced by the MAP estimator because regional phase differences bias it). They are kept for transparency; the manuscript numbers come from the scripts in the table above.

## Notes

- `frame_timing_correct` in EIDORS 3.12 resets its `method` argument to `'fft'` when called with two arguments; `exp/kt_try_ftc.m` passes three arguments so that the `'linear'` mode is really used.
- Records that store only half of the reciprocal measurement pairs (the HFOV record) are completed by reciprocity before `frame_timing_correct`, which assumes 208 measurements per frame.

## Citation

If you use this code, please cite the manuscript above (details will be updated after publication) and the source publication of every data record you use.
