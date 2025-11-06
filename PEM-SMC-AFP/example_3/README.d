# Test Case 2 (S1/S2) — Quick Guide

## Overview

A reproducible two-scenario (S1/S2) synthetic-observation case for CoLM with three streams: **NEE**, **LE**, and **RSM**.

## What each script does

* **data_prepare.m** — Orchestrates the pipeline: sets paths, creates `./data`, generates S1/S2 truth parameters, checks CoLM outputs exist, then creates synthetic observations.
* **gen_truth_params_S1S2.m** — Uses 40×2 prior bounds and a 40×1 default vector.
  S1 = defaults; S2 = Uniform(ai,bi) per parameter (seed=1234).
  Writes **`data/input_step_S1.txt`** and **`data/input_step_S2.txt`** (one value per line).
* **gen_synthetic_obs.m** — For each stream, σ = `0.05 * std(sim)`. Adds i.i.d. Gaussian noise (seed=1234) to non-NaN entries and preserves NaNs.
  Writes `obs_*` and `sigma_*` next to the model outputs.
* **CoLM (external)** — Run on your server using the S1/S2 parameter files to produce the “default” simulations `output_*_S#.txt`. This repo does **not** launch CoLM.

## Files in `./data`

**Naming pattern:** `<kind>_<VAR>_<STEP>.txt`

* `<kind>`: `input_step | output | obs | sigma`
* `<VAR>`: `NEE | LE | RSM`
* `<STEP>`: `S1 | S2`

### Meanings

* `input_step_S#.txt` — 40 truth parameters for scenario `#`, one number per line.
* `output_NEE_S#.txt` — Default CoLM NEE simulation (µmol m⁻² s⁻¹).
* `output_LE_S#.txt`  — Default CoLM LE simulation (W m⁻²).
* `output_RSM_S#.txt` — Default CoLM third stream (e.g., remote-sensing metric).
* `obs_XXX_S#.txt`    — Synthetic observations = default simulation + 5% Gaussian noise.
* `sigma_XXX_S#.txt`  — One scalar: homoscedastic 1σ used for `obs_XXX_S#.txt`.

## Sanity rules

* `bounds` is 40×2 and `defaults` has 40 values; defaults must lie within bounds.
* Within a step, LE/NEE/RSM series have the **same length**.
* Random seed fixed at **1234** → reproducible S2 parameters and noise.

## Provenance map

```
gen_truth_params_S1S2 → input_step_S1.txt, input_step_S2.txt
CoLM (external)       → output_*_S1.txt, output_*_S2.txt
gen_synthetic_obs     → obs_*_S#.txt,   sigma_*_S#.txt
```

## Minimal run checklist

1. Run `data_prepare.m` to write `input_step_S1.txt` / `input_step_S2.txt`.
2. On the server, run CoLM for S1/S2 and place the six `output_*_S#.txt` files into `./data`.
3. Re-run `data_prepare.m` (it will verify outputs and create `obs_*` and `sigma_*`).
