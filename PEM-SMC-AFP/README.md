# PEM-SMC-AFP — MATLAB Package

PEM-SMC-AFP is a tempered Sequential Monte Carlo (SMC) sampler with  
**A**daptive tempering, **F**lexible move scheduling, and **P**arallelism.  
It targets **bounded**, **moderate- to high-dimensional** log-densities (e.g., Bayesian posteriors)
and combines adaptive annealing with genetic and random-walk moves.

---

## 1. Overview

PEM-SMC-AFP implements an annealed / tempered SMC sampler in MATLAB.  
At each stage it performs

> **reweight → (policy-triggered) resample → move operators**

with the temperature increment chosen **adaptively** from a **relative CESS band**
(or from a fixed schedule if adaptive tempering is turned off).

Particle rejuvenation combines three kernels:

1. **ARM** — Adaptive Random-Walk Metropolis (with fold/reflect boundaries)  
2. **XOVER** — two-parent single-point **crossover**  
3. **DE–MH** — Differential-Evolution Metropolis–Hastings

All proposals enforce **box constraints** by **folding** (`fold` or `reflect`).
Parallel log-density evaluation is supported via MATLAB `parfor`.

You provide a **log target density** (function handle or local `target.m`); the sampler
returns the particle trajectory and diagnostics.

---

## 2. Project Layout

```text
/PEM-SMC-AFP/
  README.md

  /core/
    PEM_SMC_AFP.m          % main sampler (public API)
    install_PEM_SMC_AFP.m  % add core / postprocessing to path
    ResampSys.m            % systematic resampling
    sequencesGen.m         % optional fixed tempering schedule
    randw.m                % ARM proposals with fold/reflect
    Generatep_fold.m       % DE–MH proposals with fold/reflect
    % + other internal helpers

  /example_1/              % 2-D 20-mode Gaussian mixture
    example_1.m
    target.m

  /example_2/              % 100-D bimodal Gaussian
    example_2.m
    target.m

  /example_3/              % Synthetic CoLM calibration
    example_3.m
    target.m

  /example_4/              % Real-world CoLM calibration
    example_4.m
    target.m

  /postprocessing/
    % plotting & diagnostics utilities used by the examples
```

**Roles (brief):**

- `PEM_SMC_AFP.m` – orchestrates tempering, (re)weighting, resampling, and move operators (ARM / XOVER / DE–MH).  
- `install_PEM_SMC_AFP.m` – adds `core/` (and optionally `postprocessing/`) to the MATLAB path.  
- `example_k/` – self-contained drivers and `target.m` for each case study.  
- `postprocessing/` – optional plotting and diagnostics routines.

---

## 3. Getting Started

### 3.1 Install (once per MATLAB session)

From the project **root**:

```matlab
>> cd path/to/PEM-SMC-AFP
>> core/install_PEM_SMC_AFP
```

This adds the **core** folder (and optionally `postprocessing/`) to the MATLAB path.  
Example folders are **not** added, to avoid collisions between different `target.m` files.

To persist the path across sessions:

```matlab
>> savepath
```

### 3.2 Run the included examples

2-D 20-mode mixture:

```matlab
>> cd example_1
>> example_1
```

100-D bimodal Gaussian:

```matlab
>> cd ..
>> cd example_2
>> example_2
```

Synthetic and real-world CoLM case studies:

```matlab
>> cd ..
>> cd example_3   % or example_4
>> example_3      % or example_4
```

> **Important:** always `cd` into an `example_k/` folder so that its local `target.m` is used.

---

## 4. API

### 4.1 Basic call

```matlab
parameter_iteration = PEM_SMC_AFP(Np, S, bound)
```

### 4.2 Full call

```matlab
[parameter_iteration, out] = PEM_SMC_AFP(Np, S, bound, logpdf_handle, options)
```

- `Np` *(int)* – number of particles  
- `S`  *(int)* – stage **budget** (upper bound if adaptive tempering is used)  
- `bound` *(2 × d double)* – parameter bounds: first row = lower, second row = upper  
- `logpdf_handle` *(optional)* – function handle returning scalar log-density  
  (if omitted, `target(theta)` from the current folder is used)  
- `options` *(optional struct)* – algorithm settings (see header of `PEM_SMC_AFP.m`)

Return values:

- `parameter_iteration` *(Np × d × S_used)* – particle states per stage  
- `out` – diagnostics struct (temperature schedule, ESS / rCESS, resampling flags,
  acceptance rates for ARM / XOVER / DE–MH, log-evidence, timing, etc.)

Final posterior particles:

```matlab
theta_post = parameter_iteration(:,:,end);
```

---

## 5. Options (high-level sketch)

Detailed field names and defaults are documented at the top of `core/PEM_SMC_AFP.m`.  
Here is a minimal example:

```matlab
opts = struct();

% Tempering: rCESS-band adaptive schedule
opts.Tempering.Adaptive     = true;
opts.Tempering.TargetBand   = [0.7 0.9];
opts.Tempering.MinDeltaBeta = 1e-4;
opts.Tempering.GrowthFactor = 1.5;
opts.Tempering.BandDecay    = true;   % optional

% Resampling policy
opts.Resample.Policy   = 'ESS';       % 'ESS' | 'Periodic' | 'Never'
opts.Resample.ESSalpha = 0.8;
opts.Resample.Method   = 'systematic';

% Move scheduling
opts.Moves.Sequence = {'ARM','XOVER','DEMH'};
opts.Moves.Repeats  = [1 1 1];

% Parallelism
opts.Parallel.Enabled    = true;
opts.Parallel.NumWorkers = [];        % [] = MATLAB decides
```

Typical usage:

```matlab
[traj, out] = PEM_SMC_AFP(Np, S, bound, @(th) target(th, data), opts);
```

---

## 6. Custom Targets

### Local `target.m`

Place a `target.m` in an example folder:

```matlab
function logL = target(theta, data)
% theta : 1×d parameter vector
% data  : optional struct
% logL  : scalar log density log π(theta)
% ...
end
```

Run from that folder:

```matlab
>> cd my_example
>> my_example   % calls PEM_SMC_AFP internally
```

Or call the sampler directly with a function handle, e.g.:

```matlab
>> core/install_PEM_SMC_AFP
>> cd my_example
>> logpdf = @(th) target(th, data);
>> [traj, out] = PEM_SMC_AFP(2000, 800, bounds, logpdf);
```

---

## 7. Practical Notes

- Start with `Np ≈ max(40*d, 1000)` and an rCESS band like `[0.7 0.9]`.  
- Aim for ARM / DE–MH acceptance rates roughly in the **0.2–0.4** range.  
- For stability, `Resample.Policy='ESS'` with `ESSalpha` in `[0.7, 0.9]`
  usually works well.  
- Ensure reproducibility with, e.g., `rng(1234, 'twister');`.

Troubleshooting hints:

```matlab
which target -all      % check which target.m is in use
opts.Parallel.Enabled = false;  % disable parallelism if needed
```

---

## 8. Citation

If this package contributes to your work, please cite the associated  
PEM-SMC-AFP paper / preprint describing the algorithm and case studies.

Code © Xu Cong (Lanzhou University).  
Released to support transparent, reproducible Bayesian computation.
