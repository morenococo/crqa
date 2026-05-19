# crqa (cross-recurrence quantification analysis)

# crqa 2.1.0 (2026-05-14, final)

A major performance and modernisation update. All numerical outputs are
backward-compatible **except for `RR` when `tw > 0` or `side != "both"`**
(see "Behavioural changes" below).

## Behavioural fix: global normalisation in windowed functions (2026-05-19)

* **`wincrqa()`, `windowdrp()`, and `piecewiseRQA()`** now apply
  `normalize` to the **full** input series once, before any windowing
  or block decomposition. In previous versions `normalize` was applied
  independently inside each `crqa()` call, meaning every window had its
  own unit-interval or z-score reference — making windows
  non-comparable and causing the effective threshold to shift across
  windows (GitHub issue #19, reported by cmicek1).

  `rescale` is intentionally unchanged: it rescales the **distance
  matrix**, whose statistics (mean, max, min) are inherently local to
  each window. Pre-computing a global rescale factor would require
  materialising the full N×N distance matrix; users who need global
  rescaling should pre-compute the factor externally and set
  `rescale = 0`.

  This is a **behaviour change** when `normalize > 0`. Results for
  `normalize = 0` (the default) are identical to previous versions.

## Performance (Stage 3c additions, 2026-05-18)

* **OpenMP parallelism in the fused C++ kernel** (`src/crqa_fused.cpp`).
  The O(N·M) distance-and-threshold loop now runs in parallel across
  post-transpose columns using `#pragma omp parallel for` with static
  scheduling. Thread-local `(ii, jj)` vectors are concatenated in thread
  order after the parallel region, which preserves the column-major
  layout downstream code expects. Results are **bit-for-bit identical**
  to the serial reference at any thread count: the parallel region only
  performs floating-point comparisons (`d ≤ radius`), never reductions,
  so threshold decisions do not depend on thread count. Verified across
  an 18-case validation battery (all rescale modes, both metric
  families, Theiler windows, side masks, mdcrqa, both `rr_denom`
  conventions). The companion `crqa_rescale_stat()` is deliberately
  kept serial: it does compute an FP sum, and keeping it serial
  guarantees the rescale factor (and hence the effective radius) is
  identical to the v2.0.7 reference.

  Measured wall-time speedup at N = 20 000 (4 physical cores, Ryzen
  5800H, Linux):

  | RR regime | Serial | 4 threads | Speedup |
  |---|---|---|---|
  | Sparse (~0.3 %, typical) | 3.22 s | 1.67 s | 1.93× |
  | Dense  (~2.5 %)         | 8.31 s | 6.75 s | 1.23× |

  Sparse cases (the typical CRQA regime) scale near-linearly with
  cores. The dense case is bottlenecked by the sequential
  diagonal-line sort (O(nnz log nnz)) — a future optimisation target.
  Expect ~3–4× on 8-core systems and N ≥ 20 000.

* **Build infrastructure.** New `src/Makevars` and `src/Makevars.win`
  invoke `$(SHLIB_OPENMP_CXXFLAGS)` — R's portable OpenMP macro.
  `DESCRIPTION` adds `SystemRequirements: GNU make`. Wrapped in
  `#ifdef _OPENMP` guards, the kernel falls back to single-threaded
  execution wherever OpenMP is unavailable (e.g. some macOS source
  builds with Apple clang lacking libomp). No platform requires
  OpenMP to build, install, or run the package.

* **Thread-count control.** crqa respects the standard `OMP_NUM_THREADS`
  environment variable. By default the OpenMP runtime uses all
  available cores. To override, set `OMP_NUM_THREADS=N` in the shell
  *before* launching R — `Sys.setenv()` from inside R is too late
  because libgomp caches its thread count at first parallel region.
  When combining the wrapper `workers` argument with OpenMP, set
  `OMP_NUM_THREADS=1` in the parent session to avoid CPU
  over-subscription — `wincrqa(workers = 4)` on a 4-core machine with
  default OpenMP would otherwise spawn up to 16 threads.

## Performance (Stage 3b additions, 2026-05-14)

* **Fortran backend retired.** `src/jspd.f90` and `src/init.c` deleted;
  `spdiags()` reimplemented as ~30 lines of vectorised pure R. Package
  becomes `NeedsCompilation: yes` only for the new C++ kernel (not
  Fortran). Removes the Windows Rtools 4.5 / GCC 14 DLL crash that
  affected v2.1.0-dev builds.

* **Fused Rcpp inner loop** (`src/crqa_fused.cpp`). When `metric` is
  `"euclidean"`, `"maximum"` or `"manhattan"` (the common cases), the
  entire `cdist → dm → threshold → sparseMatrix → theiler → line_stats`
  chain is replaced by a single C++ pass that never materialises the
  N×N distance matrix or recurrence plot. Memory is O(N + nnz) instead
  of O(N²); runtime is 5–15× faster than the Stage 2 path at N ≥ 2000.
  Other metrics still go through the legacy `cdist`-based path (exact
  same outputs).

  Benchmarks at embed = 3, Gaussian input:

  | N | fused | legacy (Stage 2) | speedup | RAM ratio |
  |---|---|---|---|---|
  | 2 000 | 0.044 s | 0.653 s | 14.9× | 1.5× |
  | 4 000 | 0.169 s | 2.34 s | 13.9× | 2.7× |
  | 8 000 | 0.77 s | (1.3 GB) | — | ≥ 7× |
  | 16 000 | 4.15 s | OOM | — | — |

* **In-place rescale.** When `rescale > 0`, the rescaling factor is
  now folded into `radius` before entering the fused kernel (`dm / s <=
  r` ⟺ `dm <= r * s`), so no second N×N copy of `dm` is allocated.
  Halves peak RAM at `rescale > 0` vs the Stage 2 path.

## Performance

* **Stage 1 — outer-loop parallelism.** `wincrqa()`, `windowdrp()` and `piecewiseRQA()` now accept a `workers` argument (default `max(1L, future::availableCores() - 1L)`) and dispatch their per-window / per-block computations via `future` + `furrr`. Setting `workers = 1` reproduces the previous serial behaviour. Each worker still pays the full per-window RAM cost — see the documentation for usage notes.

* **Stage 2 — `line_stats()` replaces the `spdiags` + `tt` chain.** Inside `crqa()`, the dense `B` matrix from `spdiags()` plus the `tt()` vertical scan are no longer materialised. A single pass over the sparse recurrence indices computes `diaglines`, `lam`, `TT` and `max_vertlength` in `O(k log k)` time instead of `O(N^2)`. Empirical scaling exponent dropped from 2.26 to 1.79 on the Roessler benchmark; per-call wall time at N=10000 dropped from ~22s to ~4s; the maximum N reachable within a 30-second budget grew from ~10,500 to ~25,000.

## Behavioural changes

* **`RR` denominator now excludes Theiler-blanked and side-blanked cells.** A new helper `theiler_exclusion(m, n, w)` computes the exact number of cells in the Theiler band of width `w` for an arbitrary `m x n` matrix (works for rectangular RPs). For `tw = 0` and `side = "both"`, `RR` is identical to v2.0.7; otherwise it is larger than before by a factor of `(v1l * v2l) / region`, because Theiler/side-blanked cells no longer inflate the denominator. Concept adapted from pjbruna's community PR, generalised to rectangular matrices, and with the silent `tw = 0 -> tw = 1` coercion removed.

* **The `whiteline` argument is now ignored inside `crqa()`.** In all prior versions, `tt()` was called with this argument but its white-line output was never returned in the results list. The behaviour is therefore unchanged for users; `whiteline` remains in the function signature for backward compatibility but no longer affects timings. Callers who need white-line statistics can still invoke `tt()` directly.

## New functions

* `theiler_exclusion(m, n = m, w = 1)` — analytic cell count for the Theiler band (square and rectangular RPs).
* `line_stats(S, mindiagline, minvertline)` — single-pass diagonal and vertical line scan on a sparse recurrence matrix. Used internally by `crqa()`; exposed for users who want the same statistics from a pre-computed RP.
* `aRQA(ts1, ts2, delay, embed, radius, mindiagline, normalize)` — Stage 3a: approximative RQA following Schultz, Spiegel, Marwan & Albayrak (2015, Phys. Lett. A 379:997-1011) and Spiegel, Schultz & Marwan (2016). Computes RR, DET and L via phase-space histogram binning without forming the recurrence matrix. Scaling: ~O(N) at the time series sizes used here — measured 6 s at N=200,000 vs. exact `crqa()` OOMing at N=12,000.
* `rosslerattractor(numsteps, dt, a, b, c)` — companion to `lorenzattractor()`. Used by the package's benchmark/validation scripts and by users wanting a second canonical chaotic test system.

## New method

* `method = "aRQA"` in `crqa()` dispatches to the approximative path. Parameters honoured: `delay`, `embed`, `radius`, `mindiagline`, `normalize`. Parameters silently ignored on this path (because the algorithm never materialises the recurrence matrix): `tw`, `side`, `whiteline`, `minvertline`, `metric`, `recpt`, `rescale`. Returns the standard `crqa()` output structure with `RR`, `DET`, `L` populated; `NRLINE`, `maxL`, `ENTR`, `rENTR`, `LAM`, `TT`, `max_vertlength`, `catH`, `RP` are set to `NA` (full line-length distribution is not available from the histogram alone). Approximation error vs. exact computation: ~2 pp on DET for stochastic data, ~10 pp for weakly deterministic, larger for strongly deterministic systems with most recurrences on long diagonals.

## Validation

A 182-case test sweep verifies that `DET`, `NRLINE`, `maxL`, `L`, `ENTR`, `rENTR`, `LAM`, `TT`, and `max_vertlength` are bit-for-bit identical between v2.0.7 and v2.1.0 across all combinations of `radius`, `embed`, `delay`, `tw`, `side`, `rescale`, square/rectangular matrices, and the `rqa`/`crqa` methods. `RR` matches exactly for the 21 cases with `tw = 0` and `side = "both"`, and changes by the predicted factor `(v1l * v2l) / region` in the other 161 cases.

## Planned in the next release(s)

* Approximative RQA (`method = "aRQA"`) following Schultz et al. 2015 (Phys. Lett. A 379) and Spiegel et al. 2016, targeting `N >= 1e6`.
* Border-effect corrections (`dibo`, `kelo`, `censi`, window masking) for the diagonal-line entropy bias documented by Kraemer & Marwan 2019 (Phys. Lett. A 383).
* Lacunarity (Braun et al. 2021, Nonlinear Dynamics) as a new RQA quantifier.
* Quantile-based recurrence-threshold helpers in `optimizeParam()`.
* Vectorised parameter sweeps (inspired by the API of AccRQA) for grid exploration of `(tau, embed, threshold)`.
* Edit-distance recurrence for event-like / language data (Suzuki, Hirata & Aihara 2010, *Int. J. Bifurcat. Chaos* 20; Hirata & Aihara 2015, *Chaos* 25). Relevant for the linguistic-data community using `crqa`.

# crqa 2.0.6

* Rewritten jspd.f from 77 to Fortran compiler 90/95 (jspd.f90)

* Substituted the plotRP function with a more modern and customizable plot_rp function which is build on top of the ggplot2() ecosystem

* Added a method check to crqa() such that it will stop, rather than simply default, if one of the methods available is not correctly specified

* When method "mdcrqa" is used within crqa(), time-series will be z-scored using the scale() function as to keep each time-series independently 

* Removed the plot3D option from lorenzattractor.R and so the dependency to this package.

* Fixed the mdFnn function such that the first estimate for the embedding parameter is contingent on the data structure.  Also, relaxed the constraint that the size of the data matrix must be bigger 1, so that the function can be used with unidimensional time series.

* In optimizeParam(), when estimating false.nearest we incorporate the identified optimal delay and (if it exists) the user-specified maximum embedding dimension.

# crqa 2.0.5

* Amended plotRP() to handle better visualization 

* Checked calculation of categorical entropy (catH) and adjusted threshold to < .01 to compute it 

* Fixed issues with introducing RPs as argument for wincrqa 

# crqa 2.0.4

* Mostly fixed a couple of minor issues with wincrqa()

  * Added a warning in wincrqa() when users set delay and embedding dimension bigger than window size.
  
  * Added the calculation of nr. of windows for method `rqa`, i.e., unidimensional recurrence, which was missing

* Removed from `crqa_package.Rd` information about package version as it was conflicting with DESCRIPTION

* Added proper article CITATION to the package

# crqa 2.0.3

* Fixed minor bug on drpfromts() when method was `mdcrqa`. The argument `datatype` was forcing the input data to be in vector form, and should instead be left as a matrix. Added if statement to check (line 25)

# crqa 2.0.1

* Fixed bug on wincrqa() and windowdrp() to run windowed RQA on multidimensional data.
    * Simplified output of wincrqa() now returning a dataframe
    * Included new import from FSA() package to use diags() function. 

# crqa 2.0.0

* Major features

  * Extension of recurrence analysis to multidimensional time-series data, and significant update of computational procedures in `crqa()`
    * Added arguments: `method`, `metric` and `datatype`
    * Improved method for phase-space reconstruction (line 207:237)
    * Included computation of categorical entropy (line 411:420)
    * Removal of argument `checkl`
    
  * Simplified structure of functions and better division between core and ancillary functions:
    * `crqa_helpers` contains several functions previously exported (e.g., `theiler` or `tt`) that are now only accessed internally by the `crqa()` package.
    
  * New functions `mdDelay()` and `mdFnn()` to estimate Average Mutual Information and False Nearest Neighbours of multidimensional time-series.

  * Experimental `piecewiseRQA()` function created to better handle the computational load of large time-series.

  * `optimizeParam()` now works also with multidimensional time-series. 
    
* Minor features

  * Deprecated functions: `CTcrqa()`, `runcrqa()`, `calcphi()`, `takephi()`

  * On `crqa()` 
    * improved error checking and warning messages (line 133:170)
    * added two more ways of rescaling the distance matrix (line 252:272)
    
  * On `plotRP()` added arguments to improve the plotting of Recurrence Plots 

  * `drpdfromts()` is now called `drpfromts()` and it has been rewritten to align with the new version of `crqa()`.

  * A convenience function called `numerify` in `crqa_helpers` is automatically called when a user inputs categorical series (i.e., it contains either characters or factors) and this function is used to recode the levels of such time-series into numerical codes (to run crqa). A warning is send to the user when `crqa()`.

  * `windowdrp()` has been rewritten to align with the new version of `crqa()`.  

  * `wincrqa()` has been rewritten to align with the new version of `crqa()` and better names for the output were provided.  

# crqa 1.0.9

* On `drpdfromts()` removed a left over constant used for testing (line: 42)

# crqa 1.0.8

* Improvements and bug fixes

  * On `drpdfromts()` fixed initialisation of dimensions for empty RP (line: 52)

# crqa 1.0.7

* New functions

  * `plotRP()` convenience function based on the standard `plot()` to visualize a Recurrence Plot

* Improvements and bug fixes

  * On `crqa()` added a few more checks (`stop`) if the data inputted did not comply with the function, and send a warning message.

  * `drpdfromts()` entirely rewritten around `crqa()` to better deal with continuous valued time-series.

  * `runcrqa()` fixed to fit with the revised functions: `drpdfromts()` and `windowdrp()`, 

  * `tt()` fixed `rBind` (line 23 and 91) which was deprecated from the `Matrix()` package.

  * `windowdrp()`entirely rewritten around the new version of `drpdfromts()` 

  * `wincrqa()` adjusted indexing of windows (line: 40:41)

# crqa 1.0.6

* New functions

  * `ami()` externalised from `optimizeParam()`

  * `lorenzattractor()` simulates and plots 3D data from a Lorenz Attractor

* Improvements and bug fixes

  * On `crqa()` include a `stop` (line: 95:100) if time-series were shorter than their phase space reconstructed portaits

  * On `optimizeParam()` included argument `typeami` to set the type of `ami()` desired (either, minimum dip or maximum lag)

# crqa 1.0.5

* Improvements and bug fixes

  * On `checkts()`. added argument `pad` (line: 38:67), which gives the option, in case of series of different length to extend the shortest sequence either with mean value (if the variable is in a continuous scale) or with a random label not present in either series (if the variable is categorical).

  * On `crqa()`: 
    * Added default values for arguments in the call of the function. 
    * Included the argument `side` to select the region of the Recurrence Plot to extract measures on
    * Included the argument `checkl`, a wrapper to call the function `checkts()` directly inside this function. 

  * On `optimizeParam()`:
    * Added arguments `min.rec` and `max.rec` to the call. 
    * Added argument within `par` (`fnnpercent`) to estimate False Nearest Neighbours based on a percentage reduction with respect to first dimension (line: 199:241)
    * Simplified code throughout and markedly improved the estimation of the radius (line: 299:360)
    
  * On `runcrqa()` included argument `pad` in the call to work with the revised version of `checkts()`

  * On `wincrqa()` added calculation of TREND (line: 64:87) 

# crqa 1.0.4

  * crqa 1.0.3 did not pass CRAN check because of `DESCRIPTION` (`Depends` field) and was resubmitted as 1.0.4.

# crqa 1.0.3

* New function

  * Added `theiler()` function to choose the separation between values on the time series when specifying a delay reconstruction vector, i.e., the Theiler window.

* Improvements and bug fixes

  * On `crqa()` added theiler window (`theiler()`, line: 170:175)

  * On `optimizeParam()` improved calculation of average mutual information (line 56:119), and added estimation of radius within user specified expected recurrence values (line 237:284)

  * On `runcrqa()` added missing arguments when calling `wincrqa()` (line 104-110:112) 

  * On `tt()` simplified calculation of laminarity (line 65)

  * On `wincrqa()` added missing arguments in the function call to exploit better functionality in `crqa()`

# crqa 1.0.2

  * On `tt()` added in line comments as header of function.

# crqa 1.0.1

* Improvements and bug fixes

  * On `crqa()`: Implementation of embedding dimensions and phase space reconstruction, added in line comments in the code.

# crqa 1.0

First version of the package featuring the following original functions:

  * `calcphi`
  * `checkts`
  * `crqa`
  * `CTcrqa`
  * `drpdfromts`
  * `optimizeParam`
  * `runcrqa`
  * `simts`
  * `spdiags`
  * `takephi`
  * `tt`
  * `wincrqa`
  * `windowdrp`