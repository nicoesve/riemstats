# riemstats 0.2.0.9000 (development)

## `riem_anova()`: p-value construction and tail orientation

Two defects found by operation `package-test-hardening` (T6, 2026-08-20) and
fixed on 2026-08-27.

### The tail is now chosen by the statistic (behaviour change for `pillais_trace`)

`riem_anova()` counted permutations in the left tail unconditionally. That is
correct for `log_wilks_lambda` (Wilks' Lambda is small when the groups differ)
and backwards for `pillais_trace` (`tr((T - W) T^-1)` is large when they
differ), which is exported and documented as a `stat_fun` for the same call.
Used that way the function returned `1 - p`, and it failed silently: a clearly
separated pair of groups produced `p = 1.000` rather than an error.

The statistics shipped with the package now declare their orientation through a
`riemstats_tail` attribute, which `riem_anova()` reads. A new `tail` argument
overrides it, and a user-supplied statistic that declares no orientation is
still treated as left-tailed.

**Results computed with `stat_fun = pillais_trace` on earlier versions are
wrong and should be recomputed.** Results computed with the default
`log_wilks_lambda` are unaffected.

### The p-value estimator is now selectable, and the default is unchanged

`riem_anova()` returned `#{v < t_obs} / nperm`, omitting the observed labelling
from the reference set. Under the null the rank of the observed statistic among
the `nperm + 1` exchangeable values is uniform, so as coded
`P(p <= alpha) = (floor(alpha * nperm) + 1) / (nperm + 1)`. The gap is
negligible at `nperm = 999` and material at small `nperm`: at `nperm = 49` a
nominal 0.05 test has size 0.060.

Both estimators are now reachable through a new `p_method` argument:

- `p_method = "plugin"` (**the default**) reports `k / nperm`, counting
  strictly more extreme permutations. This is the construction the function has
  always used. It is an *unbiased* estimator of the permutation p-value, but it
  is not level-preserving and can return exactly 0.
- `p_method = "exact"` reports `(1 + k) / (nperm + 1)`, treating the observed
  labelling as one of the `nperm + 1` equally likely labellings. It is *valid*
  -- its size never exceeds the nominal level -- at the cost of a small upward
  bias, and it can never return 0.

**The default is deliberately unchanged.** Results published from earlier
versions reproduce by the default call. `"exact"` is recommended for new work,
and the default is expected to change in a future major version.

### Documentation correction

`?riem_anova` described the opposite tail from the implementation and claimed
the procedure was exact. No exactness claim is available under `"plugin"`.

# riemstats 0.2.0

## Major Changes

### Permutation Test Implementation

- **BREAKING CHANGE**: Replaced parametric bootstrap with permutation test in `riem_anova()`
  - Function parameter renamed: `den` → `nperm` (default changed from 5 to 1000)
  - Old function `one_bootstrap()` removed and replaced with `one_permutation()`
  - Users must update existing code to use the new parameter name

**Benefits of permutation test approach:**
- Faster computation (no synthetic data generation required)
- More stable p-values (no parameter estimation variability)
- Exact test under null hypothesis (when group labels are exchangeable)
- Simpler implementation (direct shuffling of observations)

**Motivation:**
- Asymptotic tests showed inflation issues
- Bootstrap was computationally expensive and under-conservative
- Permutation test is more natural for comparing sub-populations

See GitHub issue #61 and PR #62 for detailed discussion and implementation.

## Documentation

- Updated `riem_anova()` documentation with permutation test description
- Added `one_permutation()` function documentation
- Updated vignette with working permutation test examples
- Reduced vignette permutation count to 100 for faster builds

## Testing

- Updated test suite for permutation approach
- All 126 tests passing
- Fixed CSample constructor issues in tests

## Note on Validation

This release is pending comprehensive validation through simulation studies (EXP-MP-062).
Initial testing shows correct behavior, but full Type I error and power analysis is ongoing.

# riemstats 0.1.0

## Dependencies
- Depends on R (>= 4.3.0)
- Imports the following packages:
  - Matrix
  - methods
  - expm
  - R6
  - purrr
  - MASS
  - furrr
- Suggests the following packages for testing and documentation:
  - testthat (>= 3.0.0)
  - knitr
  - rmarkdown

## Miscellaneous
- Added a `LICENSE` file with the MIT license.
- Maintainer: Nicolas Escobar <nescoba@iu.edu>
