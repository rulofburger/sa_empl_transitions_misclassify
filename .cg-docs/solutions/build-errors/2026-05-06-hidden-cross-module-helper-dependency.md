---
date: 2026-05-06
title: "Hidden cross-module helper dependency via source order"
category: "build-errors"
language: "R"
tags: [source-order, module-design, helper-functions, EM-algorithm, shared-utils]
root-cause: "Helper functions defined at the bottom of one E-step file were called by two other E-step files with no explicit import, relying on source order to inject them into .GlobalEnv"
severity: "P1"
---

# Hidden cross-module helper dependency via source order

## Problem

Three E-step modules (`estep_covariates.R`, `estep_fmm.R`, `estep_inconsistency.R`)
all needed the same two helper functions (`.log_markov_trans_indiv` and
`.log_misclass_wave_ext`). The helpers were defined at the **bottom** of
`estep_covariates.R`, which happened to be sourced first.

`estep_fmm.R` and `estep_inconsistency.R` called these helpers with no
declaration, trusting that `estep_covariates.R` had already been sourced.

**Symptoms**:
- Code worked fine when sourced in full via `source_all.R`.
- Individual test files that sourced only `estep_fmm.R` or `estep_inconsistency.R`
  silently failed at runtime with `could not find function ".log_markov_trans_indiv"`.
- Review automation (`helper-source.R`) would break if extension source order
  was ever changed.
- The static code analysis could not find the definition because it was in a
  different file with no `@import` or `source()` directive.

## Root Cause

In R's `source(local=FALSE)` pattern, all sourced files share `.GlobalEnv`.
Any function defined by an earlier `source()` call is visible to later calls.
This creates an invisible runtime dependency on source order. The dependency
was not documented and not enforced by any mechanism.

The original author may have intended these as "private" helpers within the
module, but placing them at the bottom of the file rather than in a dedicated
helpers file made them invisible to other modules at the static level while
still being injected at runtime.

## Solution

Extract all shared helpers into a dedicated file (`helpers_ext.R`) and source it
**first**, before any module that depends on it:

```r
# source_all.R — correct order
source(file.path(.ext_r, "helpers_ext.R"),             local = FALSE)  # FIRST
source(file.path(.ext_r, "compute_inconsistencies.R"), local = FALSE)
source(file.path(.ext_r, "prepare_covariates.R"),      local = FALSE)
source(file.path(.ext_r, "estep_covariates.R"),        local = FALSE)
source(file.path(.ext_r, "estep_fmm.R"),               local = FALSE)
# ...
```

The same `source(helpers_ext.R)` call must appear in every entry point:
`source_all.R`, `helper-source.R`, and any standalone script.

**Anti-pattern** (what was done):
```r
# estep_covariates.R — bottom of file
.log_markov_trans_indiv <- function(...) { ... }  # used by estep_fmm.R too!
.log_misclass_wave_ext  <- function(...) { ... }  # used by estep_inconsistency.R too!
```

**Correct pattern**:
```r
# helpers_ext.R — standalone file, sourced first
.log_markov_trans_indiv <- function(...) { ... }
.log_misclass_wave_ext  <- function(...) { ... }
```

## Prevention

- **Rule**: Any helper function called by more than one module **must** live in
  a dedicated helpers file. Never define cross-module helpers at the bottom of
  one of the modules that uses them.
- **Naming convention**: Shared helper files are named `helpers_<module_prefix>.R`
  and listed first in `source_all.R` with a comment `# MUST precede all modules
  that depend on it`.
- **Verification**: When adding a helper to a module file, grep for its name
  across all `.R` files. If it appears in more than one file, extract it.

## Related

- [2026-05-05-validate-false-pattern-em-hot-path.md](../performance-issues/2026-05-05-validate-false-pattern-em-hot-path.md) — related
  pattern for safe fast paths in source-based modules
- [2026-05-06-cross-layer-function-placement.md](./2026-05-06-cross-layer-function-placement.md) —
  related: functions placed in the wrong *layer* (upward dependency), not just
  the wrong file
