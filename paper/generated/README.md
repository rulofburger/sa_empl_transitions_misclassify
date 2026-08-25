# Reproducible paper artifacts

This directory contains the compact, paper-facing outputs for revised Tables
1--3 and their associated appendix material: LaTeX tables, result summaries
and analytical diagnostics. Large fitted model objects, bootstrap draws, raw
data, caches and timestamped run outputs remain ignored.

After re-estimating the affected models and running `build_tables.R`, refresh
this snapshot from the repository root with:

```sh
Rscript --vanilla scripts/snapshot_paper_artifacts.R
```

The snapshot mirrors each file's working path below `paper/generated/` so its
provenance remains clear. `manifest.csv` records each source path, snapshot
path, size and MD5 checksum. Review the snapshot diff together with the code
diff before committing revised results.

The tracked project `.Rprofile` activates the environment recorded in
`renv.lock` on a fresh checkout. The machine-specific package library under
`renv/library/` remains ignored by `renv/.gitignore`.
