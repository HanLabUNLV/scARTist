# Known issues (pre-release)

## Open
- Reference genomes / indexes (e.g. the RSEM reference) are user-supplied, not
  built by the workflow.
- LICENSE undecided.

## Verified
- Runs end-to-end: a full 80-cell run completed (exit status 0, 80 genome BAMs),
  plus a 1-cell smoke run, using the real containers `resources/art.sif` and
  `resources/r.sif`.
- `resources/r.sif` (built from `containers/r_stack.def`) self-tests OK - the
  vendored patched scDesign2 fork loads with `fold_change` present (required for
  DE/DTU injection) - and is wired into the R rules via `container:`.
- ID/sequence misalignment (earlier concern) fixed and validated: sequences are
  attached by a name-based keyed join; round-trip recoverability 0.05 -> 0.995.
- Bundled example (`resources/example/`, ~50 genes) runs end-to-end via
  `config/config.example.yaml` (13/13 steps, 8 genome BAMs) with a mini RSEM
  reference - a fast, self-contained test for new users.
