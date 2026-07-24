# Example / smoke-test data (TODO)

Bundle a tiny synthetic input here so the workflow can be dry-run and
smoke-tested without the full references:

- a small `real_count_matrix` (a few hundred transcripts across ~10 cells, with
  valid pipe-delimited GENCODE headers)
- a matching `gencode_transcripts_tab` (sequences for those transcripts)
- a tiny `mapped_events` file
- a small RSEM reference built from the same handful of transcripts

Then add an `config/config.example.yaml` pointing here and wire it into CI so
`snakemake -n` and a `counts_only` run are exercised on every push.
