# Reference inputs

All of these are configured under `inputs:` (and `container:`) in
`config/config.yaml`. Place them under `resources/` or point the config at
wherever they live. Large files are git-ignored by design.

| Config key | File (example) | What it is | How to obtain |
|---|---|---|---|
| `real_count_matrix` | `Bcell.NumReads.filtered.rowsum55.protein.txt` | Real single-cell isoform/transcript count matrix. Row names are pipe-delimited GENCODE transcript headers (`ENST...|ENSG...|...|length|...`); columns are cells. The empirical data scDesign2 is fit to. | Your own scRNA-seq isoform quantification (e.g. per-cell RSEM/Salmon on full-length data). |
| `mapped_events` | `Mapped.EventsToExons.txt` | rMATS event -> exon -> transcript map. Used only to restrict differential-splicing genes to rMATS-detectable ones. | Product of the separate rMATS/GrASE project; supplied here as a fixed input. |
| `gencode_transcripts_tab` | `gencode.v34.transcripts.tab` | 2-column tab file: transcript header, sequence. | `seqkit fx2tab gencode.v34.transcripts.fa > gencode.v34.transcripts.tab` |
| `qual_profile_r1` / `qual_profile_r2` | `HiSeq2500L125R1.txt` / `...R2.txt` | Illumina quality profiles for art_modern. | Ship with art_modern (`data/Illumina_profiles/`) or build from real FASTQ. |
| `rsem_reference_prefix` | `rsem_reference/Homo_sapiens.GRCh38.v34` | RSEM reference prefix used by `rsem-tbam2gbam`. | `rsem-prepare-reference --gtf gencode.v34.annotation.gtf genome.fa <prefix>` (see below). |
| `container.sif` | `art.sif` | Singularity image bundling art_modern, samtools, rsem, seqkit. | Build from `containers/art_modern.def`. |

## Building the RSEM reference

The original project shipped a prebuilt `rsem_reference/` but did not record the
command. Regenerate it once with (run inside the container):

```bash
singularity exec resources/art.sif rsem-prepare-reference \
    --gtf gencode.v34.annotation.gtf \
    --num-threads 8 \
    GRCh38.primary_assembly.genome.fa \
    resources/rsem_reference/Homo_sapiens.GRCh38.v34
```

Use the **same GENCODE release** (v34) as the transcript FASTA above.

## Consistency notes

- The transcript IDs in `real_count_matrix`, `gencode_transcripts_tab`, and the
  RSEM reference must all come from the same GENCODE release.
- `model.read_len` (effective-length calculation, default 50) is independent of
  `reads.read_len` (the simulated read length, default 125). Set each to match
  your intent.
- `model.biomart_host` hits a dated Ensembl archive. For reproducible reruns set
  `use_cached_biomart: true` after the first run to reuse the cached
  gene->chromosome table.
