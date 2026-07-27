# scARTist bundled example dataset

A tiny, self-contained subset of the real B-cell inputs (**50 genes, 233
transcripts**, chosen to be multi-isoform and rMATS-detectable so the DTU
injection has candidates) for a fast end-to-end test of the whole pipeline:
scDesign2 fit -> count simulation -> art_modern reads -> genome BAM.

## Run it

From the repo root, with `art.sif` and `r.sif` present in `resources/`:

```bash
snakemake --cores 8 --use-singularity --singularity-args "--bind $PWD" \
    --configfile config/config.example.yaml -p
```

Produces `results/reads/rd_1000000_cells_4/{null,deds}/Bcell_genome.{1..4}.bam`
(4 cells x 2 arms, one shallow read depth). Paths in the config are relative to
the repo root, so this works wherever the repo is cloned.

## Files (all runtime inputs point here via `config/config.example.yaml`)

| File | Role |
|---|---|
| `Bcell.NumReads.example.txt` | subset real count matrix (233 transcripts x 80 cells) scDesign2 is fit to |
| `gencode.v34.transcripts.example.tab` | transcript sequences for those transcripts |
| `Mapped.EventsToExons.example.txt` | rMATS event map, subset to the example genes |
| `HiSeq2500L125R1.txt` / `...R2.txt` | Illumina quality profiles (bundled as-is) |
| `rsem_ref/example.*` | **mini RSEM reference** for those transcripts (transcript->genome coords, for `rsem-tbam2gbam`) |
| `example.gtf` | subset GTF the mini RSEM reference was built from (provenance / rebuild) |
| `selected_genes.txt`, `selected_tx.txt` | the selected gene / transcript IDs |

## Provenance

Subset from the real inputs under the data dir; the mini RSEM reference was built
with `rsem-prepare-reference --gtf example.gtf <GRCh38.primary_assembly.genome.fa>`
(GENCODE v34). The genome FASTA is a build-time input only and is **not** bundled;
only the small mini-reference is.
