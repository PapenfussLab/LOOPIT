
LOOPIT detects extrachromosomal DNA (ecDNA) from long- and short-read whole-genome sequencing.

## Requirements

- Nextflow and Apptainer/Singularity.
- The LOOPIT container). Download it from the Sylabs Cloud Library at `library://uretajc/loopit/loopit_env:1.0`:

  ```bash
  apptainer pull loopit_env.sif library://uretajc/loopit/loopit_env:1.0
  ```

  Or build it yourself from the definition file:

  ```bash
  apptainer build loopit.sif loopit.def
  ```


## Samplesheet

CSV with a header row.

| Column | Required | Description |
|---|---|---|
| sample | yes | Sample name |
| mode | yes | `SR` (short read) or `LR` (long read) |
| tumour_bed | yes | Tumour mosdepth BED(.gz) |
| vcf | yes | SV VCF (tested on GRIDSS and Sniffles) |
| normal_bed | no | Normal mosdepth BED(.gz). Leave empty for tumour-only |
| cnv | no | External segment-level CN table (e.g. FACETS), used only for the plot track; not used for CN calling. Read with `fread`, so CSV or TSV (header required). Needs columns `chrom`, `start`, `end` and `cn`; a `tcn.em` column (FACETS total copy number) is renamed to `cn` automatically, so a FACETS segment table can be used directly. Chromosome names may be `chr7` or `7`; `X`/`23` and `Y`/`24` are recognised  |
| purity | no | Tumour purity (default 1.0) |
| ploidy | no | Tumour ploidy (default 2.0) |
| sex | no | `M` or `F`. Used only as a cross-check. Sex is inferred from chrY/autosome depth |



## Usage


```bash
nextflow run main.nf \
    -profile slurm,apptainer \
    --input assets/samplesheet.csv \
    --outdir results \
    --reference_dir /path/to/reference \
    -resume
```



## Parameters

| Parameter | Default | Description |
|---|---|---|
| min-cn | 2,4 | Minimum segment CN for graph nodes |
| min-width | 10000 | Minimum segment width (bp) |
| max-width | 1e7 | Maximum segment width (bp) |
| min-qual-bp | 50 | Minimum breakpoint quality |
| min-qual-vaf | mode-specific | Minimum SV VAF |
| min-sv-len | 300 | Minimum SV length (bp) |
| ecdna-threshold | 2.5 | Offset added to ploidy for the ecDNA CN cutoff |
| chromosomes | NA | Restrict to chromosomes (NA = all) |

## Outputs

Per sample, under `<outdir>/<sample>/`:

| File | Content |
|---|---|
| `<sample>_cycles.txt` | Segments in each reported cycle, with strand |
| `<sample>_decomposed_cycles.csv` | Per-cycle CN, length and ecDNA call |
| `<sample>_filtered_graph.rds` | Graph at the best filter combination |
| cycle plots | One per ecDNA cycle |

