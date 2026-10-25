# goda-forbes-2025: RNA-seq quantification and alignment

<!-- version: 2 | date: 2026-10-25 -->

Snakemake workflow that downloads 12 paired-end RNA-seq runs (K562, 2x150,
NovaSeq X) from SRA, quantifies transcripts with salmon, and aligns to the
genome with STAR, on UNC's Longleaf cluster.

## Data

Two GEO series, 6 samples each, listed in `config/samples.tsv`
(built from the SRA Run Tables):

| GEO series | Comparison                 |
|------------|----------------------------|
| GSE279951  | Control vs PCBP1 knockdown |
| GSE279953  | Control vs Dox             |

Each series has its own `Control` group; use the `gse` column to tell them apart.

## What the workflow does

1. **Reference**: downloads GENCODE human release 50 (GRCh38 primary assembly
   genome, transcripts FASTA, primary-assembly GTF) with an MD5 check.
2. **SRA**: `prefetch` each run, then `fasterq-dump` + `pigz` to gzipped FASTQ.
   FASTQs are deleted once salmon and STAR are done with them
   (`keep_fastq: false`); the `.sra` files are kept in `sra/`.
3. **salmon** (2.8.0): decoy-aware index (transcripts + genome as decoys,
   `--gencode`), then selective-alignment quantification with
   `--gcBias --numGibbsSamples 30`.
4. **STAR** (2.7.11b): genome index with `sjdbOverhang 149`; coordinate-sorted
   BAMs with one read group per run, `--quantMode GeneCounts`, indexed with samtools.
5. **MultiQC**: one report over salmon and STAR output.

## Software

- Snakemake 9.27.1 + SLURM executor plugin: conda env at
  `/proj/milovelab/envs/snakemake`
- salmon 2.8.0: conda env at `/proj/milovelab/envs/salmon` (the prebuilt
  release binary needs glibc 2.35; Longleaf/RHEL 9 has 2.34)
- STAR, sra-tools, samtools, pigz, MultiQC: Longleaf modules

## Running

From a Longleaf login node, inside tmux (`tmux new -s rnaseq`):

```bash
module load star sratoolkit samtools pigz multiqc
module load anaconda
conda activate /proj/milovelab/envs/snakemake
snakemake -n                                         # dry run
snakemake --profile profiles/slurm --local-cores 4   # one SLURM job per rule
```

Downloads (GENCODE, SRA `prefetch`) run on the login node; everything else is
submitted to SLURM. Rerunning the same command resumes where it stopped.
`run_workflow.sh` is a fallback that runs everything inside a single SLURM
allocation, without the plugin.

Monitoring:

```bash
grep "steps" $(ls -t .snakemake/log/*.log | head -1) | tail -1   # x of 68 steps
squeue --me -o "%.10i %.3t %.9M %k"                               # jobs by rule
```

Per-job SLURM logs are in `.snakemake/slurm_logs/`, per-rule logs in `logs/`.

## Layout

```
Snakefile                  workflow rules
config/config.yaml         settings (GENCODE release, tool paths, salmon/STAR options)
config/samples.tsv         sample sheet with metadata (sample, run, layout, gse, treatment, ...)
profiles/slurm/            SLURM profile (partition, default resources)
run_workflow.sh            fallback single-allocation runner

not in git (see .gitignore):
ref/gencode.v50/           reference files and salmon/STAR indexes
sra/  fastq/  tmp/         downloads and scratch
salmon/<sample>/quant.sf   transcript quantifications (import with tximeta)
star/<sample>/             BAM + index, ReadsPerGene, SJ.out.tab, Log.final.out
multiqc/multiqc_report.html
logs/  .snakemake/
```
