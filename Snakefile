# =============================================================================
# Snakefile: SRA download -> salmon (decoy-aware selective alignment)
#                         -> STAR (genome BAM, coordinate sorted)
# Reference: GENCODE human (primary assembly)
#
# version: 3
# date:    2026-10-25
# changes: v3 MultiQC over salmon + STAR; STAR Log.final.out declared
#          v2 optional temp() FASTQs, longer runtimes for ~30-50 Gb runs
# =============================================================================

import pandas as pd

configfile: "config/config.yaml"

T = config["tools"]

# ---- sample sheet: one row per SRA run; runs are grouped by `sample` ---------
runs = (
    pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    .assign(layout=lambda d: d["layout"].str.upper())
    .set_index("run", drop=False)
)
assert runs["layout"].isin(["PE", "SE"]).all(), "layout must be PE or SE"
assert runs["run"].is_unique, "run accessions must be unique"
assert runs.groupby("sample")["layout"].nunique().eq(1).all(), \
    "all runs within a sample must share the same layout"

SAMPLES = sorted(runs["sample"].unique())


def sample_runs(s):
    return runs.loc[runs["sample"] == s, "run"].tolist()


def sample_layout(s):
    return runs.loc[runs["sample"] == s, "layout"].iloc[0]


def fastqs(s):
    """Named lists of FASTQ paths for a sample: r1 (and r2 if paired)."""
    rr = sample_runs(s)
    if sample_layout(s) == "PE":
        return {"r1": [f"fastq/pe/{r}_1.fastq.gz" for r in rr],
                "r2": [f"fastq/pe/{r}_2.fastq.gz" for r in rr]}
    return {"r1": [f"fastq/se/{r}.fastq.gz" for r in rr]}


# ---- reference paths ---------------------------------------------------------
REL = str(config["gencode"]["release"])
GC_URL = f"{config['gencode']['base_url']}/release_{REL}"
REF = f"ref/gencode.v{REL}"
GENOME = f"{REF}/GRCh38.primary_assembly.genome.fa"
TXOME = f"{REF}/gencode.v{REL}.transcripts.fa"
GTF = f"{REF}/gencode.v{REL}.primary_assembly.annotation.gtf"
OH = str(config["star"]["sjdb_overhang"])
SALMON_IDX = f"{REF}/salmon_index_decoy"
STAR_IDX = f"{REF}/star_index_oh{OH}"
TMP = config["tmpdir"]

# FASTQs are deleted once salmon and STAR both finish unless keep_fastq: true
# (the .sra files in sra/ are kept, so re-dumping needs no re-download)
fq_out = (lambda p: p) if config.get("keep_fastq", False) else temp

wildcard_constraints:
    run=r"[SDE]RR\d+",
    sample=r"[^/]+",

# network steps run on the host running snakemake (compute nodes often lack
# internet); set local_downloads: false if your compute nodes can reach out
if config.get("local_downloads", True):
    localrules: all, download_gencode, prefetch
else:
    localrules: all


rule all:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        expand("star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
               sample=SAMPLES),
        "multiqc/multiqc_report.html",


# =============================================================================
# Reference
# =============================================================================

rule download_gencode:
    output:
        genome=GENOME + ".gz",
        tx=TXOME + ".gz",
        gtf=GTF + ".gz",
    params:
        url=GC_URL,
        dir=REF,
    log:
        f"logs/download_gencode.v{REL}.log",
    shell:
        r"""
        for f in {output.genome} {output.tx} {output.gtf}; do
            curl -fsSL --retry 3 -o $f {params.url}/$(basename $f)
        done
        curl -fsSL --retry 3 -o {params.dir}/MD5SUMS {params.url}/MD5SUMS
        ( cd {params.dir} && \
          grep -F -e $(basename {output.genome}) -e $(basename {output.tx}) \
                  -e $(basename {output.gtf}) MD5SUMS | md5sum -c - ) > {log} 2>&1
        """


rule gunzip_ref:
    input:
        "{path}.gz",
    output:
        temp("{path}"),
    wildcard_constraints:
        path=r"ref/.*\.(fa|gtf)",
    shell:
        "gzip -dc {input} > {output}"


rule salmon_gentrome:
    input:
        tx=TXOME + ".gz",
        genome=GENOME + ".gz",
    output:
        gentrome=temp(f"{REF}/salmon_gentrome.fa.gz"),
        decoys=f"{REF}/salmon_decoys.txt",
    shell:
        r"""
        gzip -dc {input.genome} | grep '^>' | cut -d ' ' -f 1 | sed 's/^>//' \
            > {output.decoys}
        # decoys (genome) must come after transcripts
        cat {input.tx} {input.genome} > {output.gentrome}
        """


rule salmon_index:
    input:
        gentrome=f"{REF}/salmon_gentrome.fa.gz",
        decoys=f"{REF}/salmon_decoys.txt",
    output:
        directory(SALMON_IDX),
    params:
        extra=config["salmon"]["index_extra"],
    threads: 16
    resources:
        mem_mb=64000,
        runtime=480,
    log:
        f"logs/salmon_index.v{REL}.log",
    shell:
        r"""
        {T[salmon]} index -t {input.gentrome} -d {input.decoys} \
            -i {output} -p {threads} {params.extra} > {log} 2>&1
        """


rule star_index:
    input:
        fa=GENOME,
        gtf=GTF,
    output:
        directory(STAR_IDX),
    params:
        oh=OH,
        tmp=TMP,
    threads: 16
    resources:
        mem_mb=64000,
        runtime=360,
    log:
        f"logs/star_index.v{REL}.log",
    shell:
        r"""
        mkdir -p {output} {params.tmp}
        {T[STAR]} --runMode genomeGenerate --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.oh} \
            --outFileNamePrefix {output}/ \
            --outTmpDir {params.tmp}/star_genome_$$ > {log} 2>&1
        """


# =============================================================================
# SRA download
# =============================================================================

rule prefetch:
    output:
        directory("sra/{run}"),
    params:
        max_size=config["sra"]["max_size"],
    log:
        "logs/prefetch/{run}.log",
    shell:
        r"""
        {T[prefetch]} {wildcards.run} --max-size {params.max_size} \
            -O sra > {log} 2>&1
        """


rule fasterq_dump_pe:
    input:
        "sra/{run}",
    output:
        r1=fq_out("fastq/pe/{run}_1.fastq.gz"),
        r2=fq_out("fastq/pe/{run}_2.fastq.gz"),
    params:
        tmp=TMP,
    threads: 8
    resources:
        mem_mb=16000,
        runtime=720,
    log:
        "logs/fasterq_dump/{run}.log",
    shell:
        r"""
        mkdir -p {params.tmp}
        out=$(mktemp -d {params.tmp}/fqd_{wildcards.run}_XXXX)
        {T[fasterq_dump]} {input} --split-3 --threads {threads} \
            --temp $out --outdir $out > {log} 2>&1
        {T[pigz]} -p {threads} -c $out/{wildcards.run}_1.fastq > {output.r1}
        {T[pigz]} -p {threads} -c $out/{wildcards.run}_2.fastq > {output.r2}
        rm -rf $out
        """


rule fasterq_dump_se:
    input:
        "sra/{run}",
    output:
        fq_out("fastq/se/{run}.fastq.gz"),
    params:
        tmp=TMP,
    threads: 8
    resources:
        mem_mb=16000,
        runtime=720,
    log:
        "logs/fasterq_dump/{run}.log",
    shell:
        r"""
        mkdir -p {params.tmp}
        out=$(mktemp -d {params.tmp}/fqd_{wildcards.run}_XXXX)
        {T[fasterq_dump]} {input} --split-3 --threads {threads} \
            --temp $out --outdir $out > {log} 2>&1
        {T[pigz]} -p {threads} -c $out/{wildcards.run}.fastq > {output}
        rm -rf $out
        """


# =============================================================================
# salmon
# =============================================================================

def salmon_reads(wc):
    fq = fastqs(wc.sample)
    if "r2" in fq:
        return f"-1 {' '.join(fq['r1'])} -2 {' '.join(fq['r2'])}"
    return f"-r {' '.join(fq['r1'])}"


rule salmon_quant:
    input:
        unpack(lambda wc: fastqs(wc.sample)),
        index=SALMON_IDX,
    output:
        "salmon/{sample}/quant.sf",
    params:
        reads=salmon_reads,
        outdir="salmon/{sample}",
        libtype=config["salmon"]["libtype"],
        extra=config["salmon"]["quant_extra"],
    threads: 8
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        "logs/salmon_quant/{sample}.log",
    shell:
        r"""
        {T[salmon]} quant -i {input.index} -l {params.libtype} \
            {params.reads} -p {threads} {params.extra} \
            -o {params.outdir} > {log} 2>&1
        """


# =============================================================================
# STAR
# =============================================================================

def star_reads(wc):
    fq = fastqs(wc.sample)
    s = ",".join(fq["r1"])
    if "r2" in fq:
        s += " " + ",".join(fq["r2"])
    return s


def star_rg(wc):
    # one read group per SRA run, separated by " , " as STAR expects
    return " , ".join(f"ID:{r} SM:{wc.sample}" for r in sample_runs(wc.sample))


rule star_align:
    input:
        unpack(lambda wc: fastqs(wc.sample)),
        index=STAR_IDX,
    output:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        counts="star/{sample}/{sample}.ReadsPerGene.out.tab",
        sj="star/{sample}/{sample}.SJ.out.tab",
        log_final="star/{sample}/{sample}.Log.final.out",
    params:
        reads=star_reads,
        rg=star_rg,
        prefix="star/{sample}/{sample}.",
        tmp=TMP,
        sort_ram=config["star"]["limit_bam_sort_ram"],
        extra=config["star"]["extra"],
    threads: 16
    resources:
        mem_mb=48000,
        runtime=600,
    log:
        "logs/star_align/{sample}.log",
    shell:
        r"""
        mkdir -p {params.tmp}
        {T[STAR]} --runThreadN {threads} --genomeDir {input.index} \
            --readFilesIn {params.reads} --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmp}/star_{wildcards.sample}_$$ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattrRGline {params.rg} \
            --outSAMunmapped Within \
            --quantMode GeneCounts \
            --limitBAMsortRAM {params.sort_ram} \
            {params.extra} > {log} 2>&1
        """


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    threads: 4
    shell:
        "{T[samtools]} index -@ {threads} {input}"


# =============================================================================
# MultiQC: salmon (mapping rate, library type, fragment lengths) and
# STAR (alignment stats, ReadsPerGene strandedness)
# =============================================================================

rule multiqc:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        expand("star/{sample}/{sample}.Log.final.out", sample=SAMPLES),
        expand("star/{sample}/{sample}.ReadsPerGene.out.tab", sample=SAMPLES),
    output:
        "multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        r"""
        {T[multiqc]} --force -o multiqc -n multiqc_report.html \
            salmon star > {log} 2>&1
        """
