# Data preparation

Before using Nanocompore, sequencing data have to be basecalled (Dorado), aligned to a transcriptome reference (Minimap2), and resquiggled. Nanocompore supports multiple different resquigglers: Nanopolish/f5c eventalign, Uncalled4 align or Remora's API.

We're planning to release a Nextflow pipeline that would automate all those steps. However, the manual steps for performing the preprocessing are described below. The guide is based on [https://doi.org/10.1002/cpz1.683](https://doi.org/10.1002/cpz1.683). Note that this procedure should be repeated for all samples.

### Reads basecalling

Firstly, raw pod5 reads have to be basecalled with a recent version of ONT basecaller.

Example with Dorado v0.8.0:

```bash
dorado basecaller --emit-moves {dorado_model} {pod5} > {basecalled_bam}
```

### Transcriptome alignment

Basecalled reads have to be aligned to a reference. For dRNA-Seq, reads should be aligned to a reference **transcriptome (not genome)** in a non-spliced fashion. For example, one can download reference transcriptome fasta files directly from Gencode for [human](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz) and [mouse](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz).

Bam files have to be filtered to remove any reads that would be unmapped, secondary and supplementary as well as reads mapped on the reverse strand (SAM flag 2324). We also recommend to discard reads with a very low alignment score (MAPQ<10). Finally, reads have then to be sorted and indexed.

There are subtle differences in the process depending on the choice of resquiggler. Specifically, Uncalled4 requires extracting additional tags from the bam when creating a fasta file with the reads and then adding the `-y` option to Minimap2 to pass them to the aligned bam. Remora additonally requires adding the `--MD` flag to Minimap2. Examples for each resquiggler are given below.

#### Example with [Minimap2 v2.28](https://github.com/lh3/minimap2) for f5c's eventalign

Extract a fasta file from the bam:

```bash
samtools fasta {basecalled_bam} > {reads_fasta}

minimap2 -ax map-ont -L {transcriptome_fasta} {reads_fasta} \
  | samtools view -bh -F 2324 -q 10 \
  | samtools sort -O bam > {aligned_reads_bam}

samtools index {aligned_reads_bam}
```

#### Example with [Minimap2 v2.28](https://github.com/lh3/minimap2) for Uncalled4's align

Extract a fasta file from the bam, adding the tags required by Uncalled4:

```bash
samtools fasta -T "mv,ts,pi,sp,ns" {basecalled_bam} > {reads_fasta}

# Note, we're adding the -y flag to pass the tags to the output bam
minimap2 -y -ax map-ont -L {transcriptome_fasta} {reads_fasta} \
  | samtools view -bh -F 2324 -q 10 \
  | samtools sort -O bam > {aligned_reads_bam}

samtools index {aligned_reads_bam}
```

#### Example with [Minimap2 v2.28](https://github.com/lh3/minimap2) for Remora

```bash
samtools fasta -T "mv,ts,pi,sp,ns" {basecalled_bam} > {reads_fasta}

# Note, we're adding the -y and --MD tags
minimap2 -y -ax map-ont -L {transcriptome_fasta} {reads_fasta} --MD \
  | samtools view -bh -F 2324 -q 10 \
  | samtools sort -O bam > {aligned_reads_bam}

samtools index {aligned_reads_bam}
```

### Resquiggling

Nanocompore supports using signal alignments from various sources. Currently supported ones are:
- Eventalin TSVs as produced by [Nanopolish](https://github.com/jts/nanopolish) or [f5c](https://github.com/hasindu2008/f5c). Since f5c provides optimized implementations of the eventalign command from Nanopolish, most users would be advised to use f5c.
- BAM files produced by [Uncalled4](https://github.com/skovaka/uncalled4).
- [Remora](https://github.com/nanoporetech/remora) is a toolkit by ONT for training models for modification prediction. However, it also provides functionality for signal analysis, including capabilities to refine the signal alignments from the "moves" table produced by Dorado.

#### f5c
Since f5c doesn't support pod5, but uses the community-driven `slow5/blow5` formats, first you'll need to convert the pod5 to blow5 using [blue-crab](https://github.com/Psy-Fer/blue-crab).

Make sure that you have blue-crab installed. Typically, this should involve running:
```bash
pip install blue-crab
```

Convert the pod5 to blow5:
```bash
blue-crab p2s -o {blow5} {pod5}
```
Next, index the reads:
```bash
f5c index --slow5 {blow5} {reads_fasta}
```

Finally, we can resquiggle the reads. Nanocompore doesn't use the raw TSV file produced by the eventalign command for the comparison of the two conditions directly, but first preprocesses the TSV to obtain a collapsed version (with one row per k-mer). Since the raw TSV files are huge, we recommend that you directly pipe the output of f5c's eventalign to Nanocompore's eventalign_collapse subocommand.
```bash
# For RNA004 data
f5c eventalign --bam {aligned_reads_bam} \
               --genome {transcriptome_fasta} \
               --reads {reads_fasta} \
               --slow5 {blow5} \
               --scale-events \
               --print-read-name \
               --secondary=yes \
               --min-mapq 0 \
               --samples \
               --summary {summary_file} \
               --rna \
               --pore rna004 \
               --threads {num_threads}
  | nanocompore eventalign_collapse --ref {transcriptome_fasta} -o {resquiggled_db} --threads {num_threads}
```
Note that the eventalign_collapse creates rather big temporary files. By default it would use the current working directory. If you want to change that, you can add `--tmp path/to/tmp_dir`.


#### Uncalled

```bash
# For RNA004
uncalled4 align --ref {transcriptome_fasta} \
                --reads {pod5} \
                --bam-in {aligned_reads_bam} \
                --bam-out {resquiggled_bam} \
                --kit SQK-RNA004 \
                -p {num_threads}
```
This will produce a bam file with the signal alignment data that you can directly pass to Nanocompore for the analysis.


#### Remora

[Remora](https://github.com/nanoporetech/remora) provides a Python API that can be used for resquiggling. Nanocompore provides a subcommand that will invoke Remora for aligning the signal. The command will produce an SQLite database with the signal data that you can then provide to Nanocompore for the analysis.

```bash
nanocompore remora_resquiggle --ref {transcriptome_fasta} \
                              --pod5 {pod5} \
                              --bam {aligned_reads_bam} \
                              --kit RNA004 \
                              --output {resquiggled_db}
```

