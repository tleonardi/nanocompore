# Data preparation

Before using nanocompore, sequencing data have to be basecalled (Albacore or Guppy), aligned on a transcriptome reference and resquiggled with Nanopolish.

To simplify the data preprocessing we wrote a Nextflow pipeline which automatises all these steps as well as extra quality control steps: https://github.com/tleonardi/nanocompore_pipeline

### Reads basecalling

Firstly, raw fast5 reads have to be basecalled with a recent version of ONT basecaller. Basecalled fast5 files are not required for the rest of the analysis, only the raw fast5 and the basecalled fastq.

Example with [Guppy v2.3.5](https://community.nanoporetech.com/downloads)

```bash
guppy_basecaller -i {raw_fast5_dir} -s {dest_dir} --flowcell {flowcell_id} --kit {Kit_id} -r --calib_detect --enable_trimming true --trim_strategy rna --reverse_sequence true
```

Then the output fastq files should be concatenated in a single one.

```bash
cat {dir_to guppy output}/*.fastq > {basecalled_fastq}
```

### Transcriptome alignment

Basecalled reads have to be aligned to a reference. For dRNA-Seq, reads should be aligned to a reference **transcriptome (not genome)** in a non-spliced fashion. For example, one can download reference transcriptome fasta files directly from Gencode for [human](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz) and [mouse](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.transcripts.fa.gz).

 Bam files have to be filtered to remove any reads that would be unmapped, secondary and supplementary as well as reads mapped on the reverse strand (SAM flag 2324). We also recommend to discard reads with a very low alignment score (MAPQ<10). Finally, reads have then to be sorted and indexed.

Example with [Minimap2 v2.16](https://github.com/lh3/minimap2)

```bash
minimap2 -ax map-ont -L {transcriptome_fasta} {basecalled_fastq} | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > {aligned_reads_bam}

samtools index {aligned_reads_bam}
```

### Read indexing and resquiggling with

Nanopolish is required to realign raw signal to the expected reference sequence. Reads have to be indexed first with nanopolish index, realigned with nanopolish eventalign and finally the data has to be collapsed per kmer and indexed by NanopolishComp Eventalign_collapse.

Example with [Nanopolish v0.10.1](https://github.com/jts/nanopolish) and [NanopolishComp v0.4.3](https://github.com/a-slide/NanopolishComp)

```bash
nanopolish index -s {sequencing_summary.txt} -d {raw_fast5_dir} {basecalled_fastq}

nanopolish eventalign --reads {basecalled_fastq} --bam {aligned_reads_bam} --genome {transcriptome_fasta} --samples --print-read-names --scale-events --samples > {eventalign_reads_tsv}

NanopolishComp Eventalign_collapse -i {eventalign_reads_tsv} -o {eventalign_collapsed_reads_tsv}

```
