![Nanocompore](pictures/Nanocompore_logo.png)

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

---

**Nanocompore identifies differences in ONT nanopore sequencing raw signal corresponding to RNA modifications by comparing 2 samples**

Nanocompore compares 2 ONT nanopore direct RNA sequencing datasets from different experimental conditions expected to have a significant impact on RNA modifications. It is recommended to have at least 2 replicates per condition. For example one can use a control condition with a significantly reduced number of modifications such as a cell line for which a modification writing enzyme was knocked-down or knocked-out. Alternatively, on a smaller scale transcripts of interests could be synthesized in-vitro.

---

# Installation

Ideally, before installation, create a clean **python3.5+** virtual environment to deploy the package. **Python 2 is not supported**. For example you can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):
```
mkvirtualenv nanocompore -p python3.6
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):
```
conda create -n nanocompore python=3.6
```
## Dependencies

Nanocompore relies on a the following robustly maintained third party python libraries:
* numpy >= 1.14.0
* scipy >= 1.1.0
* tqdm >= 4.23.4
* pyfaidx >= 0.5.4.1
* matplotlib >= 2.2.2
* seaborn >= 0.9.0
* pandas >= 0.23.3
* statsmodels >= 0.9.0
* scikit-learn >= 0.20
* bedparse >= 0.1.2

The correct versions of packages are installed together with the software when using pip.

## Option 1: Direct installation with pip from github (recommended)

* To install the package with an ssh key
```
pip3 install git+ssh://git@github.com/tleonardi/nanocompore.git
```

* To install the package with https/ssh
```
pip3 install git+https://github.com/tleonardi/nanocompore.git
```

## Option 2: Clone the repository and install locally in develop mode

With this option, the package will be locally installed in “editable” or “develop” mode. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to modify the code and/or participate to the development of the package (see [contribution guidelines](https://github.com/tleonardi/nanocompore/blob/master/CONTRIBUTING.md)).

`git clone https://github.com/tleonardi/nanocompore.git` or bleeding edge `git clone --branch devel https://github.com/tleonardi/nanocompore.git`

`cd nanocompore`

[`chmod u+x setup.py`]

`pip3 install -e ./`

# Data preparation

Before using nanocompore, sequencing data have to be basecalled (Albacore or Guppy), aligned on a transcriptome reference and resquiggled with Nanopolish.

To simplify the data preprocessing we wrote a Nextflow pipeline that automatise all these steps as well as extra quality control steps: https://github.com/tleonardi/nanocompore_pipeline

### Reads basecalling

Firstly, raw fast5 reads have to be basecalled with a recent version of ONT basecaller. Basecalled fast5 files are not required for the rest of the analysis, only the raw fast5 and the basecalled fastq.

Example with [Guppy v2.3.5](https://community.nanoporetech.com/downloads)
```
guppy_basecaller -i {raw_fast5_dir} -s {dest_dir} --flowcell {flowcell_id} --kit {Kit_id} -r --calib_detect --enable_trimming true --trim_strategy rna --reverse_sequence true
```
Then the output fastq files should be concatenated in a single one.
```
cat {dir_to guppy output}/*.fastq > {dir_to guppy output}/reads.fastq
```

### Transcriptome alignment

Basecalled reads have to be aligned to a reference. For dRNA-Seq, reads should be aligned to a reference **transcriptome (not genome)** in a non-spliced fashion. For example, one can download reference transcriptome fasta files directly from Gencode for [human](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz) and [mouse](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.transcripts.fa.gz).

 Bam files have to be filtered to remove any reads that would be unmapped, secondary and supplementary as well as reads mapped on the reverse strand (SAM flag 2324). We also recommend to discard reads with a very low alignment score (MAPQ<10). Reads have then to be sorted and indexed.

Example with [Minimap2 v2.16](https://github.com/lh3/minimap2)
```
minimap2 -ax map-ont -L {transcriptome_fasta} {reads_fastq} | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > {reads_bam}
samtools index minimap.filt.sort.bam minimap.filt.sort.bam.bai
```

### Read indexing and resquiggling with

Nanopolish is required to realign raw signal to the expected reference sequence. Reads have to be indexed first with nanopolish index, realigned with nanopolish eventalign and finally the data has to be collapsed per kmer and indexed by NanopolishComp Eventalign_collapse.

Example with [Nanopolish v0.10.1](https://github.com/jts/nanopolish) and NanopolishComp v0.4.3

```
nanopolish index -s ${albacore_results}/sequencing_summary.txt -d 'raw_data' ${albacore_results}/workspace/*.fastq
nanopolish eventalign -t ${cpus_each} --reads ${albacore_results}/workspace/*.fastq --bam ${bam_file} --genome ${transcriptome_fasta} --samples --print-read-names --scale-events

nanopolish index -d ./raw/ ./basecall/workspace/reads.fastq

nanopolish eventalign --reads ./basecall/workspace/reads.fastq --bam ./alignment/reads.bam --genome ./reference/transcriptome.fa --samples > ./eventalign/reads.tsv

NanopolishComp Eventalign_collapse -i ./eventalign/reads.tsv -o ./eventalign/reads_collapsed.tsv

```

### Kmer level event align data collapsing with [NanopolishComp Eventalign_collapse](https://github.com/a-slide/NanopolishComp) v

```
NanopolishComp Eventalign_collapse -t ${cpus_each} -o reads_collapsed.tsv
```

# Nanocompore usage

####################################### TO DO

## Alternative and complementary packages

Here is a non-exhaustive list of alternative/complementary packages that could also be used to analyse RNA methylation from nanopore sequencing datasets:

* [Tombo](https://github.com/nanoporetech/tombo)
* ...

## Reporting bugs, feature request and contributing to development

Thanks for considering contributing to our_package!

Please see [contribution guidelines](https://github.com/tleonardi/nanocompore/blob/master/CONTRIBUTING.md) as well as the [code of conduct](https://github.com/tleonardi/nanocompore/blob/master/CODE_OF_CONDUCT.md) for more information.


## Authors

* Adrien Leger - aleg {at} ebi.ac.uk
* Tommaso Leonardi - tom {at} tleo.io

## Acknowledgements
