# Compare samples with SampComp

## Quick start

Examples command line

```bash
nanocompore sampcomp \
    --file_list1 ./data/S1_R1.tsv,./data/S1_R2.tsv \
    --file_list2 ./data/S2_R1.tsv,./data/S2_R2.tsv \
    --label1 S1 \
    --label2 S2 \
    --fasta ./reference/ref.fa \
    --outpath ./results/ \
```

Example API call

```python3
# Import package
from nanocompore.SampComp import SampComp

# Init the object
s = SampComp(
    eventalign_fn_dict = {
        'S1':{'rep1':'./data/S1_R1.tsv', 'rep2':'./data/S1_R2.tsv'},
        'S2':{'rep1':'./data/S2_R1.tsv', 'rep2':'./data/S2_R2.tsv'}},
    outpath = "./results/",
    fasta_fn = "./reference/ref.fa")

# Run the analysis
s()
```

## Description of main options

`SampComp` provides a very flexible analysis framework with a few mandatory options and many optional parameters. The full documentation command line documentation is provided at the bottom of this page or can be obtained via the help option `nanocompore sampcomp --help`. Immediately below is the detailed explanation of the main options:

#### Sample files

`SampComp` requires sample files obtained with `NanopolishComp EventalignCollapse` as explained before (see [data preparation](data_preparation.md)) for both the control and the experimental conditions. 2 conditions are expected, and at least 2 sample replicates per conditions are highly recommended. If `SampComp` is called through the CLI the files can be provides using either relevant command options or a YAML file. If using the Python API, one can pass either a python dictionary or a YAML file.

!!! info "YAML file option (CLI or API)"
    This option allows to pass a YAML formatted file indicating the sample condition labels and paths to data files with the option `--sample_yaml` for the CLI or directly to `eventalign_fn_dict` for the API. The file should be formatted as follow:

    ```yaml
    WT:
        rep1:   path/to/sample1/rep1/data
        rep2:   path/to/sample1/rep2/data

    KO:
        rep1:   path/to/sample2/rep1/data
        rep2:   path/to/sample2/rep2/data
    ```

!!! info "Command line option (CLI only)"
	This option requires to provide 1 comma separated list of files per condition using `--file_list1` and `--file_list2` arguments as well as the labels for condition each conditions using `--label1` and `--label2`.

!!! info "Python dictionary (API only)"
    This option allows to pass a multi-level python dictionary containing the sample condition labels and paths to data files. The dictionary should be formatted as follow:

    ```python
    eventalign_fn_dict = {
        "WT":  {"rep1":"path/to/sample1/rep1/data", "rep2":"path/to/sample1/rep2/data"},
        "KO": {"rep1":"path/to/sample2/rep1/data", "rep2":"path/to/sample2/rep2/data"}
        }
    ```
#### Transcriptome reference FASTA file

A transcriptome FASTA reference file is required to extract the corresponding kmer sequence when ambigous or missing in the eventalign file. The reference has to be the same as the one used at the mapping step. (CLI: `--fasta`, API: `fasta_fn`)

#### Output database path / directory

For the CLI, one has to provide a path to a directory where the program will output the result DBM database as well as additional results files and logs (`--outpath`). For the API, only the path to the file where to write the DBM database is required (`output_db_fn`).

#### Genome annotation BED file

Optionally, a BED file containing the genome annotations corresponding to the transcriptome fasta file can be provided. If this file is given, Nanocompore will also convert the transcript coordinates into the genome space (CLI: `--bed`, API: `bed_fn`)

#### Statistical options

`SampComp` implements several statistical methods to evaluate the difference between the 2 conditions (`comparison_method`).

* Gaussian Mixture Model = GMM (default)
* Kolmogorov–Smirnov test = KS (default)
* Mann–Whitney U test = MW
* T-test = TT

In addition, it is also possible to specify the number of adjacent positions to take into account for the pvalue calculation (`sequence_context`) as well as the weights to give to adjacent position, using either an "uniform" or a "harmonic" distribution (`sequence_context_weights`).

#### Coverage options

The default coverage threshold for `SampComp` to perform a statistical test is 50 reads in each replicates. This is quite conservative and can be modified if needed (`min_coverage`). In addition, to reduce the computational burden it is possible to randomly down-sample the number of reads for high coverage references (`downsample_high_coverage`).

#### Manually exclude or include references (API only)

The API allows to specify references to be included or excluded from the analysis (`select_ref_id` and
`exclude_ref_id`). This can be useful to analyse a specific transcripts set only or to run a small test before analysing the entire dataset.

## SampComp full command line help

```text
usage: nanocompore sampcomp [-h]
	--file_list1 /path/to/Condition1_rep1,/path/to/Codition1_rep2
	--file_list2 /path/to/Condition2_rep1,/path/to/Codition2_rep2
    --label1 Condition1
    --label2 Condition2
    --fasta FASTA
    --outpath OUTPATH
    [--bed BED]
    [--overwrite]
	[--force_logit]
	[--max_invalid_kmers_freq MAX_INVALID_KMERS_FREQ]
	[--min_coverage MIN_COVERAGE]
	[--downsample_high_coverage DOWNSAMPLE_HIGH_COVERAGE]
	[--comparison_methods COMPARISON_METHODS]
	[--sequence_context {0,1,2,3,4}]
	[--sequence_context_weights {uniform,harmonic}]
	[--pvalue_thr PVALUE_THR] [--nthreads NTHREADS]
	[--log_level {warning,info,debug}]

Compare 2 samples and find significant signal

optional arguments:
  -h, --help show this help message and exit

Input options:
  --file_list1 /path/to/Condition1_rep1,/path/to/Codition1_rep2 Comma separated list of NanopolishComp files for label 1 (required)
  --file_list2 /path/to/Condition2_rep1,/path/to/Codition2_rep2	Comma separated list of NanopolishComp files for label 2 (required)
  --label1 Condition1	Label for files in --file_list1 (default: Condition1)
  --label2 Condition2   Label for files in --file_list2 (default: Condition2)
  --fasta FASTA, -f FASTA Fasta file used for mapping (required)
  --bed BED	BED file with annotation of transcriptome used for mapping (optional)

Output options:
  --outpath OUTPATH, -o OUTPATH	Path to the output folder (required)
  --overwrite	Use --outpath even if it exists already (default:False)
  --force_logit	Use logistic regression testing even if all conditions have replicates (default: False)

Transcript filtering options:
  --max_invalid_kmers_freq MAX_INVALID_KMERS_FREQ	Max fequency of invalid kmers (default: 0.1)
  --min_coverage MIN_COVERAGE	Minimum coverage required in each condition to do the comparison (default: 50)
  --downsample_high_coverage DOWNSAMPLE_HIGH_COVERAGE	Used for debug: transcripts with high covergage will be downsampled (default: None)

Statistical testing options:
  --comparison_methods COMPARISON_METHODS	Comma separated list of comparison methods. Valid methods are: GMM,KS,TT,MW. (default: GMM,KS)
  --sequence_context {0,1,2,3,4}	Sequence context for combining p-values (default: 2)
  --sequence_context_weights {uniform,harmonic}	Type of weights to use for combining p-values
  --pvalue_thr PVALUE_THR	Adjusted p-value threshold for reporting significant sites (default: 0.05)

Other options:
  --nthreads NTHREADS, -n NTHREADS	Number of threads (default: 3)
  --log_level {warning,info,debug}	log level (default: info)
```
