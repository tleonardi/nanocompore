# Using the sample comparison command

Description

## Quick start

Examples command line

```
```

## Main options

`nanocompore sampcomp` requires several input files:

* **A list of paths to sample files** obtained with `NanopolishComp EventalignCollapse` as explained before (see [data preparation](data_preparation.md)), for both the control and the experimental conditions. Files can be provides though command line options or using a YAML file.

!!! info "Command line option"
	This option requires to provide 1 list of files per condition using `--file_list1` and `--file_list2` arguments as well as the labels for condition each conditions using `--label1` and `--label2`.

!!! info "YAML file option"
    This option allows to pass a YAML formatted file indicating the sample condition labels and paths to data files with the option `--sample_yaml`. The file should be formatted as follow:

    ```yaml
    Wild_Type:
        rep1:   path/to/sample1/rep1/data
        rep2:   path/to/sample1/rep2/data

    Knocked_out:
        rep1:   path/to/sample2/rep1/data
        rep2:   path/to/sample2/rep2/data
    ```

* The same **transcriptome FASTA reference file** used previously at the alignment step.

* A path to the **output directory** where to write the results and logs file.

* A **Bed file** containing the genome annotations corresponding to the transcriptome fasta file (Optional). If this file is given, Nanocompore will convert the transcript coordinates in the genome space....


## SampComp full help

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
  --bed BED	BED file with annotation of transcriptome used for	mapping (optional)

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
