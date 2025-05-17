# How to use Nanocompore

Once all samples have been preprocessed as explained in the [Data preparation]( e/data_preparation) page, we can proceed to perform the comparison of the two conditions with Nanocompore.

### Creating a configuration file

Nanocompore uses a YAML configuration file that contains all parameters for a given analysis. We can generate a new configuration using the **template** subcommand:

```bash
nanocompore template analysis.yaml
```

This will create a file called "analysis.yaml" in the current working directory. The file lists all parameters for Nanocompore documented with comments explaining their function and possible values. Here's the top part of the configuration file which describes the required fields:

```yaml
# === REQUIRED FIELDS ===

# Specify all the conditions and samples to bo analysed.
# Depending on the resquiggler used, Nanocompore expect different
# fields to be defined for all samples:
# - Uncalled4:
#   - Each sample should have the "bam" field set with a path
#     to a bam file containing the Uncalled4 resquiggling data.
# - Eventalign or Remora:
#   - Each sample should have the "db" field set with a path
#     to an SQLite database file produced by the subcommands
#     "eventalign_collapse" or "remora_resquiggle".
data:
  # Only two conditions are supported. There's
  # no limit to the number of samples per condition.
  # You can change the labels appropriately, but
  # note that all labels should be unique.
  condition1:
    sample1:
      bam: "/path/to/uncalled4.bam"
      db: "/path/to/collapsed_eventalign.sqlite"
    sample2:
      bam: "/path/to/uncalled4.bam"
      db: "/path/to/collapsed_eventalign.sqlite"
    # you can list as many samples as you want
  condition2:
    sample3:
      bam: "/path/to/uncalled4.bam"
      db: "/path/to/collapsed_eventalign.sqlite"
    sample4:
      bam: "/path/to/uncalled4.bam"
      db: "/path/to/collapsed_eventalign.sqlite"

# Which of the two conditions is the one that's
# depleted of modifications.
depleted_condition: condition1

# Reference genome/transcriptome in fasta format
fasta: "/path/to/reference.fa"

# The kit that was used for sequencing.
# Possible values: RNA002, RNA004
kit: "RNA004"

# Specify which resquiggler you have used for your data.
# Possible options are "eventalign", "uncalled4", and "remora".
resquiggler: 'eventalign'
```
The `data` field should list all samples, preprocessed with the same resquiggler. The `fasta` parameter should be set to the transcriptome reference (a FASTA file) that was used for the alignment of the reads.

Below, the file will list various optional parameters that the user can tweak. The most important ones are discussed below.

#### Output configuration

Nanocompore stores all results obtaned during the analysis in an SQLite database. After the analysis is done it will perform a postprocessing step in which it will do multiple-test correction (by default using the Benjamini-Hochberg method) and export the test results to an easy to read TSV file. One can request additional information to be exported using the configuration. The most important parameters are:

##### Map the results to genomic coordinates

The results produced by Nanocompore will contain only transcriptomic coordinates by default (i.e. the transcript reference name and the position along the transcript). Nanocompore can optionally perform the mapping to genomic coordinates. To do this, simply add the following to the configuration:

```yaml
gtf: 'path/to/annotation.gtf'
```

Where the annotation should be compatible with the reference FASTA that is used, i.e. it should provide information for all references in the FASTA.

As a result, the columns `chr`, `strand`, and `genomicPos` in the results TSV file will be populated. Note that the `genomicPos` would use a 0-based indexing (the first base on the chromosome has index 0).

##### Shift statistics

The shift statistics TSV gives summary statistics (mean, median, standard deviation) at the position level for the signal measurements (current intensity and dwell time) for the two conditions (see [Outputs](/output)). To enable the exporting of this file just add the folling to the configuration:

```yaml
export_shift_stats: True
```

#### Performance tuning

The main way to control the performance and resource usage of Nanocompore is through the `devices` parameter.

```yaml
# Control how many worker processes to run and which
# computing device they should run on.
# The parameter should be a map of devices and number of worker processes.
# The accepted device values are: "cpu", "cuda", "cuda:N"
# For example, to use 4 worker processes on the CPU and
# 8 worker processes on each of the two available GPUs
# one can set the following:
devices:
  "cpu": 4
  "cuda:0": 8
  "cuda:1": 8
```
The `devices` parameter specifies how many worker processes Nanocompore should run. Typically, the more parallelization we use, the faster the analysis will go, as long as this matches the computational resources available. However, increasing the number of workers means that more memory would be used, potentially leading to out-of-memory errors. Some guidelines to consider:

- Running the analysis on GPU devices is significantly faster than using CPUs.
- If using CPUs, a reasonable starting point is setting the parameter to `"cpu": N-1` where N is the number of available CPU cores. This would allow each worker to utilize one core while leaving one core for the main orchestrator process.
- If using GPUs, the user is advised to try experimenting with different number of workers to tune the performance according to their GPU model. In particular, one should increase the worker count to achieve high GPU utilization, while making sure that the GPU memory consumption is low enough to avoid out-of-memory errors. Each worker process would have a copy of the CUDA context, which takes a few hundred MB of memory. On top of that, the worker will load to the GPU memory the data for the transcript that it currently is processing. Monitoring the run logs can show you if out of memory errors occur. In the case of such error, Nanocompore will add the transcript back to the queue to retry it later (up to 3 times), so an occasional out-of-memory is not a reason to be very worried. However, if they are too frequent it means that the GPU is overloaded which would lead to lower performance and potential loss of results.
- For example, for a GPU model NVIDIA Tesla P100-PCIE, which has 12GB of memory we typically use 8 workers per GPU device.
- You can mix devices to improve the performance. However, bear in mind that GPU workers would also use the CPU to some degree. Hence, it's not recommended to use all CPU cores with CPU workers if you also run GPU workers. E.g. if you have 32 cores and you run 16 GPU workers, you can try running 15 CPU workers on top, but it's not wise to try running 32 CPU workers.

#### Continue a previous run

By default, if the output directory for the run already exists, Nanocompore will stop and output an error message. You can instruct Nanocompore how to proceed in such cases by changing the parameter `result_exists_strategy`: 

```yaml
# What to do if the output directory already exists.
# Possible options are:
# - stop (default): stop the execution with an error.
# - overwrite: delete the existing directory and start from scratch.
# - continue: try to continue from where the previous execution left off.
result_exists_strategy: "continue"
```
If you set the parameter to "continue", Nanocompore will try to reuse the database that was created in the previous run and won't repeat the analysis for any transcripts that have already been analysed and written to the database.

#### Other parameters

For a full list of parameters and a description of their function, inspect a template configuration created by the `nanocompore template` command.

### Run the analysis

Once you have created a configuration, you can run the analysis on your data with:

```bash
nanocompore run analysis.yaml
```

Nanocompore should immediately create a directory "nanocompore_output" in the current working directory (unless you have set the `outpath` parameter in the configuration file to a different location). All working files, logs, and results will be placed in that directory.

