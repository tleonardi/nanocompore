
# Nanocompore SampComp API demo 

## Import the SampComp module


```python
from nanocompore.SampComp import SampComp
```

## Initialise and call SampComp

#### Using a Python dictionary to specify the location of the eventalign files


```python
# Init the object
s = SampComp (
    eventalign_fn_dict = {
        'Modified': {'rep1':'./sample_files/modified_rep_1.tsv', 'rep2':'./sample_files/modified_rep_2.tsv'},
        'Unmodified': {'rep1':'./sample_files/unmodified_rep_1.tsv', 'rep2':'./sample_files/unmodified_rep_2.tsv'}},
    output_db_fn = "./results/out.db",
    fasta_fn = "./reference/ref.fa")

# Run the analysis
res = s ()
```

    Initialise SampComp and checks options
    Initialise Whitelist and checks options
    Read eventalign index files
    	References found in index: 5
    Filter out references with low coverage
    	References remaining after reference coverage filtering: 5
    Start data processing
    100%|██████████| 5/5 [00:01<00:00,  3.07 Processed References/s]


#### Using a YAML file instead to specify the files location


```python
# Init the object
s = SampComp (
    eventalign_fn_dict = "./samples.yaml",
    output_db_fn = "./results/out.db",
    fasta_fn = "./reference/ref.fa")

# Run the analysis
res = s ()
```

    Initialise SampComp and checks options
    Initialise Whitelist and checks options
    Read eventalign index files
    	References found in index: 5
    Filter out references with low coverage
    	References remaining after reference coverage filtering: 5
    Start data processing
    100%|██████████| 5/5 [00:10<00:00,  2.28s/ Processed References]


#### Tweaking statistical options


```python
# Init the object
s = SampComp (
    eventalign_fn_dict = "./samples.yaml",
    output_db_fn = "./results/out.db",
    fasta_fn = "./reference/ref.fa",
    comparison_method=["GMM", "MW", "KS"],
    sequence_context=2,
    sequence_context_weights='harmonic')

# Run the analysis
res = s ()
```
