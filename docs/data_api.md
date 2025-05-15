# Data API

The Data API provides functionality to easily read the preprocessed signal data that Nanocompore uses for the analysis. This can be used for custom plotting or other purposes. In general, you just need to load the configuration file used for the analysis via the `load_config` function and then you can query the data with `get_references`, `get_reads`, and `get_pos`.

For example:

```python
>>> from nanocompore.api import load_config, get_pos

# Load the YAML configuration file to a Config object.
>>> config = load_config('analysis.yaml')

# Get the signal data for a given position:
>>> ref_id = 'ENST00000464651.1|ENSG00000166136.16|OTTHUMG00000019346.4|OTTHUMT00000051221.1|NDUFB8-204|NDUFB8|390|retained_intron|'
>>> get_pos(config, ref_id, 243)
    condition sample                                  read  intensity  dwell
0          WT   WT_2  a6f3e188-6288-4215-acdc-fe28beba411f    -1624.0   27.0
1          WT   WT_2  09923db6-eccc-497f-8621-8adeea9b1bfb     4072.0   20.0
2          WT   WT_2  f65926cc-bf13-4396-ba92-7f2f690b71d9    -2571.0    5.0
3          WT   WT_2  aebabd0a-5260-41c4-b38b-1ebb117dc0fb      586.0   16.0
4          WT   WT_2  994256e9-afab-4b54-94ff-cc37ae4cbe08     5229.0   16.0
..        ...    ...                                   ...        ...    ...
383        WT   WT_1  79df3c74-a4c6-4335-93c5-a0ca7e3aec78    -1067.0   25.0
384        WT   WT_1  f7dad9c6-d3d9-4501-85fb-6c6246a03719    -2225.0   56.0
385        WT   WT_1  8653efdc-943f-48f8-b6f1-174cc4bb1ad5     2837.0   12.0
386        WT   WT_1  fdf524f0-5bb5-45fc-a783-7e3a592eb149      462.0   30.0
387        WT   WT_1  b05c004e-5f58-4bfa-896b-ce28b4225ab2     -469.0   27.0

[388 rows x 5 columns]

```

## Reference

::: nanocompore.api

