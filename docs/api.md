# Data API

The Data API provides functionality to easily read the preprocessed signal data that Nanocompore uses for the analysis. This can be used for custom plotting or other purposes. In general, you just need to load the configuration file used for the analysis via the `load_config` function and then you can query the data with `get_references`, `get_reads`, and `get_pos`.

::: nanocompore.api

