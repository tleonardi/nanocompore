# Plotting API

The Plotting API provides built-in plotting functions. For basic usage, you can generate plots dierctly via Nanocompore's command line interface ([Plotting guide](/plotting)). However, if you want to use the plotting functions programatically or want to customize them beyond the basic parameters provided by the command line interface, you can use the API. Nanocompore uses [Seaborn](https://seaborn.pydata.org), so all functions produce [Matplotlib](https://matplotlib.org) figures that you can modify.

For example:

```python
>>> from nanocompore.api import load_config
>>> from nanocompore.plotting import plot_coverage

# Load the YAML configuration file to a Config object.
>>> config = load_config('analysis.yaml')

# Get a coverage figure for a given transcript.
>>> ref_id = 'ENST00000464651.1|ENSG00000166136.16|OTTHUMG00000019346.4|OTTHUMT00000051221.1|NDUFB8-204|NDUFB8|390|retained_intron|'
>>> fig = plot_coverage(config, ref_id)

# We can now manipulate the figure to customize it
# beyond the parametrization of provided by the API.
# E.g. we can add a title:
>>> fig.axes[0].set_title('Coverage of NDUF8B-204')

# Then we update the layout and save the figure:
>>> fig.tight_layout()
>>> fig.savefig('NDUF8B_coverage.png')
```

## Reference

::: nanocompore.plotting

