import argparse
import math

from typing import Union

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
import torch

from gmm_gpu.gmm import GMM
from matplotlib.figure import Figure, SubFigure
from pyfaidx import Fasta

from nanocompore.api import get_pos, get_reads
from nanocompore.common import UNCALLED4
from nanocompore.common import INTENSITY_POS
from nanocompore.common import DWELL_POS
from nanocompore.common import decode_kmer
from nanocompore.comparisons import nanstd
from nanocompore.config import Config
from nanocompore.database import ResultsDB
from nanocompore.run import Uncalled4Worker, GenericWorker
from nanocompore.transcript import Transcript


TEST_REF = "ENST00000464651.1|ENSG00000166136.16|OTTHUMG00000019346.4|OTTHUMT00000051221.1|NDUFB8-204|NDUFB8|390|retained_intron|"

TEST_PVALUE_COLUMNS = {
    'GMM': ['GMM_chi2_pvalue'],
    'KS': ['KS_intensity_pvalue', 'KS_dwell_pvalue'],
    'TT': ['TT_intensity_pvalue', 'TT_dwell_pvalue'],
    'MW': ['MW_intensity_pvalue', 'MW_dwell_pvalue'],
    'AUTO': ['auto_pvalue'],
}

def plot_pvalues(config: Config,
                 reference: str,
                 start: Union[int, None]=None,
                 end: Union[int, None]=None,
                 kind: str='lineplot',
                 threshold: Union[float, None]=0.01,
                 figsize: tuple[int, int]=(30, 10),
                 tests: Union[list[str], None]=None,
                 palette: str='Dark2') -> Union[Figure, SubFigure, None]:
    """
    Plot the p-values from the statistical tests
    performed in a Nanocompore run.

    Parameters
    ----------
    config : Config
        The configuration object for the run.
    reference : str
        Transcript reference.
    start : Union[int, None]
        Start of the region that will be plotted.
    end : Union[int, None]
        End of the region that will be plotted.
    kind : str
        Kind of plot to make. The available options
        are: lineplot, barplot
    threshold : Union[float, None]
        If set, it will indicate the p-value threshold
        as a dashed horizontal line.
    figsize : tuple[int, int]
        Size of the figure.
    tests : Union[list[str], None]
        List of tests to plot. The available options
        are: GMM, KS, TT, MW.
        If set to None (default) it would plot all
        tests that were listed in the configuration.
    palette : str
        Color palette to use.

    Returns
    -------
    Union[Figure, SubFigure, None]
       The resulting figure that will contain the
       p-values for the specified tests in the
       provided region.
    """
    db = ResultsDB(config)
    cols = []
    if tests is None:
        tests = config.get_comparison_methods()
    for test in tests:
        if test not in TEST_PVALUE_COLUMNS:
            raise KeyError(f"Test {test} not supported.")
        for col in TEST_PVALUE_COLUMNS[test]:
            cols.append(col)
    pvals = db.get_columns_for_ref(cols + ['kmer'], reference)
    pvals = pd.melt(pvals,
                    id_vars=['pos', 'kmer'],
                    value_vars=[c for c in pvals.columns if c != 'pos'],
                    var_name='test',
                    value_name='pval')
    pvals['pval'] = -np.log10(pvals.pval)
    pvals['kmer'] = pvals.kmer.apply(lambda k: decode_kmer(k, config.get_kit().len))
    if start is not None and end is not None:
        pvals = pvals[pvals.pos.between(start, end)]
    elif start is not None:
        pvals = pvals[pvals.pos >= start]
    elif end is not None:
        pvals = pvals[pvals.pos <= end]

    max_pos = pvals.pos.max()
    min_pos = pvals.pos.min()
    is_detailed_plot = max_pos - min_pos <= 50

    fig, ax = plt.subplots(figsize=figsize)

    if threshold is not None:
        ax.axhline(y=-np.log10(threshold), color='grey', linestyle='--', label=f'p-value = {threshold}')

    if kind == 'lineplot' and is_detailed_plot:
        sns.pointplot(pvals, x='pos', y='pval', hue='test', ax=ax)
    elif kind == 'lineplot':
        sns.lineplot(pvals, x='pos', y='pval', hue='test', ax=ax)
    elif kind == 'barplot':
        sns.barplot(pvals, x='pos', y='pval', hue='test', ax=ax)
    else:
        raise ValueError(f'Plot kind {kind} not supported.')

    ax.set_xlabel("Reference position")
    ax.set_ylabel("-log10(pvalue)")

    if is_detailed_plot:
        ax2 = ax.twiny()
        ax2.spines['bottom'].set_position(('axes', -0.18))
        ax2.spines['bottom'].set_visible(False)
        kmer_df = pvals.loc[:, ['pos', 'kmer']].drop_duplicates()
        ax2.set_xlim(min_pos - 0.5, max_pos + 0.5)
        ax2.set_xticks(kmer_df.pos)
        ax2.set_xticklabels(kmer_df.kmer, rotation=60)


    # plt.legend()
    # print(ax2.get_legend_handles_labels())
    # plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()

    return fig


def plot_signal(config: Config,
                reference: str,
                start: Union[int, None]=None,
                end: Union[int, None]=None,
                kind: str='violinplot',
                figsize: tuple[int, int]=(30, 10),
                split_samples: bool=False,
                markersize: int=2,
                palette: str='Dark2') -> Union[Figure, SubFigure, None]:
    """
    Plot the raw signal values (intensity and dwell time)
    for the given reference and region. Note that this
    will plot all reads from the input files without
    applying all the filtering and downsampling that
    Nanocompore does.

    Parameters
    ----------
    config : Config
        The configuration object for the run.
    reference : str
        Transcript reference.
    start : Union[int, None]
        Start of the region that will be plotted.
    end : Union[int, None]
        End of the region that will be plotted.
    kind : str
        Kind of plot to make. The available options
        are: lineplot, barplot
    figsize : tuple[int, int]
        Size of the figure.
    split_samples : bool
        If True, all samples would be plotted separately.
        By default it's False and the samples are grouped
        by condition.
    markersize : int
        Size of the points (only used for the swarmplot).
    palette : str
        Color palette to use.

    Returns
    -------
    Union[Figure, SubFigure, None]
       The resulting figure that will contain the
       intensity and log-dwell-time plots for the
       specified reference and region.
    """
    data, _, samples, conditions = get_reads(config, reference)
    positions = np.arange(data.shape[1])

    data = data[:, start:end, :]
    positions = positions[start:end]

    nreads = data.shape[0]
    group_var = samples if split_samples else conditions
    group_var_name = 'Sample' if split_samples else 'Condition'
    long_data = []
    for r in range(nreads):
        for p, pos in enumerate(positions):
            long_data.append((group_var[r],
                              pos,
                              data[r, p, INTENSITY_POS],
                              np.log10(data[r, p, DWELL_POS])))
    df = pd.DataFrame(long_data, columns=[group_var_name, 'pos', 'intensity', 'dwell'])

    fig, ax = plt.subplots(2, 1, figsize=figsize)

    if kind == 'violinplot':
        sns.violinplot(df, x='pos', y='intensity', hue=group_var_name, ax=ax[0], inner='quart', split=not split_samples)
        sns.violinplot(df, x='pos', y='dwell', hue=group_var_name, ax=ax[1], inner='quart', split=not split_samples, legend=False)
    elif kind == 'swarmplot':
        sns.swarmplot(df, x='pos', y='intensity', hue=group_var_name, ax=ax[0], size=markersize, dodge=True)
        sns.swarmplot(df, x='pos', y='dwell', hue=group_var_name, ax=ax[1], size=markersize, dodge=True, legend=False)
    elif kind == 'boxenplot':
        sns.boxenplot(df, x='pos', y='intensity', hue=group_var_name, ax=ax[0])
        sns.boxenplot(df, x='pos', y='dwell', hue=group_var_name, ax=ax[1], legend=False)
    else:
        raise ValueError(f'Plot kind {kind} not supported.')

    ax[0].set_xlabel(None)
    ax[0].set_ylabel("intensity")
    ax[1].set_xlabel("Reference position")
    ax[1].set_ylabel("log10(dwell time)")

    sns.move_legend(ax[0], "upper left", bbox_to_anchor=(1, 1), markerscale=math.ceil(10/markersize))

    plt.tight_layout()

    return fig


def plot_coverage(config: Config,
                  reference: str,
                  start: Union[int, None]=None,
                  end: Union[int, None]=None,
                  figsize: tuple[int, int]=(30, 10),
                  split_samples: bool=False,
                  palette: str='Dark2') -> Union[Figure, SubFigure, None]:
    """
    Plot the read coverage over a reference for all samples analysed.
    Note that this would plot the input coverage before applying
    any filtering that Nanocompore does before the comparison.

    Parameters
    ----------
    config : Config
        The configuration object for the run.
    reference : str
        Transcript reference.
    start : Union[int, None]
        Start of the region that will be plotted.
    end : Union[int, None]
        End of the region that will be plotted.
    figsize : tuple[int, int]
        Size of the figure.
    split_samples : bool
        If True, all samples would be plotted separately.
        By default it's False and the samples are grouped
        by condition.
    palette : str
        Color palette to use.

    Returns
    -------
    Union[Figure, SubFigure, None]
        Figure with the coverage plot.

    """
    data, reads, samples, conditions = get_reads(config, reference)
    positions = np.arange(data.shape[1])[start:end]
    data = data[:, start:end, 0]

    if split_samples:
        group_var = 'sample'
        groups = np.array(samples)
    else:
        group_var = 'condition'
        groups = np.array(conditions)

    covs = {}
    all_group_labels = np.unique(groups)
    for group in all_group_labels:
        covs[group] = np.sum(~np.isnan(data[groups == group, :]), axis=0)
    covs['pos'] = positions
    covs = pd.DataFrame(covs)
    covs = pd.melt(covs,
                   id_vars=['pos'],
                   value_vars=all_group_labels,
                   var_name=group_var,
                   value_name='cov')

    fig, ax = plt.subplots(figsize=figsize)

    sns.lineplot(covs, x='pos', y='cov', hue=group_var, palette=palette, ax=ax)
    ax.axhline(y=config.get_min_coverage(), linestyle='--', color='grey', label='Minimum coverage')

    plt.legend()
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()

    ax.set_xlabel('Reference position')
    ax.set_ylabel('Coverage')

    return fig


def plot_position(config: Config,
                  reference: str,
                  position: int,
                  figsize: tuple[int, int]=(10, 10),
                  point_size: int=20,
                  xlim: Union[tuple[Union[int, None], Union[int, None]], None]=(None, None),
                  ylim: Union[tuple[Union[int, None], Union[int, None]], None]=(None, None),
                  show_kde: bool=True,
                  kde_levels: int=10,
                  palette: Union[str, None]='Dark2') -> Union[Figure, SubFigure, None]:
    """
    Plot the signal data for a given position.

    Parameters
    ----------
    config : Config
        The configuration object for the run.
    reference : str
        Transcript reference.
    position : int
        0-based index on the reference indicating the position.
    figsize : tuple[int, int]
        Size of the figure.
    point_size : int
        Size of the points.
    xlim : Union[tuple[Union[int, None], Union[int, None]], None]
        Limits of the x-axis.
    ylim : Union[tuple[Union[int, None], Union[int, None]], None]
        Limits of the y-axis.
    kde : bool
        Whether to show the KDEs.
    kde_levels : int
        How many levels of the KDEs to show.
    palette : Union[str, None]
        Palette that will be used to determine the
        colors of both points and the GMMs. If set,
        it will override point_palette and gmm_palette.

    Returns
    -------
    Union[Figure, SubFigure, None]
       The resulting figure that will contain a
       2D plot with the observations as points
       and the gaussians obtained from the GMM.
    """
    fig, ax = plt.subplots(figsize=figsize)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    df = get_pos(config, reference, position)
    sns.scatterplot(df,
                    x='dwell',
                    y='intensity',
                    hue='condition',
                    s=point_size,
                    palette=palette,
                    ax=ax)
    if show_kde:
        sns.kdeplot(df,
                    x='dwell',
                    y='intensity',
                    levels=kde_levels,
                    hue="condition",
                    palette=palette,
                    ax=ax)

    ax.set_xlabel('Dwell time')
    ax.set_ylabel('Intensity')

    return fig


def plot_gmm(config: Config,
             reference: str,
             position: int,
             figsize: tuple[int, int]=(10, 10),
             point_size: int=20,
             xlim: Union[tuple[Union[int, None], Union[int, None]], None]=(None, None),
             ylim: Union[tuple[Union[int, None], Union[int, None]], None]=(None, None),
             gmm_levels: int=4,
             palette: Union[str, None]='Dark2',
             point_palette: Union[str, None]=None,
             gmm_palette: Union[str, None]=None) -> Union[Figure, SubFigure, None]:
    """
    Plot the GMM fitted by Nanocompore for a given position.

    Parameters
    ----------
    config : Config
        The configuration object for the run.
    reference : str
        Transcript reference.
    position : int
        0-based index on the reference indicating the position.
    figsize : tuple[int, int]
        Size of the figure.
    point_size : int
        Size of the points.
    xlim : Union[tuple[Union[int, None], Union[int, None]], None]
        Limits of the x-axis.
    ylim : Union[tuple[Union[int, None], Union[int, None]], None]
        Limits of the y-axis.
    gmm_levels : int
        How many levels of the GMM to show.
    palette : Union[str, None]
        Palette that will be used to determine the
        colors of both points and the GMMs. If set,
        it will override point_palette and gmm_palette.
    point_palette : Union[str, None]
        Palette to use for the points.
    gmm_palette : Union[str, None]
        Palette to use for the GMMs.

    Returns
    -------
    Union[Figure, SubFigure, None]
       The resulting figure that will contain a
       2D plot with the observations as points
       and the gaussians obtained from the GMM.
    """
    fasta_fh = Fasta(config.get_fasta_ref())
    ref_seq = str(fasta_fh[reference])

    transcript = Transcript(1, reference, ref_seq)

    if config.get_resquiggler() == UNCALLED4:
        worker_class = Uncalled4Worker
    else:
        worker_class = GenericWorker
    worker = worker_class(1, None, None, None, None, 'cpu', config)
    data, samples, conditions = worker._read_data(transcript)
    prepared_data = worker._prepare_data(data, samples, conditions)
    data, samples, conditions, _ = prepared_data
    data = data[[position], :, :]
    max_reads = worker._conf.get_downsample_high_coverage()
    data, samples, conditions = worker._downsample(
            data, samples, conditions, max_reads)

    std = nanstd(data, 1).unsqueeze(1)
    outliers = (((data - data.nanmean(1, keepdim=True)) / std).abs() > 3).any(2)
    data[outliers] = np.nan

    # Standardize the data
    std = nanstd(data, 1)
    data = (data - data.nanmean(1).unsqueeze(1)) / std.unsqueeze(1)

    gmm = GMM(n_components=2,
              device='cpu',
              random_seed=42,
              dtype=torch.float32)
    gmm.fit(data)

    # Means is a list with the means for each component.
    # The shape of each is (Points, Dims). We have a single point.
    c1_mean = gmm.means[0][0]
    c2_mean = gmm.means[1][0]

    # Covs is a list with the cov matrices for the components.
    # The shape of each is (Points, Dims, Dims).
    c1_cov = gmm.covs[0][0]
    c2_cov = gmm.covs[1][0]

    x1, y1 = np.random.multivariate_normal(c1_mean, c1_cov, 1000).T
    x2, y2 = np.random.multivariate_normal(c2_mean, c2_cov, 1000).T
    sampled_gaussians = pd.DataFrame(
        {'x': np.concatenate([x1, x2]),
         'y': np.concatenate([y1, y2]),
         'cluster': np.concatenate([np.full((1000,), 'Cluster 1'),
                                    np.full((1000,), 'Cluster 2')])})

    df = pd.DataFrame(data[0], columns=['intensity', 'dwell'])
    cond_labels = config.get_condition_labels()
    df['condition'] = [cond_labels[c] for c in conditions]

    if palette is not None:
        palette = {cond_labels[0]: plt.get_cmap(palette).colors[0],
                   cond_labels[1]: plt.get_cmap(palette).colors[1],
                   'Cluster 1': plt.get_cmap(palette).colors[2],
                   'Cluster 2': plt.get_cmap(palette).colors[3]}
    else:
        palette = {cond_labels[0]: plt.get_cmap(point_palette).colors[0],
                   cond_labels[1]: plt.get_cmap(point_palette).colors[1],
                   'Cluster 1': plt.get_cmap(gmm_palette).colors[0],
                   'Cluster 2': plt.get_cmap(gmm_palette).colors[1]}

    fig, ax = plt.subplots(figsize=figsize)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if xlim is None and ylim is None:
        lims = (min(df.dwell.min(),
                    df.intensity.min(),
                    sampled_gaussians.x.min(),
                    sampled_gaussians.y.min()),
                max(df.dwell.max(),
                    df.intensity.max(),
                    sampled_gaussians.x.max(),
                    sampled_gaussians.y.max()))
        ax.set_xlim(*lims)
        ax.set_ylim(*lims)

    sns.kdeplot(sampled_gaussians,
                x='x',
                y='y',
                levels=gmm_levels,
                hue="cluster",
                palette=palette,
                ax=ax)
    sns.scatterplot(df,
                    x='dwell',
                    y='intensity',
                    hue='condition',
                    s=point_size,
                    palette=palette,
                    legend=False,
                    ax=ax)

    handles = [mpatches.Patch(facecolor=col, label=label)
               for label, col in palette.items()]
    plt.legend(handles=handles)

    ax.set_xlabel('Standardized log10(dwell)')
    ax.set_ylabel('Standardized intensity')
    return fig

