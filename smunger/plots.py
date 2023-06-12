"""Visualization functions."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from typing import List

from .constant import ColName, chrom_len


def get_qq_df(indf: pd.DataFrame) -> pd.DataFrame:
    """Convert GWAS P-values to quantiles."""
    q_pos = np.concatenate(
        [np.arange(99) / len(indf[ColName.P]), np.logspace(-np.log10(len(indf[ColName.P])) + 2, 0, 100)]
    )
    q_data = mquantiles(indf[ColName.P], prob=q_pos, alphap=0, betap=1, limit=(0, 1))
    qq_df = pd.DataFrame()
    qq_df['observed'] = -np.log10(q_data[1:])
    qq_df['expected'] = -np.log10(q_pos[1:])
    return qq_df


def qqplot(qq_df: pd.DataFrame, ax, **kwargs):
    """Plot quantile-quantile plot."""
    if 'expected' not in qq_df.columns or 'observed' not in qq_df.columns:
        qq_df = get_qq_df(qq_df)
    ax.scatter(qq_df['expected'], qq_df['observed'], color='k', **kwargs)
    ax.set_xlabel('Expected')
    ax.set_ylabel('Observed')
    ax.set_title('Quantile-Quantile Plot')
    ax.plot([0, qq_df['expected'].max()], [0, qq_df['expected'].max()], color='r', linestyle='-')
    return ax


def get_manh_df(indf: pd.DataFrame, reduce_size: float = 0.3) -> pd.DataFrame:
    """Convert GWAS sumstat to dataframe for manhattan plot."""
    agg_offset = 0
    bp_offset = {}
    for chrom, offset in chrom_len.items():
        bp_offset[chrom] = agg_offset
        agg_offset += offset
    plotdf = indf[[ColName.CHR, ColName.BP, ColName.P]].copy()
    plotdf['x'] = plotdf[ColName.BP] + plotdf[ColName.CHR].map(bp_offset)
    plotdf['y'] = -np.log10(plotdf[ColName.P])
    num_x_bins = int(20000 * reduce_size)
    num_y_bins = int(plotdf['y'].max() * 10 * reduce_size)

    x_bins = np.linspace(plotdf['x'].min(), plotdf['x'].max(), num_x_bins)
    y_bins = np.linspace(plotdf['y'].min(), plotdf['y'].max(), num_y_bins)

    plotdf['x_bin'] = np.digitize(plotdf['x'], x_bins)
    plotdf['y_bin'] = np.digitize(plotdf['y'], y_bins)

    plotdf = plotdf.drop_duplicates(subset=['x_bin', 'y_bin'])
    return plotdf


def manhattan(manh_df: pd.DataFrame, ax, colors: List[str] = ['red', 'blue'], **kwargs):
    """Plot manhattan plot."""
    if 'x_bin' not in manh_df.columns or 'y_bin' not in manh_df.columns:
        plotdf = get_manh_df(manh_df)
    else:
        plotdf = manh_df.copy()
    if len(colors) < 2:
        raise ValueError('colors must have at least 2 colors')
    plotdf['color'] = plotdf[ColName.CHR].apply(lambda x: colors[x % len(colors)])
    ax.scatter(plotdf['x'], plotdf['y'], c=plotdf['color'], **kwargs)
    chrom_label = {}
    agg_offset = 0
    for k, v in chrom_len.items():
        chrom_label[str(k)] = int(v / 2) + agg_offset
        agg_offset += v
    ax.set_xticks(list(chrom_label.values()))
    ax.set_xticklabels(list(chrom_label.keys()))
    return ax


def qqman(
    indf: pd.DataFrame, figsize: tuple = (24, 6), width_ratios: List[int] = [3, 1], colors=['red', 'blue'], **kwargs
):
    """Plot QQ and Manhattan plot."""
    from matplotlib import gridspec

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    qqplot(indf, ax1, **kwargs)
    manhattan(indf, ax0, colors=colors, **kwargs)
    ax0.set_xlabel('Chromosome', fontsize=14)
    ax0.set_ylabel('-$\mathregular{log_{10}}$P', fontsize=14)
    ax1.set_xlabel('Observed -$\mathregular{log_{10}}$P', fontsize=14)
    ax1.set_ylabel('Expected -$\mathregular{log_{10}}$P', fontsize=14)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    fig.tight_layout()
