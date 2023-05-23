"""Visualization functions."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats.mstats import mquantiles


def get_qq_df(indf: pd.DataFrame) -> pd.DataFrame:
    """Convert GWAS P-values to quantiles."""
    q_pos = np.concatenate([np.arange(99) / len(indf['P']), np.logspace(-np.log10(len(indf['P'])) + 2, 0, 100)])
    q_data = mquantiles(indf['P'], prob=q_pos, alphap=0, betap=1, limit=(0, 1))
    qq_df = pd.DataFrame()
    qq_df['observed'] = -np.log10(q_data[1:])
    qq_df['expected'] = -np.log10(q_pos[1:])
    return qq_df


def qqplot(qq_df: pd.DataFrame):
    raise NotImplementedError


