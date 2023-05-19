"""Main module."""

import json
import logging
from typing import Union

import numpy as np
import pandas as pd

from smunger.constant import ColName, ColRange, ColType

logger = logging.getLogger('munger')


def make_SNPID_unique(
    sumstat: pd.DataFrame,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    ea_col: str = ColName.EA,
    nea_col: str = ColName.NEA,
) -> pd.DataFrame:
    """
    Make the SNPID unique.

    The unique SNPID is chr-bp-sorted(EA,NEA)

    Parameters
    ----------
    sumstat : pd.DataFrame
        The input summary statistics.

    Returns
    -------
    pd.DataFrame
        The summary statistics with unique SNPID.
    """
    df = sumstat.copy()
    allele_df = df[[ea_col, nea_col]].copy()
    b = allele_df.values
    b.sort(axis=1)
    allele_df[[ea_col, nea_col]] = b
    allele_df[ColName.SNPID] = (
        df[chrom_col].astype(str)
        + "-"
        + df[pos_col].astype(str)
        + "-"
        + allele_df[ea_col]
        + "-"
        + allele_df[nea_col]
    )
    if ColName.SNPID in df.columns:
        df.drop(ColName.SNPID, axis=1, inplace=True)
    df.insert(loc=0, column=ColName.SNPID, value=allele_df[ColName.SNPID].values)  # type: ignore
    return df


def extract_cols(df: pd.DataFrame, colname_map: Union[dict, str]) -> pd.DataFrame:
    """Map column names."""
    # map column names
    if isinstance(colname_map, str):
        with open(colname_map, 'r') as f:
            colname_map = json.load(f)
    colname_map = dict(colname_map)  # type: ignore
    colname_map = {k: v for k, v in colname_map.items() if k in df.columns}
    outdf = df.rename(columns=colname_map).copy()
    mapped_cols = list(colname_map.values())
    outdf = outdf[mapped_cols]  # type: ignore
    outdf = rm_col_allna(outdf)
    return outdf


def check_colnames(df: pd.DataFrame) -> pd.DataFrame:
    """Check column names, fill None if not presents."""
    outdf = df.copy()
    for col in ColName.OUTCOLS:
        if col not in outdf.columns:
            outdf[col] = None
    return outdf[ColName.OUTCOLS]


def rm_col_allna(df: pd.DataFrame) -> pd.DataFrame:
    """Remove columns that are all NA."""
    outdf = df.copy()
    outdf = outdf.replace('', None)
    for col in outdf.columns:
        if outdf[col].isnull().all():
            logger.debug(f"Remove column {col} because it is all NA.")
            outdf.drop(col, axis=1, inplace=True)
    return outdf


def get_sigdf(df: pd.DataFrame, pval: float = 5e-8) -> pd.DataFrame:
    """Get significant SNPs."""
    outdf = df.copy()
    if ColName.P in outdf.columns:
        outdf = outdf[outdf[ColName.P] < pval]
    else:
        raise ValueError("Missing P column.")
    return outdf


def munge(df: pd.DataFrame) -> pd.DataFrame:
    """Munge the summary statistics."""
    outdf = df.copy()
    outdf = rm_col_allna(outdf)
    if all(
        [
            ColName.CHR in outdf.columns,
            ColName.BP in outdf.columns,
            ColName.EA in outdf.columns,
            ColName.NEA in outdf.columns,
        ]
    ):
        outdf = munge_chr(outdf)
        outdf = munge_bp(outdf)
        outdf = munge_allele(outdf)
        outdf = make_SNPID_unique(outdf)
        if ColName.P in outdf.columns:
            outdf = munge_pvalue(outdf)
            outdf = outdf.sort_values(by=ColName.P)
        elif ColName.NEGLOGP in outdf.columns:
            outdf = munge_neglogp(outdf)
            outdf[ColName.P] = outdf[ColName.NEGLOGP].apply(lambda x: 10 ** (-x))
            outdf = outdf.sort_values(by=ColName.P)
        # TODO: use zscore to calculate pvalue, if pvalue is missing and zscore is present
        else:
            pass
        pre_n = outdf.shape[0]
        outdf = outdf.drop_duplicates(subset=ColName.SNPID, keep="first")
        outdf = outdf.sort_values(by=[ColName.CHR, ColName.BP])
        after_n = outdf.shape[0]
        logger.debug(f"Remove {pre_n - after_n} duplicated SNPs.")
    else:
        raise ValueError("Missing CHR, BP, EA or NEA column.")

    # outdf = munge_rsid(outdf)
    if ColName.BETA in outdf.columns and ColName.SE in outdf.columns:
        outdf = munge_beta(outdf)
        outdf = munge_se(outdf)
    elif ColName.OR in outdf.columns and ColName.ORSE in outdf.columns:
        outdf = munge_or(outdf)
        outdf = munge_orse(outdf)
        outdf[ColName.BETA] = outdf[ColName.OR]
        outdf[ColName.SE] = outdf[ColName.ORSE] / outdf[ColName.OR]
        outdf = munge_beta(outdf)
        outdf = munge_se(outdf)
        del outdf[ColName.OR]
        del outdf[ColName.ORSE]
    else:
        logger.warning("Missing BETA or SE column.")

    if ColName.Z in outdf.columns:
        outdf = munge_z(outdf)

    if ColName.EAF in outdf.columns:
        outdf = munge_eaf(outdf)
        outdf[ColName.MAF] = outdf[ColName.EAF]
    if ColName.MAF in outdf.columns:
        outdf = munge_maf(outdf)
    outdf = check_colnames(outdf)
    return outdf


def munge_rsid(df: pd.DataFrame) -> pd.DataFrame:
    """Munge rsID column."""
    outdf = df.copy()
    outdf[ColName.RSID] = outdf[ColName.RSID].astype(ColType.RSID)
    return outdf


def munge_chr(df: pd.DataFrame) -> pd.DataFrame:
    """Munge chromosome column."""
    pre_n = df.shape[0]
    outdf = df[df[ColName.CHR].notnull()].copy()
    outdf[ColName.CHR] = outdf[ColName.CHR].astype(str)
    outdf[ColName.CHR] = outdf[ColName.CHR].str.replace("chr", "")
    # replace X, with 23
    outdf[ColName.CHR] = outdf[ColName.CHR].replace("X", 23)
    outdf[ColName.CHR] = outdf[ColName.CHR].replace("x", 23)
    # turn chromosome column into integer
    outdf[ColName.CHR] = pd.to_numeric(outdf[ColName.CHR], errors="coerce")
    outdf = outdf[outdf[ColName.CHR].notnull()]
    outdf = outdf[(outdf[ColName.CHR] >= ColRange.CHR_MIN) & (outdf[ColName.CHR] <= ColRange.CHR_MAX)]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid chromosome.")
    outdf[ColName.CHR] = outdf[ColName.CHR].astype(ColType.CHR)
    return outdf


def munge_bp(df: pd.DataFrame) -> pd.DataFrame:
    """Munge position column."""
    pre_n = df.shape[0]
    outdf = df[df[ColName.BP].notnull()].copy()
    outdf[ColName.BP] = pd.to_numeric(outdf[ColName.BP], errors="coerce")
    outdf = outdf[outdf[ColName.BP].notnull()]
    outdf = outdf[(outdf[ColName.BP] > ColRange.BP_MIN) & (outdf[ColName.BP] < ColRange.BP_MAX)]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid position.")
    outdf[ColName.BP] = outdf[ColName.BP].astype(ColType.BP)
    return outdf


def munge_allele(df: pd.DataFrame) -> pd.DataFrame:
    """Munge allele column."""
    outdf = df.copy()
    for col in [ColName.EA, ColName.NEA]:
        pre_n = outdf.shape[0]
        outdf = outdf[outdf[col].notnull()]
        outdf[col] = outdf[col].astype(str)
        outdf[col] = outdf[col].str.upper()
        # make sure all alleles only contain one or more ACGT characters
        outdf = outdf[outdf[col].str.match(r"^[ACGT]+$")]
        after_n = outdf.shape[0]
        logger.debug(f"Remove {pre_n - after_n} rows because of invalid {col}.")
    outdf = outdf[outdf[ColName.EA] != outdf[ColName.NEA]]
    return outdf


def munge_pvalue(df: pd.DataFrame) -> pd.DataFrame:
    """Munge pvalue column."""
    outdf = df.copy()
    pre_n = outdf.shape[0]
    outdf[ColName.P] = pd.to_numeric(outdf[ColName.P], errors="coerce")
    outdf = outdf[outdf[ColName.P].notnull()]
    outdf = outdf[(outdf[ColName.P] > ColRange.P_MIN) & (outdf[ColName.P] < ColRange.P_MAX)]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid pvalue.")
    outdf[ColName.P] = outdf[ColName.P].astype(ColType.P)
    return outdf


def munge_neglogp(df: pd.DataFrame) -> pd.DataFrame:
    """Munge neglogp column."""
    outdf = df.copy()
    pre_n = outdf.shape[0]
    outdf[ColName.NEGLOGP] = pd.to_numeric(outdf[ColName.NEGLOGP], errors="coerce")
    outdf = outdf[outdf[ColName.NEGLOGP].notnull()]
    outdf = outdf[outdf[ColName.NEGLOGP] > ColRange.NEGLOGP_MIN]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid neglogp.")
    outdf[ColName.NEGLOGP] = outdf[ColName.NEGLOGP].astype(ColType.NEGLOGP)
    return outdf


def munge_beta(df: pd.DataFrame) -> pd.DataFrame:
    """Munge beta column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.BETA] = pd.to_numeric(outdf[ColName.BETA], errors="coerce")
    outdf = outdf[outdf[ColName.BETA].notnull()]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid beta.")
    outdf[ColName.BETA] = outdf[ColName.BETA].astype(ColType.BETA)
    return outdf


def munge_se(df: pd.DataFrame) -> pd.DataFrame:
    """Munge se column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.SE] = pd.to_numeric(outdf[ColName.SE], errors="coerce")
    outdf = outdf[outdf[ColName.SE].notnull()]
    outdf = outdf[outdf[ColName.SE] > ColRange.SE_MIN]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid se.")
    outdf[ColName.SE] = outdf[ColName.SE].astype(ColType.SE)
    return outdf


def munge_or(df: pd.DataFrame) -> pd.DataFrame:
    """Munge or column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.OR] = pd.to_numeric(outdf[ColName.OR], errors="coerce")
    outdf = outdf[outdf[ColName.OR].notnull()]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid or.")
    outdf[ColName.OR] = outdf[ColName.OR].astype(ColType.OR)
    return outdf


def munge_orse(df: pd.DataFrame) -> pd.DataFrame:
    """Munge orse column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.ORSE] = pd.to_numeric(outdf[ColName.ORSE], errors="coerce")
    outdf = outdf[outdf[ColName.ORSE].notnull()]
    outdf = outdf[outdf[ColName.ORSE] > ColRange.ORSE_MIN]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid orse.")
    outdf[ColName.ORSE] = outdf[ColName.ORSE].astype(ColType.ORSE)
    return outdf


def munge_z(df: pd.DataFrame) -> pd.DataFrame:
    """Munge z column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.Z] = pd.to_numeric(outdf[ColName.Z], errors="coerce")
    outdf = outdf[outdf[ColName.Z].notnull()]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid z.")
    outdf[ColName.Z] = outdf[ColName.Z].astype(ColType.Z)
    return outdf


def munge_eaf(df: pd.DataFrame) -> pd.DataFrame:
    """Munge eaf column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.EAF] = pd.to_numeric(outdf[ColName.EAF], errors="coerce")
    outdf = outdf[outdf[ColName.EAF].notnull()]
    outdf = outdf[(outdf[ColName.EAF] >= ColRange.EAF_MIN) & (outdf[ColName.EAF] <= ColRange.EAF_MAX)]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid eaf.")
    outdf[ColName.EAF] = outdf[ColName.EAF].astype(ColType.EAF)
    return outdf


def munge_maf(df: pd.DataFrame) -> pd.DataFrame:
    """Munge maf column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.MAF] = pd.to_numeric(outdf[ColName.MAF], errors="coerce")
    outdf = outdf[outdf[ColName.MAF].notnull()]
    outdf[ColName.MAF] = outdf[ColName.MAF].apply(lambda x: 1 - x if x > 0.5 else x)
    outdf = outdf[(outdf[ColName.MAF] >= ColRange.MAF_MIN) & (outdf[ColName.MAF] <= ColRange.MAF_MAX)]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid maf.")
    outdf[ColName.MAF] = outdf[ColName.MAF].astype(ColType.MAF)
    return outdf


def calculate_lambda(df: pd.DataFrame) -> float:
    """Calculate lambda."""
    from scipy import stats

    observed = np.median((df[ColName.BETA] / df[ColName.SE]) ** 2)
    lambda_gc = observed / stats.chi2.ppf(0.5, 1)
    return lambda_gc  # type: ignore


def harmonize(sumstat1, sumstat2) -> pd.DataFrame:
    """Harmonize two sumstats."""
    sumstat1 = make_SNPID_unique(sumstat1)
    sumstat2 = make_SNPID_unique(sumstat2)
    merged = pd.merge(sumstat1, sumstat2, on=ColName.SNPID, how="inner", suffixes=("_1", "_2"))
    merged[f'{ColName.BETA}_2'] = merged[f'{ColName.BETA}_2'].where(
        merged[f'{ColName.EA}_1'] == merged[f'{ColName.EA}_2'], -merged[f'{ColName.BETA}_2']
    )
    del merged[f'{ColName.EA}_2']
    del merged[f'{ColName.NEA}_2']
    del merged[f'{ColName.CHR}_2']
    del merged[f'{ColName.BP}_2']
    merged = merged.rename(
        columns={
            f'{ColName.EA}_1': ColName.EA,
            f'{ColName.NEA}_1': ColName.NEA,
            f'{ColName.CHR}_1': ColName.CHR,
            f'{ColName.BP}_1': ColName.BP,
        }
    )
    return merged
