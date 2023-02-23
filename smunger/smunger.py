"""Main module."""

import pandas as pd
import json
from typing import Union
import logging

from smunger.constant import ColName, ColType, ColRange


logger = logging.getLogger('munger')


def make_SNPID_unique(sumstat: pd.DataFrame, replace_rsIDcol: bool = False, remove_duplicates: bool = True):
    """
    Make the SNPID unique.

    The unique SNPID is chr-bp-sorted(EA,NEA)

    Parameters
    ----------
    sumstat : pd.DataFrame
        The input summary statistics.
    replace_rsIDcol : bool, optional
        Whether to replace the rsID column with the unique SNPID, by default False
    remove_duplicates : bool, optional
        Whether to remove the duplicated SNPs, keep the one with smallest P-value, by default True

    Returns
    -------
    pd.DataFrame
        The summary statistics with unique SNPID.
    """
    df = sumstat.copy()
    allele_df = df[[ColName.EA, ColName.NEA]].copy()
    b = allele_df.values
    b.sort(axis=1)
    allele_df[[ColName.EA, ColName.NEA]] = b
    allele_df[ColName.SNPID] = (
        df[ColName.CHR].astype(str)
        + "-"
        + df[ColName.BP].astype(str)
        + "-"
        + allele_df[ColName.EA]
        + "-"
        + allele_df[ColName.NEA]
    )
    if replace_rsIDcol:
        df[ColName.RSID] = allele_df[ColName.SNPID]
    else:
        df.insert(loc=0, column=ColName.SNPID, value=allele_df[ColName.SNPID].values)  # type: ignore
    if remove_duplicates:
        df.sort_values(ColName.P, inplace=True)
        if replace_rsIDcol:
            df.drop_duplicates(subset=[ColName.RSID], keep="first", inplace=True)
        else:
            df.drop_duplicates(subset=[ColName.SNPID], keep="first", inplace=True)
        df.sort_values([ColName.CHR, ColName.BP], inplace=True)
        df.reset_index(drop=True, inplace=True)
    return df


def map_colnames(df: pd.DataFrame, colname_map: Union[dict, str]) -> pd.DataFrame:
    """Map column names."""
    # map column names
    if isinstance(colname_map, str):
        with open(colname_map, 'r') as f:
            colname_map = json.load(f)
    outdf = df.rename(columns=colname_map).copy()  # type: ignore
    mapped_cols = list(colname_map.values())  # type: ignore
    outdf = outdf[mapped_cols]
    outdf = rm_col_allna(outdf)
    return outdf


def checl_colnames(df: pd.DataFrame) -> pd.DataFrame:
    """Check column names, fill None if not presents."""
    outdf = df.copy()
    for col in ColName.OUTCOLS:
        if col not in outdf.columns:
            outdf[col] = None
    return outdf


def rm_col_allna(df: pd.DataFrame) -> pd.DataFrame:
    """Remove columns that are all NA."""
    outdf = df.copy()
    for col in outdf.columns:
        if outdf[col].isnull().all():
            logger.info(f"Remove column {col} because it is all NA.")
            outdf.drop(col, axis=1, inplace=True)
    return outdf


def munge_chr(df: pd.DataFrame) -> pd.DataFrame:
    """Munge chromosome column."""
    pre_n = df.shape[0]
    outdf = df[df[ColName.CHR].notnull()].copy()
    outdf[ColName.CHR] = outdf[ColName.CHR].astype(ColType.CHR)
    outdf[ColName.CHR] = outdf[ColName.CHR].str.replace("chr", "")
    # replace X, with 23
    outdf[ColName.CHR] = outdf[ColName.CHR].replace("X", 23)
    outdf[ColName.CHR] = outdf[ColName.CHR].replace("x", 23)
    # turn chromosome column into integer
    outdf[ColName.CHR] = pd.to_numeric(outdf[ColName.CHR], errors="coerce")
    outdf = outdf[outdf[ColName.CHR].notnull()]
    outdf = outdf[(outdf[ColName.CHR] >= ColRange.CHR_MIN) & (outdf[ColName.CHR] <= ColRange.CHR_MAX)]
    after_n = outdf.shape[0]
    logger.info(f"Remove {pre_n - after_n} rows because of invalid chromosome.")
    outdf[ColName.CHR] = outdf[ColName.CHR].astype(int)
    return outdf


def munge_bp(df: pd.DataFrame) -> pd.DataFrame:
    """Munge position column."""
    pre_n = df.shape[0]
    outdf = df[df[ColName.BP].notnull()].copy()
    outdf[ColName.BP] = pd.to_numeric(outdf[ColName.BP], errors="coerce")
    outdf = outdf[outdf[ColName.BP].notnull()]
    outdf = outdf[(outdf[ColName.BP] > ColRange.BP_MIN) & (outdf[ColName.BP] < ColRange.BP_MAX)]
    after_n = outdf.shape[0]
    logger.info(f"Remove {pre_n - after_n} rows because of invalid position.")
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
        logger.info(f"Remove {pre_n - after_n} rows because of invalid {col}.")
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
    logger.info(f"Remove {pre_n - after_n} rows because of invalid pvalue.")
    outdf[ColName.P] = outdf[ColName.P].astype(ColType.P)
    return outdf


def munge_beta(df: pd.DataFrame) -> pd.DataFrame:
    """Munge beta column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.BETA] = pd.to_numeric(outdf[ColName.BETA], errors="coerce")
    outdf = outdf[outdf[ColName.BETA].notnull()]
    after_n = outdf.shape[0]
    logger.info(f"Remove {pre_n - after_n} rows because of invalid beta.")
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
    logger.info(f"Remove {pre_n - after_n} rows because of invalid se.")
    outdf[ColName.SE] = outdf[ColName.SE].astype(ColType.SE)
    return outdf


def munge_or(df: pd.DataFrame) -> pd.DataFrame:
    """Munge or column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.OR] = pd.to_numeric(outdf[ColName.OR], errors="coerce")
    outdf = outdf[outdf[ColName.OR].notnull()]
    after_n = outdf.shape[0]
    logger.info(f"Remove {pre_n - after_n} rows because of invalid or.")
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
    logger.info(f"Remove {pre_n - after_n} rows because of invalid orse.")
    outdf[ColName.ORSE] = outdf[ColName.ORSE].astype(ColType.ORSE)
    return outdf


def munge_z(df: pd.DataFrame) -> pd.DataFrame:
    """Munge z column."""
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.Z] = pd.to_numeric(outdf[ColName.Z], errors="coerce")
    outdf = outdf[outdf[ColName.Z].notnull()]
    after_n = outdf.shape[0]
    logger.info(f"Remove {pre_n - after_n} rows because of invalid z.")
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
    logger.info(f"Remove {pre_n - after_n} rows because of invalid eaf.")
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
    logger.info(f"Remove {pre_n - after_n} rows because of invalid maf.")
    outdf[ColName.MAF] = outdf[ColName.MAF].astype(ColType.MAF)
    return outdf
