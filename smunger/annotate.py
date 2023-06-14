"""Module for annotating a file with the results of a Smunger run."""

import logging
import tabix
import pandas as pd
from smunger.constant import ColName
from smunger.smunger import make_SNPID_unique

logger = logging.getLogger('annotate')


def annotate_rsid(
        indf: pd.DataFrame,
        database: str,
        rsid_col: str = ColName.RSID,
        chrom_col: str = ColName.CHR,
        pos_col: str = ColName.BP,
        ea_col: str = ColName.EA,
        nea_col: str = ColName.NEA,) -> pd.DataFrame:
    """Annotate a dataframe with rsids."""
    tb = tabix.open(database)
    chunk_df = make_SNPID_unique(indf, chrom_col, pos_col, ea_col, nea_col)
    # chunk_df = chunk_df.drop_duplicates(subset=[ColName.SNPID])
    if len(chunk_df) == 0:
        return pd.DataFrame()
    chrom = chunk_df[chrom_col].iloc[0]
    start = chunk_df[pos_col].min()
    end = chunk_df[pos_col].max()
    rsid_map = pd.DataFrame(
        columns=[ColName.CHR, ColName.BP, 'rsid', 'ref', 'alt'], data=tb.query(str(chrom), start - 1, end)
    )
    rsid_map = make_SNPID_unique(rsid_map, ColName.CHR, ColName.BP, 'ref', 'alt')
    rsid_map = rsid_map.drop_duplicates(subset=[ColName.SNPID])
    rsid_map = pd.Series(data=rsid_map['rsid'].values, index=rsid_map[ColName.SNPID].values)  # type: ignore
    chunk_df[rsid_col] = chunk_df[ColName.SNPID].map(rsid_map)
    del chunk_df[ColName.SNPID]
    return chunk_df


def annotate_rsid_file(
    infile: str,
    outfile: str,
    database: str,
    chunksize: int = 2000000,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    ea_col: str = ColName.EA,
    nea_col: str = ColName.NEA,
) -> None:
    """Annotate a file with rsids."""
    ith = 0
    for df in pd.read_csv(infile, sep='\t', chunksize=100000):
        for chrom, chr_df in df.groupby(chrom_col):
            chr_df = chr_df.sort_values(pos_col)
            for i in range(chr_df[pos_col].min(), chr_df[pos_col].max() + 1, chunksize):
                chunk_df = chr_df[(chr_df[pos_col] >= i) & (chr_df[pos_col] < i + chunksize)].copy()
                chunk_df = annotate_rsid(chunk_df, database, rsid_col, chrom_col, pos_col, ea_col, nea_col)
                logging.info(f'Processing {chrom}:{i}-{i + chunksize}, chunk No.{ith}, {len(chunk_df)} rows.')
                ith += 1
                if ith == 1:
                    chunk_df.to_csv(outfile, sep='\t', index=False, mode='w')
                else:
                    chunk_df.to_csv(outfile, sep='\t', index=False, mode='a', header=False)
