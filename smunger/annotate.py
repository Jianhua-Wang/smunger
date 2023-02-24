"""Module for annotating a file with the results of a Smunger run."""

import logging
import tabix
import pandas as pd
from smunger.constant import ColName

logger = logging.getLogger('annotate')


def annotate_rsid(
    infile: str,
    outfile: str,
    database: str,
    chunksize: int = 2000000,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
) -> None:
    """Annotate a file with rsids."""
    tb = tabix.open(database)
    ith = 0
    for df in pd.read_csv(infile, sep='\t', chunksize=100000):
        for chrom, chr_df in df.groupby(chrom_col):
            chr_df = chr_df.sort_values(pos_col)
            for i in range(chr_df[pos_col].min(), chr_df[pos_col].max(), chunksize):
                chunk_df = chr_df[(chr_df[pos_col] >= i) & (chr_df[pos_col] < i + chunksize)].copy()
                if len(chunk_df) == 0:
                    continue
                rsid_map = pd.DataFrame(
                    columns=['chr', 'bp', 'rsid', 'ref', 'alt'], data=tb.query('1', i - 1, i + chunksize)
                )
                rsid_map = rsid_map.drop_duplicates(subset=['bp'])
                rsid_map['bp'] = rsid_map['bp'].astype(int)
                rsid_map = pd.Series(data=rsid_map['rsid'].values, index=rsid_map['bp'].values)
                chunk_df[rsid_col] = chunk_df[pos_col].map(rsid_map)
                logging.info(f'Processing {chrom}:{i}-{i + chunksize}, chunk No.{ith}, {len(chunk_df)} rows.')
                ith += 1
                if ith == 1:
                    chunk_df.to_csv(outfile, sep='\t', index=False, mode='w')
                else:
                    chunk_df.to_csv(outfile, sep='\t', index=False, mode='a', header=False)
