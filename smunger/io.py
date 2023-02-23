"""Read and write data from/to files."""

import pandas as pd
from pathlib import Path
from typing import Optional
import gzip
from subprocess import run, PIPE
import logging
import time

logger = logging.getLogger('io')


def load_sumstats(
    filename: Path,
    sep: Optional[str] = None,
    nrows: Optional[int] = None,
    skiprows: int = 0,
    comment: Optional[str] = None,
    gzipped: Optional[bool] = None,
) -> pd.DataFrame:
    """Load summary statistics from a file."""
    # determine whether the file is gzipped
    if gzipped is None:
        gzipped = filename.suffix.endswith('gz')

    # read the first line of the file to determine the separator
    if sep is None:
        if gzipped:
            with gzip.open(filename, 'rt') as f:
                line = f.readline()
        else:
            with open(filename, 'rt') as f:
                line = f.readline()
        if '\t' in line:
            sep = '\t'
        elif ',' in line:
            sep = ','
        else:
            sep = ' '

    # determine the separator, automatically if not specified
    return pd.read_csv(
        filename, sep=sep, nrows=nrows, skiprows=skiprows, comment=comment, compression='gzip' if gzipped else None
    )


def save_sumstats(sumstats: pd.DataFrame, filename: Path, tabix: bool = True):
    """Save summary statistics to a file."""
    # save the summary statistics to a file
    if isinstance(filename, str):
        filename = Path(filename)
    if filename.exists():
        logger.warning(f'File {filename} already exists. Overwriting.')

    if filename.suffix == '.gz':
        filename = filename.with_suffix('')
    sumstats.to_csv(filename, sep='\t', index=False, header=True)

    # compress the file
    run(['bgzip', '-f', filename], stdout=PIPE, stderr=PIPE, check=True)
    # index the file
    if tabix:
        time.sleep(10)
        run(['tabix', '-f', '-p', 'vcf', filename.with_suffix('.gz')], stdout=PIPE, stderr=PIPE, check=True)
