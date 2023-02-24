"""Read and write data from/to files."""

import gzip
import logging
import shutil
from pathlib import Path
from subprocess import PIPE, run
from typing import Optional

import pandas as pd

from .constant import ColName

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
    logger.info(f'File {filename} is gzipped: {gzipped}')
    logger.info(f'Separator is {sep}')
    logger.info(f'loading data from {filename}')
    # determine the separator, automatically if not specified
    return pd.read_csv(
        filename, sep=sep, nrows=nrows, skiprows=skiprows, comment=comment, compression='gzip' if gzipped else None
    )


def check_tool(tool: str) -> str:
    """
    Check if tool is installed.

    Parameters
    ----------
    tool : str
        The tool name.

    Returns
    -------
    str
        The path of the tool if it is installed, otherwise exit the program.
    """
    tool_path = shutil.which(tool)
    if tool_path:
        return tool_path
    else:
        logger.error(f"{tool} is not installed. Please install it first and make sure it is in your PATH.")
        raise ValueError(f"{tool} is not installed. Please install it first and make sure it is in your PATH.")


def save_sumstats(sumstats: pd.DataFrame, filename: Path, build_index: bool = True):
    """Save summary statistics to a file."""
    # save the summary statistics to a file
    if isinstance(filename, str):
        filename = Path(filename)
    if filename.exists():
        logger.warning(f'File {filename} already exists. Overwriting.')

    if filename.suffix == '.gz':
        filename = filename.with_suffix('')
    sumstats = sumstats.sort_values(by=[ColName.CHR, ColName.BP])
    logger.info(f'Saving summary statistics to {filename}')
    sumstats.to_csv(filename, sep='\t', index=False, header=True)

    # compress the file
    compress(filename)
    # index the file
    if build_index:
        index(str(filename) + '.gz')


def compress(filename: Path):
    """Compress a file with bgzip."""
    bgzip = check_tool('bgzip')
    logger.info(f'Compressing {filename} with bgzip')
    run([bgzip, '-f', str(filename)], stdout=PIPE, stderr=PIPE, check=True)


def index(filename: str, start: int = 1, end: int = 2, skip: int = 1):
    """Index a file with tabix."""
    tabix = check_tool('tabix')
    logger.info(f'Indexing {filename} with tabix')
    run(
        [tabix, '-f', '-S', str(skip), '-s', str(start), '-b', str(end), '-e', str(end), filename],
        stdout=PIPE,
        stderr=PIPE,
        check=True,
    )
