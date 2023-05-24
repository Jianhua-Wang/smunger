"""Read and write data from/to files."""

import gzip
import logging
import shutil
import os
from pathlib import Path
from subprocess import PIPE, run
from typing import Optional

import pandas as pd
import tabix
from .constant import ColName
from .smunger import munge

logger = logging.getLogger('io')


def load_sumstats(
    filename: str,
    sep: Optional[str] = None,
    nrows: Optional[int] = None,
    skiprows: int = 0,
    comment: Optional[str] = None,
    gzipped: Optional[bool] = None,
) -> pd.DataFrame:
    """Load summary statistics from a file."""
    # determine whether the file is gzipped
    if gzipped is None:
        gzipped = filename.endswith('gz')

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


def check_header(filename) -> bool:
    """Check if the header of a file contains the required columns."""
    header = load_sumstats(filename, nrows=5)
    if list(header.columns) != list(ColName.OUTCOLS):
        logger.error(f'Header of {filename} does not contain the required columns.')
        raise ValueError(f'Header of {filename} does not contain the required columns.')
    else:
        return True


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


def save_sumstats(sumstats: pd.DataFrame, filename: str, build_index: bool = True, bgzipped: bool = True):
    """Save summary statistics to a file."""
    # save the summary statistics to a file
    filename_path = Path(filename)
    if filename_path.exists():
        logger.warning(f'File {filename} already exists. Overwriting.')

    if filename_path.suffix == '.gz':
        filename_path = filename_path.with_suffix('')
    sumstats = sumstats.sort_values(by=[ColName.CHR, ColName.BP])
    logger.info(f'Saving summary statistics to {filename_path}')
    sumstats.to_csv(filename_path, sep='\t', index=False, header=True, float_format='%g')

    # compress the file
    if bgzipped:
        compress(str(filename_path))
    # index the file
    if build_index and bgzipped:
        index(str(filename_path) + '.gz')


def compress(filename: str):
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


def export_sumstats(
    filename: str,
    chrom: Optional[int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    rename_headers: Optional[dict] = None,
    out_filename: Optional[str] = None,
    bgzipped: bool = True,
) -> pd.DataFrame:
    """Export summary statistics to a file."""
    if chrom and start and end:
        logger.info(f'Loading summary statistics from {filename} for {chrom}:{start}-{end}')
        if not os.path.exists(filename + '.tbi'):
            raise FileNotFoundError(f'Index file {filename}tbi does not exist. Please index the file first.')
        else:
            tb = tabix.open(filename)
            indf = pd.DataFrame(columns=ColName.OUTCOLS, data=tb.query(str(chrom), start, end))
            indf = munge(indf)
    else:
        logger.info(f'Loading summary statistics from {filename}')
        indf = load_sumstats(filename)
    if rename_headers:
        logger.info(f'Renaming headers to {rename_headers}')
        indf = indf[list(rename_headers.keys())].copy()
        indf = indf.rename(columns=rename_headers)
    if out_filename:
        save_sumstats(indf, out_filename, build_index=False, bgzipped=bgzipped)
    return indf
