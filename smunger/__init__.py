"""Top-level package for smunger."""

import logging

from rich.logging import RichHandler

from .console import console
from .constant import ColAllowNA, ColName, ColRange, ColType
from .io import load_sumstats, save_sumstats
from .mapheader import map_colnames
from .smunger import (
    make_SNPID_unique,
    extract_cols,
    munge,
    munge_allele,
    munge_beta,
    munge_bp,
    munge_chr,
    munge_eaf,
    munge_maf,
    munge_or,
    munge_pvalue,
    munge_se,
    munge_z,
)

__author__ = """Jianhua Wang"""
__email__ = 'jianhua.mert@gmail.com'
__version__ = '0.0.6'

# Set up logging
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
