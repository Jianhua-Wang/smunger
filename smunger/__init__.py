"""Top-level package for smunger."""

import logging

from rich.logging import RichHandler

from .constant import ColName, ColType, ColAllowNA, ColRange
from .smunger import (
    make_SNPID_unique,
    munge_allele,
    munge_pvalue,
    munge_beta,
    munge_or,
    munge_z,
    munge_se,
    munge_eaf,
    munge_maf,
    munge_chr,
    munge_bp,
)

__author__ = """Jianhua Wang"""
__email__ = 'jianhua.mert@gmail.com'
__version__ = '0.0.3'

# Set up logging
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
