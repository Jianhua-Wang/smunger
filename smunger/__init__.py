"""Top-level package for smunger."""

import logging

from rich.logging import RichHandler

from .console import console
from .constant import ColAllowNA, ColName, ColRange, ColType
from .io import load_sumstats, save_sumstats, check_header, export_sumstats
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
    harmonize,
)
from .liftover import liftover, liftover_file
from .annotate import annotate_rsid, annotate_rsid_file
from .plots import qqplot, get_qq_df, manhattan, get_manh_df, qqman

__author__ = """Jianhua Wang"""
__email__ = 'jianhua.mert@gmail.com'
__version__ = '0.1.1'

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format="%(name)s - %(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)
