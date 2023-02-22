"""Top-level package for smunger."""

import logging

from rich.logging import RichHandler

__author__ = """Jianhua Wang"""
__email__ = 'jianhua.mert@gmail.com'
__version__ = '0.0.2'

# Set up logging
logging.basicConfig(
    level="NOTSET",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
