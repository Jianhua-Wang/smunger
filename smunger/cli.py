"""Console script for smunger."""

import logging
from pathlib import Path

import typer

from smunger import __version__, console

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)


@app.callback(invoke_without_command=True, no_args_is_help=True)
def main(
    version: bool = typer.Option(False, '--version', '-V', help='Show version.'),
    verbose: bool = typer.Option(False, '--verbose', '-v', help='Show verbose info.'),
):
    """smunger: munger for GWAS summary statistics."""
    console.rule("[bold blue]smunger[/bold blue]")
    console.print(f"Version: {__version__}", justify="center")
    console.print("Author: Jianhua Wang", justify="center")
    console.print("Email: jianhua.mert@gmail.com", justify="center")
    if version:
        typer.echo(f'smunger version: {__version__}')
        raise typer.Exit()
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.info('Verbose mode is on.')
    else:
        logging.getLogger().setLevel(logging.INFO)


@app.command()
def mapheader(
    infile: Path = typer.Argument(..., help='Input summary statistics.'),
    outfile: Path = typer.Argument(..., help='Output file, json.'),
    sep: str = typer.Option(None, '--sep', '-s', help='Separator of the input file.'),
    nrows: int = typer.Option(5, '--nrows', '-n', help='Number of rows to display.'),
    skiprows: int = typer.Option(0, '--skiprows', '-k', help='Number of rows to skip.'),
    comment: str = typer.Option(None, '--comment', '-c', help='Comment character.'),
    gzipped: bool = typer.Option(None, '--gzipped', '-z', help='Input file is gzipped.'),
):
    """Map column names."""
    from smunger.io import load_sumstats
    from smunger.mapheader import map_colnames
    df = load_sumstats(infile, nrows=nrows, sep=sep, skiprows=skiprows, comment=comment, gzipped=gzipped)
    map_colnames(df, outfile=outfile)


@app.command()
def munge(
    infile: Path = typer.Argument(..., help='Input summary statistics.'),
    outfile: Path = typer.Argument(..., help='Output file, json.'),
    colmap: str = typer.Argument(..., help='Column map file, json.'),
    sep: str = typer.Option(None, '--sep', '-s', help='Separator of the input file.'),
    skiprows: int = typer.Option(0, '--skiprows', '-k', help='Number of rows to skip.'),
    comment: str = typer.Option(None, '--comment', '-c', help='Comment character.'),
    gzipped: bool = typer.Option(None, '--gzipped', '-z', help='Input file is gzipped.'),
    build_index: bool = typer.Option(True, '--build-index', '-b', help='Build tabix index.'),
):
    """Munge summary statistics."""
    import smunger
    from smunger.io import load_sumstats, save_sumstats
    df = load_sumstats(infile, sep=sep, skiprows=skiprows, comment=comment, gzipped=gzipped)
    df = smunger.extract_cols(df, colname_map=colmap)
    df = smunger.munge(df)
    save_sumstats(df, outfile, build_index=build_index)


if __name__ == "__main__":
    app()
