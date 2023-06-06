"""Console script for smunger."""

import logging
from pathlib import Path
from enum import Enum
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
    infile: str = typer.Argument(..., help='Input summary statistics.'),
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
    infile: str = typer.Argument(..., help='Input summary statistics.'),
    outfile: str = typer.Argument(..., help='Output munged summary statistics.'),
    colmap: str = typer.Argument(..., help='Column map file, json.'),
    sep: str = typer.Option(None, '--sep', '-s', help='Separator of the input file.'),
    skiprows: int = typer.Option(0, '--skiprows', '-k', help='Number of rows to skip.'),
    comment: str = typer.Option(None, '--comment', '-c', help='Comment character.'),
    gzipped: bool = typer.Option(None, '--gzipped', '-z', help='Input file is gzipped.'),
    build_index: bool = typer.Option(True, '--build-index', '-b', help='Build tabix index.'),
    sigsnps: str = typer.Option(None, '--sigsnps', '-S', help='save significant SNPs to file.'),
    sigsnps_pval: float = typer.Option(5e-8, '--sigsnps-pval', '-P', help='p-value threshold for significant SNPs.'),
    report: str = typer.Option(None, '--report', '-R', help='save report to file.'),
):
    """Munge summary statistics."""
    import smunger
    from smunger.io import load_sumstats, save_sumstats

    df = load_sumstats(infile, sep=sep, skiprows=skiprows, comment=comment, gzipped=gzipped)
    pre_nrow = len(df)
    df = smunger.extract_cols(df, colname_map=colmap)
    df = smunger.munge(df)
    after_nrow = len(df)
    save_sumstats(df, outfile, build_index=build_index)
    from smunger.smunger import get_sigdf

    df_sig = get_sigdf(df, pval=sigsnps_pval)
    if sigsnps:
        if len(df_sig) > 0:
            save_sumstats(df_sig, sigsnps, build_index=False)
        else:
            console.print('[bold red]No significant SNPs found.[/bold red]')
    non_null_cols = df.columns[df.notnull().any()].tolist()
    report_json = {
        'in_rows': pre_nrow,
        'out_rows': after_nrow,
        'sigsnps': len(df_sig),
        'non_null_cols': non_null_cols,
    }
    if report:
        import json

        with open(report, 'w') as f:
            json.dump(report_json, f, indent=4)


@app.command()
def extract(
    infile: str = typer.Argument(..., help='Input summary statistics.'),
    outfile: str = typer.Argument(..., help='Output munged summary statistics.'),
    rename_headers: str = typer.Option(None, '--rename-headers', '-r', help='Rename headers, json.'),
    chrom: int = typer.Option(None, '--chrom', '-c', help='Chromosome.'),
    start: int = typer.Option(None, '--start', '-s', help='Start position.'),
    end: int = typer.Option(None, '--end', '-e', help='End position.'),
    no_bgzip: bool = typer.Option(True, '--no-bgzip', '-z', help='Do not bgzip output file.'),
):
    """Extract columns."""
    from smunger.io import export_sumstats
    import json
    if rename_headers:
        headers_map = json.loads(rename_headers)
    else:
        headers_map = None
    export_sumstats(
        filename=infile,
        out_filename=outfile,
        rename_headers=headers_map,
        chrom=chrom,
        start=start,
        end=end,
        bgzipped=no_bgzip,
    )


class Build(str, Enum):
    """Genome Builds."""

    hg17 = 'hg17'
    hg18 = 'hg18'
    hg19 = 'hg19'
    hg38 = 'hg38'


@app.command()
def liftover(
    infile: str = typer.Argument(..., help='Input summary statistics.'),
    outfile: str = typer.Argument(..., help='Output munged summary statistics.'),
    inbuild: Build = typer.Option(..., '--inbuild', '-i', help='Input build.'),
    outbuild: Build = typer.Option(..., '--outbuild', '-o', help='Output build.'),
    chromcol: str = typer.Option('CHR', '--chromcol', '-C', help='chromosome column.'),
    poscol: str = typer.Option('BP', '--poscol', '-p', help='position column.'),
):
    """Liftover summary statistics."""
    from smunger.liftover import liftover_file

    liftover_file(infile, outfile, inbuild, outbuild, chromcol, poscol)


@app.command()
def annorsid(
    infile: str = typer.Argument(..., help='Input summary statistics.'),
    outfile: str = typer.Argument(..., help='Output munged summary statistics.'),
    database: str = typer.Option(..., '--database', '-d', help='Database.'),
    chunksize: int = typer.Option(2000000, '--chunksize', '-c', help='Chunk size.'),
    rsidcol: str = typer.Option('rsID', '--rsidcol', '-r', help='rsid column.'),
    chromcol: str = typer.Option('CHR', '--chromcol', '-C', help='chromosome column.'),
    poscol: str = typer.Option('BP', '--poscol', '-p', help='position column.'),
    eacol: str = typer.Option('EA', '--eacol', '-e', help='effect allele column.'),
    neacol: str = typer.Option('NEA', '--neacol', '-n', help='non-effect allele column.'),
):
    """Annotate rsid."""
    from smunger.annotate import annotate_rsid_file

    annotate_rsid_file(infile, outfile, database, chunksize, rsidcol, chromcol, poscol, eacol, neacol)


if __name__ == "__main__":
    app()
