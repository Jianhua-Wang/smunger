"""Console script for smunger."""

import typer

app = typer.Typer()

def main():
    """Main entrypoint."""
    typer.echo("smunger")
    typer.echo("=" * len("smunger"))
    typer.echo("munger for GWAS summary statistics")


if __name__ == "__main__":
   app(main)