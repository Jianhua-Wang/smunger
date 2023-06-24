"""Liftover summary statistics from one genome build to another."""

import logging
from typing import Tuple

import pandas as pd
import requests  # type: ignore
from liftover import get_lifter

from smunger.constant import ColName
from smunger.smunger import munge_bp, munge_chr

logger = logging.getLogger('liftover')


def liftover_singlesnp(inbuild: str, outbuild: str, chrom: int, pos: int) -> Tuple[int, int]:
    """Liftover a single SNP from one genome build to another."""
    lo = get_lifter(inbuild, outbuild)
    out = lo.query(chrom, pos)
    if len(out) > 0:
        return int(out[0][0][3:]), out[0][1]
    else:
        return 0, 0


def guess_build_by_singleSNP(rsid: str, pos: int) -> str:
    """Guess the genome build of a single SNP."""
    build_guess = ''
    try:
        a = requests.get(f'https://myvariant.info/v1/variant/{rsid}')
        a = a.json()
        chrom = a['chrom']
        hg19 = int(a['hg19']['end'])

        coordinates = {'hg19': hg19}
        for build in ['hg17', 'hg18', 'hg38']:
            coordinates[build] = liftover_singlesnp('hg19', build, chrom, hg19)[1]

        if coordinates['hg19'] == coordinates['hg17'] == coordinates['hg18'] == coordinates['hg38'] == pos:
            build_guess = ''
        elif coordinates['hg19'] == pos:
            build_guess = 'hg19'
        elif coordinates['hg17'] == pos:
            build_guess = 'hg17'
        elif coordinates['hg18'] == pos:
            build_guess = 'hg18'
        elif coordinates['hg38'] == pos:
            build_guess = 'hg38'
        else:
            build_guess = ''
    except Exception:
        pass
    finally:
        return build_guess


def guess_genome_build(df: pd.DataFrame) -> str:
    """Guess the genome build of a summary statistics file."""
    build_guess = ''
    if ColName.RSID in df.columns:
        df = df[df[ColName.RSID].notnull()]
        df = df[df[ColName.RSID].str.startswith('rs')]
        if len(df) > 0:
            guess_df = df.sample(5)
            all_guess = []
            for i in range(len(guess_df)):
                guess = guess_build_by_singleSNP(guess_df.iloc[i][ColName.RSID], guess_df.iloc[i][ColName.BP])
                if guess != '':
                    all_guess.append(guess)
            if len(all_guess) > 0:
                build_guess = max(set(all_guess), key=all_guess.count)
    return build_guess


def liftover(
    df: pd.DataFrame,
    inbuild: str,
    outbuild: str,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
) -> pd.DataFrame:
    """Liftover summary statistics from one genome build to another."""
    lo = get_lifter(inbuild, outbuild)
    df = df.rename(columns={chrom_col: ColName.CHR, pos_col: ColName.BP})
    df = munge_chr(df)
    df = munge_bp(df)
    df[ColName.CHR] = 'chr' + df[ColName.CHR].astype(str)
    df[ColName.CHR] = df[ColName.CHR].str.replace('chr23', 'chrX')
    df[outbuild] = df[[ColName.CHR, ColName.BP]].apply(lambda x: lo.query(x[0], x[1]), axis=1)
    df[ColName.CHR] = df[outbuild].apply(lambda x: x[0][0] if len(x) > 0 else 0)
    df[ColName.BP] = df[outbuild].apply(lambda x: x[0][1] if len(x) > 0 else 0)
    df.drop(outbuild, axis=1, inplace=True)
    df = munge_chr(df)
    df = munge_bp(df)
    df = df.rename(columns={ColName.CHR: chrom_col, ColName.BP: pos_col})
    return df


def liftover_file(
    infile: str,
    outfile: str,
    inbuild: str,
    outbuild: str,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
) -> None:
    """Liftover summary statistics from one genome build to another."""
    chunksize = 100000
    logger.info(f'liftover {infile}...')
    ith = 0
    for df in pd.read_csv(infile, sep='\t', chunksize=chunksize):
        logger.info(f'processing chunk {ith}...')
        df = liftover(df, inbuild, outbuild, chrom_col, pos_col)
        ith += 1
        if ith == 1:
            df.to_csv(outfile, sep='\t', index=False)
        else:
            df.to_csv(outfile, sep='\t', index=False, mode='a', header=False)
