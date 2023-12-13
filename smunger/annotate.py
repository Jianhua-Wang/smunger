"""Module for annotating a file with the results of a Smunger run."""

import logging
import tabix
import pandas as pd
from subprocess import check_output
from io import StringIO
from smunger.constant import ColName
from smunger.smunger import make_SNPID_unique

logger = logging.getLogger("annotate")


def annotate_rsid(
    indf: pd.DataFrame,
    database: str,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    ea_col: str = ColName.EA,
    nea_col: str = ColName.NEA,
) -> pd.DataFrame:
    """Annotate a dataframe with rsids."""
    tb = tabix.open(database)
    chunk_df = make_SNPID_unique(indf, chrom_col, pos_col, ea_col, nea_col)
    # chunk_df = chunk_df.drop_duplicates(subset=[ColName.SNPID])
    if len(chunk_df) == 0:
        return pd.DataFrame()
    chrom = chunk_df[chrom_col].iloc[0]
    start = chunk_df[pos_col].min()
    end = chunk_df[pos_col].max()
    rsid_map = pd.DataFrame(
        columns=[ColName.CHR, ColName.BP, "rsid", "ref", "alt"],
        data=tb.query(str(chrom), start - 1, end),
    )
    rsid_map = make_SNPID_unique(rsid_map, ColName.CHR, ColName.BP, "ref", "alt")
    rsid_map = rsid_map.drop_duplicates(subset=[ColName.SNPID])
    rsid_map = pd.Series(data=rsid_map["rsid"].values, index=rsid_map[ColName.SNPID].values)  # type: ignore
    chunk_df[rsid_col] = chunk_df[ColName.SNPID].map(rsid_map)
    del chunk_df[ColName.SNPID]
    return chunk_df


def annotate_rsid_file(
    infile: str,
    outfile: str,
    database: str,
    chunksize: int = 2000000,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    ea_col: str = ColName.EA,
    nea_col: str = ColName.NEA,
) -> None:
    """Annotate a file with rsids."""
    ith = 0
    for df in pd.read_csv(infile, sep="\t", chunksize=100000):
        for chrom, chr_df in df.groupby(chrom_col):
            chr_df = chr_df.sort_values(pos_col)
            for i in range(chr_df[pos_col].min(), chr_df[pos_col].max() + 1, chunksize):
                chunk_df = chr_df[(chr_df[pos_col] >= i) & (chr_df[pos_col] < i + chunksize)].copy()
                chunk_df = annotate_rsid(
                    chunk_df, database, rsid_col, chrom_col, pos_col, ea_col, nea_col
                )
                logger.info(
                    f"Processing {chrom}:{i}-{i + chunksize}, chunk No.{ith}, {len(chunk_df)} rows."
                )
                ith += 1
                if ith == 1:
                    chunk_df.to_csv(outfile, sep="\t", index=False, mode="w")
                else:
                    chunk_df.to_csv(outfile, sep="\t", index=False, mode="a", header=False)


def update_rsid(
    indf: pd.DataFrame,
    database: str,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    ea_col: str = ColName.EA,
    nea_col: str = ColName.NEA,
) -> pd.DataFrame:
    """Update rsids in a dataframe."""
    # tb = tabix.open(database)
    pass


def annotate_alleles_from_rsid_pos(
    indf: pd.DataFrame,
    database: str,
    rsid_col: str,
    chrom_col: str,
    pos_col: str,
    a1_col: str,
    a2_col: str,
) -> pd.DataFrame:
    """Annotate a dataframe with rsids."""
    chunk_df = indf.copy()
    tb = tabix.open(database)
    chunk_df = chunk_df[
        (chunk_df[rsid_col].notnull())
        & (chunk_df[chrom_col].notnull())
        & (chunk_df[pos_col].notnull())
    ]
    if chunk_df.shape[0] == 0:
        logger.warning("No rsids found in the input file.")
        return pd.DataFrame()
    chunk_df["index"] = chunk_df.index
    chunk_df[rsid_col] = chunk_df[rsid_col].astype(str)
    chrom = chunk_df[chrom_col].iloc[0]
    start = chunk_df[pos_col].min()
    end = chunk_df[pos_col].max()
    query_res = pd.DataFrame(
        columns=["chr_dbsnp", "bp_dbsnp", "rsid_dbsnp", "ref_dbsnp", "alt_dbsnp"],
        data=tb.query(str(chrom), start - 1, end),
    )
    query_res = query_res.astype(
        {"chr_dbsnp": str, "bp_dbsnp": int, "rsid_dbsnp": str, "ref_dbsnp": str, "alt_dbsnp": str}
    )
    chunk_df = chunk_df.merge(query_res, how="left", left_on=rsid_col, right_on="rsid_dbsnp")
    n_wrong_chrom = len(
        chunk_df[
            (chunk_df[chrom_col].astype(str) != chunk_df["chr_dbsnp"])
            & (chunk_df["chr_dbsnp"].notnull())
            & (chunk_df[chrom_col].notnull())
        ]["rsid_dbsnp"].unique()
    )
    if n_wrong_chrom > 0:
        logger.warning(f"found {n_wrong_chrom} rsids with different chrom")
    n_wrong_pos = len(
        chunk_df[
            (chunk_df[pos_col] != chunk_df["bp_dbsnp"])
            & (chunk_df["bp_dbsnp"].notnull())
            & (chunk_df[pos_col].notnull())
        ]["rsid_dbsnp"].unique()
    )
    if n_wrong_pos > 0:
        logger.warning(f"found {n_wrong_pos} rsids with different pos")
    chunk_df.loc[chunk_df[a1_col].isnull() & chunk_df[a2_col].isnull(), a1_col] = chunk_df[
        "ref_dbsnp"
    ]
    chunk_df.loc[chunk_df[a1_col].isnull() & chunk_df[a2_col].isnull(), a2_col] = chunk_df[
        "alt_dbsnp"
    ]
    chunk_df.loc[
        (chunk_df[a1_col].isnull() & chunk_df[a2_col].notnull())
        & (
            (chunk_df[a2_col] == chunk_df["ref_dbsnp"])
            | (chunk_df[a2_col] == chunk_df["alt_dbsnp"])
        ),
        a1_col,
    ] = chunk_df["ref_dbsnp"].where(
        chunk_df[a2_col] == chunk_df["alt_dbsnp"], chunk_df["alt_dbsnp"]
    )
    chunk_df.loc[
        (chunk_df[a1_col].notnull() & chunk_df[a2_col].isnull())
        & (
            (chunk_df[a1_col] == chunk_df["ref_dbsnp"])
            | (chunk_df[a1_col] == chunk_df["alt_dbsnp"])
        ),
        a2_col,
    ] = chunk_df["alt_dbsnp"].where(
        chunk_df[a1_col] == chunk_df["ref_dbsnp"], chunk_df["ref_dbsnp"]
    )
    chunk_df = chunk_df.drop_duplicates(subset=["index"])
    for col in ["chr_dbsnp", "bp_dbsnp", "rsid_dbsnp", "ref_dbsnp", "alt_dbsnp", "index"]:
        del chunk_df[col]
    return chunk_df


def annotate_alleles_from_rsid_pos_file(
    infile: str,
    outfile: str,
    database: str,
    chunksize: int = 2000000,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    a1_col: str = ColName.EA,
    a2_col: str = ColName.NEA,
) -> None:
    """Annotate a file with rsids."""
    logger.warning("Annotate alleles from rsID and coordinates is not recommended!!!")
    ith = 0
    for df in pd.read_csv(infile, sep="\t", chunksize=100000):
        for chrom, chr_df in df.groupby(chrom_col):
            chr_df = chr_df.sort_values(pos_col)
            start = chr_df[pos_col].min()
            end = start + chunksize
            while not pd.isna(start):
                chunk_df = chr_df[(chr_df[pos_col] >= start) & (chr_df[pos_col] < end)].copy()
                chunk_df = annotate_alleles_from_rsid_pos(
                    chunk_df, database, rsid_col, chrom_col, pos_col, a1_col, a2_col
                )
                logger.info(
                    f"Processing {chrom}:{start}-{end}, chunk No.{ith}, {len(chunk_df)} rows."
                )
                ith += 1
                if ith == 1:
                    chunk_df.to_csv(outfile, sep="\t", index=False, mode="w")
                else:
                    chunk_df.to_csv(outfile, sep="\t", index=False, mode="a", header=False)
                start = chr_df[chr_df[pos_col] > end][pos_col].min()
                end = start + chunksize


def annotate_pos_alleles_from_rsid(
    indf: pd.DataFrame,
    database: str,
    rsid_col: str,
    chrom_col: str,
    pos_col: str,
    a1_col: str,
    a2_col: str,
    chunksize: int = 1000000,
    overwrite: bool = False,
    remove_failed: bool = False,
) -> pd.DataFrame:
    """Annotate a dataframe with rsids."""
    chunk_df = indf.copy()
    tb = tabix.open(database)
    chunk_df[rsid_col] = chunk_df[rsid_col].astype(str)
    chunk_df = chunk_df[
        (chunk_df[rsid_col].str.startswith("rs")) & (chunk_df[rsid_col].str.len() >= 3)
    ]
    chunk_df["rsid_fake_pos"] = "1" + chunk_df[rsid_col].str[4:]
    chunk_df["rsid_fake_pos"] = pd.to_numeric(chunk_df["rsid_fake_pos"], errors="coerce")
    chunk_df = chunk_df[chunk_df["rsid_fake_pos"].notnull()]
    chunk_df["rsid_fake_chr"] = chunk_df[rsid_col].str[2:4]
    chunk_df["rsid_fake_chr"] = pd.to_numeric(chunk_df["rsid_fake_chr"], errors="coerce")
    chunk_df = chunk_df[chunk_df["rsid_fake_chr"].notnull()]

    if len(chunk_df) == 0:
        logging.warning('No rsids start with "rs" found in the input file.')
        return pd.DataFrame()
    chunk_df["rsid_fake_pos"] = chunk_df["rsid_fake_pos"].astype(int)
    chunk_df["rsid_fake_chr"] = chunk_df["rsid_fake_chr"].astype(int)
    chunk_df = chunk_df.sort_values(["rsid_fake_chr", "rsid_fake_pos"])
    logging.info(f"{chunk_df.shape[0]} rsids found in the input file.")
    out_df = []
    ith = 0
    for rsid_fake_chr, subdf in chunk_df.groupby("rsid_fake_chr"):
        chrom = rsid_fake_chr
        start = subdf["rsid_fake_pos"].min()
        end = start + chunksize
        while not pd.isna(start):
            df = subdf[(subdf["rsid_fake_pos"] >= start) & (subdf["rsid_fake_pos"] < end)].copy()
            ith += 1
            logging.info(f"Processing chunk No.{ith}, {len(df)} rows.")
            # df.to_csv("test1.csv", index=False)
            df["index"] = df.index
            end = df["rsid_fake_pos"].max()
            logging.info(f"Querying rsID: {chrom}:{start-1}-{end} from database.")
            # query_res = check_output(
            #     f"tabix {database} {chrom}:{start - 1}-{end}", shell=True
            # )
            # query_res = query_res.decode("utf-8")
            # query_res = pd.read_csv(
            #     StringIO(query_res),
            #     sep="\t",
            #     names=[
            #         "0",
            #         "rsid_fake_pos",
            #         "chr_dbsnp",
            #         "bp_dbsnp",
            #         "ref_dbsnp",
            #         "alt_dbsnp",
            #     ],
            # )
            query_res = pd.DataFrame(
                columns=[
                    "rsid_fake_chr",
                    "rsid_fake_pos",
                    "rsid_dbsnp",
                    "chr_dbsnp",
                    "bp_dbsnp",
                    "ref_dbsnp",
                    "alt_dbsnp",
                ],
                data=tb.query(str(chrom), start - 1, end),
            )
            query_res = query_res.astype(
                {
                    "rsid_fake_chr": int,
                    "rsid_fake_pos": int,
                    "chr_dbsnp": str,
                    "bp_dbsnp": int,
                    "ref_dbsnp": str,
                    "alt_dbsnp": str,
                }
            )
            df = df.merge(query_res, how="left", on=["rsid_fake_chr", "rsid_fake_pos"])
            if overwrite:
                df[chrom_col] = df["chr_dbsnp"]
                df[pos_col] = df["bp_dbsnp"]
                df[a1_col] = df["ref_dbsnp"]
                df[a2_col] = df["alt_dbsnp"]
            else:
                df[chrom_col] = df["chr_dbsnp"].where(df[chrom_col].isnull(), df[chrom_col])
                df[pos_col] = df["bp_dbsnp"].where(df[pos_col].isnull(), df[pos_col])
                n_wrong_chrom = len(
                    df[
                        (df[chrom_col].astype(str) != df["chr_dbsnp"])
                        & (df["chr_dbsnp"].notnull())
                        & (df[chrom_col].notnull())
                    ]["rsid_dbsnp"].unique()
                )
                if n_wrong_chrom > 0:
                    logging.warning(f"found {n_wrong_chrom} rsids with different chrom")
                n_wrong_pos = len(
                    df[
                        (df[pos_col] != df["bp_dbsnp"])
                        & (df["bp_dbsnp"].notnull())
                        & (df[pos_col].notnull())
                    ]["rsid_dbsnp"].unique()
                )
                if n_wrong_pos > 0:
                    logging.warning(f"found {n_wrong_pos} rsids with different pos")
                df.loc[df[a1_col].isnull() & df[a2_col].isnull(), a1_col] = df["ref_dbsnp"]
                df.loc[df[a1_col].isnull() & df[a2_col].isnull(), a2_col] = df["alt_dbsnp"]
                df.loc[
                    (df[a1_col].isnull() & df[a2_col].notnull())
                    & ((df[a2_col] == df["ref_dbsnp"]) | (df[a2_col] == df["alt_dbsnp"])),
                    a1_col,
                ] = df["ref_dbsnp"].where(df[a2_col] == df["alt_dbsnp"], df["alt_dbsnp"])
                df.loc[
                    (df[a1_col].notnull() & df[a2_col].isnull())
                    & ((df[a1_col] == df["ref_dbsnp"]) | (df[a1_col] == df["alt_dbsnp"])),
                    a2_col,
                ] = df["alt_dbsnp"].where(df[a1_col] == df["ref_dbsnp"], df["ref_dbsnp"])
            df = df.drop_duplicates(subset=["index"])
            for col in [
                "rsid_fake_pos",
                "rsid_fake_chr",
                "rsid_dbsnp",
                "chr_dbsnp",
                "bp_dbsnp",
                "ref_dbsnp",
                "alt_dbsnp",
                "index",
            ]:
                del df[col]
            out_df.append(df)
            start = subdf[subdf["rsid_fake_pos"] > end]["rsid_fake_pos"].min()
            end = start + chunksize
    out_df = pd.concat(out_df, ignore_index=True)
    out_df[pos_col] = pd.to_numeric(out_df[pos_col], errors="coerce")
    out_df = out_df[out_df[pos_col].notnull()]
    out_df[pos_col] = out_df[pos_col].astype(int)
    out_df = out_df.sort_values([chrom_col, pos_col])
    if remove_failed:
        out_df = out_df[
            (out_df[chrom_col].notnull())
            & (out_df[pos_col].notnull())
            & (out_df[a1_col].notnull())
            & (out_df[a2_col].notnull())
        ]
    return out_df


def annotate_pos_alleles_from_rsid_file(
    infile: str,
    outfile: str,
    database: str,
    chunksize: int = 2000000,
    rsid_col: str = ColName.RSID,
    chrom_col: str = ColName.CHR,
    pos_col: str = ColName.BP,
    a1_col: str = ColName.EA,
    a2_col: str = ColName.NEA,
    overwrite: bool = False,
) -> None:
    """Annotate a file with rsids."""
    logging.info("Loading input file...")
    indf = pd.read_csv(infile, sep="\t")
    logging.info(f"{indf.shape[0]} SNPs loaded.")
    out_df = annotate_pos_alleles_from_rsid(
        indf,
        database,
        rsid_col,
        chrom_col,
        pos_col,
        a1_col,
        a2_col,
        chunksize,
        overwrite,
    )
    out_df.to_csv(outfile, sep="\t", index=False)
