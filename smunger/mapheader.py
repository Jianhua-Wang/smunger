"""Main module."""

import json
from typing import Optional, Iterable
from pathlib import Path

import pandas as pd
from rich.prompt import Confirm
from rich.table import Table

from smunger.console import console
from smunger.constant import COMMON_COLNAMES, ColName

# display first five rows of input dataframe, and guess the column names
# that are most likely to be the column names of munged dataframe
# display the guessed column map, let user to confirm
# if user confirms, then use the guessed column map
# if user does not confirm, then ask user to input the column map
# TODO:batch mode: remember the column map for each input file


def display_df(df: pd.DataFrame, colname_map: dict, nrows: int = 5):
    """Display the first five rows of a dataframe."""
    if len(df) < 5:
        nrows = len(df)
    console.print(f"Displaying first {nrows} rows of the dataframe:")
    console.print("First row is the number of the column.")
    console.print("Second row is the original column name.")
    console.print("Third row is the gussed mapped column name.")
    table = Table(show_header=True, header_style="bold magenta")
    original_colnames = []
    mapped_colnames = []
    value_colidx = []
    for i, col in enumerate(df.columns):
        table.add_column(f"No.{i}", justify="center", style="cyan", no_wrap=False)

        original_colnames.append(col)
        value_colidx.append(i)
        if col in colname_map:
            mapped_colnames.append(colname_map[col])
        else:
            mapped_colnames.append('')
        if (i + 1) % 10 == 0 or i == len(df.columns) - 1:
            table.add_row(*original_colnames, style="bold red", end_section=True)
            table.add_row(*mapped_colnames, style="bold green", end_section=True)
            for n in range(nrows):
                table.add_row(*[str(_) for _ in df.iloc[n, value_colidx].values])
            console.print(table)
            table = Table(show_header=True, header_style="bold magenta")
            original_colnames = []
            mapped_colnames = []
            value_colidx = []
    missing_colnames = set(ColName.OUTCOLS) - set(colname_map.values())
    if missing_colnames:
        console.print(f"Missing columns: {missing_colnames}")


def manul_map(df_cols: Iterable, guessed_map: dict = {}) -> dict:
    """Manual map."""
    colname_map = {}
    df_cols = list(df_cols)
    guess_orignal_out = {v: k for k, v in guessed_map.items()}
    for col in ColName.OUTCOLS:
        if col in guess_orignal_out:
            map_correct = Confirm.ask(
                f"Is {guess_orignal_out[col]} the correct column name for {col}?", default=True
            )
            if map_correct:
                colname_map[guess_orignal_out[col]] = col
                continue
        manual_map = input(f"Please input the correct column number for {col}: ")
        if manual_map == '':
            continue
        else:
            manual_map = df_cols[int(manual_map)]
            colname_map[manual_map] = col
    # use OR if BETA is not found
    if ColName.BETA not in colname_map.values() or ColName.SE not in colname_map.values():
        console.print("BETA is not found.")
        for col in [ColName.OR, ColName.ORSE, ColName.Z]:
            if col in guess_orignal_out:
                map_correct = Confirm.ask(
                    f"Is {guess_orignal_out[col]} the correct column name for {col}?", default=True
                )
                if map_correct:
                    colname_map[guess_orignal_out[col]] = col
                else:
                    manual_map = input(f"Please input the correct column number for {col}: ")
                    if manual_map == '':
                        continue
                    else:
                        manual_map = df_cols[int(manual_map)]
                        colname_map[manual_map] = col
    # use -log10(P) if P is not found
    if ColName.P not in colname_map.values():
        console.print("P is not found.")
        for col in [ColName.NEGLOGP]:  # TODO: use NEGLOGP if it is found
            if col in guess_orignal_out:
                map_correct = Confirm.ask(
                    f"Is {guess_orignal_out[col]} the correct column name for {col}?", default=True
                )
                if map_correct:
                    colname_map[guess_orignal_out[col]] = col
                else:
                    manual_map = input(f"Please input the correct column number for {col}: ")
                    if manual_map == '':
                        continue
                    else:
                        manual_map = df_cols[int(manual_map)]
                        colname_map[manual_map] = col
    return colname_map


def map_colnames(df: pd.DataFrame, outfile: Optional[Path] = None) -> dict:
    """Map column names."""
    guessed_map = guess_colnames(df.columns)
    display_df(df, guessed_map)
    use_guessed_map = Confirm.ask("Do you want to use the guessed column map?", default=True)
    if use_guessed_map:
        colname_map = guessed_map
    else:
        colname_map = manul_map(df.columns, guessed_map)
    console.print(f"Column map is: {colname_map}")
    display_df(df, colname_map)
    if outfile:
        console.print(f"Saving column map to {outfile}")
        with open(outfile, 'w') as f:
            json.dump(colname_map, f, indent=4)
    return colname_map


def guess_colnames(df_cols: Iterable[str]) -> dict:
    """Guess column names."""
    colnames = {}
    for col in df_cols:
        if col in COMMON_COLNAMES:
            colnames[col] = COMMON_COLNAMES[col]
    return colnames
