################################################################################
#                           PROJECT INFORMATION                                #
#                              Name: MUTADOCK                                  #
#                           Author: Naisarg Patel                              #
#                                                                              #
#       Copyright (C) 2026 Naisarg Patel (https://github.com/naisarg14)        #
#                                                                              #
#          Project: https://github.com/naisarg14/mutadock                      #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
################################################################################


import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Union

from .exceptions import MutationError
from .helpers import backup

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

try:
    import pandas as pd
except ImportError:
    msg = "Error with importing pandas module for mutation using mutadock.\n"
    msg += (
        "Easiest way to fix this is to install pandas using the following command:\n\n"
    )
    msg += "python -m pip install pandas\n"
    msg += "If you already have pandas installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    logger.error(msg)
    sys.exit(2)


def main() -> None:
    """CLI entry point for ``md_csv_sort``."""
    parser = argparse.ArgumentParser(
        description="This program takes as input a CSV file and sorts it according to the coloumn given.",
        epilog="Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument(
        "-i", "--input", help="CSV File to sort", metavar="CSV", required=True
    )
    parser.add_argument("-o", "--output", help="Sorted Output CSV", metavar="CSV")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-a", action="store_true", help="Sort in ascending order")
    group.add_argument("-d", action="store_true", help="Sort in descending order")

    parser.add_argument(
        "-n", "--name", metavar="NAME", help="Name of the coloumn to sort"
    )
    parser.add_argument(
        "-N", "--num", metavar="NUM", type=int, help="Number of the coloumn to sort"
    )

    args = parser.parse_args()
    if args.d:
        order = False
    else:
        order = True

    in_file = args.input

    p = Path(in_file)
    if not p.is_absolute():
        logger.info("Assuming current directory as root since path not specified.")
        full_csv_path = Path.cwd() / in_file
    else:
        full_csv_path = p.resolve()

    if not full_csv_path.is_file():
        sys.exit(
            "No such file found in current directory, enter full path for other directories."
        )

    if full_csv_path.suffix != ".csv":
        sys.exit("Given file is not a CSV file, input should be a CSV file.")

    sort_csv(
        in_file=str(full_csv_path),
        out_file=args.output,
        col_num=args.num,
        col_name=args.name,
        order=order,
    )


def sort_csv(
    in_file: str,
    out_file: Optional[str] = None,
    col_num: Optional[Union[str, int]] = None,
    col_name: Optional[str] = None,
    order: bool = True,
) -> str:
    """Sort a CSV by a chosen column and renumber the ``sr`` serial column.

    The ``sr`` column is dropped and re-inserted as a sequential index after
    sorting.  If neither *col_num* nor *col_name* is provided, the available
    column names are logged and a :class:`MutationError` is raised.

    Args:
        in_file: Path to the input CSV (the ``.csv`` suffix is optional).
        out_file: Destination path.  Defaults to ``<in_file>_sorted.csv``.
        col_num: Zero-based column index to sort by.
        col_name: Column header to sort by (case-insensitive; takes priority
            over *col_num* when both are supplied).
        order: ``True`` for ascending order (default), ``False`` for descending.

    Returns:
        Path to the written sorted CSV.

    Raises:
        FileNotFoundError: If *in_file* does not exist.
        MutationError: If neither *col_num* nor *col_name* resolves to a column.
    """
    in_file = in_file.removesuffix(".csv")
    if col_num is not None:
        col_num = int(col_num)
    if not out_file:
        out_file = f"{in_file.removesuffix('.csv')}_sorted.csv"
    backup(out_file)
    try:
        df = pd.read_csv(f"{in_file}.csv")
    except FileNotFoundError:
        logger.error(
            "No such file found in current directory, enter full path for other directories."
        )
        raise

    if col_name:
        header = list(df.columns.values)
        for i in range(len(header)):
            if header[i].lower() == col_name.lower():
                col_num = i
                break

    if col_num is None:
        count = 0
        for i in df.columns.values:
            logger.info(f"{count}: {i}")
            count += 1
        raise MutationError(
            "No column specified to sort by; provide col_num or col_name "
            "(see logged column list above)."
        )
    df = df.sort_values(by=df.columns[col_num], ascending=order)
    df = df.iloc[:, 1:]
    df.insert(0, "sr", range(1, 1 + df.shape[0]))
    df.to_csv(out_file, index=False)

    return out_file


if __name__ == "__main__":
    main()
