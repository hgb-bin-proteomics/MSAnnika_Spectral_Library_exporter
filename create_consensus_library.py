#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.7"
# dependencies = [
#   "pandas",
#   "tqdm",
# ]
# ///

# CONSENSUS LIBRARY EXPORTER
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


# version tracking
__version = "0.0.1"
__date = "2025-08-06"

# import packages
import argparse
import os

import pandas as pd
from tqdm import tqdm

from typing import Any
from typing import Dict
from typing import List

def get_unique_precursor(row: pd.Series) -> str:
    return f"{str(row['ModifiedPeptide']).strip()}.{int(row['PrecursorCharge'])}"

def create_consensus_library(filename: str, save_output: bool = True) -> pd.DataFrame:
    spec_lib = pd.read_csv(filename, low_memory=False)

def main(argv = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "file",
                        help = "Name/Path of the Spectronaut result file to process.",
                        type = str,
                        nargs = 1)
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args(argv)

    filename = args.file[0]
    consensus_library = create_consensus_library(filename)
    consensus_library.to_csv(filename + ".consensus.csv", index=False)

    return 0

if __name__ == "__main__":

    print(main())
