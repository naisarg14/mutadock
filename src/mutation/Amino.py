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

import logging
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


aa_dict = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def main() -> None:
    logger.info(
        "This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's mutation module."
    )


def check_3(aa: str) -> bool:
    aa = str(aa).strip().upper()
    if aa in aa_dict:
        return True
    return False


def check_1(aa: str) -> bool:
    aa = str(aa).strip().upper()
    for a in aa_dict:
        if aa_dict[a] == aa:
            return True
    return False


def get_1(three: str) -> Optional[str]:
    three = str(three).strip().upper()
    if check_1(three):
        return three
    elif three in aa_dict:
        return aa_dict[three]
    else:
        return None


def get_3(one: str) -> Optional[str]:
    one = str(one).strip().upper()
    if check_3(one):
        return one
    for aa in aa_dict:
        if one == aa_dict[aa]:
            return aa
    return None


def get_dict() -> dict[str, str]:
    return aa_dict


if __name__ == "__main__":
    main()
