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
"""MUTADOCK reporting subpackage.

Public data contract for the report generators (figures / HTML / PPTX).
"""

from __future__ import annotations

from mutadock.report.data import RunData, discover_run, mutation_label, summary_stats
from mutadock.report.report import generate_report

__all__ = [
    "RunData",
    "discover_run",
    "generate_report",
    "mutation_label",
    "summary_stats",
]
