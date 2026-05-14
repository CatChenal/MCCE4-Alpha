#!/usr/bin/env python

"""
Holds objects and libraries available to all in MCCE_bin/.
"""

import argparse
import os
from pathlib import Path
import sys

__version__ = "4.0.0"


CLONE_PATH = Path(__file__).parent.parent.parent
CLONE = Path(__file__).parent.parent.parent.name

# for cli parsers: give correct repo:
CLI_EPILOG = f"\nReport issues & feature requests here: https://github.com/GunnerLab/{CLONE}/issues\n"


class Common:
    """Shared functions"""

    def print_class_elements(self) -> None:
        """
        Print all available variable elements of a class
        """
        for key, value in self.__dict__.items():
            if key[0] != "_":
                print("%15s: %30s    %-15s" % (key, str(value), str(type(value)).strip("<>")[6:].strip("'")))


__all__ = ["argparse",
           "Common",
           "os",
           "Path",
           "sys",
           "mcce"
           "runprm",
           "pdbio",
           "geom"
          ]
