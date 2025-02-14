#!/usr/bin/env python

"""
Holds objects and libraries available to all in MCCE_bin/.
"""

import argparse
import os
from pathlib import Path
import sys

__version__ = "Initial stage of developing"

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
