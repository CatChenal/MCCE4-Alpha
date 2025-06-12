#!/usr/bin/env python

"""
Module: args.py

Convert command line arguments into classes.
"""

import argparse

class PrecheckOptions:
    def __init__(self) -> None:
        """Get the command line arguments for tool precheck."""

        helpmsg = "Precheck input structure files and report findings"
        parser = argparse.ArgumentParser(description=helpmsg)
        parser.add_argument("fnames", nargs="+", default=[], metavar="filename", help="specify input file names in pdb format")
        args = parser.parse_args()
        self.fnames = args.fnames


class Step1Options:
    def __init__(self) -> dict:
        """Get the command line arguments for program step1."""
        from .mcce import _dist_folder

        helpmsg = "Run mcce step 1, process structure and convert to MCCE format."
        parser = argparse.ArgumentParser(description=helpmsg)
        parser.add_argument("--noter", default=False, help="Do not label terminal residues (for making ftpl).", action="store_true")
        parser.add_argument("--deleteh", default=False, help="Delete all H atoms.", action="store_true")
        parser.add_argument("--dry", default=False, help="Delete all water molecules.", action="store_true")
        parser.add_argument("-sas_cut", metavar="cutoff", help="Fraction exposure threshold to be cut. Default is 0.05.", default=0.05, type=float)
        parser.add_argument("-ftpl_folder", metavar="ftpl_folder", help="Load ftpl files from alternative folder.", default="")
        parser.add_argument("-load_runprm", nargs="+", default=[], metavar="file", help="Load additional run.prm files in this order.")
        parser.add_argument("-load_options", default="", metavar="file", help="Command options can be loaded from a file, one option per line")
        parser.add_argument("inpdb", metavar="inpdb", help="Input pdb file, default is prot.pdb", default="prot.pdb", nargs="?")
        parser.add_argument("--debug", action="store_true", default=False, help="Enable debug mode")
        self.args = parser.parse_args()


class Step2Options:
    def __init__(self) -> dict:
        """Get the command line arguments for program step2."""
        from .mcce import _head1


        helpmsg = "Run mcce step 2, make conformers, and prepare for step 3 to make energy lookup table."
        parser = argparse.ArgumentParser(description=helpmsg)
        parser.add_argument("-d", metavar="epsilon", help="Dielectric constant for electrostatic interaction.", default=4.0, type=float)
        parser.add_argument("-level", metavar="level", help="Heavy atom conformer making level.", default=0, type=int)
        parser.add_argument("-ftpl_folder", metavar="ftpl_folder", help="Load ftpl files from alternative folder.", default="")
        parser.add_argument("--rot_specific", help="Use %s to overwrite heavy atom conformer making level." % _head1, default=False, action="store_true")
        parser.add_argument("--ga ", help="Run Genetic Algorithm to generate conformers.", default=False, action="store_true")
        parser.add_argument("-load_runprm", nargs="+", default=[], metavar="file", help="Load additional run.prm files in this order.")
        parser.add_argument("-load_options", default="", metavar="file", help="Command options can be loaded from a file, one option per line")
        parser.add_argument("--debug", action="store_true", default=False, help="Enable debug mode")
        self.args = parser.parse_args()


class vdw_tpl2ftp_Options:
    def __init__(self) -> dict:
        """Get the command line arguments for program vdw_tpl2ftp."""

        helpmsg = "Create RADIUS parameters from 00always_needed.tpl."
        parser = argparse.ArgumentParser(description=helpmsg)
        parser.add_argument("-i", default=False, help="Write records in place. Recommend to run without -i first.", action="store_true")
        parser.add_argument("param_folder", metavar="param_folder", help="The parameter folder that contains 00always_needed.tpl and ftpl files", nargs=1)
        self.args = parser.parse_args()

class detect_link_Options:
    def __init__(self) -> dict:
        """Get the command line arguments for program detect_link."""

        helpmsg = "Detect non-peptide bond by rule file ligand_detect_rules.ftpl."
        parser = argparse.ArgumentParser(description=helpmsg)
        parser.add_argument("-ftpl_folder", help="The folder that contains rule.", default="")
        parser.add_argument("inpdb", metavar="inpdb", help="Input pdb file in PDB format", default="prot.pdb", nargs="?")
        self.args = parser.parse_args()


class VDWOptions:
    def __init__(self) -> dict:
        """Get the command line arguments for vdw."""

        helpmsg = "Calculate breakdown vdw energy terms between conformers."
        parser = argparse.ArgumentParser(description=helpmsg)
        parser.add_argument("-ftpl_folder", metavar="ftpl_folder", help="Load ftpl files from alternative folder.", default="")
        parser.add_argument("-cutoff", metavar="cutoff", help="Display values above the cutoff", default=-0.0001, type=float)
        parser.add_argument("confs", metavar="confID", help="Specify a confID for vdw calculation", default="", nargs=1)
        self.args = parser.parse_args()


