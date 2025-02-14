#!/usr/bin/env python

"""
This module serves the same function as run.prm file in Stable-MCCE.

These are ideas to be implemented in the 1st iteration:

* The data structure uses key - value pairs, the same way as Stable-MCCE's
  run.prm
* It uses RunPrm class elements to hold key-value entries
* Method to convert the elements to a dictionary
* Load default values from a default run.prm file
* Can be overwritten by (in this order):
    * run.prm files
    * command arguments file
    * command arguments directly
    * direct access of elements within Python
* An auxilary command line options processing module will help converting
  options from command line to the values in dictionary
* Method to display all the elements
* Method to dump these values to a file as part of log
* Allow loading values from secondary file (1st for customized default, second
  for frequently customized values)
"""

import os
import datetime
from .mcce import logging
from .mcce import _dist_folder
from .mcce import _runprm_record_fname


# run.prm filepaths relative to the distribution / package root:
runprm_default = "runprms/run.prm.default"
runprm_deprecated = "runprms/run.prm.deprecated"


BY_CMD = "# Set by command option"


class RunPrmEntry:
    """This class defines the content of a run.prm entry with
    attribues: value, description, set_by.
    """

    def __init__(self, value="", description="", set_by="") -> None:
        """
        :param value: value of the run.prm entry
        :type value: varies, could be int, float, str, bool
        :param description: one-line description of this entry
        :type description: str
        """
        self.value = value
        self.description = description
        self.set_by = set_by


class RunPrm:
    """This class stores mcce run parameters. The parameters are given default
    values, and can be altered by command line options and direct access within
    python.

    Global values are given in __init__().
    Values can be altered by command line options by method - update_by_cmd()
    Values can be altered by reading from a file - update_from_file()
    Values can be altered by Python script by direct access.
    Last set value method is recorded.
    """

    def __init__(self) -> None:
        """Define runprm default values."""
        # from .mcce import _dist_folder
        self._dist_folder = _dist_folder

        runprm_default_file = "%s/%s" % (self._dist_folder, runprm_default)
        if os.path.isfile(runprm_default_file):
            logging.info("Load default run.prm at %s ..." % runprm_default_file)
            self.update_by_files([runprm_default_file])
        else:
            logging.warning("Default parameter %s doesn't exist" % runprm_default_file)

    def update_by_files(self, fnames):
        """Update runprm key value pairs by files"""
        for fname in fnames:
            with open(fname) as fh:
                for line in fh:
                    entry_str = line.strip().split("#")[0]
                    fields = entry_str.split()
                    if not len(fields) > 1:
                        continue
                    key_str = fields[-1]
                    if key_str[0] == "(" and key_str[-1] == ")":
                        key = key_str.strip("()").strip()
                        value = fields[0]
                        description = ""
                        if len(fields) > 2:
                            description = " ".join(fields[1:-1])
                        set_by = "# Set by %s" % fname
                        entry = RunPrmEntry(
                            value=value, description=description, set_by=set_by
                        )
                        setattr(self, key, entry)

    def warn_deprecated(self):
        runprm_deprecated_file = "%s/%s" % (self._dist_folder, runprm_deprecated)
        # load deprecated warning messages
        deprecated_warning = {}
        with open(runprm_deprecated_file) as fh:
            for line in fh:
                entry_str = line.strip().split("#")[0]
                fields = entry_str.split()
                if not len(fields) > 1:
                    continue
                key_str = fields[-1]
                if key_str[0] == "(" and key_str[-1] == ")":
                    key = key_str.strip("()").strip()
                    value = " ".join(fields[:-1])
                    deprecated_warning[key] = value

        # check loaded entries for deprecated warning messages
        attributes = vars(self)
        # for key, entry in attributes.items():
        for key in attributes:
            if not key.startswith("_"):
                if key in deprecated_warning:
                    message = "%s: %s" % (key, deprecated_warning[key])
                    logging.warning(message)

    def _process_one_option(self, fields):
        key = fields[0]
        if key == "-load_runprm":
            # must be the first to execute so that other options can overwrite the values
            prm_files = fields[1:]
            self.update_by_files(prm_files)
        elif key == "--noter":
            self.TERMINALS.value = "f"
            self.TERMINALS.set_by = BY_CMD
        elif key == "--deleteh":
            self.IGNORE_INPUT_H.value = "t"
            self.IGNORE_INPUT_H.set_by = BY_CMD
        elif key == "--dry":
            self.H2O_SASCUTOFF.value = "-0.01"
            self.H2O_SASCUTOFF.set_by = BY_CMD
        elif key == "-sas_cut":
            self.H2O_SASCUTOFF.value = fields[1]
        elif key == "-load_options":
            pass

    def update_by_cmdfile(self, fname):
        """Update runprm key value pairs by command line options from a file"""
        logging.info("Loading run options from file %s" % fname)
        with open(fname) as fh:
            for line in fh:
                option_line = line.split("#")[0].strip()
                fields = [x.strip() for x in option_line.split()]
                if len(fields) > 0 and fields[0] != "-load_options":
                    # exclude -load_options file to ensure loading file recursively
                    self._process_one_option(fields)

    def update_by_cmd(self, args):
        """Update runprm key value pairs by command line options
        Ignored options (processed separately): --load_runprm, --load_options
        """
        logging.info("Loading run options from command line.")
        # update if set and different:
        if hasattr(args, "noter") and args.noter and self.TERMINALS.value == "t":
            self.TERMINALS.value = "f"
            self.TERMINALS.set_by = BY_CMD
        if (
            hasattr(args, "deleteh")
            and args.deleteh
            and self.IGNORE_INPUT_H.value == "f"
        ):
            self.IGNORE_INPUT_H.value = "t"
            self.IGNORE_INPUT_H.set_by = BY_CMD
        if hasattr(args, "dry") and args.dry:
            self.H2O_SASCUTOFF.value = "-0.01"
            self.H2O_SASCUTOFF.set_by = BY_CMD
        if hasattr(args, "sas_cut") and args.sas_cut:
            self.H2O_SASCUTOFF.value = args.sas_cut
            self.H2O_SASCUTOFF.set_by = BY_CMD
        if hasattr(args, "ftpl_folder") and args.ftpl_folder:
            self.FTPL_FOLDER.value = args.ftpl_folder
            self.FTPL_FOLDER.set_by = BY_CMD
        if hasattr(args, "inpdb") and args.inpdb:
            self.INPDB.value = args.inpdb
            self.INPDB.set_by = BY_CMD
        if hasattr(args, "level") and args.level:
            if args.level == 1:
                logging.info("Rotamer making level 1 is set by command line.")
                self.PACK.value = "f"
                self.PACK.set_by = BY_CMD
                self.SWING.value = "f"
                self.SWING.set_by = BY_CMD
            elif args.level == 2:
                logging.info("Rotamer making level 2 is set by command line.")
                self.PACK.value = "f"
                self.PACK.set_by = BY_CMD
                self.SWING.value = "t"
                self.SWING.set_by = BY_CMD
            elif args.level == 3:
                logging.info("Rotamer making level 3 is set by command line.")
                self.PACK.value = "t"
                self.PACK.set_by = BY_CMD
                self.SWING.value = "f"
                self.SWING.set_by = BY_CMD
            else:
                runprm_default_file = "%s/%s" % (self._dist_folder, runprm_default)
                print("Using settings in run.prm.default" % runprm_default_file)

        

    def dump(self, comment=""):
        # For now use file run.prm.record in working directory to record,
        # better to move the file to a single place
        # from .mcce import _runprm_record_fname
        attributes = vars(self)
        current_time = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
        with open(_runprm_record_fname, "w") as fh:
            fh.write(f"# Recorded on {current_time}\n{comment}\n")
            for key, entry in attributes.items():
                if not key.startswith("_"):
                    line = "%-10s %-50s %20s %s\n" % (
                        entry.value,
                        entry.description,
                        "(" + key + ")",
                        entry.set_by,
                    )
                    fh.write(line)
        logging.info("Run parameters are recorded in file %s" % _runprm_record_fname)

    def amend_by_cmd(self, cmd_options):
        """Load additional run.prm entries based on the instruction in command line options
        The code here determine the order of runprm entries are loaded.
        """
        # Load runprm from various souces
        # Source 1: default run.prm
        # Always load default values, tolerate non-existent default run.prm but will give a warning.
        # runprm = RunPrm()       # create a runprm object and load default values
        # runprm.record(comment="Recorded by step1")  write to record file with optional comment
        # runprm.warn_deprecated()

        # Source 2: additional run.prm including presets
        if hasattr(cmd_options.args, "load_runprm"):
            prm_files = cmd_options.args.load_runprm
            self.update_by_files(prm_files)

        # Source 3: Command option file except -load_options which is ignored when in a file
        if hasattr(cmd_options.args, "load_options"):
            cmdfile = cmd_options.args.load_options
            if os.path.exists(cmdfile):
                self.update_by_cmdfile(cmdfile)
            elif cmdfile:  # file name specified but not exist
                logging.warning("Command option file %s not found." % cmdfile)

        # Source 4: Command options
        self.update_by_cmd(cmd_options.args)

        # # Source 5: Direct access
        # runprm.H2O_SASCUTOFF.value = "0.10"
        # runprm.H2O_SASCUTOFF.set_by = "# Set by python statement"
        # runprm.record()
