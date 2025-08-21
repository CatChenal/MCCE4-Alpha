#!/usr/bin/env python
"""Module: job_setup.py

Contains functions to prepare a user's benchmarking folder using user-provided options
(from cli args if cli is used).

Functions:
----------
* setup_pdbs_folder(bench_dir:str) -> None:
    Replicate current setup.
    - Create a copy of BENCH_PDBS (packaged data) in user_pdbs_folder = `bench_dir`/runs,
      or in user_pdbs_folder = `./runs` if called from within `bench_dir`;
    - Soft-link the relevant pdb as "prot.pdb";
    - Copy the "queue book" and default script files (BENCH.BENCH_Q_BOOK, BENCH.DEFAULT_JOB_SH, respectively)
      in `user_pdbs_folder`;
    - Copy ancillary files BENCH.BENCH_WT, BENCH.BENCH_PROTS `bench_dir`.

* delete_sentinel(bench_dir:str, sentinel_file:str) -> None:
    Part of the job preparation for each new script.
    Delete sentinel_file from 'bench_dir'/runs subfolders.

* write_run_script(job_name, steps_options_dict)
    Beta Phase : job_name = "default_run" (or soft link to 'default_run.sh' if different).
    Write a shell script in user_job_folder similar to RUN_SH_DEFAULTS.

    Current default template: (BENCH.DEFAULT_JOB_SH, "default_run.sh"):
     ```
     #!/bin/bash

     step1.py prot.pdb --dry
     step2.py
     step3.py
     step4.py --xts

     sleep 10
     ```
"""

from argparse import Namespace
import logging
import os
import pandas as pd
from pathlib import Path
import shutil
import sys
from typing import Union
from mcce4.mcce_benchmark import BenchResources, BENCH_DATA, Q_BOOK, DEFAULT_JOB, DEFAULT_JOB_SH
from mcce4.mcce_benchmark import SUB1, RUNS_DIR, N_PDBS, N_BATCH
from mcce4.mcce_benchmark import audit
from mcce4.mcce_benchmark.io_utils import make_executable, Pathok


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


CREATED_PRM_MSG = """
# custom.prm created by bench app to store user key(s) in -u options.
# Note: The -u option in the bench app is used for relinking customized name.txt or extra.tpl files ONLY;
#       Additional options must be given in a customized run.prm file passed to -load_runprm.
#
"""


def get_pkdb_list(dataset:str,
                  return_df: bool = False,
                  n_pdbs: int = None) -> Union[pd.DataFrame, Path]:
    """Write non-excluded 'PDBIDS' file to enable selection by file.
    If n_pdbs is not None, creates 'PDBIDS_<n_pdbs>_smallest'.
    - n_pdbs: To save & return a list of n smaller pdbids.
    File is sorted ascendingly by number of residues, Res#.
    """
    df = pd.read_csv(BenchResources(dataset).BENCH_PROTS, sep="\t")
    cols = "PDB Res# Method Function Biounits Resolution".split()
    df = df[cols]
    # filter out commented rows
    msk = df.PDB.str.startswith("#")
    df = df[~msk]
    df.sort_values(by="Res#", ascending=True, inplace=True)

    if n_pdbs is None or n_pdbs >= N_PDBS:
        out_fp = Path(f"PDBIDS_{dataset}")
    else:
        out_fp = Path(f"PDBIDS_{dataset}_{n_pdbs}_smallest")
        df = df[:n_pdbs]

    out_fp.write_text(df.to_string(index=False) + "\n")
    logger.info(f"Saved pkDB pdbs list in: {out_fp!s}")

    if return_df:
        return df
    return out_fp


def create_next_commands(args: Namespace):
    """Create scripts 'launch.sh' and 'analyse.sh' with the pre-populated commands.
    The launch command uses the scheduler (bench_setup launch).
    """
    if not hasattr(args, "n_batch"):
        setattr(args, "n_batch", N_BATCH)
    n_batch = args.n_batch

    cmd = ("#!/bin/bash\n\n"
           "# command to launch runs with the scheduler:\n\n"
           "bench_setup launch -bench_dir . -job_name "
           f"{args.job_name} -n_batch {n_batch} -sentinel_file {args.sentinel_file}\n\n"
    )
    
    launch_sh = Path(args.bench_dir).joinpath("launch.sh")
    if launch_sh.exists():
        launch_sh.unlink()
    launch_sh.write_text(cmd)
    make_executable(launch_sh)
    logger.info("Use launch.sh to submit the jobs to the scheduler.")

    cmd = ("#!/bin/bash\n\n"
           "# command to run the analysis when the runs have completed.\n\n"
           f"bench_analyze {args.subparser_name} -bench_dir . \n\n"
    )
    ana_sh = Path(args.bench_dir).joinpath("analyze.sh")
    if ana_sh.exists():
        ana_sh.unlink()
    ana_sh.write_text(cmd)
    make_executable(ana_sh)
    logger.info("Use analyze.sh to run the analysis when all jobs have completed.")

    return

def relink_files(bench_dir: Path, custom_files: dict):
    """
    custom_files (dict): key=filename local to prot folder, value=source
    """
    for d in bench_dir.joinpath(RUNS_DIR).glob("*/"):
        if not d.is_dir():
            continue

        os.chdir(d)
        for fname, src_path in custom_files.items():
                #local file:
                fp = d.joinpath(fname)
                try:
                    fp.symlink_to(src_path)
                except FileExistsError:
                    if not fp.is_symlink() or (fp.resolve().name != src_path.name):
                        fp.unlink()
                        fp.symlink_to(src_path)
                        logger.info(f"Soft-linked {src_path} as {fname} in {d.name}.")
        os.chdir("../")

    return


def setup_customized_files(args: Namespace) -> Namespace:
    """Soft-link the file paths given to the '-u' or '-load_runprm command line options into
    each of the proteins foldeers.
    1) If `args.u` contains other keys besides EXTRA and RENAME_RULES, a custom runprm is created
       to pass the additional keys.
    2) `args.u` is reset to its default value ("").
    3) `args.load_runprm` is reset to the linked filename in the protein folders: 'run.prm.user'.
    
    * Expected, but not checked: Entries that must appear in customized runprm (args.load_runprm)
      if the corresponding key is given in args.u:
        name.txt     (RENAME_RULES)
        extra.tpl    (EXTRA)

    Arguments:
      - args (argparse.Namespace) :: arguments parsed from the cli.
    Returns:
      - the input args (argparse.Namespace) with the reset -u option and updated -load_runprm option
        to match the linked files in the prot folders.
    """
    # case 0: function was called without prior checks:
    if not args.load_runprm and not args.u:
        logger.info("No customized files to soft-link.")
        return args
    
    bdir = Path(args.bench_dir)
    user_prm = args.load_runprm != ""
    created_prm = False

    # other cases: user prm or -u data (or both)
    if user_prm:
        customprm = bdir.joinpath(args.load_runprm)
        if not customprm.exists():
            logger.info(f"Customized run.prm not found: {customprm!s}")
            sys.exit(1)

        # simplest case, customized run.prm, but no u option:
        if not args.u:
            relink_files(bdir, {"run.prm.user": customprm})
            args.load_runprm = "run.prm.user"
            return args

    if args.u:
        u_info = {}
        # read keys, values into dict:
        for field in args.u.split(","):
            try:
                key, value = field.split("=")
                u_info[key] = value
            except ValueError:
                logger.critical(f"Each component format is 'KEY=VALUE'. Unrecognized: {field}.")
                sys.exit(1)

        to_relink = {}   # will only hold the EXTRA or RENAME_RULES data
        custom_found = False

        k = "EXTRA" 
        if u_info.get(k) is not None:
            custom_found = True
            # save the info:
            to_relink["extra.tpl"] = bdir.joinpath(u_info[k])
            # remove
            _ = u_info.pop(k)

        k = "RENAME_RULES"
        if u_info.get(k) is not None:
            custom_found = True
            # save the info:
            to_relink["name.txt"] = bdir.joinpath(u_info[k])
            # remove
            _ = u_info.pop(k)

        if len(u_info):   # not empty: there were more than 2 keys in the -u option
            if not user_prm:
                # create a user run.prm file:
                customprm = bdir.joinpath("custom.prm")
                customprm.touch()
                customprm.write_text(CREATED_PRM_MSG)
                to_relink["run.prm.user"] = customprm
                user_prm = True
                created_prm = True
                args.load_runprm = "run.prm.user"

            # update with additional entries:
            with open(customprm, "a") as fp:
                for k,v in u_info.items():
                    fp.write(f"{v:30s}    ({k})\n")
    
        if custom_found:
            if not user_prm:
                # create a user run.prm file:
                customprm = bdir.joinpath("custom.prm")
                customprm.touch()
                customprm.write_text(CREATED_PRM_MSG)
                # update with customized files info:
                with open(customprm, "a") as fp:
                    k = "name.txt"
                    if to_relink.get(k) is not None:
                        fp.write(f"{k:30s}    (RENAME_RULES)\n")
                    k = "extra.tpl"
                    if to_relink.get(k) is not None:
                        fp.write(f"{k:30s}    (EXTRA)\n")

                to_relink["run.prm.user"] = customprm
                args.load_runprm = "run.prm.user"
            else:
                # assume custom prm properly setup (has customized filenames local to prot folders)
                to_relink["run.prm.user"] = customprm
                args.load_runprm = "run.prm.user"  # reset to filename local to prot folder
        
    args.u = ""  # reset to default

    # finally perform the relinking:
    relink_files(bdir, to_relink)
    logger.info("All files relinked into the proteins run folders.")

    return args


def setup_ancillaries(bench_dir: Path, subcmd: str=SUB1):
    """Install a copy of the book.txt and default script files in <bench_dir>/runs."""
    runs_dir = bench_dir.joinpath(RUNS_DIR)
    if not runs_dir.exists():
        sys.exit("The runs subfolder should have been setup prior to installing the book and default script files!")
    bench = BenchResources(subcmd)

    # copy script file:
    dest = runs_dir.joinpath(DEFAULT_JOB_SH)
    if not dest.exists():
        shutil.copy(bench.DEFAULT_JOB_SH, dest)
        logger.info(f"Setup default script: {dest!s}")

    book_fp = runs_dir.joinpath(Q_BOOK)
    if not book_fp.exists():
        audit.rewrite_book_file(book_fp, subcmd=subcmd)
        logger.info(f"Setup bookkeeping file: {book_fp!s}")

    return


def setup_user_runs(args: Namespace) -> None:
    """
    - Create subfolders for the pdbs found in args.pdbs_list, which is a file or dir path,
      in <bench_dir>/runs;
    - Soft-link the relevant pdb as "prot.pdb";
    - Create a "queue book" and default script files in <bench_dir>/runs;
    """
    bench_dir = Path(args.bench_dir)
    curr = Path.cwd()
    in_benchmarks = curr.name == bench_dir.name
    if in_benchmarks:
        logger.info(f"Call from within {bench_dir!s}, not re-created.")
    else:
        if not bench_dir.exists():
            bench_dir.mkdir()

    pdbs_in_file = False
    pdbs_lst = []
    # args.pdbs_list: from file or dir
    p = Path(args.pdbs_list)
    if p.is_dir():
        pdbs_lst = list(p.glob("*.pdb"))
        if not pdbs_lst:
            logger.error(f"No pdbs in {p}.")
            sys.exit(1)
    else:
        pdbs_in_file = True
        with open(p) as f:
            for lin in f:
                line = lin.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                pdb_fp = Path(line)
                if pdb_fp.exists():
                    if pdb_fp.is_symlink():
                        logger.error(f"Cannot use a linked file as pdb source: {pdb_fp!s}")
                        continue
                    pdbs_lst.append(pdb_fp)

        if not pdbs_lst:
            logger.error(f"None of the pdbs in {p!s} were found.")
            sys.exit(1)

    runs_dir = bench_dir.joinpath(RUNS_DIR)
    if not runs_dir.exists():
        runs_dir.mkdir()

    for i, fp in enumerate(pdbs_lst):
        if pdbs_in_file:
            # use parent folder name as prot folder:
            pname = fp.parent.name
        else:
            pname = fp.stem
        # create pdb dir:
        pd = runs_dir.joinpath(pname.upper())
        if not pd.is_dir():
            pd.mkdir()

        fp_dest = pd.joinpath(fp.name)
        if not fp_dest.exists():
            shutil.copy(fp, fp_dest, follow_symlinks=False)

        # cd to avoid links with long names:
        os.chdir(pd)
        prot = Path("prot.pdb")
        try:
            prot.symlink_to(fp.name)
        except FileExistsError:
            if not prot.is_symlink() or (prot.resolve().name != fp.name):
                prot.unlink()
                prot.symlink_to(fp.name)
                logger.info(f"Reset soft-linked pdb to prot.pdb for {pd.name}")
        os.chdir("../")  # pd.parent: runs_dir)

    os.chdir(curr)
    
    setup_ancillaries(bench_dir, subcmd=args.subparser_name)
    logger.info(f"The data setup in {runs_dir!s} went beautifully!")

    return


def list_user_pdbids(pdbids_file: str) -> list:
    """Extract the PDB from user-supplied file, which
    is presumed to be a modification of the file created by
    mcce_benchmark.cli.get_pkdb_list. The PDB ids are assumed to
    be in the first column, named 'PDB'.
    """
    return pd.read_fwf(Path(pdbids_file)).PDB.to_list()


def setup_pkdb_runs(bench_dir: str, n_pdbs: int, pdbids_file: str, subcmd: str) -> None:
    """
    Replicate current setup.
    - Create a copy of BENCH_PDBS (packaged data) in <bench_dir>/runs, or a subset
      of size (1, n_pdbs) if n_pdbs < N_PDBS.
    - Soft-link the relevant pdb as "prot.pdb";
    - Copy the "queue book" and default script files (BENCH.BENCH_Q_BOOK, BENCH.DEFAULT_JOB_SH)
      in <bench_dir>/runs;
    """
    bench_dir = Path(bench_dir)
    curr = Path.cwd()
    in_benchmarks = curr.name == bench_dir.name
    if in_benchmarks:
        logger.info(f"Call from within {bench_dir}, not re-created.")
    else:
        if not bench_dir.exists():
            logger.debug(f"Creating {str(bench_dir)}.")
            bench_dir.mkdir()

    if pdbids_file:
        logger.info(f"Loading pdbids from file {pdbids_file}.")
        PDBIDS = list_user_pdbids(pdbids_file)
        if not PDBIDS:
            logger.error("Empty pdbids list.")
            sys.exit(1)
    elif n_pdbs < N_PDBS:
        # create pdbids file with n_pdbs smallest
        pdbids_file = get_pkdb_list(subcmd, n_pdbs=n_pdbs)
        PDBIDS = list_user_pdbids(pdbids_file)
        if not PDBIDS:
            logger.error("Empty pdbids list.")
            sys.exit(1)

    runs_dir = bench_dir.joinpath(RUNS_DIR)
    if not runs_dir.exists():
        runs_dir.mkdir()

    BENCH = BenchResources(subcmd)
    # if subcmd == SUB1:
    #     valid, invalid = audit.list_all_valid_pdbs()
    # else:
    valid = [f"{d.name}/{d.name.lower()}.pdb"
                for d in list(BENCH.BENCH_PDBS.glob("./*"))
                if (d.is_dir() and d.name in PDBIDS)]

    for i, v in enumerate(valid):
        if i == n_pdbs:
            break
        # v :: PDBID/pdbid.pdb
        p = runs_dir.joinpath(v)
        d = p.parent
        #if d.name not in PDBIDS:
        #    continue
        if not d.is_dir():
            logger.info(f"Creating {d!s}.")
            d.mkdir()

        if not p.exists():
            shutil.copy(BENCH.BENCH_PDBS.joinpath(v), p)

        # if subcmd == SUB1:
        #     # also copy full if prot is multi;
        #     # needed for name validation in audit.valid_pdb:
        #     if p.name.startswith(f"{d.name.lower()}_"):
        #         if not d.joinpath(f"{d.name.lower()}.pdb.full").exists():
        #             try:
        #                 shutil.copy(
        #                     BENCH.BENCH_PDBS.joinpath(d.name, f"{d.name.lower()}.pdb.full"),
        #                     d,
        #                 )
        #                 # logger.info(f"Copied .pdb.full for {d.name}")
        #             except Exception as e:
        #                 logger.exception(f".pdb.full not found for {d.name}?", e)
        #                 raise

        # cd to avoid links with long names:
        os.chdir(d)

        prot = Path("prot.pdb")
        try:
            prot.symlink_to(p.name)
        except FileExistsError:
            if not prot.is_symlink() or (prot.resolve().name != p.name):
                prot.unlink()
                prot.symlink_to(p.name)
                logger.info(f"Reset soft-linked pdb to prot.pdb for {d.name}")

        os.chdir("../")  # d.parent)

    os.chdir(curr)

    setup_ancillaries(bench_dir)
    logger.info(f"The data setup in {runs_dir!s} went beautifully!")

    return


def delete_sentinel(bench_dir: str, sentinel_file: str) -> None:
    """Delete sentinel_file from 'bench_dir'/runs/ subfolders."""

    bench_dir = Path(bench_dir)
    fl = list(bench_dir.joinpath(RUNS_DIR).glob("./*/" + sentinel_file))
    for f in fl:
        f.unlink()
    logger.info(f"{len(fl)} {sentinel_file!r} file(s) deleted.")

    return


def get_default_script(pdb_dir: str) -> Path:
    """Re-install BENCH_DATA/DEFAULT_JOB_SH in pdb_dir if not found.
    Return its path.
    """
    pdb_dir = Path(pdb_dir)
    sh_path = pdb_dir.joinpath(DEFAULT_JOB_SH)
    if not sh_path.exists():
        shutil.copy(BENCH_DATA.joinpath(DEFAULT_JOB_SH), sh_path)
        logger.info(f"Re-installed {DEFAULT_JOB_SH}")

    return sh_path


def write_default_run_script(bench_dir: str, job_name: str = DEFAULT_JOB) -> None:
    """
    To use when cli args are all default.
    If job_name is different from "default_run", the default script is soft-linked to it
    as <job_name>.sh
    """
    bench_dir = Path(bench_dir)
    curr = Path.cwd()
    in_benchmarks = curr.name == bench_dir.name
    if in_benchmarks:
        bench_dir = curr

    user_pdbs = Pathok(bench_dir.joinpath(RUNS_DIR))

    # reinstall the default script if not found:
    default_sh = get_default_script(user_pdbs)

    sh_name = f"{job_name}.sh"
    if job_name == DEFAULT_JOB:
        sh_path = default_sh
    else:
        # soft-link default_sh to sh_name
        if not in_benchmarks:
            os.chdir(bench_dir)

        os.chdir(user_pdbs)

        sh_path = Path(sh_name)
        try:
            sh_path.symlink_to(default_sh)
        except FileExistsError:
            sh_path.unlink()
            sh_path.symlink_to(default_sh)

        logger.info(f"Soft-linked {default_sh} as {sh_name}")

        # reset path:
        sh_path = user_pdbs.joinpath(sh_name)

    logger.info(f"Script contents:\n{sh_path.read_text()}\n")
    os.chdir(curr)

    return
