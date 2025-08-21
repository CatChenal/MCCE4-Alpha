#!/usr/bin/env python

"""
Module: batch_submit.py

Main functions:
--------------
* read_book_entries(book:str = Q_BOOK) -> list:
    Read book file data using ENTRY class.
    Return a list of entry instances.

* get_running_jobs_dirs(job_name:str) -> list:
    Query shell for user's processes with job_name.
    Return a list of runs/ sub-directories where the jobs are running.

* batch_run(job_name:str, n_batch:int = N_BATCH, sentinel_file:str = "pK.out") -> None:
    Update Q_BOOK according to user's running jobs' statuses.
    Launch new jobs inside runs subfolders until the number of
    job equals n_batch.
    To be run in /runs subfolder, which is where Q_BOOK resides.

* launch_job(bench_dir:str),
             job_name:str = None,
             n_batch:int = N_BATCH,
             sentinel_file:str = "pK.out") -> None:
    Go to bench_dir/runs directory & call batch_run.

* launch_cli(argv=None)
    Entry point function.

Q book status codes:
     "i": initial setup (not submitted)
     "r": running
     "c": completed - was running, disapeared from job queue, sentinel_file generated
     "e": error - was running, disapeared from job queue and no sentinel_file
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
from collections import Counter
import logging
import os
from pathlib import Path
import sys
from typing import Union

from mcce4.mcce_benchmark import DEFAULT_JOB, Q_BOOK, FILES, RUNS_DIR, N_BATCH, USER, ENTRY_POINTS
from mcce4.mcce_benchmark.cli_utils import arg_valid_dirpath
import mcce4.mcce_benchmark.io_utils as iou
from mcce4.mcce_benchmark import mcce_env as mcenv


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


CLI_NAME = ENTRY_POINTS["launch"]  # as per pyproject.toml entry point


class ENTRY:
    def __init__(self):
        self.name = ""
        self.state = "i"

    def __str__(self):
        return f"{self.name:6s} {self.state:1s}"


class ENTRIES:
    """Collection of ENTRY instances to keep track
    of entry types count.
    """

    def __init__(self, bench_dir: Path):
        self.bench_dir = bench_dir
        # attributes are set by update_entries & update_status:
        self.all = None
        self.N = 0
        self.unsubmitted = 0
        self.running = 0
        self.completed = 0
        self.failed = 0

    def get_book_fp(self) -> Path:
        if Path.cwd().name == RUNS_DIR:
            book_fp = Path(Q_BOOK)
        elif Path.cwd().name == self.bench_dir.name:
            book_fp = self.bench_dir.joinpath(RUNS_DIR, Q_BOOK)
        else:
            book_fp = Path(RUNS_DIR).joinpath(Q_BOOK)
        return book_fp

    def update_entries(self):
        entries = []
        first_time = self.all is None

        book = self.get_book_fp()
        with open(book) as bk:
            for line in bk:
                rawtxt = line.strip().split("#")[0]
                # skip blank lines:
                if rawtxt == "":
                    continue
                fields = rawtxt.split()
                entry = ENTRY()
                entry.name = fields[0]
                if len(fields) > 1:
                    entry.state = fields[1].lower()
                entries.append(entry)

        self.all = entries
        if first_time:
            self.N = len(entries)
            self.unsubmitted = self.N

        return

    def update_counts(self):
        """Get count of each status type."""
        counts = Counter([e.state for e in self.all])
        self.unsubmitted = counts[" "] + counts["i"]
        self.running = counts["r"]
        self.completed = counts["c"]
        self.failed = counts["e"]

        return

    def tot_check(self) -> bool:
        return self.N == (
            self.unsubmitted + self.running + self.completed + self.failed
        )

    def __str__(self):
        return f"""
- Unsubmitted: {self.unsubmitted}
-     Running: {self.running}
-   Completed: {self.completed}
-      Failed: {self.failed}
-       Total: {self.N}
"""

    def final_counts(self):
        return f"""
- Completed: {self.completed}
-    Failed: {self.failed}
-     Total: {self.N}
"""


def read_book_entries(book: str = Q_BOOK) -> list:
    """Read book file data using ENTRY class.
    Return a list of entry instances.
    """
    entries = []
    with open(book) as bk:
        for line in bk:
            rawtxt = line.strip().split("#")[0]
            # skip blank lines:
            if rawtxt == "":
                continue
            fields = rawtxt.split()
            entry = ENTRY()
            entry.name = fields[0]
            if len(fields) > 1:
                entry.state = fields[1].lower()
            entries.append(entry)

    return entries


def get_running_jobs_dirs(job_name: str) -> list:
    """
    Query shell for user's processes with job_name.
    Return a list of runs/ sub-directories where the jobs are running.
    """
    # get the process IDs that match job_name from the user's running processes
    data = iou.subprocess_run(f"pgrep -u {USER} {job_name}.sh")
    if data is iou.subprocess.CalledProcessError:
        logger.error("Error with pgrep cmd.")
        raise data

    dirs = []
    for uid in data.stdout.splitlines():
        # get the current working directory of a process
        out = iou.subprocess_run(f"pwdx {uid}")
        if out is iou.subprocess.CalledProcessError:
            logger.error("Error with pwdx cmd.")
            raise out

        d = Path(out.stdout.split(":")[1].strip())
        dirs.append(d.name)

    return dirs


def batch_run(args: Union[dict, Namespace]) -> None:
    """
    Update Q_BOOK according to user's running jobs' states.
    Launch new jobs inside the runs subfolders until the number of
    job equals n_batch.
    To be run in /runs folder, which is where Q_BOOK resides.

    Args:
    args.job_name (str): Name of the job and script to use in /runs folder.
    args.n_batch (int, 10): Number of jobs/processes to maintain.
    args.sentinel_file (str, "pK.out"): File whose existence signals a completed job;
      When running all 4 MCCE steps (default), this file is 'pK.out', while
      when running only the first 2, this file is 'step2_out.pdb'.
    """
    if isinstance(args, dict):
        args = Namespace(**args)
    job_name = args.job_name
    job_script = f"{job_name}.sh"

    logger.info("Read book entries")
    # list of entry instances from Q_BOOK:
    entries = read_book_entries()
    running_jobs = get_running_jobs_dirs(job_name)
    n_jobs = len(running_jobs)
    logger.info(f"Running jobs: {n_jobs}")

    new_entries = []
    for entry in entries:
        if entry.state in [" ", "i"]:  # unsubmitted
            n_jobs += 1
            if n_jobs <= args.n_batch:
                iou.subprocess.Popen(
                    f"../{job_script}",
                    cwd=f"./{entry.name}",
                    close_fds=True,
                    stdout=open(f"./{entry.name}/run.log", "a"),
                    stderr=open(f"./{entry.name}/err.log", "a"),
                )
                entry.state = "r"
                logger.info(f"Running: {entry.name}")

        elif entry.state == "r":
            if entry.name not in running_jobs:
                # was running => completed or error
                entry.state = "e"
                sentin_fp = Path(entry.name).joinpath(args.sentinel_file)
                if sentin_fp.exists():
                    entry.state = "c"
                logger.info(f"Changed {entry.name}: 'r' -> {entry.state!r}")

        new_entries.append(entry)

    # update q-book
    with open(Q_BOOK, "w") as bk:
        bk.writelines([f"{e}\n" for e in new_entries])

    return


def launch_job(args: Namespace) -> None:
    """
    Go to bench_dir/runs directory & call batch_run.

    Args:
    Options from the command line:
      [required] bench_dir (Path, None): Path of the folder containing the /runs folder.
      job_name (str, None): Name of the job and script to use in /runs folder.
      n_batch (int, 10): Number of jobs/processes to maintain.
      sentinel_file (str, "pK.out"): File whose existence signals a completed step;
          When running all 4 MCCE steps (default), this file is 'pK.out', while
          when running only the first 2, this file is 'step2_out.pdb'.
    """
    bench_dir = Path(args.bench_dir)
    if Path.cwd().name != bench_dir.name:
        os.chdir(bench_dir)
    
    if args.job_name != DEFAULT_JOB:
        sh_name = f"{args.job_name}.sh"
        sh_path = bench_dir.joinpath(RUNS_DIR, sh_name)
        sentinel = iou.get_sentinel(sh_path)
        if sentinel != args.sentinel_file:
            # update current args:
            args.sentinel_file = sentinel

            # maybe update pickled:
            saved_senti = iou.get_setup_args_vals(bench_dir, "sentinel_file")
            if saved_senti:
                if saved_senti[0] != sentinel:
                    args_pk = iou.from_pickle(args.bench_dir.joinpath(FILES.CLI_ARGS_PKL.value))
                    setattr(args_pk, "prev_sentinel", sentinel)
                    args_pk.sentinel_file = sentinel
                    iou.to_pickle(args_pk, args.bench_dir.joinpath(FILES.CLI_ARGS_PKL.value))

                    # append change to cli arg txt file:
                    argstxt_fp = bench_dir.joinpath(FILES.CLI_ARGS_TXT.value)
                    with open(argstxt_fp, "a") as fp:
                        fp.write("\n" + "."*20 + f"\nCommand line options at submission:\n{vars(args)}\n")
    
    pbes = iou.get_pbe_solver_with_env(bench_dir)
    if pbes is not None:
        if pbes == "ZAP":
            # ZAP is currently the only solver with specific envir:
            oepath = mcenv.check_openeye_env_activated()
            if oepath is None:
                sys.exit("Activate a conda environment for OpenEye software.")
            _ = mcenv.get_openeye_license()

    os.chdir(RUNS_DIR)

    batch_run(args)

    os.chdir("../")

    return


def batch_parser():
    """Command line arguments parser for batch_submit.launch_job."""
    parser = ArgumentParser(
        prog=CLI_NAME,
        description="Cli for launching one batch of runs.",
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-bench_dir",
        required=True,
        type=arg_valid_dirpath,
        help="""The user's choice of directory for setting up the benchmarking run(s); this is where
        the /runs folder reside. The directory is created if it does not exists unless this cli is
        called within that directory.
        """,
    )
    parser.add_argument(
        "-job_name",
        type=str,
        default=DEFAULT_JOB,
        help="""The descriptive name, without spaces, for the current job; 
        Required for non-default jobs. Used to identify the shell script in 'bench_dir' that launches
        the MCCE simulation in 'bench_dir/runs/' subfolders; default: %(default)s.
        """,
    )
    parser.add_argument(
        "-n_batch",
        type=int,
        default=N_BATCH,
        help="The number of jobs to keep launching; default: %(default)s.",
    )
    parser.add_argument(
        "-sentinel_file",
        type=str,
        default="pK.out",
        help="File whose existence signals a completed step; default: %(default)s.",
    )

    return parser


def launch_cli(argv=None):
    """
    Command line interface for MCCE benchmarking entry point 'bench_batch'.
    Launch one batch of size args.n_batch (repeatedly if used by scheduling.py).
    """
    launch_parser = batch_parser()
    args = launch_parser.parse_args(argv)

    bench_dir = Path(args.bench_dir)
    if Path.cwd().name != bench_dir.name:
        os.chdir(bench_dir)

    entries = ENTRIES(bench_dir)
    # initial loading:
    entries.update_entries()
    entries.update_counts()
    book_fp = entries.get_book_fp()

    launch_job(args)

    entries.update_entries()
    entries.update_counts()
    pct = iou.pct_completed(book_fp)
    if pct == 1.0:
        logger.info(f"All jobs processed. Final counts:\n{entries.final_counts()}")
    else:
        logger.info(f"Percentage of jobs processed: {pct:.1%}")
        logger.info(
            f"{'If using bench_batch at the command line'.upper()}: Repeat until pct = 100%."
        )
    return


if __name__ == "__main__":
    launch_cli(sys.argv[1:])
