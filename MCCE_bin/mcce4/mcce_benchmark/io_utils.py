#!/usr/bin/env python

"""
Module: io_utils

Module with generic (enough) functions related to loading and saving files.
"""
import argparse
import json
import logging
from pathlib import Path
import pickle
from typing import Any, List, Tuple, Union
import pandas as pd
from mcce4.io_utils import subprocess, subprocess_run
from mcce4.mcce_benchmark import FILES


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def Pathok(
    pathname: str, check_fn: str = "exists", raise_err=True
) -> Union[Path, bool]:
    """Return path if check passed, else raise error.
    check_fn: one of 'exists', 'is_dir', 'is_file'.
    if raise_err=False, return False instead of err.
    """

    pathname = Path(pathname)
    if check_fn not in ["exists", "is_dir", "is_file"]:
        check_fn = "exists"

    # failure msg:
    if check_fn == "exists":
        msg = f"Path not found: {pathname}"
    elif check_fn == "is_dir":
        msg = f"Directory not found: {pathname}"
    elif check_fn == "is_file":
        msg = f"Path is not a file: {pathname}"

    if not pathname.__getattribute__(check_fn)():
        if not raise_err:
            return False
        raise FileNotFoundError(msg)

    return pathname


def pct_completed(book_fpath: str) -> Union[float, None]:
    """Return the pct of runs that are completed or finished with error."""
    book_fp = Pathok(book_fpath)
    # 2 cmds: n processed; total count => 2 lines, e.g.:
    # 1
    # 6
    cmd = f"grep '[ce]' {book_fp} | wc -l; cat {book_fp} | wc -l"
    logger.info(f"{cmd = }")
    data = subprocess_run(cmd)

    if isinstance(data, subprocess.SubprocessError):
        logger.error("Error fetching pct completed.")
        return None

    if not data.stdout.strip():
        logger.error("No data from book file.")
        return None

    logger.info(f"{data.stdout = }")

    out = data.stdout.splitlines()
    logger.info(f"out = data.stdout.splitlines() :: {out = }")
    try:
        n = int(out[0])
        logger.info(f"{n = }")
    except ValueError:
        logger.error("No completed count from book file.")
        return None
    try:
        tot = int(out[1])
        logger.info(f"{tot = }")
    except ValueError:
        logger.error("No line count from book file.")
        return None

    pct = n / tot
    logger.info(f"pct processed: {pct:.1%}")

    return pct


def make_executable(sh_path: str) -> None:
    """Alternative to os.chmod(sh_path, stat.S_IXUSR): permission denied."""

    sh_path = Pathok(sh_path)
    cmd = f"chmod +x {str(sh_path)}"

    p = subprocess_run(cmd, capture_output=False, check=True)
    if isinstance(p, subprocess.CalledProcessError):
        logger.exception(f"Error in subprocess cmd 'chmod +x':\nException: {p}")
        raise

    return


def df_to_txt(df: pd.DataFrame, fpath: str, replace: bool = True):
    """Save a pandas.DataFrame as a text file."""
    fp = Path(fpath)
    if fp.exists():
        if not replace:
            logger.error(f"File already exists: {fp}; {replace = }")
            return None
        fp.unlink()
    fp.write_text(df.to_string(index=False) + "\n")


def txt_to_df(fp: Path, header=None) -> pd.DataFrame:
    """Load a fixed-width file into a pandas.DataFrame."""
    return pd.read_fwf(fp, header=header)


# aliases
df2txt = df_to_txt
txt2df = txt_to_df


def get_file_header(fp: str) -> str:
    hdr = ""
    with open(fp) as f:
        for i, line in enumerate(f):
            if i > 0:
                break
            hdr = line

    return hdr


def get_book_dirs_for_status(book_fpath: str, status: str = "c") -> list:
    """Return a list of folder names from book_fp, the Q_BOOK file path,
    if their status codes match a 'completion status', i.e. completed ('c', default),
    or errorneous ('e').
    """

    status = status.lower()
    if not status or status not in ["c", "e"]:
        logger.error("Invalid 'status'; choices are 'c' or 'e'.")
        raise ValueError("Invalid 'status'; choices are 'c' or 'e'.")

    book_fp = Pathok(book_fpath)
    book_dirs = []
    with open(book_fp) as book:
        for line in book:
            # select the portion preceding any appended comment
            rawtxt = line.strip().split("#")[0]
            fields = rawtxt.split()
            if len(fields) == 2:
                if fields[1].lower() == status:
                    book_dirs.append(fields[0])

    return book_dirs


def to_pickle(obj: Any, fp: str):
    pickle.dump(obj, open(fp, "wb"))
    return


def from_pickle(fp: str) -> Any:
    obj = pickle.load(open(fp, "rb"))
    return obj


def to_json(obj: Any, filepath: Path):
    """Serialize obj as a JSON formatted stream to filepath."""
    if filepath.suffix != ".json":
        filepath = filepath.parent.joinpath(filepath.stem + ".json")
        logger.warning(f"File extension replaced; file is: {str(filepath)!r}")
    with open(filepath, "w") as fh:
        json.dump(obj, fh)
    return


def from_json(fp: Path) -> Any:
    """Deserialize a json file given its filepath."""
    with open(fp) as fh:
        out = json.load(fh)
    return out


def get_setup_args_vals(dirpath: Path, option: Union[str, List[str]]) -> list:
    """Retrieve one or more values from the pickled bench_setup cli args.
    Args:
     - dirpath (Path): benchmark dir path
     - option (str, list): A single cli option, if type is string, or multiple options to get
       from the file.
    """
    out = []
    is_string = isinstance(option, str)

    # get val from saved setup args:
    setup_args_fp = Path(dirpath).joinpath(FILES.CLI_ARGS_PKL.value)
    if not setup_args_fp.exists():
        logger.warning("Saved setup args file not found.")
        return out

    # unpickle the Namespace instance:
    setup_d = vars(from_pickle(setup_args_fp))
    if is_string:
        out = [setup_d.get(option)]
    else:
        out = [setup_d.get(opt) for opt in option]

    return out


def get_pbe_solver_with_env(dirpath: Path) -> Union[None, str]:
    """Get the pdb solver string from the cli options.
    Return either 'NGPB' or 'ZAP' as these solvers need
    special environment, or None for other cases.
    Notes:
    `get_pbe_solver_with_env` returns the PBE solvers KNOWN to need a special environment;
    currently, only NGPB and ZAP need special environment.
    `get_pbe_solver_with_env` returns only the PBE solvers KNOWN to need a special environment.
    """
    vals = get_setup_args_vals(dirpath, "s")
    if not vals:
        # error? nothing to do
        return None
    if vals[0] is None:
        return None

    pbe_solver = vals[0].upper()
    if pbe_solver in ["NGPB", "ZAP"]:
        return pbe_solver
    else:
        return None


def matched_pks_to_df(matched_fp: str) -> pd.DataFrame:
    """
    Load MATCHED_PKAS file into a pandas DataFrame.
    PRE: file MATCHED_PKAS file created via mcce_benchmark.pkanalysis.
    Columns: resid; <set1 pk>; <set2 pk>. The last 2 columns are
    dynamically named.
    """
    fp = Path(matched_fp)
    fname = FILES.MATCHED_PKAS_TXT.value
    if fp.name != fname:
        logger.error(f"Only {fname} is a valid file name.")
        raise ValueError(f"Only {fname} is a valid file name.")

    if not fp.exists():
        logger.error(f"Not found: {fp}; run pkanalysis to create.")
        raise FileNotFoundError(f"Not found: {fp}; run pkanalysis to create.")

    return txt2df(fp, header=0)


# 1 :: aka quick
levels2names = {1: "isosteric", 2: "medium", 3: "comprehensive"}


def get_level_tup(level: Union[int, Tuple[int, int]] = None) -> Tuple:
    if level is None:
        return (None, "unknown")

    if isinstance(level, int):
        lev = levels2names.get(level)
        return (None, "unknown") if lev is None else (level, lev)

    if isinstance(level, tuple):
        l1, l2 = level
        lev1 = levels2names.get(l1)
        lt1 = (None, "unknown") if lev1 is None else (l1, lev1)
        lev2 = levels2names.get(l2)
        lt2 = (None, "unknown") if lev2 is None else (l2, lev2)
        return lt1, lt2


def sentinel_from_args(args: Union[argparse.Namespace, dict]) -> Union[str, None]:
    if isinstance(args, dict):
        args = argparse.Namespace(**args)

    if args.s4_norun:
        if args.s3_norun:
            if args.s2_norun:
                if args.s1_norun:
                    # cannot have args.s1_norun as there is no sentinel_file
                    # associated with "run no steps":
                    sentinel = None
                else:
                    sentinel = "step1_out.pdb"
            else:
                sentinel = "step2_out.pdb"
        else:
            sentinel = "head3.lst"
    else:
        sentinel = "pK.out"

    return sentinel


def sh_steps(script_fp: Path) -> dict:
    """Recover the values for the --sx_norun options for each steps in the run script.
    Only mcce steps 1 through 4 are considered.
    """
    sh_lines = [
    line.strip()
    for line in script_fp.read_text().splitlines()[1:-1]
    if line.strip()
    ]
    sh_d = {}  # will need to have 4 keys
    seen = set()
    for line in sh_lines:
        try:
            str1, str2 = line.split("#")
            commented = True
        except ValueError:
            # No hash in line:
            str2 = line
            commented = False
        finally:
            idx = str2.find("step")
            if idx == -1:  # not a step line
                continue

            s = str2[idx : idx + 5][-1]
            if s not in "1234":
                continue
            seen.add(s)
            if "norun" in line or commented:
                sh_d[f"s{s}_norun"] = True
            else:
                sh_d[f"s{s}_norun"] = False

    if len(seen) != 4:
        missing = set("1234").difference(seen)
        for m in missing:
            print("missing, m:", m)
            sh_d[f"s{m}_norun"] = True

    return sh_d


def get_sentinel(script_fp: Path) -> Union[str, None]:
    """Return the name of the sentinel_file from the custom script path.
    """
    runsh_steps = sh_steps(script_fp)
    sentinel = sentinel_from_args(runsh_steps)

    return sentinel


def update_pkl_args_from_sh(args: argparse.Namespace, new_values: dict):
    """Update pickled args with new values for these keys:
     -sentinel_file
     -si_norun
    """
    # update pickled:
    pkl_fp = args.bench_dir.joinpath(FILES.CLI_ARGS_PKL.value)
    args_pkl = from_pickle(pkl_fp)
    args_d = vars(args_pkl)
    for k in new_values:
        if k == "sentinel_file":
            if args_d.get("prev_sentinel") is None:
                args_d["prev_sentinel"] = args_d[k]
            else:
                args_d["prev_sentinel"] = args_d[k]
            args_d[k] = new_values.get(k)
        else:
            if args_d.get(k) is None:
                args_d[k] = new_values.get(k)
            else:
                args_d[k] = new_values.get(k)

    to_pickle(argparse.Namespace(**args_d), pkl_fp)

    return
