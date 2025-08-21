#!/usr/bin/env python
"""
Module: scheduling.py

For automating the crontab creation for scheduling batch_submit every minute.
"""
from argparse import Namespace
import logging
from pathlib import Path
import subprocess
from typing import Union
from mcce4.mcce_benchmark import USER_MCCE, LAUNCHJOB
from mcce4.mcce_benchmark import io_utils as iou
from mcce4.mcce_benchmark import mcce_env as mcenv


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def create_single_crontab(args: Namespace, debug: bool = False) -> Union[None, str]:
    """
    Create a crontab entry.
    If debug: return crontab_text w/o creating the crontab.
    """
    bdir = str(args.bench_dir)
    tools_dir = str(Path(LAUNCHJOB).parent)  # MCCE_bin/

    # get_pbe_solver_with_env returns the PBE solvers KNOWN to need a special environment
    pbes = iou.get_pbe_solver_with_env(bdir)
    if pbes is None or pbes != "ZAP":
        PATH = "PATH={}:{}:/usr/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin/env:/sbin:/bin\n"
        ct_text = PATH.format(USER_MCCE, tools_dir)

        SCHED = "* * * * * {} -bench_dir {} -job_name {} -n_batch {} -sentinel_file {}"
        ct_text = ct_text + SCHED.format(
            LAUNCHJOB,
            bdir,
            args.job_name,
            args.n_batch,
            args.sentinel_file,
        )
    else:
        # zap
        conda_exe, conda_env_bin, conda_sh, oe_licence, oe_env = mcenv.get_zap_env()
        conda_path = str(Path(conda_exe).parent)
        zap_license = f"OE_LICENSE={oe_licence}\n"

        PATH = "PATH={}:{}:{}:{}:/usr/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin\n"
        ct_text = zap_license + PATH.format(
            USER_MCCE, tools_dir, conda_path, conda_env_bin
        )

        SCHED = "* * * * * . {} && conda activate {} && {} -bench_dir {} -job_name {} -n_batch {} -sentinel_file {}"
        ct_text = ct_text + SCHED.format(
            conda_sh,
            oe_env,  # for activate
            LAUNCHJOB,
            bdir,
            args.job_name,
            args.n_batch,
            args.sentinel_file,
        )

    crontab_txt = f"{ct_text} >> {bdir}/cron.log 2>&1\n"
    if not debug:
        cron_in = subprocess.Popen(
            ["crontab", "-l"],
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        cur_crontab, _ = cron_in.communicate()
        if len(cur_crontab) and "PATH" in cur_crontab:
            crontab_txt = cur_crontab + "\n" + crontab_txt
        cron_out = subprocess.Popen(["crontab", "-"], stdin=subprocess.PIPE)
        cron_out.communicate(input=bytes(crontab_txt, "utf-8"))

        logger.info("Created crontab:\n```\n%s\n```" % crontab_txt)
        return

    return crontab_txt


def schedule_job(launch_args: Namespace) -> None:
    """Create a contab entry for batch_submit.py with args from the
    `bench_setup launch` command.
    """
    create_single_crontab(launch_args)
    logger.info("Scheduled batch submission with crontab every minute.")

    return
