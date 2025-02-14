

import argparse
from pathlib import Path
import shlex
from mcce4.mcce_benchmark import cli_utils as cliu
from mcce4.mcce_benchmark.cli import bench_parser


def test_bench_parser_fromfile():
    opt_fp = Path("bench_setup_options.txt")
    cliu.write_default_setup_options_file(opt_fp)
    print(opt_fp.read_text())

    argv = shlex.split("pkdb_pdbs -bench_dir tests @bench_setup_options.txt")

    parser = bench_parser()
    args = parser.parse_args(argv)
    print(vars(args))
