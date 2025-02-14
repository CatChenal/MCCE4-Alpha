#!/usr/bin/env python

from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory
import pytest
from mcce4.mcce_benchmark import RUNS_DIR
from mcce4.mcce_benchmark.custom_sh import (
    all_default_opts,
    all_opts_are_defaults, 
    cli_args_to_dict,
    populate_custom_template,
    ScriptChoices,
    write_run_script_from_template,
)


class TestPopulateCustomTemplate:

    def test_return_empty_dictionary_if_no_mcce_steps_args_present(self):
        # ns with w/o step args
        sh_args = Namespace(
            subparser_name="subparser",
            bench_dir="bench",
            n_pdbs=10,
            pdbs_list=["pdb1"],
            sentinel_file="sentinel",
            job_name="job",
            launch=True,
            func="func",
            help="help",
        )

        result = cli_args_to_dict(sh_args)
        assert result == {}

    def test_valid_input(self):

        job_args = Namespace(dry=False)
        #step2.py {conf_making_level}{d}{s2_norun}{e}{u}
        #{pound}step3.py {d}{c}{s}{t}{p}{ftpl}{salt}{refresh}{vdw}{fly}{debug}{l}
        #step4.py --xts {titr_type}{i}{interval}{n}{ms}{s4_norun}{e}{u}   
        expected_output = """#!/bin/bash

step1.py prot.pdb 
step2.py 
step3.py 
step4.py --xts 

sleep 10
"""
        result = populate_custom_template(job_args)

        assert result == expected_output

    def test_empty_namespace_input(self):
        with pytest.raises(ValueError):
            job_args = Namespace()
            populate_custom_template(job_args)

    def test_none_input(self):
        with pytest.raises(TypeError):
            job_args = None
            populate_custom_template(job_args)

    def test_special_case_arguments(self):
        job_args = Namespace(s3_norun=True)
        # note: trailing spaces in all step lines
        expected_output = """#!/bin/bash

step1.py prot.pdb 
step2.py 
# step3.py 
step4.py --xts 

sleep 10
"""
        result = populate_custom_template(job_args)
        assert result == expected_output


class TestWriteRunScriptFromTemplate:

    def test_all_opts_are_defaults_with_default_args(self):

        # Create a Namespace with default options
        default_args = Namespace(dry=False, noter=False, s1_norun=False, s2_norun=False, s3_norun=False, s4_norun=False)

        # Update all_default_opts to match the default_args
        all_default_opts.update(cli_args_to_dict(default_args))

        # Test
        assert all_opts_are_defaults(default_args) == True

    def test_all_opts_are_defaults_with_empty_namespace(self):

        # Create an empty Namespace object
        empty_args = Namespace()
    
        # Test
        assert all_opts_are_defaults(empty_args) == True

    def test_custom_script_creation(self):
        """The function successfully writes a custom shell script in
        <bench_dir>/runs/ when script_template is CUSTOM and job_args is not None.
        """

        tmp_benchdir = TemporaryDirectory(prefix="TMP_")
        tmp_runsdir = Path(tmp_benchdir.name).joinpath(RUNS_DIR)
        tmp_runsdir.mkdir()
        job_name = "custom"

        script_template = ScriptChoices.CUSTOM
        job_args = Namespace(dry=False)
        write_run_script_from_template(
            tmp_benchdir.name, job_name, script_template, job_args
        )

        expected_path = tmp_runsdir.joinpath("custom.sh")
        assert expected_path.exists()

        expected_content = """#!/bin/bash

step1.py prot.pdb 
step2.py 
step3.py 
step4.py --xts 

sleep 10
"""
        with open(expected_path, "r") as fh:
            assert fh.read() == expected_content

        tmp_benchdir.cleanup()
        assert not Path(tmp_benchdir.name).exists()
        assert not tmp_runsdir.exists()
