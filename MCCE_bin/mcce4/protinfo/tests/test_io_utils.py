#!/usr/bin/env python

import shutil
from mcce4 import os, Path
from mcce4.protinfo.io_utils import get_rcsb_pdb


tmp_dir = Path.cwd().joinpath("tmp_pinfo")


class TestGetRcsbPdb:
    tmp_path = tmp_dir

    @classmethod
    def setup_class(cls):
        """setup any state specific to the execution of the given class"""
        if not cls.tmp_path.exists():
            cls.tmp_path.mkdir()
        cls.here = Path.cwd()

    @classmethod
    def teardown_class(cls):
        """teardown any state that was previously setup with a call to
        setup_class.
        """
        for dp in cls.tmp_path.iterdir():
            shutil.rmtree(dp)
        cls.tmp_path.rmdir()

    def test_valid_pdbid(self):
        """Given a valid pdbid, the function should download the pdb
        file containing the biological assembly from rcsb.org and
        return a Path object."""

        pdbid = "4lzt"
        d = self.tmp_path / pdbid
        if not d.is_dir():
            d.mkdir()
        os.chdir(d)
        result = get_rcsb_pdb(pdbid)
        os.chdir(self.here)

        assert isinstance(result, Path)

    def test_invalid_pdbid(self):
        """Given an invalid pdbid, the function should return a
        tuple: Tuple[None, str]."""

        pdbid = "xyz"
        d = self.tmp_path / pdbid
        if not d.is_dir():
            d.mkdir()
        os.chdir(d)

        result = None
        result = get_rcsb_pdb(pdbid)

        os.chdir(self.here)

        assert result[0] is None
        assert (
            result[1] == "Error: Could neither download the bio assembly or pdb file."
        )

    def test_overwrite_existing_file(self):
        """Given a pdbid that already exists in the current directory,
        the function should overwrite the existing file with the new download.
        """

        pdbid = "1ans"
        d = self.tmp_path / pdbid
        if not d.is_dir():
            d.mkdir()

        existing_file = d.joinpath(f"{pdbid}.pdb")  # .resolve()
        existing_file.touch()
        assert existing_file.exists()

        os.chdir(d)
        result = get_rcsb_pdb(pdbid)
        os.chdir(self.here)

        assert isinstance(result, Path)
        assert result == existing_file
