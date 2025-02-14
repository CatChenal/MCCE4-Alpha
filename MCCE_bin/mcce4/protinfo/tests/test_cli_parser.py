#!/usr/bin/env python

from mcce4.protinfo.cli import pi_parser
import pytest


class TestPiParser:
    # Verify that the parser correctly initializes with the program name 'ProtInfo'
    def test_parser_initialization(self):
        parser = pi_parser()
        assert (
            parser.prog == "ProtInfo"
        ), "The program name should be initialized as 'ProtInfo'"

    def test_fetch_only(self):
        parser = pi_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--fetch"])

    # Test with empty string input for 'pdb' to see if it is handled or raises an error
    def test_empty_pdb_input(self):
        parser = pi_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    # Test with input strings significantly longer than typical use cases
    def test_long_input_strings_handling(self):
        parser = pi_parser()
        long_input = "A" * 100 + ".pdb"
        args = parser.parse_args([long_input])
        assert (
            args.pdb == long_input
        ), "The parser should handle extremely long input strings for 'pdb' correctly"
