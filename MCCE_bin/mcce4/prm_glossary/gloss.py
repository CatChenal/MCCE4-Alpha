#!/usr/bin/env python

"""
Module: gloss.py

Holds classes to create and query a glossary of MCCE parameters.

CHANGELOG:
* 11-08-24:
  The run.prm file used is MCCE4/runprms/run.prm.default
"""

import argparse
from collections import defaultdict
from functools import cache
from pathlib import Path
import sys
from typing import Any, List, Tuple, Union


runprm = Path(__file__).parent.joinpath("gloss.run.prm.full")


class TrieNode:
    def __init__(self):
        self.children = {}
        self.isEndOfWord = False
        self.data = None  # Container for custom data


class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word: str, data: Any = None):
        node = self.root
        for char in word:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.isEndOfWord = True
        node.data = data

    def searchPrefix(self, prefix):
        node = self.root
        for char in prefix:
            if char in node.children:
                node = node.children[char]
            else:
                return None
        return node

    def _findAllWords(self, node, prefix, words: List[Tuple]):
        if node.isEndOfWord:
            words.append((prefix, node.data))
        for char, nextNode in node.children.items():
            self._findAllWords(nextNode, prefix + char, words)

    def StartsWith(self, prefix):
        node = self.searchPrefix(prefix)
        words = []
        if node:
            self._findAllWords(node, prefix, words)
        if not words:
            print(f"Not found: {prefix}")
            return None
        return words


class EnvGlossary:
    def __init__(self):
        self.glossary = defaultdict(list)
        self.query = None

    def _setup_search(self):
        self.query = Trie()
        for k, v in self.glossary.items():
            self.query.insert(k, {"Description": v[0], "Default": v[1]})
            self.query.insert(self.glossary[k][0], {"Key": k, "Default": v[1]})
        return

    @cache
    def load_runprm(self, filename: str = runprm):
        """Format in run.prm:  default value - explanation - (KEY)  # maybe a comment 
        """
        with open(filename) as fh:
            for line in fh:
                if not line:
                    continue
                entry_str = line.strip().split("#")[0]
                if not entry_str:
                    continue
                try:
                    line_start, key_fld = entry_str.rsplit(maxsplit=1)
                except ValueError:
                    continue
                if key_fld[0] == "(" and key_fld[-1] == ")":
                    key = key_fld.strip("()").strip()
                    def_val, desc = line_start.split(maxsplit=1)
                    self.glossary[key].extend([desc, def_val])

        self._setup_search()
        return

    def Search(self, prefix, list_result: bool = True) -> Union[list, None]:
        result = self.query.StartsWith(prefix)
        if result is None:
            return
        else:
            if list_result:
                if prefix:
                    print(f"Result of query: {prefix!r}\n")
                for pref, d in result:
                    detail = ""
                    for k, v in d.items():
                        detail = detail + f"  {k}: {v}\n"
                    print(f"{pref}\n{detail}")
            else:
                return result

    def __str__(self) -> str:
        if not self.glossary:
            return "Attribute 'glossary' is not set."
        return "\n".join(
            f"{k}\n   Description: {v[0]}\n   Default: {v[1]}" for k, v in self.glossary.items()
        )


def cli_parser():
    p = argparse.ArgumentParser(
        description="""
    Query a glossary of MCCE parameters (based on run.prm.full).
    Notes:
        - Search is case-sensitive!
        - To query for a specific step, the query should be 'StepN', e.g. Step1.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Post an issue for all errors and feature requests at:
    https://github.com/GunnerLab/MCCE4/issues
    """,
    )
    p.add_argument(
        "query",
        type=str,
        nargs="?",
        default="",
        help="Query the glossary with the given prefix string; Note: search is case-sensitive!.",
    )
    p.add_argument(
        "--print",
        default=False,
        action="store_true",
        help="Print the entire glossary in usual query result format.",
    )

    return p


def gloss_cli(argv=None):
    """Cli 'main' function: creates a searchable glossary of MCCE parameters."""
    parser = cli_parser()
    args = parser.parse_args(argv)

    envgloss = EnvGlossary()
    envgloss.load_runprm()

    if args.query:
        envgloss.Search(args.query)
    elif args.print:
        print(envgloss)


if __name__ == "__main__":
    gloss_cli(sys.argv)
