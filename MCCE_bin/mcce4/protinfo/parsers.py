#!/usr/bin/env python

"""
Module: parsers.py

Functions and classes to process information from the input pdb and
from mcce step1 run1.log.
For this app purpose, only the output from step1 is considered, which
is handled by the class RunLog1.
"""

from argparse import Namespace
from collections import defaultdict
from dataclasses import dataclass
import logging
from pathlib import Path
from pprint import pformat
from time import sleep
from typing import Tuple, Union

import pandas as pd

from mcce4 import pdbio
from mcce4.protinfo import MCCE4,RPT, RUN1_LOG, USER_MCCE
from mcce4.protinfo import run
from mcce4.protinfo.io_utils import ENV, get_path_keys, retry


logger = logging.getLogger(__name__)


WARN_MULTI_MODELS = ("MCCE cannot handle multi-model proteins; Model 1 is parsed by default "
                     "to obtain basic information, but MCCE4 pdb loader can only handle single "
                     "model pdbs."
                     )
WARN_MALFORMED_PDB = ("MCCE could not parse the pdb into at least one model, possibly "
                      "due to a missing MODEL line."
                      )
MSG_KEEP_H2O = ("NOTE: Include the '--wet' option at the command line to keep buried "
                "waters and cofactors. Alternatively, change the water SAS cutoff to a "
                "non-zero, positive number using the command line 'u' option:\n"
                "  > protinfo 1fat.pdb -u H2O_SASCUTOFF=0.05"
                )

# kept from when biopython was used to output buried res;
# may be useful if/when acc.files are parsed.
#BURIED_THRESH = 0.05  # mcce default; res with sasa < this are buried.
#BURIED_THR_MSG = f"(using default mcce SASA threshold of {BURIED_THRESH:.0%}):\n"
# This ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7817970/
# has benchmarked a SASA threshold of 20%.


# pdbio parser ........................................
def info_input_prot(pdb: Path) -> dict:
    """Return information about 'pdb' from mcce4.pdbio."""
    structure = pdbio.Structure()
    # load with MCCE4 renaming rules file
    name_rules = MCCE4.joinpath("name.txt")
    structure.load_pdb(pdb, str(name_rules))

    info_d = structure.get_prerun_dict()
    if info_d["PDB.Structure"].get("Models") is not None:
        if info_d["PDB.Structure"]["Models"] == "0":
            info_d["PDB.Structure"]["Malformed PDB"] = WARN_MALFORMED_PDB
        elif info_d["PDB.Structure"]["Models"] != "1":
            info_d["PDB.Structure"]["MultiModels"] = WARN_MULTI_MODELS

    return info_d


# run1.log parser ........................................
def get_pubchem_compound_link(compound_id: str) -> str:
    """Return the link of the PubChem subtance tab for compound_id.
    The link is prepended with ": " as it will follow compound_id
    in the report line; this saves checking for empty str.
    """
    if compound_id:
        url_fstr = ": https://pubchem.ncbi.nlm.nih.gov/#query={}&tab=substance"
        if compound_id.startswith("_"):
            return url_fstr.format(compound_id[1:])
        else:
            return url_fstr.format(compound_id)
    else:
        return ""


def extract_content_between_tags(
    text: str, tag1: str, tag2: str = "   Done"
) -> Union[str, None]:
    """Extracts the content between two string tags in a text.
    Args:
      text: A text.
      tag1: The first tag.
      tag2: The second tag, default: "   Done".

    Returns:
      A string containing the content between the two tags if both found;
      A string containing the content from tag1 if tag2 not found;
      None if tag1 is not found.
    """
    try:
        pos1 = text.find(tag1)
        if pos1 == -1:
            return None
    except AttributeError:
        return None

    start_pos = pos1 + len(tag1)
    end_pos = text.find(tag2, start_pos)
    if end_pos != -1:
        return text[start_pos:end_pos]
    else:
        return text[start_pos:]


@dataclass
class LogSection:
    """Dataclass to store information and line transformations for
    each 'processing section' returned in run1.log after step1 has run.
    Each processing section starts with a header line (hdr) and ends
    with a common '   Done' line.

    Attributes:
      idx (int, start=1): Index of the section in order of appearance
      hdr (str): header line
      rpt_hdr (str): Corresponding header in the report
      line_start (Union[None, str]): Remove the substring from the start of the line
      skip_lines (Union[None, list, tuple]): If list: lines to skip; if tuple:
                                             skip a line if substr in line.
      debuglog (bool): True if a line in the given section mentions 'debug.log',
                       the debug.log file will then be parsed.
    """
    idx: int
    hdr: str
    rpt_hdr: str
    line_start: Union[None, str] = None
    skip_lines: Union[None, list, tuple] = None
    debuglog: bool = False
    get_key: str = None

    def has_debuglog(self, line: str, calling_key: int = 5) -> bool:
        """Set debuglog proprerty to True if line ends with
        'saved in debug.log.'
        calling_key is meant to be the index key of the calling dict.
        The only known case where debug.log is mentioned is in 'block 5'
        of step1.py: 'free cofactor stripping' (see runlog_headers list).
        '"""
        # condition to be removed if debug.log found in other steps:
        if calling_key != 5:
            return

        if not self.debuglog:
            self.debuglog = line.endswith("saved in debug.log.")
        return


def get_log1_specs(pdb: Path) -> dict:
    """Return a dict of LogSection classes for processing step1 sections in run1.log.
    """
    # list of run1.log headers returned by mcce step1:
    runlog1_headers = [
        "   Rename residue and atom names...",
        "   Identify NTR and CTR...",
        "   Label backbone, sidechain and altLoc conformers...",
        "   Load pdb lines into data structure...",
        # 5: dynamic header
        "   Strip free cofactors with SAS >  {: .0%}...",
        "   Check missing heavy atoms and complete altLoc conformers...",
        # 7: dynamic header
        "   Find distance clash (<{:.3f})...",
        "   Make connectivity network ...",
    ]

    all = defaultdict(dict)
    for i, hdr in enumerate(runlog1_headers, start=1):
        if i == 1:
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Renamed",
                line_start="   Renaming ",
            )
        elif i == 2:
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Termini",
                line_start="      Labeling ",
            )
        elif i == 3:
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Labeling",
                line_start="      Labeling ",
                skip_lines=(
                    "Creating temporary parameter file for unrecognized",
                    "Trying labeling again",
                    "Try delete this entry and run MCCE again",
                    "Error! premcce_confname()",
                    # "STOP",
                    "is already loaded somewhere else.",
                ),
            )
        elif i == 4:
            # keep as is until error found
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Load Structure",
            )
        elif i == 5:
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Free Cofactors",
                skip_lines=(
                    "free cofactors were stripped off in this round",
                    "saved in debug.log.",
                ),
                get_key="H2O_SASCUTOFF",
            )
        elif i == 6:
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Missing Heavy Atoms",
                line_start="   Missing heavy atom  ",
                skip_lines=["   Missing heavy atoms detected."],
            )
        elif i == 7:
            all[i] = LogSection(
                i, hdr, rpt_hdr="Distance Clashes", get_key="CLASH_DISTANCE"
            )
        elif i == 8:
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Connectivity",
            )
        else:
            # unknown
            all[i] = LogSection(
                i,
                hdr,
                rpt_hdr="Other",
            )

    return dict(all)


class RunLog1Parser:
    """A class to parse mcce run1.log into sections pertaining to step1, and
    process each one of them into a simplified output.
    """
    def __init__(self, pdb: Path) -> None:
        self.pdb = pdb
        self.pdbid = self.pdb.stem
        self.s1_dir = self.pdb.parent
        # id of block with debug.log mentions, if any:
        self.check_debuglog_idx = [5]
        self.blocks_specs = get_log1_specs(pdb)
        self.txt_blocks = None
        self.runprm = None  # self.get_runprm()
        self.dry_opt = None  # float(self.runprm["H2O_SASCUTOFF"]) == -0.01

    def get_runprm(self) -> dict:
        env = ENV(self.s1_dir)

        return env.runprm

    @retry()   # 6 times: default
    def get_log_contents(self) -> str:
        """Read entire contents of run1.log file."""
        text = ""
        log_fp = self.s1_dir.joinpath(RUN1_LOG)
        if not log_fp.exists:
            sleep(1)
            raise ValueError("Step1 log not yet created.")

        text = log_fp.read_text()
        if not text:
            sleep(1)
            raise ValueError("Step1 log is empty.")

        last_line = "Step 1 Done."
        last_found = text.find(last_line)
        if last_found == -1:
            logger.error("Step1 log missing completion statement: run ended with error.")
        
        return text

    def get_debuglog_species(self) -> list:
        """Summarize info in debug.log.
        Note:
        Parsing of get_connect12 warnings is likely not exhaustive;
        Are all going to debug.log?
        Also found in triplicate:
        '''
            Error! get_connect12(): connectivity of atom " CA  ASPBK A0513" is not complete
                   get_connect12(): atom  CB  in the same residue is not found
        '''
        """
        dbgl_props_vdw = []
        dbgl_props_tor = []
        dbgl_empty_conect = []
        empty_slot = "   Warning! get_connect12(): An empty ligand connectivity slot found for atom"

        fp = self.s1_dir.joinpath("debug.log")
        for line in fp.read_text().splitlines():
            if not line:
                continue
            if line.endswith("not put in the connectivity list "):
                # only pertains to free cofactors? skip
                continue
            if line.startswith(("    Error! get_connect12(): conn",
                                "           get_connect12(): atom")):
                continue

            if not line.startswith("   Warning"):
                if line.startswith("TORSION"):
                    dbgl_props_tor.append(line.split())
                else:
                    dbgl_props_vdw.append(line.split())
            else:
                if line.startswith(empty_slot):
                    dbgl_empty_conect.append(
                        line.removeprefix(empty_slot).strip().split(" in residue ")
                    )
                else:
                    # TODO: handle "Error! get_connect12()..."
                    continue

        # populate output list:  # TODO? output dicts
        if dbgl_props_vdw:
            out = [
                "Species and properties with assigned default values in debug.log:\n"
            ]
            df = pd.DataFrame(dbgl_props_vdw)
            for k in df[1].unique():
                out.append(f"{k}: {list(df[df[1] == k][0].unique())}\n")
        if dbgl_props_tor:
            df = pd.DataFrame(dbgl_props_tor)
            for k in df[1].unique():
                out.append(f"{k}: {list(df[df[1] == k][0].unique())}\n")
        if dbgl_empty_conect:
            df = pd.DataFrame(dbgl_empty_conect)
            df_uniq = df[0].unique()
            if df_uniq.shape[0]:
                out.append("Empty connection slot(s):\n")
                for k in df_uniq:
                    out.append(f"{k.strip()}: {list(df[df[0] == k][1].unique())}\n")

        return out

    def process_content_block(self, content: list, log_section: LogSection) -> list:
        out = []
        skip = log_section.skip_lines is not None
        change = log_section.line_start is not None
        newtpl = None
        tpl_mismatch = None
        tpl_mismatch_atoms = None
        tpl_err = "   Error! The following atoms of residue "
        
        # section that lists new.tpl creation:
        if log_section.idx == 3:
            newtpl = ""

        for line in content:
            if not line:
                continue

            if log_section.idx == 3:
                if line.startswith("   Error! premcce_confname()"):
                    # add conf name & link:
                    conf = line.rsplit(maxsplit=1)[1]
                    newtpl += f"{conf}{get_pubchem_compound_link(conf)}; "

                elif line.startswith(tpl_err):
                    if tpl_mismatch is None:
                        tpl_mismatch = defaultdict(list)
                    if tpl_mismatch_atoms is None:
                        tpl_mismatch_atoms = defaultdict(set)

                    #   Error! The following atoms of residue ASN A  65 can not be loaded to conformer type ASN01
                    res_info, tpl_conf = line.removeprefix(tpl_err).split(
                        " can not be loaded to conformer type "
                    )
                    _, resloc = res_info.split(maxsplit=1)
                    resloc = resloc.strip().replace("  ", " ")
                    tpl_mismatch[tpl_conf.strip()].append(resloc)
                    tpl_mismatch_atoms[tpl_conf.strip()]
                    continue

                elif line.startswith("          ") and tpl_mismatch is not None:
                    #          DE22
                    last_key = list(tpl_mismatch_atoms)[-1]
                    tpl_mismatch_atoms[last_key].add(line.strip())
                    continue
                elif line.startswith("   STOP:  "):
                    # remove multiple spaces:
                    line = " ".join(line.split())
                else:
                    continue

            if log_section.idx == 5:
                # flag if 'debug.log' found in line:
                if not log_section.debuglog: 
                    log_section.has_debuglog(line)

                if line.startswith("   Total deleted cofactors"):
                    if int(line.rsplit(maxsplit=1)[1][:-1]) != 0:
                        line = line.strip()
                    else:
                        continue

            if skip:
                if isinstance(log_section.skip_lines, tuple):
                    found = False
                    for t in log_section.skip_lines:
                        found = found or (t in line)
                    if found:
                        continue
                else:
                    if line in log_section.skip_lines:
                        continue

            if change:
                # remove common start:
                if line.startswith(log_section.line_start):
                    line = line.removeprefix(log_section.line_start)

            out.append(line)

        if log_section.idx == 5:
            if self.dry_opt:
                out.insert(0, MSG_KEEP_H2O)

        # check if new tpl confs:
        if log_section.idx == 3:
            if newtpl:
                out.append("Generic topology file created for")
                out.append(newtpl)

            if tpl_mismatch:
                d = get_path_keys(self.pdb)
                out.append("Unloadable topology")
                fmt = "Unmatched {} topology for these atoms: {}, in these residues: {}.\n"
                for k in tpl_mismatch:
                    atms = ", ".join(a for a in tpl_mismatch_atoms[k])
                    reslocs = ", ".join(lr for lr in tpl_mismatch[k])
                    out.append(fmt.format(k, atms, reslocs))

                out.append(
                    (
                        "Likely cause: the renaming file is missing entries for these species, resulting in "
                        f"unloadable topology files;\n(renaming file: {d['renaming file']}; topologies {d['topologies']}.\n"
                    )
                )

        return out

    def get_blocks(self, text: str) -> dict:
        """Extract 'processing blocks' from contents of run1.log file
        passed into text argument.
        """
        block_txt = {}
        for k in self.blocks_specs:
            lhdr = self.blocks_specs[k]
            rpt_k = lhdr.rpt_hdr
            if k in [5, 7]:
                # dynamic headers
                h = lhdr.hdr
                lhdr.hdr = h.format(float(self.runprm[lhdr.get_key]))

            content = extract_content_between_tags(text, lhdr.hdr)
            if content is None:
                continue
            else:
                content = content.splitlines()

            if k == 1:
                content = sorted(content)

            if (lhdr.line_start is not None) or (lhdr.skip_lines is not None):
                content = self.process_content_block(content, lhdr)

            block_txt[rpt_k] = [line for line in content if line.strip()]

        # process termini; group res into NTR, CTR
        b2_hdr = self.blocks_specs[2].rpt_hdr
        if block_txt.get(b2_hdr) is not None:
            if block_txt[b2_hdr]:
                termi = defaultdict(list)
                for line in block_txt[b2_hdr]:
                    try:
                        i = line.index('"', 3) + 1
                        termi[line[-3:]].append(line[:i])
                    except ValueError:
                        continue

                block_txt[b2_hdr] = []
                for k in termi:
                    block_txt[b2_hdr].append((k, termi[k]))

        if self.blocks_specs[5].debuglog:
            # add extra line for each info found:
            rk = self.blocks_specs[5].rpt_hdr
            block_txt[rk].extend(self.get_debuglog_species())

        # collapse dist clashes block 7:
        if block_txt.get("Distance Clashes") is not None:
            if block_txt["Distance Clashes"]:
                new7 = []
                new7.append("Clashes found")
                for d in block_txt["Distance Clashes"]:
                    new7.append(d.strip())
                new7.append("end_clash")  # tag for formatting section

                block_txt["Distance Clashes"] = new7

        self.txt_blocks = block_txt

        return


def filter_heavy_atm_section(pdb: Path, s1_info_d: dict) -> dict:
    """Process the 'Missing Heavy Atoms' section to remove
    lines for missing backbone atoms of terminal residues.
    """
    # termini values are [2-tuples]
    termi = s1_info_d["MCCE.Step1"].get("Termini")
    heavy = s1_info_d["MCCE.Step1"].get("Missing Heavy Atoms")
    if heavy is None or termi is None:
        return s1_info_d

    if len(heavy) > 1:
        _ = heavy.pop(-1)
        # == Ignore warning messages if they are in the terminal res

    hvy_lst = []
    for line in heavy:
        conf, res = line.split(" in ")
        is_bkb = conf.rsplit(maxsplit=1)[1].endswith("BK")
        if is_bkb and (res in T[1] for T in termi):
            continue
        hvy_lst.append(line)

    # update dict
    s1_info_d["MCCE.Step1"]["Missing Heavy Atoms"] = hvy_lst

    return s1_info_d


def info_s1_log(pdb: Path) -> dict:
    s1log = RunLog1Parser(pdb)
    s1log.runprm = s1log.get_runprm()
    s1log.dry_opt = float(s1log.runprm["H2O_SASCUTOFF"]) == -0.01

    s1_text = s1log.get_log_contents()
    s1log.get_blocks(s1_text)

    # set the section data with dict & cleanup heavy atoms section:
    dout = {}
    dout = {"MCCE.Step1": s1log.txt_blocks}
    dout = filter_heavy_atm_section(pdb, dout)

    return dout



def hetero_count_per_chain(atomlines: list) -> Union[dict, None]:
    """Given a list of pdb coordinates lines, parse the hetero atoms to
    return a dictionary holding the count of free cofactors and waters,
    per chain (key).
    """
    # get unique hetero species per chain:
    hetero_set = defaultdict(set)
    for line in atomlines:
        if not line.startswith("HETATM"):
            continue
        # name = line[17:20]; chn = line[21]; seq = int(line[22:26])
        hetero_set[line[21]].add((line[17:20], int(line[22:26])))
    if not hetero_set:
        return None
    # get species count per chain:
    tots_per_chain = defaultdict(dict)
    for c in hetero_set:
        cntr = defaultdict(int)
        for val in hetero_set[c]:
            cntr[val[0]] += 1
        tots_per_chain[c]=  dict(cntr)
        
    return dict(tots_per_chain)


def get_cofactors_change(pdb_heteros: dict, step1_heteros: dict) -> list:
    """Return a list of changes in cofactors counts between the starting structure
    dictionary entry key prot_d["PDB.Structure"]["Model 1 Free Cofactors & Waters"] passed in pdb_heteros
    and step1_d["MCCE.Step1"]["Free Cofactors"] passed in step1_heteros.
    """
    # check what was removed:
    diff = []
    for chn in pdb_heteros:
        for k in pdb_heteros[chn]:
            if step1_heteros.get(chn) is None:
                continue
            if step1_heteros[chn].get(k) is None:
                diff.append(f"Removed all {pdb_heteros[chn][k]} {k} in {chn}.")
            else:
                # get non-zero count difference
                k_rem = step1_heteros[chn][k]
                d = pdb_heteros[chn][k] - k_rem
                if d:
                    diff.append(f"Removed {d} {k} in {chn}; {k_rem} remaining.")

    return diff


def update_s1_dict_cofactors_change(pdb: Path, prot_d: dict, step1_d: dict):
    """Update step1_d["MCCE.Step1"]["Free Cofactors"] list with cofactor changes if any.
    """
    if prot_d["PDB.Structure"].get("Model 1 Free Cofactors & Waters") is None:
        return

    pdb_heteros = prot_d["PDB.Structure"]["Model 1 Free Cofactors & Waters"]

    # get cofactors count from step1 pdb:
    s1_pdb = pdb.parent.joinpath("step1_out.pdb")
    s1_heteros = hetero_count_per_chain(s1_pdb.read_text().splitlines())
    if s1_heteros is None:
        return

    diff = get_cofactors_change(pdb_heteros, s1_heteros)
    if diff:
        step1_d_heteros = []
        idx = 0
        # initial list
        step1_d_heteros = step1_d["MCCE.Step1"]["Free Cofactors"].copy()
        if step1_d_heteros[idx].startswith("Total deleted cofactors"):
            idx = 1
        changes = " ".join(h1 for h1 in diff)
        if step1_d_heteros[idx] != changes:
            step1_d_heteros.insert(idx, changes)
            step1_d["MCCE.Step1"]["Free Cofactors"] = step1_d_heteros

    return


def collect_info(pdb: Path, args: Namespace) -> Tuple[dict, Union[dict, None]]:
    """Return at least one dict holding info from mcce4.pdbio.
    The second dict is None when step1 cannot be run, otherwise
    it contains info from the parsed run1.log file.
    Args:
      pdb (Path): File path of validated pdb.
      args (argparse.Namespace): Cli arguments (for creating step1 script).
    Returns:
      A 2-tuple of dicts when step1 can run, else (dict1, None).
    """
    step1_d = None

    # structural info:
    prot_d = info_input_prot(pdb)

    DO_STEP1 = (USER_MCCE is not None
                and prot_d["PDB.Structure"].get("MultiModels") is None
                and prot_d["PDB.Structure"].get("UNUSABLE") is None
    )

    if DO_STEP1:
        result = run.do_step1(pdb, args)
        if result is None:  # no error message
            step1_d = info_s1_log(pdb)
            # update cofactors section changes if any:
            update_s1_dict_cofactors_change(pdb, prot_d, step1_d)
        else:
            step1_d = {"MCCE.Step1": f"Error in run script 's1.sh': {result}."}

    return prot_d, step1_d


def write_report(pdb: Path, prot_d: dict, s1_d: Union[dict, None]):
    """Write the prerun report for pdb in its parent folder using the
    information passed in the parsers dicts, prot_d and step1_d.
    Args:
      pdb (Path): the pdb filepath; expected: path of pdb in prerun folder.
      prot_d (dict): The dictionary of sections from mcce4.pdbio.
      s1_d ([dict, None]): The dictionary of sections from step1 log parser.
    """
    name = prot_d.pop("Name")
    if s1_d is None:
        dict_lst = [prot_d]
    else:
        dict_lst = [prot_d, s1_d]

    # save report in the call folder (e.g. ../prerun/.)
    rpt_fp = pdb.parent.parent.joinpath(f"{pdb.stem}_{RPT}")
    with open(rpt_fp, "w") as rpt:
        rpt.write(f"---\n# {name}\n")

        for i, subd in enumerate(dict_lst):
            # h2: section hdrs, PDB.Structure or MCCE.Step1
            h2 = list(subd.keys())[0]
            rpt.write(f"## {h2}\n")

            for k in subd[h2]:
                if not subd[h2][k]:
                    continue

                if i == 0:
                    if isinstance(subd[h2][k], (str, int, float)):
                        rpt.write(f"### {k}: {subd[h2][k]}\n")
                        continue

                    rpt.write(f"### {k}:\n")
                    if isinstance(subd[h2][k], list):
                        for val in subd[h2][k]:
                            rpt.write(f"  - {val}\n")
                    elif isinstance(subd[h2][k], dict):
                        for kk, vv in subd[h2][k].items():
                            if isinstance(vv, (dict, tuple)):
                                if isinstance(vv, dict):
                                    wid = max(len(str(list(vv.values()))), 100)
                                else:
                                    wid = max(len(str(vv)), 100)
                                out = pformat(vv, sort_dicts=True, compact=True, width=wid)[1:-1]
                                rpt.write(f"  - {kk}:\n {out}\n")
                            else:
                                rpt.write(f"  - {kk}: {vv}\n")
                else: 
                    rpt.write(f"### {k}:\n")
                    for val in subd[h2][k]:
                        if k == "Distance Clashes":
                            if isinstance(val, str) and val.startswith("Clashes"):
                                rpt.write(f"<details><summary>{val}</summary>\n\n")
                            elif isinstance(val, str) and val.endswith("end_clash"):
                                rpt.write("\n</details>\n")
                            else:
                                rpt.write(f"- {val}\n")

                        elif k == "Labeling":
                            if val.startswith("Generic") or val.startswith(
                                "Unloadable"
                            ):
                                rpt.write(
                                    f"<strong><font color='red'>{val}:</font></strong>  \n"
                                )
                            elif val.startswith("Likely"):
                                val = val.replace(
                                    "Likely cause",
                                    "<strong><font color='red'>Likely cause</font></strong>",
                                )
                                rpt.write(f"{val}  \n")
                            else:
                                rpt.write(f"{val}\n")
                        else:
                            if isinstance(val, tuple):
                                ter, lst = val
                                rpt.write(
                                    f" - <strong>{ter}</strong>: {', '.join(lst)}\n"
                                )
                            elif isinstance(val, list):
                                ter, lst = val
                                rpt.write(
                                    f" - <strong>{ter}</strong>: {', '.join(x for x in lst)}\n"
                                )
                            elif isinstance(val, dict):
                                for kk, vv in val.items():
                                    rpt.write(f"  - {kk}: {vv}\n")
                            else:
                                rpt.write(f"  - {val}\n")

                rpt.write("\n")

    return
