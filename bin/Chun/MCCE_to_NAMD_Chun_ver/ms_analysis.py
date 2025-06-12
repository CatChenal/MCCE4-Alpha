#!/usr/bin/env python

from collections import defaultdict
import operator
from pathlib import Path
import sys
from typing import List, Tuple, Union
import zlib
import numpy as np
import pandas as pd

#--The is a new MSana on Aug 15
# TODO: Fix functions so that they work via import;

ph2Kcal = 1.364
Kcal2kT = 1.688


class Conformer:
    """Minimal Conformer class for use in microstate analysis.
    Attributes: iconf, confid, ires, resid, crg.
    """
    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.ires = 0
        self.resid = ""
        self.crg = 0.0

    def load_from_head3lst(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3] + self.confid[5:11]
        self.crg = float(fields[4])

    def __str__(self):
        return f"{self.confid} ({self.iconf}), {self.resid} ({self.ires}), {self.crg}"


class Conformers:
    """Holds a collection of Conformer objects.
    Has methods and attributes pertaining to that collection."""

    def __init__(self, head3_path: str):
        """Given the path to 'head3.lst', process it to populate these attributes:
          - conformers (list): List of Conformer objects
          - N (int): size of conformers list
        """
        self.h3_fp = Path(head3_path)
        self.conformers = None
        self.N = None

        self.load_data()

    def load_data(self) -> list:
        """Populate the class attributes."""
        if not self.h3_fp.exists():
            #print(f"File not found: {self.h3_fp}")
            return
        with open(self.h3_fp) as fp:
            lines = fp.readlines()
        lines.pop(0)
        self.conformers = []
        for line in lines:
            conf = Conformer()
            conf.load_from_head3lst(line)
            self.conformers.append(conf)

        self.N = len(self.conformers)

        return


try:
    # conformers will be an empty list if module is imported
    # or called outside of a MCCE output folder
    # or if head3.lst is missing or corrupted.
    Confs = Conformers("head3.lst")
    conformers = Confs.conformers
except FileNotFoundError:
    conformers = []


def get_confs_collection(head3_fp: Path) -> Conformers:
    """Returns:
       An instance of the Conformer collections class.
    """
    Confs = Conformers(head3_fp)

    return Confs


def read_conformers(head3_path: str):
    """Legacy function. See Conformers.load_data."""
    conformers = []
    lines = open(head3_path).readlines()
    lines.pop(0)
    for line in lines:
        conf = Conformer()
        conf.load_from_head3lst(line)
        conformers.append(conf)

    return conformers


class Microstate:
    """Sortable class for mcce conformer microstates."""
    def __init__(self, state: list, E: float, count: int):
        self.state = state
        self.E = E
        self.count = count

    def __str__(self):
        return f"Microstate(\n\tcount={self.count:,},\n\tE={self.E:,},\n\tstate={self.state}\n)"

    def _check_operand(self, other):
        """Fails on missing attribute."""

        if not (hasattr(other, "state") and hasattr(other, "E") and hasattr(other, "count")):
            return NotImplemented("Comparison with non Microstate object.")
        return

    def __eq__(self, other):
        self._check_operand(other)
        return (self.state, self.E, self.count) == (
            other.state,
            other.E,
            other.count,
        )

    def __lt__(self, other):
        self._check_operand(other)
        return (self.state, self.E, self.count) < (
            other.state,
            other.E,
            other.count,
        )


class Charge_Microstate:
    """
    Sortable class for charge microstates.
    Notes:
      For usage symmetry with Microstate class, the energy accessor E exists,
      but here: Charge_Microstate.E == Charge_Microstate.average_E. Thus, accessing
      E via Charge_Microstate.E, returns the average energy of a charge microstate.
      The Charge_Microstate.crg_stateid is implemented as compressed bytes as in the
      ms_analysis.py Microstate class in Junjun Mao's demo, but the 'key' that
      is encoded is a 2-tuple: (resid, crg) to facilitate further processing.
    """
    def __init__(self, crg_state: list, total_E: float, count: int):
        # crg_state is a list of (resid, crg) tuples:
        self.crg_stateid = zlib.compress(" ".join([x for x in crg_state]).encode())
        self.average_E = self.E = 0  # .E -> average E
        self.total_E = total_E
        self.count = count

    def state(self):
        #recover (resid, crg) tuple:
        return [i for i in zlib.decompress(self.crg_stateid).decode().split()]

    def __str__(self):
        s = f"Charge_Microstate(\n\tcount = {self.count:,},\n"
        s = s + f"\tE = {self.E:,.2f},\n\tstate = {self.state()}\n"

        return s

    def _check_operand(self, other):
        """Fails on missing attribute."""
        if not (hasattr(other, "crg_stateid") and hasattr(other, "E") and hasattr(other, "count")):
            return NotImplemented("Comparison with non Charge_Microstate object.")
        return

    def __eq__(self, other):
        self._check_operand(other)
        return (self.crg_stateid, self.E, self.count) == (
            other.crg_stateid,
            other.E,
            other.count,
        )

    def __lt__(self, other):
        self._check_operand(other)
        return (self.crg_stateid, self.E, self.count) < (
            other.crg_stateid,
            other.E,
            other.count,
        )


class MSout:
    def __init__(self, fname):
        self.T = 273.15
        self.pH = 7.0
        self.Eh = 0.0
        self.N_ms = 0
        self.N_uniq = 0
        self.lowest_E = 0.0
        self.highest_E = 0.0
        self.average_E = 0.0
        #self.fixed_crg = 0.0  unused
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.fixed_iconfs = []
        self.free_residues = []  # free residues, referred by conformer indices, iconf
        self.iconf2ires = {}     # from conformer index to free residue index
        self.microstates = {}    # dict of Microstate objects

        self.load_msout(fname)

    def load_msout(self, fname):
        """ Process the 'msout file' to populate these attributes:
          - fixed_iconfs (list)
          - T, pH, Eh (float)
          - fixed_iconfs (list)
          - free_residues (list)
          - iconf2ires (dict)
          - N_ms, N_uniq (int)
          - microstates (dict)
          - lowest_E, average_E, highest_E (float)
          """

        with open(fname) as fh:
            lines = fh.readlines()

        # Get a valid line
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        fields = line.split(",")
        for field in fields:
            key, value = field.split(":")
            key = key.strip().upper()
            value = float(value)
            if key == "T":
                self.T = value
            elif key == "PH":
                self.pH = value
            elif key == "EH":
                self.Eh = value

        # second line, confirm this is from Monte Carlo sampleing
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        key, value = line.split(":")
        if key.strip() != "METHOD" or value.strip() != "MONTERUNS":
            print("This file %s is not a valid microstate file" % fname)
            sys.exit(-1)

        # Third line, fixed conformer indicies
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        _, iconfs = line.split(":")
        self.fixed_iconfs = [int(i) for i in iconfs.split()]

        # 4th line, free residues
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        _, residues_str = line.split(":")
        residues = residues_str.split(";")
        self.free_residues = []
        for f in residues:
            if f.strip():
                self.free_residues.append([int(i) for i in f.split()])
        for i_res in range(len(self.free_residues)):
            for iconf in self.free_residues[i_res]:
                self.iconf2ires[iconf] = i_res

        # find the next MC record
        found_mc = False
        newmc = False
        self.N_ms = 0

        for line in lines:
            if line.find("MC:") == 0:  # ms starts
                found_mc = True
                newmc = True
                continue
            elif newmc:
                f1, f2 = line.split(":")
                current_state = [int(c) for c in f2.split()]
                newmc = False
                continue
            elif found_mc:
                fields = line.split(",")
                if len(fields) >= 3:
                    state_e = float(fields[0])
                    count = int(fields[1])
                    flipped = [int(c) for c in fields[2].split()]

                    for ic in flipped:
                        ir = self.iconf2ires[ic]
                        current_state[ir] = ic

                    ms = Microstate(list(current_state), state_e, count)
                    key = ",".join(["%d" % i for i in ms.state])
                    if key in self.microstates:
                        self.microstates[key].count += ms.count
                    else:
                        self.microstates[key] = ms

        # find N_ms, lowest, highest, averge E
        self.N_ms = 0
        E_sum = 0.0
        self.lowest_E = next(iter(self.microstates.values())).E
        self.highest_E = next(iter(self.microstates.values())).E
        msvalues = self.microstates.values()
        self.N_uniq = len(msvalues)
        for ms in msvalues:
            self.N_ms += ms.count
            E_sum += ms.E * ms.count
            if self.lowest_E > ms.E:
                self.lowest_E = ms.E
            if self.highest_E < ms.E:
                self.highest_E = ms.E
        self.average_E = E_sum / self.N_ms

    def get_sampled_ms(
        self,
        size: int,
        kind: str = "deterministic",
        seed: Union[None, int] = None,
    ) -> list:
        """
        Implement a sampling of MSout.microstates depending on `kind`.
        Args:
            size (int): sample size
            kind (str, 'deterministic'): Sampling kind: one of ['deterministic', 'random'].
                 If 'deterministic', the microstates in ms_list are sampled at regular intervals
                 otherwise, the sampling is random. Case insensitive.
            seed (int, None): For testing purposes, fixes random sampling.
        Returns:
            A list of lists: [[selection index, selected microstate], ...]
        """

        if not len(self.microstates):
            print("The microstates dict is empty.")
            return []

        kind = kind.lower()
        if kind not in ["deterministic", "random"]:
            raise ValueError(f"Values for `kind` are 'deterministic' or 'random'; Given: {kind}")

        ms_sampled = []
        ms_list = list(self.microstates.values())
        counts = ms_counts(ms_list)  # total number of ms
        sampled_cumsum = np.cumsum([mc.count for mc in ms_list])

        if kind == "deterministic":
            sampled_ms_indices = np.arange(size, counts - size, counts / size, dtype=int)
        else:
            rng = np.random.default_rng(seed=seed)
            sampled_ms_indices = rng.integers(low=0, high=counts, size=size, endpoint=True)

        for i, c in enumerate(sampled_ms_indices):
            ms_sel_index = np.where((sampled_cumsum - c) > 0)[0][0]
            ms_sampled.append([ms_sel_index, ms_list[ms_sel_index]])

        return ms_sampled

    def sort_microstates(self, sort_by: str = "E", sort_reverse: bool = False) -> list:
        """Return the list of microstates sorted by one of these attributes: ["count", "E"],
        and in reverse order (descending) if sort_reverse is True.
        Args:
        microstates (list): list of Microstate instances;
        sort_by (str, "E"): Attribute as sort key;
        sort_reverse (bool, False): Sort order: ascending if False (default), else descending.
        Return None if 'sort_by' is not recognized.
        """

        if sort_by not in ["count", "E"]:
            print("'sort_by' must be a valid microstate attribute; choices: ['count', 'E']")
            return None

        return sorted(list(self.microstates.values()), key=operator.attrgetter(sort_by), reverse=sort_reverse)


class Charge_Microstates:
    """
    Holds a collection of charge microstates. Has methods over the collection.
    Attributes:
      - charge_microstates (list): List of ms_analysis.Charge_Microstate objects
                                   processed from ms_analysis.MSOut.microstates.
      - orig_ms_by_crg_stateid (dict): Grouping of original conformer microstates
                                   for each charge state in charge_microstates.
      - crg_group_stats (dict): Values: (low_state, lowE), (hi_state, hi_E), average.
    """
    def __init__(self, microstates: dict, conformers: list):
        """
        microstates (dict): MSout.microstates
        conformers (list): conformers, i.e. Confs.conformers
        """
        self.charge_microstates = []
        self.orig_ms_by_crg_stateid = {}
        self.crg_group_stats = {}

        # populate charge_microstates & orig_ms_by_crg_stateid:
        self.ms_to_charge_ms(microstates, conformers)
        # populate crg_group_stats:
        self.get_crg_state_energy_stats()

    def ms_to_charge_ms(self, microstates: dict, conformers: list) -> Union[list, dict]:
        """Process the conformer microstates collection to
        obtain a list of charge microstate objects.
        """
        # 1. populate dict charge_ms_by_id:
        charge_ms_by_id = {}
        # crgms_microstate gathers all ms objects in a charge state:
        crgms_microstate = defaultdict(list)

        for ms in list(microstates.values()):
            current_crg_state = [f"{conformers[ic].resid}|{round(conformers[ic].crg)}"
                                 for ic in ms.state]
            crg_ms = Charge_Microstate(current_crg_state, ms.E * ms.count, ms.count)
            crg_id = crg_ms.crg_stateid  # compressed bytes for key
            crgms_microstate[crg_id].append(ms)

            if crg_id in charge_ms_by_id:
                charge_ms_by_id[crg_id].count += crg_ms.count
                charge_ms_by_id[crg_id].total_E += crg_ms.total_E
            else:
                charge_ms_by_id[crg_id] = crg_ms
        self.orig_ms_by_crg_stateid = dict(crgms_microstate)

        # 2. complete population of class attribute w/aver E
        for k in charge_ms_by_id:
            crg_ms = charge_ms_by_id[k]
            crg_ms.average_E = crg_ms.E = crg_ms.total_E / crg_ms.count
            self.charge_microstates.append(crg_ms)

        return

    @staticmethod
    def ms_energy_stats(conf_microstates: list) -> tuple:
        """
        Returns a tuple that includes the state list that correspond
        to lower & highest E:
          [m_state_lo, lowest_E], [m_state_hi, highest_E], average_E;
        The lowest and highest energies and corresponding state (iconf list),
        along with the average energy of the listed conformer microstates.
        """
        ms = next(iter(conf_microstates))
        # initialize:
        lowest_E, highest_E = ms.E, ms.E
        m_state_lo, m_state_hi = None, None
        N_ms = 0
        total_E = 0.0
        for ms in conf_microstates:
            if lowest_E > ms.E:
                lowest_E = ms.E
                m_state_lo = ms.state
            elif highest_E < ms.E:
                highest_E = ms.E
                m_state_hi = ms.state
            N_ms += ms.count
            total_E += ms.E * ms.count
        average_E = total_E / N_ms

        return [m_state_lo, lowest_E], [m_state_hi, highest_E], average_E

    def get_crg_state_energy_stats(self):
        """Get the lowest, aver, highest energy of a collection of microstates with
        the same charge state.
        """
        if not self.orig_ms_by_crg_stateid:
            print("Empty self.orig_ms_by_crg_stateid")
            return None
        for stateid in self.orig_ms_by_crg_stateid:
            # tup:: [m_state_lo, lowest_E], [m_state_hi, highest_E], average_E
            tup = self.ms_energy_stats(self.orig_ms_by_crg_stateid[stateid])
            self.crg_group_stats[stateid] = tup

        return

    def get_topN_lowestE_crgms(self, n_top: int = 5) -> list:
        """Return list of [[m_state_lo, lowest_E], count]."""
        lst = []
        for i, cms in enumerate(self.charge_microstates):
            # first item is [m_state_lo, lowest_E], 0 index:
            lst.append([self.crg_group_stats[cms.crg_stateid][0], cms.count])

        # sort by lowest_E, keep topN:
        return sorted(lst, key=lambda x: x[0][1])[:n_top]


def free_residues_df(free_res: list, conformers: list, colname: str = "FreeRes") -> pd.DataFrame:
    """Return the free residues' ids in a pandas DataFrame."""

    free_residues = [conformers[res[0]].resid for res in free_res]

    return pd.DataFrame(free_residues, columns=[colname])


def fixed_res_crg(conformers: list,
                  fixed_iconfs: list,
                  res_of_interest: list = None,
                  return_df: bool = False,
                  rm_trailing_underscore: bool = False
                  ) -> Tuple[float, Union[dict, pd.DataFrame]]:
    """
    Args:
      fixed_iconfs (list): List of fixed conformers
      conformers (list): List of Conformers objects.
      res_of_interest (list, None): List of resid for filetering.
      return_df (bool, False): If True, the second item of the output tuple
                               will be a pandas.DataFrame, else a dict.
      rm_trailing_underscore (bool, False): If True, output resid w/o trailing "_".

    Returns:
      A 2-tuple:
      The net charge contributed by the fixed residues in `fixed_iconfs`;
      A dictionary: key=conf.resid, value=conf.iconf, int(conf.crg).
    """
    fixed_net_charge = 0.
    dd = defaultdict(float)
    for conf in conformers:
        if conf.iconf in fixed_iconfs:
            fixed_net_charge += conf.crg
            resid = conf.resid if not rm_trailing_underscore else conf.resid[:-1]
            dd[resid] = conf.iconf, int(conf.crg)
    if res_of_interest:
        dd = {k: dd[k] for k in dd if k[:3] in res_of_interest}
    if return_df:
        fixed_res_crg_df = pd.DataFrame([[k, v[0], v[1]] for k, v in dd.items()],
                                        columns=["residues", "iconf", "crg"])
        return fixed_net_charge, fixed_res_crg_df

    return fixed_net_charge, dict(dd)


def ms_counts(microstates: Union[dict, list]) -> int:
    """Sum the microstates count attribute."""
    if isinstance(microstates, dict):
        return sum(ms.count for ms in microstates.values())
    else:
        return sum(ms.count for ms in microstates)


def ms_charge(ms: Microstate):
    """Compute microstate charge"""
    crg = 0.0
    for ic in ms.state:
        crg += conformers[ic].crg
    return crg


def get_ms_crg(ms: Microstate, conformers: list):
    """Compute microstate charge.
    Alternate version of `ms_charge` for use when
    conformers is not a global variable (case when
    module is imported).
    """
    crg = 0.0
    for ic in ms.state:
        crg += conformers[ic].crg
    return crg


def groupms_byenergy(microstates: list, ticks: List[float]) -> list:
    """
    Group the microstates' energies into bands provided in `ticks`.
    Args:
      microstates (list): List of microstates
      ticks (list(float)): List of energies.
    """

    N = len(ticks)
    ticks.sort()
    ticks.append(1.0e100)  # add a big number as the last boundary
    resulted_bands = [[] for i in range(N)]

    for ms in microstates:
        it = -1
        for itick in range(N):
            if ticks[itick] <= ms.E < ticks[itick + 1]:
                it = itick
                break
        if it >= 0:
            resulted_bands[it].append(ms)

    return resulted_bands


def groupms_byiconf(microstates: list, iconfs: list) -> tuple:
    """
    Divide the microstates by the conformers indices provided in `iconfs`
    into 2 groups: the first contains one of the given conformers, the
    second one contains none of the listed conformers.
    Args:
      microstates (list): List of microstates
      iconfs (list): List of conformer indices.
    Return:
      A 2-tuple: (Microstates with any of `iconfs`,
                  microstates with none)
    """

    ingroup = []
    outgroup = []
    for ms in microstates:
        contain = False
        for ic in iconfs:
            if ic in ms.state:
                ingroup.append(ms)
                contain = True
                break
        if not contain:
            outgroup.append(ms)

    return ingroup, outgroup


def groupms_byconfid(microstates: list, confids: list) -> tuple:
    """
    Divide the microstates by the conformers ids provided in `confids`
    into 2 groups: the first contains ALL of the given conformers, the
    second one does not; contains only some but not all.
    Note: An ID is a match if it is a substring of the conformer name.
    Args:
      microstates (list): List of microstates
      confids (list): List of conformer ids.
    Return:
      A 2-tuple: Microstates with all of `confids`, microstates with some or none.
    """

    ingroup = []
    outgroup = []
    for ms in microstates:
        contain = True
        names = [conformers[ic].confid for ic in ms.state]
        for confid in confids:
            innames = False
            for name in names:
                if confid in name:
                    innames = True
                    break
            contain = contain and innames
        if contain:
            ingroup.append(ms)
        else:
            outgroup.append(ms)

    return ingroup, outgroup


def ms_energy_stat(microstates: list) -> tuple:
    """
    Return the lowest, average, and highest energies of the listed
    microstates.
    """

    ms = next(iter(microstates))
    lowest_E = highest_E = ms.E
    N_ms = 0
    total_E = 0.0
    for ms in microstates:
        if lowest_E > ms.E:
            lowest_E = ms.E
        elif highest_E < ms.E:
            highest_E = ms.E
        N_ms += ms.count
        total_E += ms.E * ms.count

    average_E = total_E / N_ms

    return lowest_E, average_E, highest_E


def ms_convert2occ(microstates: list) -> dict:
    """
    Given a list of microstates, convert to conformer occupancy
    for conformers that appear at least once in the microstates.
    Return:
      A dict: {ms.state: occ}
    """
    # TODO: use Counter

    occurance = {}  # dict of conformer occurance
    occ = {}
    N_ms = 0
    #!!- Chun - Jul 24: need to account for the case where the input list is empty
    #! my adjustment: added an if-else condition
    if microstates:
        for ms in microstates:
            N_ms += ms.count
            for ic in ms.state:
                if ic in occurance:
                    occurance[ic] += ms.count
                else:
                    occurance[ic] = ms.count

        for key in occurance:
            occ[key] = occurance[key] / N_ms

    return occ


def ms_convert2sumcrg(microstates: list, free_res: list) -> list:
    """
    Given a list of microstates, convert to net charge of each free residue.
    """
    # FIX: dependence on global conformers variable

    iconf2ires = {}
    for i_res in range(len(free_res)):
        for iconf in free_res[i_res]:
            iconf2ires[iconf] = i_res

    charges_total = [0.0 for i in range(len(free_res))]
    N_ms = 0
    for ms in microstates:
        N_ms += ms.count
        for ic in ms.state:
            ir = iconf2ires[ic]
            #!--should have a separte input for conformers below
            charges_total[ir] += conformers[ic].crg * ms.count

    charges = [x / N_ms for x in charges_total]

    return charges


def e2occ(energies: list) -> float:
    """Given a list of energy values in unit Kacl/mol,
    calculate the occupancy by Boltzmann Distribution.
    """

    e = np.array(energies)
    e = e - min(e)
    Pi_raw = np.exp(-Kcal2kT * e)
    Pi_sum = sum(Pi_raw)
    Pi_norm = Pi_raw / Pi_sum

    return Pi_norm


def bhata_distance(prob1: list, prob2: list) -> float:
    """Bhattacharyya distance between 2 probability distributions."""

    d_max = 10000.0  # Max possible value
    #!--Chun [Jul 24] unable to handle empty list which can happen
    #in some runs [added if else conditions below]
    if prob1 and sum(prob1) != 0:
        p1 = np.array(prob1) / sum(prob1)
    else:
        p1 = np.array([0])
    if prob2 and sum(prob2) != 0:
        p2 = np.array(prob2) / sum(prob2)
    else:
        p2 = np.array([0])
    if len(p1) != len(p2):
        d = d_max
    else:
        bc = sum(np.sqrt(p1 * p2))
        if bc <= np.exp(-d_max):
            d = d_max
        else:
            d = -np.log(bc)

    return d


def whatchanged_conf(msgroup1: list, msgroup2: list) -> dict:
    "Given two group of microstates, calculate what changed at conformer level."

    occ1 = ms_convert2occ(msgroup1)
    occ2 = ms_convert2occ(msgroup2)

    #!!-Chun bug fix below [Jul 24]: should be occ2.keys() instead of occ2.key()
    #all_keys = list(set(occ1.keys()) | set(occ2.key()))
    all_keys = list(set(occ1.keys()) | set(occ2.keys()))
    all_keys.sort()
    diff_occ = {}
    for key in all_keys:
        if key in occ1:
            p1 = occ1[key]
        else:
            p1 = 0.0
        if key in occ2:
            p2 = occ2[key]
        else:
            p2 = 0.0
        diff_occ[key] = p2 - p1

    return diff_occ


def whatchanged_res(msgroup1: list, msgroup2: list, free_res: list) -> list:
    "Return a list of Bhattacharyya distance of free residues."

    occ1 = ms_convert2occ(msgroup1)
    occ2 = ms_convert2occ(msgroup2)

    bhd = []
    for res in free_res:
        p1 = []
        p2 = []
        for ic in res:
            if ic in occ1:
                p1.append(occ1[ic])
            else:
                p1.append(0.0)
            if ic in occ2:
                p2.append(occ2[ic])
            else:
                p2.append(0.0)
        #!!-Chun bug fix below [Jul 24]: the distance is computed between 2 distributions
        # but sometimes, p2 can be null
        bhd.append(bhata_distance(p1, p2))

    return bhd


def example(msout_file:str):
    msout = MSout(msout_file)

    n_bands = 20
    e_step = (msout.highest_E - msout.lowest_E)/n_bands
    ticks = [msout.lowest_E + e_step*(i) for i in range(n_bands)]
    ms_in_bands = groupms_byenergy(msout.microstates.values(), ticks)
    print([len(band) for band in ms_in_bands])
    netural, charged = groupms_byiconf(msout.microstates.values(), [12, 13, 14, 15])
    lo_E, av_E, hi_E = ms_energy_stat(msout.microstates.values())
    print(lo_E, av_E, hi_E)

    # charge over energy bands
    e_step = (msout.highest_E - msout.lowest_E) / n_bands
    ticks = [msout.lowest_E + e_step*(i+1) for i in range(n_bands - 1)]
    ms_in_bands = groupms_byenergy(msout.microstates.values(), ticks)
    for band in ms_in_bands:
        band_total_crg = 0.0
        for ms in band:
            band_total_crg += ms_charge(ms)
        print(band_total_crg/ms_counts(band))

    netural, charged = groupms_byiconf(msout.microstates.values(), [12, 13, 14, 15])
    diff_occ = whatchanged_conf(netural, charged)
    for key in diff_occ:
        print("%3d, %s: %6.3f" % (key, conformers[key].confid, diff_occ[key]))

    diff_bhd = whatchanged_res(netural, charged, msout.free_residues)
    for ir in range(len(msout.free_residues)):
         print("%s: %6.4f" % (conformers[msout.free_residues[ir][0]].resid, diff_bhd[ir]))
    charges = ms_convert2sumcrg(msout.microstates.values(), msout.free_residues)
    for ir in range(len(msout.free_residues)):
         print("%s: %6.4f" % (conformers[msout.free_residues[ir][0]].resid, charges[ir]))

    microstates = list(msout.microstates.values())
    glu35_charged, _ = groupms_byconfid(microstates, ["GLU-1A0035"])
    print(len(microstates))
    print(len(glu35_charged))

    return


if __name__ == "__main__":
    example("ms_out/pH4eH0ms.txt")

