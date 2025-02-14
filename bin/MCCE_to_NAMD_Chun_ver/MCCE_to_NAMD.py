from MCCE_to_NAMD_Chun_ver.base.data import Info
from pathlib import Path
import pandas as pd
from collections import defaultdict
import os
import subprocess
import sys

#import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
import math
import pickle
import MCCE_to_NAMD_Chun_ver.ms_analysis as msa


class MCCE_pdb_to_VMD_pdb(object):
    '''
    Create an object to clean PDB from MCCE into a format that VMD can process for NAMD.
    This object will also store relevent information for post-processing
    towards psfgen.
    '''
    def __init__(
            self,
            in_mcce_pdb_path,
        ):
        self.in_mcce_pdb_path = in_mcce_pdb_path
        if os.path.exists(self.in_mcce_pdb_path):
            with open(self.in_mcce_pdb_path, 'r') as fin:
                self.in_mcce_pdb_lines = fin.readlines()
            print(f"{self.in_mcce_pdb_path} has been read.")
            self.out_VMD_pdb_path = os.path.join(
                os.path.dirname(self.in_mcce_pdb_path),
                f"VMD_{os.path.basename(self.in_mcce_pdb_path)}"
            )
        else:
            print(f"{self.in_mcce_pdb_path} is not found.")
            sys.exit()


    def mcce_pdb_to_vmd_pdb(self):
        if self.out_VMD_pdb_path:
            lines_tuples = [
                self.__mcce_line_to_vmd_line__(line) for line in self.in_mcce_pdb_lines
            ]
            self.out_VMD_pdb_lines = [elm[3] for elm in lines_tuples]
            self.chain_ids = [elm[0] for elm in lines_tuples]
            self.res_ids = [elm[1] for elm in lines_tuples]
            self.conformer_ids = [elm[2] for elm in lines_tuples]
            self.unique_chain_ids = set(self.chain_ids)
        pass


    def __mcce_line_to_vmd_line__(
            self,
            line,
            part_1_slice_ind: int=21,
            part_2_slice_ind: int=30,
            part_3_slice_ind: int=62,
            end_slice_ind: int=-16,
            regex_chain=r'([A-Z]+)([0-9]+)_([0-9]+)',
            regex_charge=r'(-*[0-9]+).([0-9]+)',
            pad_size: int=3,
        ):
        line = line.strip()
        part_1 = line[:part_1_slice_ind]
        chain_part = line[part_1_slice_ind:part_2_slice_ind]
        part_2 = line[part_2_slice_ind:part_3_slice_ind]
        charge_part = line[part_3_slice_ind:end_slice_ind]
        match_chain = re.search(regex_chain, chain_part)
        if match_chain:
            chain_id = self.__fill_chain_id__(chain_id=match_chain[1])
            res_id = self.__get_resid_id__(res_id=match_chain[2])
            if not res_id:
                return None
        else:
            return None
        match_charge = re.search(regex_charge, charge_part)
        if match_charge:
            charge = self.__fill_charge__(
                charge_1=match_charge[1], charge_2=match_charge[2]
            )
        else:
            return None
        out_line = part_1 + chain_id + resid_catch + ' ' * pad_size + part_2 + charge
        return (match_chain[1], match_chain[2], match_chain[3], out_line)


    def __fill_chain_id__(self, chain_id, vmd_width: int=1,):
        chain_id = str(chain_id)
        if len(chain_id) >= vmd_width:
            return chain_id[:vmd_width]
        else:
            return = " " * (vmd_width - len(chain_id)) + chain_id


    def __get_resid_id__(self, res_id, vmd_width: int=4,):
        for ind, elm in enumerate(res_id):
            if elm != '0':
                return " " * (vmd_width - len(res_id[ind:])) + res_id[ind:]
        return None


    def __fill_charge__(self, charge_1, charge_2, vmd_width: int=3,):
        return " " * (vmd_width - len(charge_1)) + charge_1 + '.' + charge_2




class PDBtoMSAna(object):
    '''
    Task: Perform pH dependent MCCE an dsubsequent MS analyses to prepare
    info. for MSAnaToNAMD
    '''
    def __init__(
            self,
            in_pdb_path,
            to_mcce_path,
            ms_base_dir,
            pbs_solver_choice: str='delphi'
        ):
        '''
        Task: Initialize the object

        Args:
            in_pdb_path:
                - the path to the input pdb

            to_mcce_path:
                - the path to the mcce executible, which is at the same place as other executibles,
                namely, step1.py .. etc,

            ms_base_dir:
                - the directory to hold outputs for MS analysis work.
        '''
        #--Check for input
        if os.path.exists(in_pdb_path):
            print("Launching calculations for MCCE to NAMD.")
            self.in_pdb_path = in_pdb_path
            self.in_pdb_name = os.path.basename(self.in_pdb_path)
            self.to_mcce_path = to_mcce_path
            self.to_mcce_dir = os.path.dirname(self.to_mcce_path) if os.path.dirname(self.to_mcce_path) \
            else os.getcwd()
            self.alt_mcce_dir = self.to_mcce_dir.replace("bin", "MCCE_bin")
            self.ms_base_dir = ms_base_dir if not str(ms_base_dir) == "." or str(ms_base_dir) == "./"\
            else os.getcwd()
            self.pbs_solver_choice = pbs_solver_choice
        else:
            print(f"{in_pdb_path} not found. Operations aborted.")
            return None


    def run_mcce(
            self,
            pH_value: float=7.0,
            run_mcce: bool=False,
        ):
        '''
        Task: run MCCE calculations or generate a bash script for MCCE calculations

        Args:
            pH_value:
                - the pH-value to be accessed in MCCE calculations.

            run_mcce:
                - if False, a bash script for running MCCE calculations will be generated.
                This script will only be run if "run_mcce" is set to be True.
        '''
        self.pH_value = pH_value
        self.out_sh = f'''#! /bin/bash

## stop on errors
set -e

export PATH={self.to_mcce_dir}:{self.alt_mcce_dir}:$PATH

in_pdb_path="{self.in_pdb_path}"
out_pdb_path="{self.ms_base_dir}/prot.pdb"
work_dir="{self.ms_base_dir}"
pH_value={self.pH_value}
pbs_solver_choice="{self.pbs_solver_choice}"

cp $in_pdb_path $out_pdb_path
cd $work_dir

step1.py --dry prot.pdb
step2.py
step3.py
step3.py -s $pbs_solver_choice
step4.py -i 1 -n $pH_value --ms

cp head3.lst ms_out/.
cp pK.out ms_out/.
cp respair.lst ms_out/.
'''
        self.out_sh_path = os.path.join(self.ms_base_dir, f'MCCE_run_at_ph_{self.pH_value}.sh')
        #--Should make the output directory here.
        Path(self.ms_base_dir).mkdir(parents=True, exist_ok=True)
        with open(self.out_sh_path, 'w') as fout:
            fout.write(self.out_sh)
        print(f"The script for running MCCE calculations has been created at: {self.out_sh_path}")

        if run_mcce:
            print(f"Running {self.out_sh_path}.")
            subprocess.call(['bash', self.out_sh_path])
            print("MCCE calculations complete.")


    def run_MSana(
            self,
            residue_of_interest='all',
            num_of_ms_picked: int=5,
            out_pickle: bool=False,
        ):
        '''
        Task: Run the essential MS analysis steps to extract information for
        MCCE to NAMD run later on.
        '''
        self.num_of_ms_picked = num_of_ms_picked
        self.in_ms_path = os.path.join(self.ms_base_dir, 'ms_out', f"pH{int(self.pH_value)}eH0ms.txt")
        self.in_head3_path = os.path.join(self.ms_base_dir, 'ms_out', "head3.lst")
        if not os.path.exists(self.in_ms_path) or not os.path.exists(self.in_head3_path):
            print(f'''Either {self.in_ms_path} or {self.in_head3_path} is missing.
Please run MCCE calculations.
            ''')
            return None
        else:
            self.mc = msa.MSout(self.in_ms_path)

            if out_pickle:
                self.out_pickle = os.path.join(self.ms_base_dir, f"{self.in_pdb_name}_pH{self.pH_value}.pickle")
                with open(self.out_pickle, 'wb') as fout:
                    pickle.dump(self.mc, fout, protocol=pickle.HIGHEST_PROTOCOL)

            if residue_of_interest == 'all':
                self.interested_res = list(Info().charmm_dict_for_MCCE.keys())
            else:
                self.interested_res = residue_of_interest
            print(f"Residues to be included in MS analysis report are: {self.interested_res}")

            #--!!Old
            self.mc.conformers = msa.read_conformers(self.in_head3_path)
            #--!!New Aug 15
            self.conformers = mas.get_confs_collection(self.in_head3_path)
            #--!!Old
            #crg_orig_list = self.__convert_ms_crg__(self.__get_ms_orig_list__(), self.__get_id_vs_charge__())
            #--!!New Aug 15
            self.charge_ms = msa.Charge_Microstates(
                microstates=mc.microstates,
                conformers=self.conformers.conformers
            )
            #--!!Old
            #charge_ms_file = self.__find_unique_crgms_count_order__(crg_orig_list)
            #--!!New Aug 15
            self.charge_ms_file = self.__get_charge_ms_file__(self.charge_ms)
            #--!!Old
            #all_crg_count_res = self.__conca_crg_ms_pandas__(
            #    charge_ms_file[0],
            #    charge_ms_file[1],
            #    charge_ms_file[2],
            #    self.__get_free_residues__(),
            #    self.__get_background_charge__(),
            #    self.interested_res
            #)
            #--!!New
            self.free_residues = msa.free_residues_df(
                self.mc.free_residues, self.conformers.conformers, colname='Residue',
            )
            self.background_charge_info = msa.fixed_res_crg(
                conformers=self.conformers.conformers,
                fixed_iconfs=self.mc.fixed_iconfs,
                res_of_interest=self.interested_res,
                return_df=False,
            )
            all_crg_count_res = self.__conca_crg_ms_pandas__(
                self.charge_ms_file[0],
                self.charge_ms_file[1],
                self.charge_ms_file[2],
                self.free_residues,
                self.background_charge_info[0],
                self.interested_res,
            )
            self.ms_result = all_crg_count_res.head(self.num_of_ms_picked)
            self.out_csv_for_namd = os.path.join(self.ms_base_dir, f"{self.in_pdb_name}_pH{self.pH_value}_MSAna.csv")
            self.ms_result.to_csv(self.out_csv_for_namd, index=False)
            print(f"The CSV file needed for converting MS analysis results to PSF for NAMD run has been created at: {self.out_csv_for_namd} .")
            return True


    def __get_charge_ms_file__(self, charge_ms):
        '''
        Args:
            charge_ms:
                - Charge microstates object parsed from MSA
        '''
        crgs_ = []
        counts_ = []
        orders_ = []
        for ind, elm in enumerate(charge_ms.get_topN_lowestE_crgms()):
            #--the list of unique charges
            energy = elm[0][1]
            count = elm[1]
            conf_ids = elm[0][0]
            conf_crgs = [self.conformers.conformers[id].crg for id in conf_ids]
            crgs_.append(conf_crgs)
            counts_.append(count)
            orders_.append(ind + 1)

        return [crgs_, counts_, orders_, ]


    def __get_ms_orig_list__(self):
        #--!!Old depreciated
        ms_orig_list = [
            [ms.E, ms.count, ms.state] for ms in list((self.mc.microstates.values()))
        ]
        return sorted(ms_orig_list, key = lambda x: x[0])


    def __get_free_residues__(self):
        #--!!Old depreciated
        free_residues = []
        for res in self.mc.free_residues:
            try:
                #free_residues.append(msa.conformers[res[0]].resid)
                free_residues.append(self.mc.conformers[res[0]].resid)
            except:
                print(res[0])
        return pd.DataFrame(free_residues, columns = ["Residue"])


    def __get_id_vs_charge__(self):
        id_vs_charge = {}
        for conf in self.mc.conformers:
            id_vs_charge[conf.iconf] = conf.crg
        return id_vs_charge


    def __get_background_charge__(self):
        fixed_res_crg_dict = {}
        for conf in self.mc.conformers:
            if conf.iconf in self.mc.fixed_iconfs:
                try:
                    #initializing the dict
                    if conf.resid not in fixed_res_crg_dict:
                        fixed_res_crg_dict[conf.resid] = conf.crg
                except:
                    print("Error in ms file! Fixed residues are duplicating")

        return sum(fixed_res_crg_dict.values())


    def __convert_ms_crg__(self, l, d):
        #A recursive function
        crg_lst = [
            [
                y[0], y[1], [self.__convert_ms_crg__(x, d) if isinstance(x, list) else d.get(x, x) for x in y[2]]
            ] for y in l
        ]
        return crg_lst


    def __find_unique_crgms_count_order__(
            self,
            crg_list_ms,
            begin_energy = None,
            end_energy = None
        ):
        #--!!Old Depreciated
        """
        Args:
            - crg_list_ms: the charge list of microstates.
            * the crg file must have been sorted with increasing order.

            - begin_energy: min. enegy used for filtering charge id based on the ms' energy.
            - end_energy: min. enegy used for filtering charge id based on the ms' energy.
            * supplied together with begin_energy

        Return:
            - a unique_crg_state_order that gives the order of unique charge state based on energy.
            * Lowest energy charge state will give the order 1 and then second unique charge
            state will give the order 2. This order is based on unique charge ms order.
        """
        if not begin_energy and not end_energy:
            print("All energy microstates are selected.")
            begin_energy = crg_list_ms[0][0]
            eng_energy = crg_list_ms[-1][0]
        elif begin_energy and end_energy:
            crg_list_ms = [
                [x[0], x[1], x[2]]\
                for x in crg_list_ms \
                if x[0] >= begin_energy and x[0] <= end_energy
            ]
        else:
            #!! error here in the original code
            sys.exit('Give the lower and the upper energy bound.')

        #unique charge as key and energy, count and order
        crg_all_count = {}
        unique_crg_state_order = 1
        for x, array in enumerate(crg_list_ms):
            #adding new items
            if tuple(array[2]) not in crg_all_count.keys():
                crg_all_count[(tuple(array[2]))] = [
                    array[1], [array[0]], [unique_crg_state_order]
                ]
                unique_crg_state_order += 1
            else:
                #in the case of a key has already been generated
                #update values, stored as a list, associated with that key
                crg_all_count[(tuple(array[2]))][0] += array[1]
                #update th min energy by:
                #1st, search for the min value in the updated list
                min_energy = min(min(
                    crg_all_count[(tuple(array[2]))][1]
                ), array[0])
                max_energy = max(max(
                    crg_all_count[(tuple(array[2]))][1]
                ), array[0])
                # clear the energy list and append min. and max. energy to it
                crg_all_count[(tuple(array[2]))][1].clear()
                crg_all_count[(tuple(array[2]))][1].append(min_energy)
                crg_all_count[(tuple(array[2]))][1].append(max_energy)

        # make a list of count, unique charge microstate, energy difference and order.
        all_crg_ms_unique = []
        all_count = []
        energy_diff_all = []
        unique_crg_state_order = []
        for u, v in crg_all_count.items():
            all_crg_ms_unique.append(list(u))
            all_count.append(v[0])
            unique_crg_state_order.append(v[2][0])
            if len(v[1]) == 2:
                energy_diff_all.append(round(v[1][1] - v[1][0], 6))
            elif len(v[1]) == 1:
                energy_diff_all.append(0)
            else:
                sys.exit("There is error while creating unique charge state.")

        print(f"Total number of state: {len(crg_list_ms)}")
        print(f"Total number of unique charge ms: {len(all_crg_ms_unique)}")
        return all_crg_ms_unique, all_count, unique_crg_state_order, energy_diff_all


    def __conca_crg_ms_pandas__(
            self,
            unique_crg_ms_list,
            ms_count,
            ms_order,
            free_residues,
            background_charge,
            residue_interest_list
        ):
        unique_crg_ms_list_pd = pd.DataFrame(unique_crg_ms_list).T
        ms_count_pd = pd.DataFrame(ms_count, columns=["Count"]).T
        ms_order_pd = pd.DataFrame(ms_order, columns=["Order"]).T
        crg_ms_count_pd = pd.concat(
            [unique_crg_ms_list_pd, ms_count_pd, ms_order_pd]
        )
        crg_count_res_1 = pd.concat(
            [free_residues, crg_ms_count_pd], axis=1
        )
        crg_count_res_1.loc["Count", "Residue"] = 'Count'
        crg_count_res_1.loc["Order", "Residue"] = 'Order'
        all_crg_count_res = crg_count_res_1.set_index("Residue")
        # sort based on the count
        all_crg_count_res = all_crg_count_res.sort_values(
            by = "Count", axis = 1, ascending=False
        )
        all_crg_count_res.columns = range(all_crg_count_res.shape[1])
        all_crg_count_res = all_crg_count_res.T.set_index("Order")
        all_crg_count_res["Occupancy"] = round(
            all_crg_count_res["Count"] / sum(all_crg_count_res["Count"]), 3
        )
        all_crg_count_res['Sum_crg_protein'] = \
        all_crg_count_res.iloc[:, :-2].sum(axis=1) + background_charge
        crg_count_res = all_crg_count_res.copy()
        for i in all_crg_count_res.columns:
            if i[:3] not in residue_interest_list\
            and i != "Occupancy" \
            and i != "Count" \
            and i != "Sum_crg_protein":
                crg_count_res.drop([i], axis=1, inplace=True)

        return crg_count_res



class AggMCCEResults(object):
    '''
    Task: Aggregrate results from multiple MCCE runs into individual cvss
    for later use for NAMD
    '''
    def __init__(
        self,
        in_ms_csvs: list,
        out_dir: Path=os.getcwd(),
        number_of_chosen_cases: int=5
    ):
        '''
        Task: Initiaize the object.

        Args:
            in_ms_csvs:
                - A list of paths to input csvs for MCCE MS results.

            out_idr:
                - The output directory to deposit aggregrate result.

            number_of_chosen_cases:
                - Number of unique aggregated popular MS output for later NAMD use
        '''
        self.in_ms_csvs = in_ms_csvs
        if out_dir == '.' or out_dir == './':
            out_dir = os.getcwd()
        self.out_dir = out_dir
        self.number_of_chosen_cases = number_of_chosen_cases


    def run(self, write_result: bool=True):
        a = self.__modify_and_join_mcce_results__(in_csv_paths=self.in_ms_csvs)
        b = self.__reconstruct_dfs_from_agg_df__(a)[:self.number_of_chosen_cases]
        out_paths = []
        if write_result:
            for ind, _ in enumerate(b):
                out_path = os.path.join(self.out_dir, f"Agg_MCCE_run_{ind}.csv")
                df = b[ind]
                df.to_csv(out_path, index=False)
                print(f"Case {ind} of the aggregrate result has been writtedn to {out_path}.")
                out_paths.append(out_path)
        return b, out_paths



    def __gen_key__(self, row, columns):
        return '%'.join([str(elm) for elm in row.to_list()[:-3]])  + '|' + '%'.join(columns[:-3])


    def __modify_and_join_mcce_results__(self, in_csv_paths: list):
        dfs = []
        for in_csv_path in in_csv_paths:
            df = pd.read_csv(in_csv_path)
            cols = list(df.columns)
            df['Key'] = df.apply(lambda row: self.__gen_key__(row, columns=cols), axis=1)
            df_mod = df[['Key', 'Count', 'Occupancy', 'Sum_crg_protein']]
            dfs.append(df_mod)
        df_full = pd.concat(dfs)
        a = df_full.groupby('Key').aggregate(['mean'])
        b = a.sort_values(by=[('Count', 'mean')], ascending=False)
        return b


    def __reconstruct_resid_info__(self, resid_key):
        resid_crgs = [float(elm) for elm in resid_key.split("|")[0].split("%")]
        resid_cols = resid_key.split("|")[1].split("%")
        return resid_crgs, resid_cols


    def __reconstruct_df__(self, mcce_entry):
        resid_crgs, resid_cols = self.__reconstruct_resid_info__(mcce_entry[0])
        df_dict = {}
        for val, key in zip(resid_crgs, resid_cols):
            df_dict[key] = [val]
        df_dict['Count'] = [mcce_entry[1][0]]
        df_dict['Occupancy'] = [mcce_entry[1][1]]
        df_dict['Sum_crg_protein'] = [mcce_entry[1][2]]
        return pd.DataFrame(df_dict)


    def __reconstruct_dfs_from_agg_df__(self, agg_df):
        out_dfs = []
        indices = agg_df.index
        for ind, key in enumerate(indices):
            c = [key, agg_df.iloc[ind].to_list()]
            d = self.__reconstruct_df__(c)
            out_dfs.append(d)
        return out_dfs



class MSAnaToNAMD(object):
    '''
    Task: Convert a microstat analysis output from MCCE to
    structural files ready for NAMD runs.
    '''
    def __init__(
            self,
            in_ms_csv: Path=None,
            in_pdb: Path=None,
            in_topologies: list=None,
            protein_only: bool=True,
            default_topo_dir: Path='./param_Chun'
        ):
        '''
        Task: Initialize the object

        Args:
            in_ms_csv:
                - the input path for the CSV that holds MCCE microstate analyses results for coming NAMD run.

            in_pdb:
                - the input PDB file for MCCE calculation.

            in_topologies:
                - a list of paths pointing to topology files to b sourced.
        '''
        default_param_files = os.listdir(default_topo_dir)
        default_topo_files = [elm for elm in default_param_files if ".prm" not in elm]

        self.in_ms_csv = in_ms_csv
        self.in_ms_df = pd.read_csv(in_ms_csv)
        self.Info = Info()
        self.in_pdb = in_pdb
        self.topologies = default_topo_files + in_topologies if in_topologies\
        else default_topo_files
        self.protein_only = protein_only
        self.out_dir = os.path.dirname(self.in_ms_csv)


    def get_example(self):
        '''
        Task: Provide an example on the format for the input CSV.
        '''
        example_dict = {
            'GLUA0007_': [0.0, 0.0, 0.0, 0.0, 0.0],
            'HISA0015_': [1.0, 1.0, 1.0, 1.0, 1.0],
            'ASPA0018_': [0.0, 0.0, 0.0, -1.0, 0.0],
            'TYRA0020_': [0.0, 0.0, 0.0, 0.0, 0.0],
            'CTRA0129_': [0.0, 0.0, 0.0, 0.0, 0.0],
            'Count': [1129301.0, 190004.0, 53465.0, 36175.0, 36096.0],
            'Occupancy': [0.753, 0.127, 0.036, 0.024, 0.024],
            'Sum_crg_protein': [19.0, 18.0, 18.0, 18.0, 18.0],
        }
        self.example = pd.DataFrame(example_dict)
        print(self.example.head())


    def parse_info_from_ms_ana(self, number_of_chosen_cases=None):
        '''
        Task: Parse information from the input MMCE calculations and convert them
        into actionable lines for users and psfgen in VMD.

        Args:
            number_of_chosen_cases
                - The number of cases to be chosen from MCCE microstate analyses results.
        '''
        #--1st find residues involved, residues types involved, and their charge states present
        self.__get_chain_and_residue_types__()
        self.__get_per_type_charge_states__()
        #--2nd select a subset of data as being instructed
        if not number_of_chosen_cases or number_of_chosen_cases >= len(self.in_ms_df):
            self.sel_ms_df = self.in_ms_df.copy(deep=True)
            number_of_chosen_cases = len(self.sel_ms_df)
        elif number_of_chosen_cases > 0:
            self.sel_ms_df = self.in_ms_df.head(number_of_chosen_cases).copy(deep=True)
        else:
            print('The "number_of_chosen_cases" must be set to be larger than 0.')
        pass
        #--3rd collect messages for each case present in self.sel_ms_df
        actions_dict = defaultdict(lambda: [])
        for ind in range(0, number_of_chosen_cases, 1):
            actions = self.__get_message_and_psfgen_lines_for_a_row__(
                self.sel_ms_df.loc[[ind]]
            )
            actions_dict[f'MCCE_case_{str(ind)}'] = actions
        #--4th group actions by chain id for psfgen processing
        self.actions_by_chains = self.__sort_actions_by_chains__(actions_dict)
        print('Information from the input MCCE calculations has been parsed and can be accessed by the variable "self.actions_by_chains".')


    def print_action_messages(self):
        for case_key in self.actions_by_chains.keys():
            print(f'Actions to be taken for {case_key} are .......')
            for chain_key in self.actions_by_chains[case_key].keys():
                print(f'Actions for chain {chain_key} of {case_key}:')
                for ind, line in enumerate(self.actions_by_chains[case_key][chain_key]):
                    print(f'Action {ind}: {line}')
            print('''
Next case ....
''')


    def write_psfgen_file(self, out_dir: Path):
        '''
        Task: generate a compact TCL script for psfgen as to build the NAMD system from
        MCCE calculations.

        Args:
            out_dir:
                - the output directory to deposit psfgen files generated.
        '''
        self.out_dir = out_dir
        Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        #if not os.path.exists(out_dir):
        #    os.makedirs(out_dir)
        #--generate that for 1 case
        for case_key in self.actions_by_chains.keys():
            self.out_path = os.path.join(self.out_dir, f'{case_key}_mcce_to_namd.tcl')
            out_lines = self.__write_psfgen_file_for_one_single_case__(case_key)
            with open(self.out_path, 'w') as fout:
                fout.write(out_lines)
            print(f'The TCL scritp for converting MMCE calculations to NAMD for {case_key} has been deposited at {self.out_path}.')


    def __write_psfgen_file_for_one_single_case__(self, case_key):
        '''
        Task: Write out the psfgen file needed to build the system for 1 single case
        '''
        #-1st: split protein into chains
        split_pdb_by_chain_lines = f'mol new {self.in_pdb}\n\n'
        pdb_prefix = os.path.basename(self.in_pdb).replace('.pdb', '').strip()
        #pdb_dir = os.path.dirname(self.in_pdb)
        pdb_dir = self.out_dir
        for chain in self.chains:
            if self.protein_only:
                split_pdb_by_chain_lines += f'set sel [atomselect top "chain {chain} and protein"]\n$sel writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}.pdb\n\n'
            else:
                split_pdb_by_chain_lines += f'set sel [atomselect top "chain {chain}"]\n$sel writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}.pdb\n\n'

        #-2nd: rename PDB by chain
        rename_lines = ''
        for chain in self.chains:
            rename_lines += f'mol new {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}.pdb\n'
            for line in self.actions_by_chains[case_key][chain]:
                psfgen_priority = line[1]
                if psfgen_priority == 1:
                    rename_lines += line[2] + '\n'
            rename_lines += f'''set sel [atomselect top "all"]
$sel writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_renamed.pdb
\n'''

        #-3rd: build psf & patches for each chain
        psf_lines = ''
        for chain in self.chains:
            psf_lines += f'''package require psfgen
resetpsf\n\n'''
            if not self.topologies:
                pass
            else:
                for topo in self.topologies:
                    psf_lines += f'topology {topo}\n'
            psf_lines += '\n'
            psf_lines += f'segment SEG{chain} ' + '{\n' +\
            f'  pdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_renamed.pdb\n'
            for line in self.actions_by_chains[case_key][chain]:
                psfgen_priority = line[1]
                if psfgen_priority == 2 or psfgen_priority == 3:
                    psf_lines += '  ' + line[2] +'\n'
            psf_lines += '}\n\n'
            for line in self.actions_by_chains[case_key][chain]:
                psfgen_priority = line[1]
                if psfgen_priority == 4:
                    psf_lines += '  ' + line[2] +'\n'
            psf_lines += f'''
pdbalias atom ILE CD1 CD
pdbalias atom HOH O OH2
pdbalias residue HOH TIP3
coordpdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_renamed.pdb SEG{chain}

guesscoord

writepsf {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.psf
writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.pdb

'''
        #-4th: merge PSF & PDB from different chains
        merge_lines = ''
        if len(self.chains) == 1:
            merge_lines += f'''mol load psf {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.psf pdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.pdb
set al [atomselect top "all"]
$al writepdb {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.pdb
$al writepsf {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.psf
'''
        elif len(self.chains) > 1:
            merge_lines += '''package require psfgen
resetpsf

'''
            for chain in self.chains:
                merge_lines += f'readpsf  {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.psf pdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.pdb\n'
            merge_lines += f'''writepsf {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.psf
writepdb {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.pdb'''

        #-5th combine all lines
        psfgen_lines = split_pdb_by_chain_lines + rename_lines + psf_lines + merge_lines
        return psfgen_lines


    def __sort_actions_by_chains__(self, actions_dict: dict=None):
        '''
        Task: Group message for the same chain into 1 group.

        !! May NOT need that as VMD psfgen should be able to build psf for PDB with multiple chians
        as long as each chain is being loaded separately.
        1) just output each chain separately [read Catherine's notes;
        2) then load the corresponding PDB for psfgen during segment creation.
        '''
        if not actions_dict:
            print(f'No action is needed/provided for the current case.')
            return None
        else:
            for case_key in actions_dict.keys():
                lines_by_chains = defaultdict(lambda: [])
                in_lines = actions_dict[case_key]
                for chain in self.chains:
                    for line in in_lines:
                        chain_id = line[0][-1]
                        if chain_id == chain:
                            lines_by_chains[chain].append(line)
                    lines_by_chains[chain] = sorted(lines_by_chains[chain], key=lambda elm: elm[1])
                actions_dict[case_key] = lines_by_chains
        return actions_dict


    def __get_chain_and_residue_types__(self):
        '''
        Task: Get the type of residues involved in MCCE calculations.
        '''
        residues_involved = self.in_ms_df.columns[:-3]
        residue_types_involved = list(set(
            [key[:3] for key in residues_involved]
        ))
        self.residues_involved = residues_involved
        self.residue_types_involved = residue_types_involved
        self.chains = list(set([res[3] for res in self.residues_involved]))


    def __get_per_type_charge_states__(self):
        '''
        Task: Get the charge states involved for each residue type
        present in MCCE calculations.
        '''
        crg_dict = defaultdict(lambda: [])
        for res in self.residues_involved:
            resname = res[:3]
            crg_dict[resname].extend(self.in_ms_df[res].to_list())
        for resname in crg_dict.keys():
            crg_dict[resname] = list(set(crg_dict[resname]))
        self.charge_states_involved = crg_dict


    def __get_message_and_psfgen_lines_for_a_row__(self, in_single_row_df):
        '''
        Task: Generate action messages and TCL lines based on a single input row
        from a larger CSV table containing many MCCE calculations. Each row consists
        of many columns, depending of the number of residues being considered.

        Args:
            in_single_row_df:
                - A single row of self.sel_ms_df being casted in a form of dataFrame,
                following the format of the input example.
        '''
        actions = []
        for res in self.residues_involved:
            resname = res[:3]
            chain_id = res[3]
            resid = int(res[4:-1])
            charge_state = in_single_row_df[res].item()
            segment_id = 'SEG' + str(chain_id)
            segment_id, psfgen_priority, psfgen_line, message_line = \
            self.__write_message_and_psfgen_lines__(
                resname, chain_id, resid, charge_state, segment_id
            )
            actions.append((segment_id, psfgen_priority, psfgen_line, message_line))
        return actions


    def __write_message_and_psfgen_lines__(
            self,
            resname,
            chain_id,
            resid,
            charge_state,
            segment_id,
            his_default: str='ND1',
        ):
        '''
        Task: Generate messages and TCL lines based on info for the input residue.
        The messages are for users to take action if they are to be done manually.
        The lines are for creating a corresponding psfgen TCL file that will build the
        needed NAMD system through VMD.
        '''
        if resname != 'HIS':
            patch_to_use = self.Info.charmm_dict_for_MCCE[resname][str(charge_state)]
        elif resname == 'HIS':
            if float(charge_state) != 0:
                patch_to_use = self.Info.charmm_dict_for_MCCE[resname][str(charge_state)]
            elif int(charge_state) == 0:
                patch_to_use = self.Info.charmm_dict_for_MCCE[resname]["0.0"][his_default]
            else:
                print("Not implemented. Resolve to default neutral HIS.")
                patch_to_use = self.Info.charmm_dict_for_MCCE[resname]["0.0"][his_default]
        else:
            print("Not implemented. Resolve to default neutral HIS.")
            patch_to_use = self.Info.charmm_dict_for_MCCE['HIS']["0.0"][his_default]

        patch_key = patch_to_use[5:].strip()
        if resname == 'CTR':
            message_line = f'Apply patch "{patch_to_use}" to the C-terminal (resid {resid}) of chain {chain_id}.'
            psfgen_line = f'last {patch_key}'
            psfgen_priority = 3
        elif resname == 'NTR':
            #--This function does not work yet because NTR is not yet included in our dict.
            #message_line = f'Apply patch "{patch_to_use}" to the N-terminal (resid {resid}) of chain {chain_id}.'
            #psfgen_line = f'first {patch_key}'
            #message_line = None
            message_line = 'Not implemented.'
            #psfgen_line = None
            psfgen_line = 'Not implemented.'
            psfgen_priority = 2
        else:
            # assuming the resname has been documented as either a regular residue or a patch
            if 'RESI' in patch_to_use:
                #--Regular residue
                message_line = f'Name/Rename residue with resid {resid} of chain {chain_id} to {patch_key}.'
                psfgen_line = f'''set sel [atomselect top "resid {resid} and chain {chain_id}"]
$sel set resname {patch_key}'''
                psfgen_priority = 1
            elif 'PRES' in patch_to_use:
                message_line = f'Apply patch ({patch_to_use}) to resid {resid} of chain {chain_id}.'
                psfgen_line = f'patch {patch_key} {segment_id}:{resid}'
                psfgen_priority = 4
            else:
                #--not implemented
                pass
        return segment_id, psfgen_priority, psfgen_line, message_line

