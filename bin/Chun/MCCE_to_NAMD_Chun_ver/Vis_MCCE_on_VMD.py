import pandas as pd
import numpy as np
import re
import pprint
import os
from pathlib import Path



class Vis_MCCE_on_VMD(object):
    '''
    Task: Initiate an object to write visualization scripts for MCCE on VMD.
    '''
    def __init__(
        self,
        in_pdb_path,
        vis_option=None,
        in_ms_ana_csv_path: Path=None,
        in_pK_out_path: Path=None,
        in_res_pair_lst_path: Path=None,
        out_dir: Path='.',
        avail_opt: list=[1, 2, 3, 4],
    ):
        '''
        Task: Initiate an object for visualizing MCCE results on VMD.

        Agrs:
            in_pdb_path:
                - The path to the PDB that the visualized MCCE run was performed on.

            vis_option,
                - The visualization option to be chosen. Different options give out different VMD scripts.

            in_ms_ana_csv_path:
                - The path to the MS analysis results, which is hosted as a CVS.
                    (By default: this is the same CVS used in MCCE_to_NAMD)

            in_pK_out_path:
                - The path of "pK.out" file produced by MCCE runs.

            in_res_pair_lst_path
                - The path of "respair.lst" produced by MCCE runs

            out_dir:
                - The directory that will host all VMD visualization scripts generated.

            avail_opt:
                - Available visualization options.
        '''
        self.in_pdb_path = in_pdb_path
        self.in_pdb_name = os.path.basename(self.in_pdb_path).replace('.pdb', '')

        if not vis_option or vis_option == 'all':
            print("Visualization options are set to all available default options.")
            self.vis_option = avail_opt
        else:
            if not isinstance(vis_option, list):
                vis_option = [vis_option]

            self.vis_option = [elm for elm in vis_option if elm in avail_opt]

            if len(self.vis_option) < len(vis_option):
                print("Only visualization options available are retained.")

        self.in_ms_ana_csv_path = in_ms_ana_csv_path
        self.in_pK_out_path = in_pK_out_path
        self.in_res_pair_lst_path = in_res_pair_lst_path
        self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.option_details = {
            1: {
                "Explanation": '''Highlight residues demonstrating multiple protonate/charge states in MCCE calculations.
The color of each highlighted residue shows the standard deviation (STD) across charge states found.
The deeper the blue color is on a highlighted residue, the higher the sampled STD is.
''',
                "Internal_command": self.__gen_vmd_vis_script_opt_1__,
                "Input_file": self.in_ms_ana_csv_path,
                "Output_file": os.path.join(self.out_dir, 'vis_option_1.tcl'),
            },
            2: {
                "Explanation": '''Show the likelihood for each highlighted residue to be in its deprotonated or protonated state.
A blue color indicates a protonated state; and a red color indicates a deprotonated state.
''',
                 "Internal_command": self.__gen_vmd_vis_script_opt_2__,
                 "Input_file": self.in_ms_ana_csv_path,
                 "Output_file": os.path.join(self.out_dir, 'vis_option_2.tcl'),
            },
            3: {
                "Explanation": '''Show the pKa of each residue sampled in MCCE runs.
The deeper the blue color of a residue, the higher its pKa is.
''',
                 "Internal_command": self.__gen_vmd_vis_script_opt_3__,
                 "Input_file": self.in_pK_out_path,
                 "Output_file": os.path.join(self.out_dir, 'vis_option_3.tcl'),
            },
            4: {
                "Explanation": '''A visualizatoin will be created for each residue sampled by MCCE.
In each of this visualization, neighboring residues of the residues highlighted through VDW representation
will be colored according to their interaction (energy) with the hightlighted residue.
Stabilizing residues (with negative itneraction energy) are shown in red;
residues with positive interaction energy are shown in blue.
''',
                "Internal_command": self.__gen_vmd_vis_script_opt_4__,
                "Input_file": self.in_res_pair_lst_path,
                "Output_file": os.path.join(self.out_dir, 'vis_option_4.tcl'),
            }
        }


    def gen_VMD_vis_script(
        self,
    ):
        for index, vis_opt in enumerate(self.vis_option):
            #---something here
            #--check if the needed input is here, then write scripts
            input_file_path = self.option_details[vis_opt]["Input_file"]
            if input_file_path:
                vis_fn = self.option_details[vis_opt]["Internal_command"]
                vis_fn()
                print(f"Visualization opton {index + 1}:")
                print(f'VMD script(s) are ready at: {self.option_details[vis_opt]["Output_file"]}')
                #pprint.pprint(self.option_details[vis_opt]["Explanation"])
                print(self.option_details[vis_opt]["Explanation"])
            else:
                print(f"{input_file_path} is missing. Generation of VMD visualization script for visualization option {vis_opt} aborts.")
            pass


    def __gen_vmd_vis_script_opt_4__(
        self,
    ):
        in_res_pair_lst_path = self.in_res_pair_lst_path
        res_pair_df = pd.read_csv(
            in_res_pair_lst_path, header=0, delim_whitespace=True
        )
        for_vmd_results = self.__get_respair_interaction__(
            in_res_pair_df=res_pair_df,
        )
        out_text = self.__gen_vmd_vis_lines_opt_4_full__(
            for_vmd_results=for_vmd_results,
            in_pdb_path=self.in_pdb_path,
            mid_pt=0.5
        )
        out_path = os.path.join(self.out_dir, 'vis_option_4.tcl')
        with open(out_path, 'w') as fout:
            fout.write(out_text)


    def __get_respair_interaction__(
        self,
        in_res_pair_df,
        res_ref_col_index: int=0,
        res_partner_col_index: int=1,
        interaction_energy_col_index: int=4,
        res_pattern: str=r"([A-Z]{3})[\+-]*([A-Za-z]+)0*([1-9][0-9]*)"
    ):
        vmd_data_commands_per_res = {}
        unique_res_info = set(in_res_pair_df.iloc[:, res_ref_col_index].to_list())
        for res in unique_res_info:
            sel_df = in_res_pair_df[in_res_pair_df.iloc[:, res_ref_col_index] == res]
            partner_names = sel_df.iloc[:, res_partner_col_index].to_list()
            int_energies = [float(elm) for elm in sel_df.iloc[:, interaction_energy_col_index].to_list()]
            vmd_data_commands = []
            for res_partner, eng in zip(partner_names, int_energies):
                res_partner_info = self.__mcce_name_to_vmd_comand__(in_mcce_name=res_partner, pattern=res_pattern)
                vmd_data_commands.append((eng, res_partner_info['suggested_vmd_command']))
            vmd_data_commands_per_res[res] = {
                "for_ref_residue": self.__mcce_name_to_vmd_comand__(
                    in_mcce_name=res, pattern=res_pattern
                ),
                "for_partners": vmd_data_commands
            }
        return vmd_data_commands_per_res


    def __gen_vmd_vis_lines_opt_4__(
        self,
        for_vmd_results,
        in_pdb_path,
        mid_pt=None,
    ):
        out_text_dict = {}
        if for_vmd_results:
            for key in for_vmd_results.keys():
                val_list = []
                per_res_vmd_results = for_vmd_results[key]["for_partners"]
                mol_name = self.in_pdb_name + key
                out_text = f'''mol new {in_pdb_path}
set al [atomselect top "all"]
$al set beta 0
$al set occupancy 0
'''
                ref_vmd_line = for_vmd_results[key]["for_ref_residue"]['suggested_vmd_command']
                out_text += f'''{ref_vmd_line}
$sel set occupancy 1
'''
                for vmd_line in per_res_vmd_results:
                    out_text += f'''{vmd_line[1]}
$sel set beta {vmd_line[0]}
'''
                    val_list.append(float(vmd_line[0]))
                val_arr = np.array(val_list)
                if not mid_pt:
                    mid_pt = (0.5 * (val_arr.min() + val_arr.max()) - val_arr.min())/(val_arr.max() - val_arr.min())
                out_text += f'''
mol modselect 0 top "noh"
mol modmaterial 0 top AOChalky
mol modcolor 0 top ColorID 8
mol modstyle 0 top NewCartoon 0.3 10.0 4.1

mol addrep top
mol modselect 1 top "(not beta = 0) and occupancy = 0"
mol modmaterial 1 top AOChalky
mol modcolor 1 top Beta
mol modstyle 1 top Licorice 0.3 12.0 12.0

mol addrep top
mol modselect 2 top "not occupancy = 0 and noh"
mol modmaterial 2 top AOChalky
mol modcolor 2 top ResType
mol modstyle 2 top VDW 1.0 12.0

color scale midpoint {mid_pt}
color scale min 0.1
color scale max 0.0
mol rename top {mol_name}
'''
                out_text_dict[key] = out_text

        return out_text_dict


    def __gen_vmd_vis_lines_opt_4_full__(
        self,
        for_vmd_results,
        in_pdb_path,
        mid_pt=None,
    ):
        out_text_dict = self.__gen_vmd_vis_lines_opt_4__(
            for_vmd_results=for_vmd_results,
            in_pdb_path=in_pdb_path,
            mid_pt=mid_pt
        )
        out_text = ''
        for key in out_text_dict.keys():
            out_text += f"{out_text_dict[key]}\nmol off top\n"

        return out_text


    def __gen_vmd_vis_script_opt_3__(
        self,
    ):
        in_pK_out_path = self.in_pK_out_path
        pK_df = pd.read_csv(
            in_pK_out_path, header=0, delim_whitespace=True
        )
        for_vmd_results = self.__get_residues_pKa__(
            in_pK_df=pK_df,
        )
        out_text = self.__gen_vmd_vis_lines_opt_3__(
            for_vmd_results=for_vmd_results,
            in_pdb_path=self.in_pdb_path,
            mid_pt=None
        )
        out_path = os.path.join(self.out_dir, 'vis_option_3.tcl')
        with open(out_path, 'w') as fout:
            fout.write(out_text)


    def __get_residues_pKa__(
        self,
        in_pK_df,
        residue_name_col_index: int=0,
        pK_col_index: int=1,
        res_pattern: str=r"([A-Z]{3})[\+-]([A-Za-z]+)0*([1-9][0-9]*)"
    ):
        vmd_data_commands = []
        res_info = in_pK_df.iloc[:, residue_name_col_index].to_list()
        pK_info = [float(elm.replace('>', '').replace('<', '')) for elm in in_pK_df.iloc[:, pK_col_index].to_list()]

        for res, pKa in zip(res_info, pK_info):
            resid_info = self.__mcce_name_to_vmd_comand__(in_mcce_name=res, pattern=res_pattern)
            vmd_data_commands.append((pKa, resid_info['suggested_vmd_command']))
        return vmd_data_commands


    def __gen_vmd_vis_lines_opt_3__(
        self,
        for_vmd_results,
        in_pdb_path,
        mid_pt=None,
    ):
        if for_vmd_results:
            mol_name = f'{self.in_pdb_name}_show_pKa'
            out_text = f'''mol new {in_pdb_path}
set al [atomselect top "all"]
$al set beta 0
'''
            val_list = []
            for vmd_line in for_vmd_results:
                out_text += f'''{vmd_line[1]}
$sel set beta {vmd_line[0]}
'''
                val_list.append(float(vmd_line[0]))
            val_arr = np.array(val_list)
            if not mid_pt:
                mid_pt = (0.5 * (val_arr.min() + val_arr.max()) - val_arr.min())/(val_arr.max() - val_arr.min())

            out_text += f'''
mol modselect 0 top "noh"
mol modmaterial 0 top AOChalky
mol modcolor 0 top ColorID 8
mol modstyle 0 top NewCartoon 0.3 10.0 4.1

mol addrep top
mol modselect 1 top "not beta = 0 and noh"
mol modmaterial 1 top AOChalky
mol modcolor 1 top Beta
mol modstyle 1 top Licorice 0.3 12.0 12.0

color scale midpoint {mid_pt}
color scale min 0.1
color scale max 0.0
mol rename top {mol_name}
'''
            return out_text
        else:
            return ''


    def __gen_vmd_vis_script_opt_2__(
        self,
    ):
        in_ms_ana_csv_path = self.in_ms_ana_csv_path
        ms_ana_df = pd.read_csv(
            in_ms_ana_csv_path
        )
        for_vmd_results = self.__compute_residues_protonation_likelihood__(in_ms_ana_df=ms_ana_df)
        out_text = self.__gen_vmd_vis_lines_opt_2__(
            for_vmd_results=for_vmd_results,
            in_pdb_path=self.in_pdb_path,
        )
        out_path = os.path.join(self.out_dir, 'vis_option_2.tcl')
        with open(out_path, 'w') as fout:
            fout.write(out_text)


    def __compute_residues_protonation_likelihood__(
        self,
        in_ms_ana_df,
        not_needed_cols: int=3
    ):
        '''
        Args:
            in_ms_ana_df:
                - The Dataframe of the input MS analysis results.

            not_needed_cols:
                - The number of end columns to be ommitted for handling "in_ms_ana_df".
        '''
        vmd_data_commands = []
        for col in in_ms_ana_df.columns[:-not_needed_cols]:
            col_vals = np.array(in_ms_ana_df[col].to_list())
            col_rate = col_vals.sum()/col_vals.shape[0]
            if col_rate != 0:
                resid_info = self.__mcce_name_to_vmd_comand__(in_mcce_name=col)
                vmd_data_commands.append((col_rate, resid_info['suggested_vmd_command']))
        return vmd_data_commands


    def __gen_vmd_vis_lines_opt_2__(
        self,
        for_vmd_results,
        in_pdb_path,
        mid_pt=0.5,
    ):
        if for_vmd_results:
            mol_name = f'{self.in_pdb_name}_show_protonation_likelihood'
            out_text = f'''mol new {in_pdb_path}
set al [atomselect top "all"]
$al set beta 0
'''
            for vmd_line in for_vmd_results:
                out_text += f'''{vmd_line[1]}
$sel set beta {vmd_line[0]}
'''
            out_text += f'''
mol modselect 0 top "noh"
mol modmaterial 0 top AOChalky
mol modcolor 0 top ColorID 8
mol modstyle 0 top NewCartoon 0.3 10.0 4.1

mol addrep top
mol modselect 1 top "not beta = 0 and noh"
mol modmaterial 1 top AOChalky
mol modcolor 1 top Beta
mol modstyle 1 top Licorice 0.3 12.0 12.0

color scale midpoint {mid_pt}
color scale min 0.1
color scale max 0.0
mol rename top {mol_name}
'''
            return out_text
        else:
            return ''


    def __gen_vmd_vis_script_opt_1__(
        self,
    ):
        in_ms_ana_csv_path = self.in_ms_ana_csv_path
        ms_ana_df = pd.read_csv(
            in_ms_ana_csv_path
        )
        for_vmd_results = self.__pick_residues_with_multiple_charge_states__(
            in_ms_ana_df=ms_ana_df
        )
        out_text = self.__gen_vmd_vis_lines_opt_1__(
            for_vmd_results=for_vmd_results,
            in_pdb_path=self.in_pdb_path
        )
        out_path = os.path.join(self.out_dir, 'vis_option_1.tcl')
        with open(out_path, 'w') as fout:
            fout.write(out_text)


    def __gen_vmd_vis_lines_opt_1__(
        self,
        for_vmd_results,
        in_pdb_path,
        mid_pt=0.0,
    ):
        '''
        Args:
            for_vmd_results:
                - Organized results that are ready to be converted into TCL commands on VMD.

            in_pdb_path:
                - The input PDB where MCCE runs were performed on.

            mid_point:
                - A parameter controlling the color scale used in VMD.
                It is used to enhance aesthetic view for users.
        '''
        if for_vmd_results:

            mol_name = f'{self.in_pdb_name}_highlight_charge_states_diversity'
            out_text = f'''mol new {in_pdb_path}
set al [atomselect top "all"]
$al set beta 0
'''
            for vmd_line in for_vmd_results:
                out_text += f'''{vmd_line[1]}
$sel set beta {vmd_line[0]}
'''
            out_text += f'''
mol modselect 0 top "noh"
mol modmaterial 0 top AOChalky
mol modcolor 0 top ColorID 8
mol modstyle 0 top NewCartoon 0.3 10.0 4.1

mol addrep top
mol modselect 1 top "not beta = 0 and noh"
mol modmaterial 1 top AOChalky
mol modcolor 1 top Beta
mol modstyle 1 top Licorice 0.3 12.0 12.0

color scale midpoint {mid_pt}
color scale min 0.1
color scale max 0.0
mol rename top {mol_name}
'''
            return out_text
        else:
            return ''


    def __pick_residues_with_multiple_charge_states__(
        self,
        in_ms_ana_df,
        not_needed_cols: int=3
    ):
        '''
        Args:
            in_ms_ana_df:
                - The dataframe onject of the input MS analysis table.

            not_needed_cols:
                - The number of end columns to be ommitted for handling "in_ms_ana_df".
        '''
        vmd_data_commands = []
        for col in in_ms_ana_df.columns[:-not_needed_cols]:
            col_vals = np.array(in_ms_ana_df[col].to_list())
            col_std = col_vals.std()
            if col_std > 0:
                resid_info = self.__mcce_name_to_vmd_comand__(in_mcce_name=col)
                vmd_data_commands.append((col_std, resid_info['suggested_vmd_command']))
        return vmd_data_commands


    def __mcce_name_to_vmd_comand__(
        self,
        in_mcce_name,
        pattern: str=r"([A-Z]{3})([A-Za-z]+)0*([1-9][0-9]*)"
    ):
        '''
        Args:
            in_mcce_name:
                - The residue name used in showing MCCE results.

            pattern:
                - A regular expression pattern used to extract information from "in_mcce_name".
        '''
        m = re.search(pattern, in_mcce_name)
        if m:
            return {
                'resname': m.group(1),
                'chain_id': m.group(2),
                'resid': m.group(3),
                'suggested_vmd_command': f'''set sel [atomselect top "chain {m.group(2)} and resname {m.group(1)} and resid {m.group(3)}"]'''
            }
        else:
            return None

