from APBS_for_MCCE_Chun_ver.base.config_data import Config
import os
from timeit import default_timer as timer
import subprocess
from pathlib import Path
import re
import numpy as np
import sys


class APBSForMCCE(object):

    def __init__(
            self,
            in_pqrs: list,
            temperatures: list=[300],
            mol_dies: list=[12.0],
            sol_dies: list=[78.54],
            pos_salt_concs: list=[0.15],
            neg_salt_concs: list=[0.15],
            pos_ion_radii: list=[2.0],
            neg_ion_radii: list=[2.0],
            pos_ion_charges: list=[1.0],
            neg_ion_charges: list=[-1.0],
            calc_xyz_dims: list=[[100, 100, 100]],
            cg_xyz_dims: list=[[100, 100, 100]],
            cg_xyz_centers: list=['mol 1'],
            fg_xyz_dims: list=[[100, 100, 100]],
            fg_xyz_centers: list=['mol 1'],
            out_file_names: list=[None],
            cg_pad: float=20,
            fg_pad: float=10,
        ):
        '''
        Task: Initiate an object with necessary information for APBS runs.

        Args:
            The inputs here have the same meaning as for that in the Config object.
            One can check out definition there.
            The difference is that, here, eevery argument is now made to be a list as an inputs
            instead of a single input value.
            This addresses the need to work on multiple objects.

            P.S. For most the above keywords, the Config object has provided default values.
            Thus, a pre-set None value was given.
        '''
        self.in_pqrs = in_pqrs
        if any([not os.path.exists(in_pqr) for in_pqr in in_pqrs]):
            print(f'Some paths from {in_pqrs} are missing. Operations aborted.')
            sys.exit()
        self.work_dirs = [os.path.dirname(os.path.abspath(str(in_pqr))) for in_pqr in in_pqrs]
        self.pqr_files = [os.path.basename(str(in_pqr)) for in_pqr in in_pqrs]
        self.num_tasks = len(self.pqr_files)
        self.cg_pad = cg_pad
        self.fg_pad = fg_pad
        xyz_dims, cg_dims, fg_dims  = self.__gen_default_grid_dim_for_apbs__(
            cg_pad = self.cg_pad,
            fg_pad = self.fg_pad,
        )
        print(xyz_dims, cg_dims, fg_dims)
        self.mol_dims = xyz_dims
        self.calc_xyz_dims = calc_xyz_dims + cg_dims[len(calc_xyz_dims):] if calc_xyz_dims else cg_dims
        #calc_xyz_dims + [[100, 100, 100]] * int(self.num_tasks - len(calc_xyz_dims))
        self.cg_xyz_dims = cg_xyz_dims + cg_dims[len(cg_xyz_dims):] if cg_xyz_dims else cg_dims
        #cg_xyz_dims + [[100, 100, 100]] * int(self.num_tasks - len(cg_xyz_dims))
        self.cg_xyz_centers = cg_xyz_centers + ['mol 1'] * int(self.num_tasks - len(cg_xyz_centers)) if cg_xyz_centers else ['mol 1'] * int(self.num_tasks)
        #cg_xyz_centers + ['mol 1'] * int(self.num_tasks - len(cg_xyz_centers))
        self.fg_xyz_dims = fg_xyz_dims + fg_dims[len(fg_xyz_dims):] if fg_xyz_dims else fg_dims
        #fg_xyz_dims + [[100, 100, 100]] * int(self.num_tasks - len(fg_xyz_dims))
        self.fg_xyz_centers = fg_xyz_centers + ['mol 1'] * int(self.num_tasks - len(fg_xyz_centers)) if fg_xyz_centers else ['mol 1'] * int(self.num_tasks)
        #fg_xyz_centers + ['mol 1'] * int(self.num_tasks - len(fg_xyz_centers))
        self.temperatures = temperatures + [300] * int(self.num_tasks - len(temperatures))
        self.mol_dies = mol_dies + [12.0] * int(self.num_tasks - len(mol_dies))
        self.sol_dies = sol_dies + [78.54] * int(self.num_tasks - len(sol_dies))
        self.pos_salt_concs = pos_salt_concs + [0.15] * int(self.num_tasks - len(pos_salt_concs))
        self.neg_salt_concs = neg_salt_concs + [0.15] * int(self.num_tasks - len(neg_salt_concs))
        self.pos_ion_radii = pos_ion_radii + [2.0] * int(self.num_tasks - len(pos_ion_radii))
        self.neg_ion_radii = neg_ion_radii + [2.0] * int(self.num_tasks - len(neg_ion_radii))
        self.pos_ion_charges = pos_ion_charges + [1.0] * int(self.num_tasks - len(pos_ion_charges))
        self.neg_ion_charges = neg_ion_charges + [-1.0] * int(self.num_tasks - len(neg_ion_charges))
        self.out_file_names = out_file_names + [None] * int(self.num_tasks - len(out_file_names))
        pass


    def gen_apbs_configs(
            self,
            compute_mode: str='lpbe',
            write_mode: str='pot',
            write_format: str='dx',
        ):
        '''
        Task: Get the configuraiton text for each case to be computed.
        '''
        self.configs = []
        self.compute_mode = compute_mode
        self.write_mode = write_mode
        self.write_format = write_format
        config_obj = Config()
        for ind, in_dir in enumerate(self.work_dirs):
            self.configs.append(
                config_obj.get_apbs_config(
                    in_pqr=os.path.join(in_dir, self.pqr_files[ind]),
                    compute_mode=self.compute_mode,
                    temperature=self.temperatures[ind],
                    mol_die=self.mol_dies[ind],
                    sol_die=self.sol_dies[ind],
                    pos_salt_conc=self.pos_salt_concs[ind],
                    neg_salt_conc=self.neg_salt_concs[ind],
                    pos_ion_radius=self.pos_ion_radii[ind],
                    neg_ion_radius=self.neg_ion_radii[ind],
                    pos_ion_charge=self.pos_ion_charges[ind],
                    neg_ion_charge=self.neg_ion_charges[ind],
                    calc_xyz_dim=self.calc_xyz_dims[ind],
                    cg_xyz_dim=self.cg_xyz_dims[ind],
                    cg_xyz_center=self.cg_xyz_centers[ind],
                    fg_xyz_dim=self.fg_xyz_dims[ind],
                    fg_xyz_center=self.fg_xyz_centers[ind],
                    write_mode=self.write_mode,
                    write_format=self.write_format,
                    out_file_name=self.out_file_names[ind],
                )
            )


    def run_apbs(
            self,
            apbs_path,
            time_run: bool=False,
            run_apbs: bool=True,
        ):
        #--Write config
        if os.path.exists(apbs_path):
            self.apbs_path = apbs_path
        elif not run_apbs:
            self.apbs_path = apbs_path
        else:
            print(f"{apbs_path} is not found and 'run_apbs' is set to True.")
            sys.exit()
        self.apbs_config_file_paths = []
        for ind, (work_dir, pqr_file) in enumerate(zip(self.work_dirs, self.pqr_files)):
            #if not work_dir or work_dir == ".":
            #    work_dir = os.getcwd()
            apbs_config_file_path = os.path.join(
                work_dir, pqr_file.replace('.pqr', '.apbs')
            )
            #print(f"Writing APBS configuration file at: {apbs_config_file_path}")
            print("Writing APBS configuration file at: {}".format(apbs_config_file_path))
            with open(apbs_config_file_path, 'w') as fout:
                fout.write(self.configs[ind])

            current_dir = os.getcwd()
            print(current_dir)

            self.apbs_config_file_paths.append(apbs_config_file_path)
            if run_apbs:
                #cmd_1 = ' '.join(['cd', work_dir])
                #print(cmd_1)
                #result_1 = subprocess.run([cmd_1], shell=True, capture_output=True, text=True)
                #result_1 = subprocess.run([cmd_1], shell=True)
                #os.chdir(work_dir)
                print(f"Navigate to {os.path.dirname(apbs_config_file_path)}")
                os.chdir(os.path.dirname(apbs_config_file_path))
                #print(result_1.stdout)
                if time_run:
                    start = timer()
                cmd_1 = f'chmod +777 {self.apbs_path}'
                #result_1 = subprocess.run([cmd_1], shell=True, capture_output=True, text=True)
                subprocess.run(cmd_1, shell=True, capture_output=True, text=True)
                #cmd_2 = ' '.join([str(self.apbs_path), os.path.basename(apbs_config_file_path)])
                #cmd_2 = ' '.join([self.apbs_path, apbs_config_file_path])
                cmd_2 = ' '.join([self.apbs_path, apbs_config_file_path])
                print(cmd_2)
                #result_2 = subprocess.run([cmd_2], shell=True, capture_output=True, text=True)
                subprocess.run(cmd_2, shell=True, capture_output=True, text=True)
                #result_2 = subprocess.run([cmd_2], shell=True)
                #print(result_2.stdout)
                if time_run:
                    elapsed_time = timer() - start
                    #print(f"The last APBS run took {elapsed_time}s.")
                    print("The last APBS run took {}s.".format(elapsed_time))
                os.chdir(current_dir)
        pass


    def write_visualizable_objects(self, apbs_output_file_paths):
        self.vis_obj = []
        self.apbs_output_file_paths = apbs_output_file_paths
        for ind, elm in enumerate(zip(self.in_pqrs, self.apbs_output_file_paths)):
            with open(elm[0], 'r') as fin:
                pqr_lines = fin.readlines()
            with open(elm[1], 'r') as fin:
                pot_lines = fin.readlines()
            self.vis_obj.append(
                VisEMfromMCCE(
                    pqr_lines=pqr_lines,
                    pot_lines=pot_lines,
                    pot_type=self.write_mode,
                    pot_format=self.write_format,
                    pot_from='APBS',
                    convert_for='VMD',
                    work_dir=os.path.dirname(os.path.abspath(elm[1])),
                    file_prefix=os.path.basename(elm[1]),
                )
            )
            self.vis_obj[ind].write_visualizable_object()


    def __gen_default_grid_dim_for_apbs__(
        self,
        cg_pad,
        fg_pad,
        regex=r'ATOM\s+\d+\s+[A-Z]{1,2}\s+[A-Z]{3}\s+\d*\s+(-*\d+.\d+)\s+(-*\d+.\d+)\s+(-*\d+.\d+)'
    ):
        #--get minmax-for each in_pqr
        xyz_dims = []
        for in_pqr in self.in_pqrs:
            #--Get data
            with open(in_pqr, 'r') as fin:
                pqr_lines = fin.readlines()
            #--Convert to array
            vals = []
            for line in pqr_lines:
                match = re.search(regex, line.strip())
                if match:
                    vals.append([float(match[1]), float(match[2]), float(match[3])])
            val_arr = np.array(vals)
            #--Get dim
            xyz_dims.append(np.max(val_arr, axis=0) - np.min(val_arr, axis=0))
            cg_dims = [elm + cg_pad for elm in xyz_dims]
            fg_dims = [elm + fg_pad for elm in xyz_dims]
            #print("hi")
        return xyz_dims, cg_dims, fg_dims



class VisEMfromMCCE(object):
    '''
    This is a auxillary class to provide functions needed.
    '''
    def __init__(
        self,
        pqr_lines,
        pot_lines,
        pot_type,
        pot_format,
        pot_from: str='APBS',
        convert_for: str='VMD',
        implemented: list=['VMD'],
        work_dir=os.getcwd(),
        file_prefix=None,
    ):
        self.pqr_lines = pqr_lines
        self.pot_lines = pot_lines
        self.pot_type = pot_type
        self.pot_format = pot_format
        self.pot_from = pot_from
        self.work_dir = work_dir
        self.file_prefix = file_prefix
        if convert_for in implemented:
            self.convert_for = convert_for
        else:
            print(f"{convert_for} has not been implemented yet.")
            sys.exit()


    def write_visualizable_object(self):
        if self.pot_from == 'APBS':
            self.__write_visualizable_object_APBS__()
        else:
            print(f"Conversions for outputs from {self.pot_from} have not been implemented.")
            sys.exit()


    def __write_visualizable_object_APBS__(self):
        if self.pot_format == 'flat' and self.pot_type == 'atompot':
            self.converted_pot_lines = self.__APBS_convert_atompot_to_pqr__()
            print(f"APBS outputs with the format {self.pot_format} and the values for {self.pot_type} are best being visualized as pqr/pdb objects.")
            if not self.file_prefix:
                out_file_name = 'apbs_atompot.pqr'
            else:
                out_file_name = self.file_prefix + '.pqr'
            self.out_path = os.path.join(self.work_dir, out_file_name)
            with open(self.out_path, 'w') as fout:
                fout.write('\n'.join(self.converted_pot_lines))
            print(f"Conevrsion to pqr object for visualization at {self.out_path} has been completed.")
            print(f"To use, please load the pqr and color the atoms by its 'charge' attribute in VMD.")
        else:
            print('The format for the input is already visualizable; or its conversaion has not been implemented.')


    def __APBS_convert_atompot_to_pqr__(
        self,
        output_width: int=7,
        decimal_places: int=4,
        part_1_end_ind: int=54,
        part_2_beg_ind: int=64,
    ):
        #--extract atompot values
        atompot_vals = []
        for line in self.pot_lines:
            if line.strip():
                if not '#' in line.strip():
                    val = str(round(float(line.strip()), decimal_places))
                    to_pqr_val = ' ' * (output_width - len(val)) + val
                    atompot_vals.append(to_pqr_val)
        #--clearn pqr lines
        cleaned_pqr_lines = []
        for line in self.pqr_lines:
            if line.strip():
                cleaned_pqr_lines.append(line.strip())
        #--add atompot vals to pqr
        out_lines = []
        for ind, line in enumerate(cleaned_pqr_lines):
            out_lines.append(
                ' '.join([
                    line[:part_1_end_ind],
                    atompot_vals[ind],
                    line[part_2_beg_ind:]
                ])
            )
        return out_lines


#---The below needs to be changed in the future since
#--we want to visualize EM maps of various form without running psfgen or APSB
class VisEMfromMCCEtoNAMD(object):
    '''
    Task: Generate visualization scripts in VMD for APBS computed EM maps
    '''
    def __init__(
        self,
        in_psfgen_path,
        vmd_path: Path='/Common/linux/bin/vmd',
        apbs_path: Path='/Common/linux/bin/apbs',
        map_dimensions: list=[200, 200, 200],
        salinity: float=0.15,
        run_psfgen: bool=False,
        run_apbs: bool=False,
    ):
        self.in_psfgen = in_psfgen_path
        self.in_dir = os.path.dirname(in_psfgen_path)
        self.vmd = vmd_path
        self.map_dimensions = map_dimensions
        self.salinity = salinity
        self.apbs_path = apbs_path
        self.run_psfgen = run_psfgen
        self.run_apbs = run_apbs
        pass


    def get_APBS_EM_map_and_Vis(
        self
    ):
        #--Step 1: generate pqr path
        print("Obtaining paths info.")
        self.__get_output_psf_path_from_psfgen_script__()
        print(self.pqr_path)
        #--Step 2: modify psfgen script
        print("Adjusting outputs from step 2.")
        self.__modify_mcce_to_namd_psf_scripts__()
        #--Step 3: Run VMD to get pqr_file
        if self.run_psfgen:
            print("Calling psfgen.")
            self.__run_psfgen__()
        else:
            print("Not calling VMD psfgen to generate structures.")
        #--Step 4: compute APBS
        print("Calling APBS.")
        apbs_obj = APBSForMCCE(
            in_pqrs=[self.pqr_path],
            pos_salt_concs=[self.salinity],
            neg_salt_concs=[self.salinity],
            calc_xyz_dims=[self.map_dimensions],
            cg_xyz_dims=[self.map_dimensions],
            fg_xyz_dims=[self.map_dimensions],
        )
        apbs_obj.gen_apbs_configs()
        if self.run_apbs:
            print("Run APBS to get the EM map.")
            apbs_obj.run_apbs(
                apbs_path=self.apbs_path,
                run_apbs=True,
                time_run=False,
            )
        else:
            apbs_file = self.pqr_path.replace('.pqr', '.dx')
            print(f"Not running APBS. Please provide the EM map as being named as {apbs_file}")
        #--Step 5: Write Vis script that load the molecule and the map
        print("Producing visualization scripts.")
        self.__run_vis_script__()


    def __run_vis_script__(
        self,
        out_name: str='vis_option_5.tcl'
    ):
        if self.pqr_path:
            self.vis_path = os.path.join(os.path.dirname(self.pqr_path), out_name)
            out_text = f'''
mol load psf {self.psf_path} pdb {self.pdb_path}
mol modselect 0 top "noh"
mol modmaterial 0 top AOChalky
mol modcolor 0 top ColorID 8
mol modstyle 0 top NewCartoon 0.3 10.0 4.1
mol new {self.pqr_path.replace(".pqr", ".dx")}
mol modcolor 0 top ColorID 0
mol modstyle 0 top Isosurface 0.5 0 0 1 1 1
mol addrep top
mol modcolor 1 top ColorID 1
mol modstyle 0 top Isosurface -0.5 0 0 1 1 1
'''
            with open(self.vis_path, 'w') as fout:
                fout.write(out_text)
        else:
            print("No psf/pqr has been generated.")
            pass


    def __run_psfgen__(
        self,
    ):
        #result_vmd = subprocess.run([self.vmd, '-dispdev text <', self.in_psfgen], shell=True, capture_output=True, text=True)
        #print([str(self.vmd), '-dispdev text <', self.in_psfgen])
        cmd = ' '.join([str(self.vmd), '-dispdev text <', self.in_psfgen])
        subprocess.run(cmd, shell=True, capture_output=True, text=True)
        #subprocess.call([self.vmd, '-dispdev text <', self.in_psfgen])
        #result_vmd = subprocess.run([cmd_2], shell=True)
        #print(result_vmd.stdout)


    def __get_output_psf_path_from_psfgen_script__(
        self,
        regex=r'\$al writepsf (\S+).psf',
    ):
        with open(self.in_psfgen, 'r') as fin:
            in_text = fin.read()
        match = re.search(regex, in_text)
        self.pqr_path = match[1] + '.pqr' if match else None
        self.pdb_path = match[1] + '.pdb' if match else None
        self.psf_path = match[1] + '.psf' if match else None


    def __modify_mcce_to_namd_psf_scripts__(
        self,
    ):
        if self.pqr_path:
            out_text = f'''$al writepqr {self.pqr_path}\n'''
        else:
            out_text = ''
        with open(self.in_psfgen, 'a') as fout:
            fout.write(out_text)




class pdb_to_pqr(object):

    def __init__(
        self,
        in_psfgen_path,
    ):
        pass
    def pdb_to_psf():
        pass
    def pdb_psf_to_pqr():
        pass
    def write_psfgen():
        pass
