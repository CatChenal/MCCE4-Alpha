import os

class Config(object):


    def __init__(self):
        self.help=''''''


    def get_apbs_config(
            self,
            in_pqr,
            compute_mode: str='lpbe',
            temperature: float=300.0,
            mol_die: float=12.0,
            sol_die: float=78.54,
            pos_salt_conc: float=0.15,
            neg_salt_conc: float=0.15,
            pos_ion_radius: float=2.0,
            neg_ion_radius: float=2.0,
            pos_ion_charge: float=1.0,
            neg_ion_charge: float=-1.0,
            calc_xyz_dim: list=[100, 100, 100],
            cg_xyz_dim: list=[100, 100, 100],
            cg_xyz_center='mol 1',
            fg_xyz_dim: list=[100, 100, 100],
            fg_xyz_center='mol 1',
            write_mode: str='pot',
            write_format: str='dx',
            out_file_name=None,
        ):
        '''
        Args:
            in_pqr:
                - The path to the input pqr

            compute_mode:
                - The computation object and method to be used in APBS:
                    (1) npbe: non-linear Poisson-Boltzman equation (gives potential);
                    (2) lpbe: linearized Poisson-Boltzman equation (gives potential);
                    (3) nrpbe: non-linear regularized Poisson-Boltzman equation (gives reaction field);
                    (4) lrpbe: linear regularized Poisson-Boltzman equation (gives reaction field);

            temperature:
                - The temperature for APBS calculations.

            mol_die:
                - The dielectric constant for biomolecules in the pqr.

            sol_die:
                - The solvent dielectric constant.

            pos_salt_conc:
                - Concentration of positive ions.

            neg_salt_conc:
                - Concentration of negative ions.

            pos_ion_radius:
                - The probing radius for positive ions.

            neg_ion_radius:
                - The probing radius for negative ions.

            pos_ion_charge:
                - Charges for positive ions.

            neg_ion_charge:
                - Charges for negative ions.

            calc_xyz_dim:
                - Number of grid points along X Y Z dimensions for APBS calculations.

            cg_xyz_dim:
                - X Y Z dimensions for the coarse grid.

            cg_xyz_center
                - Center of the coarse grid. It can be a list of x, y, z coordinates or a string of 'mol 1'.
                - Indicte using the input molecule's COM as center.

            fg_xyz_dim:
                - X Y Z diensions for the fine grid.

            fg_xyz_center:
                - Center of the fine grid.

            write_mode:
                - The object to be written out by APBS after a PB calculation:
                    (1) pot: the potential on grid points (can use multigrid & finite element) in KBT/e;
                    (2) charge: the charge distributions on multigrid;
                    (3) atompot: the potential at each atom on KBT/e;
                    (4) smol: solvent accessibility defined by molecular surface representation (unitless, from 0 to 1);
                    (5) sspl: spline-based solvent accessibility;
                    (6) vdw: van der Waals-based solvent accessibility;
                    (7) ivdw: the inflated van der Waals-based ion accessibility
                    (8) lap: the Laplacian of the potential (multigrid);
                    (9) edens: energy density;
                    (10) ndens: the total mobile ion number density for all ion species in units of M (multigrid only);
                    (11) qden: the total mobile ion charge density for all ion species in units of ec M;
                    (12) dielx/diely/dielz: the dielectric map shifted by 1/2 grid spacing in the {x, y, z}-direction

            write_format:
                - The format for writing APBS outputs:
                    (1) dx: regular output
                    (2) flat:

            out_file_name:
                - If None, a default name will be generated.
        '''
        self.work_dir = os.path.dirname(in_pqr)
        self.pqr_file = os.path.basename(in_pqr)
        self.temperature = temperature
        self.mol_die = mol_die
        self.sol_die = sol_die
        self.pos_salt_conc = pos_salt_conc
        self.neg_salt_conc = neg_salt_conc
        self.pos_ion_radius = pos_ion_radius
        self.neg_ion_radius = neg_ion_radius
        self.pos_ion_charge = pos_ion_charge
        self.neg_ion_charge = neg_ion_charge
        self.calc_xyz_dim = calc_xyz_dim
        self.cg_xyz_dim = cg_xyz_dim
        self.cg_xyz_center = cg_xyz_center
        self.fg_xyz_dim = fg_xyz_dim
        self.fg_xyz_center = fg_xyz_center
        self.compute_mode = compute_mode
        self.write_mode = write_mode
        self.write_format = write_format
        if out_file_name:
            self.out_path = os.path.join(self.work_dir, out_file_name)
        else:
            self.out_path = os.path.join(self.work_dir, self.pqr_file.replace('.pqr' , ''))

        calc_xyz_dim_in = ' '.join([str(elm) for elm in calc_xyz_dim])
        cg_xyz_dim_in = ' '.join([str(elm) for elm in cg_xyz_dim])
        if isinstance(cg_xyz_center, list):
            cg_xyz_center_in = ' '.join([str(elm) for elm in cg_xyz_center])
        else:
            cg_xyz_center_in = cg_xyz_center
        fg_xyz_dim_in = ' '.join([str(elm) for elm in fg_xyz_dim])
        if isinstance(fg_xyz_center, list):
            fg_xyz_center_in = ' '.join([str(elm) for elm in fg_xyz_center])
        else:
            fg_xyz_center_in = fg_xyz_center

        self.config_text = '''read
mol pqr {}
end
elec
mg-auto
dime {}
cglen {}
cgcent {}
fglen {}
fgcent {}
mol 1
{}
bcfl sdh
srfm smol
chgm spl2
ion {} {} {}
ion {} {} {}
pdie  {}
sdie  {}
sdens  10.0
srad  1.4
swin  0.3
temp  {}
gamma  0.105
calcenergy no
calcforce no
write {} {} {}
end
quit'''.format(
    self.pqr_file,
    calc_xyz_dim_in,
    cg_xyz_dim_in,
    cg_xyz_center_in,
    fg_xyz_dim_in,
    fg_xyz_center_in,
    self.compute_mode,
    self.pos_ion_charge,
    self.pos_salt_conc,
    self.pos_ion_radius,
    self.neg_ion_charge,
    self.neg_salt_conc,
    self.neg_ion_radius,
    self.mol_die,
    self.sol_die,
    self.temperature,
    self.write_mode,
    self.write_format,
    os.path.basename(self.out_path),
)
        return self.config_text

