

class Info(object):
    '''
    Task: Construct an object that holds information and data necessary for
    converting langauges used by MCCE to those for NAMD.
    '''

    def __init__(self):
        '''
        Task: Initialize the object as being described above.
        '''
        self.help = '''
        Hold information and data necessary for converting langauges used by MCCE to those for NAMD.

        self.charmm_dict_for_MCCE holds info. to match residue name and its charge state from MCCE to topology names in Charmm force field.

        self.charmm_dict_for_PDB holds info. to match residue names and atom types from regular PDBs to topology names in Charmm force field.
        '''
        self.charmm_dict_for_MCCE = {
            'GLU': {
                '-1.0': 'RESI GLU',
                '0.0': 'PRES GLUP',
            },
            'HIS': {
                '0.0': {
                    'ND1': 'RESI HSD',
                    'NE2': 'RESI HSE',
                },
                '1.0': 'RESI HSP',
            },
            'ASP': {
                '-1.0': 'RESI ASP',
                '0.0': 'PRES ASPP',
            },
            'TYR': {
                '0.0': 'RESI TYR'
            },
            'SER': {
                '0.0': 'RESI SER'
            },
            'THR': {
                '0.0': 'RESI THR',
            },
            'CTR': {
                '-1.0': 'PRES CTER',
                '0.0': 'PRES CT3'
            }
        }
        self.charmm_dict_for_PDB = {}

