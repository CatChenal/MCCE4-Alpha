from mcce4.pdbio import *
from mcce4.geom import *

def center_protein(self):
    """Center the protein.
    """
    atoms = []
    for res in self.protein.residue:
        for conf in res.conf:
            for atom in conf.atom:
                atoms.append(atom)

    x_lo = x_up = atoms[0].xyz[0]
    y_lo = y_up = atoms[0].xyz[1]
    z_lo = z_up = atoms[0].xyz[2]
    for atom in atoms[1:]:
        if x_lo > atom.xyz[0]: x_lo = atom.xyz[0]
        if x_up < atom.xyz[0]: x_up = atom.xyz[0]
        if y_lo > atom.xyz[1]: y_lo = atom.xyz[1]
        if y_up < atom.xyz[1]: y_up = atom.xyz[1]
        if z_lo > atom.xyz[2]: z_lo = atom.xyz[2]
        if z_up < atom.xyz[2]: z_up = atom.xyz[2]

    c = ((x_lo+x_up)/2, (y_lo+y_up)/2, (z_lo+z_up)/2)
    logging.debug("Centering the structure by shifting (%.3f, %.3f, %.3f)" % (-c[0], -c[1], -c[2]))
    op = OPERATION()
    op.move((-c[0], -c[1], -c[2]))

    # dynamic flag _applied
    for atom in atoms:
        atom._applied = False

    for atom in atoms:
        if not atom._applied:  
            # atom is a reference, as atom may be shared by conformers, so we need to check this
            atom.xyz = op.apply(atom.xyz)
            atom._applied = True

    # clean up
    for atom in atoms:
        if hasattr(atom, "_applied"):
            del atom._applied


    return