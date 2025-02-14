#!/usr/bin/env python
# vdw functions
import numpy as np
from ..geom import *
from ..constants import VDW_CUTOFF_FAR
from ..constants import VDW_CUTOFF_NEAR
from ..constants import VDW_UPLIMIT
from ..constants import VDW_SCALE14
import math

VDW_CUTOFF_FAR2 = VDW_CUTOFF_FAR * VDW_CUTOFF_FAR
VDW_CUTOFF_NEAR2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR

class Blob:
    """
    Define a sphere with center and radius to cover all atoms in a cobnformer including their vdw_r
    """
    def __init__(self, conf):
        self.center = (0.0, 0.0, 0.0)
        self.radius = 0.0
        x = 0.0
        y = 0.0
        z = 0.0
        for atom in conf.atom:
            x += atom.xyz[0]
            y += atom.xyz[1]
            z += atom.xyz[2]
        n = len(conf.atom)

        if n > 0:
            self.center = (x/n, y/n, z/n)

            r_max = max([x.r_vdw for x in conf.atom])

            d_far2 = 0.0
            for atom in conf.atom:
                d2 = ddvv(self.center, atom.xyz)
                if d2 > d_far2:
                    d_far2 = d2
            d_far = math.sqrt(d_far2)

            self.radius = d_far + r_max



def assign_vdw_param(self):
    """Assign vdw r0 and eps to atoms and make blob for conformers.
    A blob is a sphere that a conformer may have significant vdw interaction.
    """
    for res in self.protein.residue:
        for conf in res.conf:
            for atom in conf.atom:
                key = ("RADIUS", conf.confType, atom.name)
                atom.r_vdw = self.tpl.db[key].r_vdw
                atom.e_vdw = self.tpl.db[key].e_vdw
    
    # make blob based on 

def make_blob(self):
    for res in self.protein.residue:
        for conf in res.conf:
            conf.blob  = Blob(conf)


def vdw_conf(conf1, conf2, cutoff=0.001, verbose=False, display=False, dirty=False):
    vdw = 0.0

    d = dvv(conf1.blob.center, conf2.blob.center)
    #print(d, conf1.blob.radius + 6 + conf2.blob.radius)
    if d > conf1.blob.radius + 6 + conf2.blob.radius:
        if display:
            print("%s at (%.3f, %.3f, %.3f) <-> %s at (%.3f, %.3f, %.3f): d = %.3f" %
                  (conf1.confID, conf1.blob.center[0], conf1.blob.center[1], conf1.blob.center[2],
                   conf2.confID, conf2.blob.center[0], conf2.blob.center[1], conf2.blob.center[2],
                   d))

    else:
        if display and verbose:
            print("%12s %16s     %10s %8s %6s   %6s %6s %6s %6s %6s %6s" % ("ATOM1",
                                                                 "ATOM2",
                                                                 "vdw",
                                                                 "dist",
                                                                 "cnct",
                                                                 "r1",
                                                                 "e1",
                                                                 "r2",
                                                                 "e2",
                                                                            "R_sum",
                                                                            "E_par"))
        for atom1 in conf1.atom:
            for atom2 in conf2.atom:
                if dirty:
                    vdw_a2a = vdw_atom_dirty(atom1, atom2)
                else:
                    vdw_a2a = vdw_atom(atom1, atom2)
                vdw += vdw_a2a
                if verbose and abs(vdw_a2a) >= cutoff:
                    d2 = ddvv(atom1.xyz, atom2.xyz)
                    if atom1 == atom2:
                        connect = "self"
                    elif atom1 in atom2.connect12:
                        connect = "1--2"
                    elif atom1 in atom2.connect13:
                        connect = "1--3"
                    elif atom1 in atom2.connect14:
                        connect = "1--4"
                    else:
                        connect = "none"
                    if display:  # display details
                        R_sum = atom1.r_vdw+atom2.r_vdw
                        E_par = np.sqrt(atom1.e_vdw*atom2.e_vdw)
                        print("%s -> %s: %8.3f %8.3f %6s   %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f" % (atom1.atomID,
                                                                                                   atom2.atomID,
                                                                                                   vdw_a2a,
                                                                                                   np.sqrt(d2),
                                                                                                   connect,
                                                                                                   atom1.r_vdw,
                                                                                                   atom1.e_vdw,
                                                                                                   atom2.r_vdw,
                                                                                                   atom2.e_vdw,
                                                                                                   R_sum,
                                                                                                   E_par))
                    else:  # display essential information
                        print(
                            "%s -> %s: %8.3f %8.3f %6s" % (atom1.atomID, atom2.atomID, vdw_a2a, np.sqrt(d2), connect))
        if vdw >= VDW_UPLIMIT:
            vdw = 999.0

        if conf1 == conf2:
            vdw = 0.5*vdw

    return vdw

def vdw_atom(atom1, atom2):
    # A good post: https://mattermodeling.stackexchange.com/questions/4845/how-to-create-a-lookup-table-of-%CF%B5-and-%CF%83-values-for-lennard-jones-potentials
    # Parameter source: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
    #   CHARM36-jul2022.ff/ffnonbonded.itp
    #   σ and ε values are in nm and kJ/mol in this file
    # Lorentz-Berthelot combining rules:
    #   σij = 0.5*(σi + σj)
    #   ϵij = sqrt(ϵi * ϵj)
    # p_lj = ϵij[(σij/r)^12 - 2(σij/r)^6]
    # σij is the distance where LJ potential reaches minimum: -ϵij
    # r is the atom distance

    p_lj = 0.0

    # Using the equation in vdw.c. Need to work on new parameter set.
    if abs(atom1.xyz[0] - atom2.xyz[0]) < VDW_CUTOFF_FAR and \
       abs(atom1.xyz[1] - atom2.xyz[1]) < VDW_CUTOFF_FAR and \
       abs(atom1.xyz[2] - atom2.xyz[2]) < VDW_CUTOFF_FAR:

        if atom1 != atom2 and atom2 not in atom1.connect12 and atom2 not in atom1.connect13:
            d2 = ddvv(atom1.xyz, atom2.xyz)
            if d2 > VDW_CUTOFF_FAR2:
                p_lj = 0.0
            elif d2 < VDW_CUTOFF_NEAR2:
                p_lj = 999.0
            else:
                r1 = atom1.r_vdw
                e1 = atom1.e_vdw
                r2 = atom2.r_vdw
                e2 = atom2.e_vdw
                if atom2 in atom1.connect14:
                    scale = VDW_SCALE14
                else:
                    scale = 1.0

                sig_min = r1 + r2
                eps = math.sqrt(e1 * e2)

                sig_d2 = sig_min * sig_min / d2
                sig_d6 = sig_d2 * sig_d2 * sig_d2
                sig_d12 = sig_d6 * sig_d6

                p_lj = scale * (eps * sig_d12 - 2. * eps * sig_d6)
                # #print("===%s===%s===" % (atom1.atomID, atom2.atomID))
                # if (atom1.atomID == " HB2ASP0018A001" and atom2.atomID == " OD2ASP0018A001") or \
                #    (atom1.atomID == " HB2ASP0018A003" and atom2.atomID == " OD2ASP0018A003"):
                #     print("===%s -> %s: %8.3f===" % (atom1.atomID, atom2.atomID, p_lj))
                #     print("%s: r_vdw=%8.3f, e_vdw=%8.3f" % (atom1.atomID, atom1.r_vdw, atom1.e_vdw))
                #     print("%s: r_vdw=%8.3f, e_vdw=%8.3f" % (atom2.atomID, atom2.r_vdw, atom2.e_vdw))


    return p_lj

def vdw_atom_dirty(atom1, atom2):
    eps = 0.16
    r0 = 3.8
    p_lj = 0.0

    # Using the equation in vdw.c. Need to work on new parameter set.
    if abs(atom1.xyz[0] - atom2.xyz[0]) < VDW_CUTOFF_FAR and \
       abs(atom1.xyz[1] - atom2.xyz[1]) < VDW_CUTOFF_FAR and \
       abs(atom1.xyz[2] - atom2.xyz[2]) < VDW_CUTOFF_FAR:
        if atom1 != atom2 and atom2 not in atom1.connect12 and atom2 not in atom1.connect13:
            d = dvv(atom1.xyz, atom2.xyz)
            if d > VDW_CUTOFF_FAR:
                p_lj = 0.0
            elif d < VDW_CUTOFF_NEAR:
                p_lj = 999.0
            else:
                if d < r0*0.9:
                    p_lj = -8*r0*(d-r0*0.9) - eps
                elif r0*0.9 <= d < r0:
                    p_lj = -eps
                elif r0 <= d < r0*2:
                    p_lj = eps*(d-r0)/r0 - eps
                else:
                    p_lj = 0
                    
                if atom2 in atom1.connect14:
                    p_lj = p_lj * VDW_SCALE14

    return p_lj
