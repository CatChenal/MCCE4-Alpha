from mcce4.pdbio import *

def write_head1(self, fname):
    """Write head1.lst
    """

    new_lines = ["#Rotamer making site specific instruction:\n"]
    for res in self.protein.residue:
        # resID (resName, chainID, seqNum, iCode)
        line = "%3s %c%04d%c R f 00 S f 0.0 H f\n" % (res.resID[0],
                                 res.resID[1],
                                 int(res.resID[2]),
                                 res.resID[3])
        new_lines.append(line)
    
    open(fname, "w").writelines(new_lines)
    return