from ._vdw import vdw_conf

def clean_hvrot(self):
    Forbiden_vdw = float(self.prm.VDW_CUTOFF.value)  # self vdw can not exceed this value
    # remove high self energy rotamers
    backbone_confs = [res.conf[0] for res in self.protein.residue]
    #self.make_connect12() # internal atom connect12 inherited during conformer making
    self.make_connect13()
    self.make_connect14()
    self.make_blob()
    for res in self.protein.residue:
        if len(res.conf) > 2: # more than 1 side chain conformers
            for conf in res.conf[2:]:
                conf.keep = True
                if conf.history[2] == "R":
                    vdw_total = vdw_conf(conf, conf)
                    for conf_other in backbone_confs:
                        vdw_total += vdw_conf(conf, conf_other)
                        if vdw_total > Forbiden_vdw:
                            conf.keep = False
                            break
            
            res.conf = res.conf[:2] + [conf for conf in res.conf[2:] if conf.keep]
    