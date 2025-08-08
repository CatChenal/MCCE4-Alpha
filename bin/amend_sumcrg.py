#!/usr/bin/env python
"""
Module: ammend_sumcrg.py

This module amends the sum_crg.out file using information from head3.lst and fort.38.
"""

import os

def amend_sum_crg():
    """
    Ammend sum_crg.out, use head3.lst and fort.38 to compose sum_crg.out.
    """
    class ConfProperties:
        def __init__(self, conf_id, crg, ne, nH):
            self.conf_id = conf_id
            self.crg = crg
            self.ne = ne
            self.nH = nH
    
    if not os.path.exists("head3.lst"):
        print("head3.lst not found, cannot compose sum_crg.out.")
        return

    if not os.path.exists("fort.38"):
        print("fort.38 not found, cannot compose sum_crg.out.")
        return

    with open("head3.lst", "r") as head_file, open("fort.38", "r") as fort_file:
        head_lines = head_file.readlines()
        fort_lines = fort_file.readlines()

    # group head3.lst lines into residues
    residues = {}
    for line in head_lines[1:]: # skip the header line
        fields = line.split()
        if len(fields) < 17:
            continue
        conf_id = fields[1]
        res_id = conf_id[:3] + conf_id[5:11]
        crg = float(fields[4])
        ne = int(fields[7])
        nH = int(fields[8])
        conf_properties = ConfProperties(conf_id, crg, ne, nH)
        if res_id not in residues:
            residues[res_id] = []
        residues[res_id].append(conf_properties)

    # convert fort.38 lines into a dictionary
    title_line = fort_lines[0]
    fields = title_line.split()
    titration_type = fields[0]
    titration_range = fields[1:]
    titration_points = len(titration_range)
    fort38_data = {}
    for line in fort_lines[1:]:  # skip the header line
        fields = line.split()
        if len(fields) != titration_points+1:
            print(f"Skipping line in fort.38 with unexpected number of values: {line.strip()}")
            continue
        conf_id = fields[0]
        occ_values = [float(v) for v in fields[1:]]
        fort38_data[conf_id] = occ_values

    # create sum_crg.out
    sum_crg_lines = [title_line]
    total_crg = [0.0 for _ in range(titration_points)]
    total_ne = [0 for _ in range(titration_points)]
    total_nH = [0 for _ in range(titration_points)]
    for res_id, conf_list in residues.items():
        # determine the ionized state charge of this residue by detecting the first none conf crg
        first_conf_crg = None
        for conf in conf_list:
            if conf.crg > 1e-4:
                first_conf_crg = "+"
                break
            elif conf.crg < -1e-4:
                first_conf_crg = "-"
                break
        if first_conf_crg is None:
            continue  # skip residues with no ionized state

        total_res_crg = [0.0 for _ in range(titration_points)]
        total_res_ne = [0 for _ in range(titration_points)]
        total_res_nH = [0 for _ in range(titration_points)]
        for conf in conf_list:
            if conf.conf_id in fort38_data:
                conf_crg = [fort38_data[conf.conf_id][i] * conf.crg for i in range(titration_points)]
                conf_ne = [fort38_data[conf.conf_id][i] * conf.ne for i in range(titration_points)]
                conf_nH = [fort38_data[conf.conf_id][i] * conf.nH for i in range(titration_points)]
                total_res_crg = [total_res_crg[i] + conf_crg[i] for i in range(titration_points)]
                total_res_ne = [total_res_ne[i] + conf_ne[i] for i in range(titration_points)]
                total_res_nH = [total_res_nH[i] + conf_nH[i] for i in range(titration_points)]
                total_crg = [total_crg[i] + conf_crg[i] for i in range(titration_points)]
                total_ne = [total_ne[i] + conf_ne[i] for i in range(titration_points)]
                total_nH = [total_nH[i] + conf_nH[i] for i in range(titration_points)]
            else:
                print(f"Warning: {conf.conf_id} not found in fort.38, skipping.")
        #add the sign to the residue id
        res_id_with_sign = res_id[:3] + first_conf_crg + res_id[3:]
        res_line = f"{res_id_with_sign:14s} {' '.join(f'{v:5.2f}' for v in total_res_crg)}\n"
        sum_crg_lines.append(res_line)
    # add total line
    sum_crg_lines.append("-"*len(title_line) + "\n")
    total_crg_line = f"{'Net_Charge':14s} {' '.join(f'{v:5.2f}' for v in total_crg)}\n"
    total_ne_line = f"{'Electrons':14s} {' '.join(f'{v:5.2f}' for v in total_ne)}\n"
    total_nH_line = f"{'Protons':14s} {' '.join(f'{v:5.2f}' for v in total_nH)}\n"
    sum_crg_lines.append(total_crg_line)
    sum_crg_lines.append(total_nH_line)
    sum_crg_lines.append(total_ne_line)

    with open("sum_crg.out", "w") as sum_crg_file:
        sum_crg_file.writelines(sum_crg_lines)


if __name__ == "__main__":
    amend_sum_crg()