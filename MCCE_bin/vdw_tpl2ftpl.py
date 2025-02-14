#!/usr/bin/env python

"""Create RADIUS parameters from 00always_needed.tpl
"""

from mcce4.args import vdw_tpl2ftp_Options
from mcce4.pdbio import *

def read_00always(folder):
    fname = "%s/00always_needed.tpl" % folder
    vdw_db = {}
    lines = open(fname).readlines()
    for line in lines:
        if line.startswith("VDW_"):
            line = line.split("#")[0].strip()
            if len(line) >= 20:
                key1 = line[:9].strip()
                key2 = line[9:15].strip()
                key3 = line[15:19]
                value_str = float(line[20:].split()[0])
                vdw_db[(key1, key2, key3)] = float(value_str) 

    return vdw_db


def update_radius(fname, db={}, tpl=None, inplace=False):

    atoms = []
    lines = open(fname).readlines()
    for line in lines:
        end = line.find("#")
        line = line[:end]
        fields = line.split(":")

        if len(fields) == 2:
            key_string = fields[0].strip()
            keys = key_string.split(",")
            key1 = keys[0].strip().strip("\"")
            if len(keys) > 1:
                key2 = keys[1].strip().strip("\"")
            else:
                key2 = ""
            if len(keys) > 2:
                key3 = keys[2].strip().strip("\"")   # No more trip() after stripping off quotes so atom names with space are preserved
            else:
                key3 = ""
            value_string = fields[1].strip()

            # Use CONNECT to identify atoms that need RADIUS record
            if key1 == "CONNECT":
                atom = (key3, key2)
                atoms.append(atom)        
            

    # Compose new RADIUS records
    new_lines = ["# %s\n" % fname]
    for atom in atoms:
        radius = RADIUS_param("0, 0, 0")
        new_key = ("RADIUS", atom[0], atom[1])
        old_key_rad = ("VDW_RAD", atom[0], atom[1])
        old_key_eps = ("VDW_EPS", atom[0], atom[1])
        if new_key in tpl.db:
            radius = tpl.db[new_key]
            print("Initilized value for %s from ftpl" % str(new_key))
        else:
            print("No RADIUS defined for %s from ftpl, 0 assumed" % str(new_key))
        if old_key_rad in db:
            radius.r_vdw = db[old_key_rad]
            print("Updated value by %s from 00always_needed.tpl" % str(old_key_rad))
        if old_key_eps in db:
            radius.e_vdw = db[old_key_eps]
            print("Updated value by %s from 00always_needed.tpl" % str(old_key_eps))

        if not (abs(radius.r_bound) < 0.001 and abs(radius.r_vdw) < 0.001 and abs(radius.e_vdw) < 0.001):
            new_lines.append("#RADIUS, %s, \"%s\": %.3f, %.3f, %.3f\n" % \
                             (atom[0], atom[1], radius.r_bound, radius.r_vdw, radius.e_vdw))
            
    if inplace:
        # clean up
        cleanup_radius(fname)
        open(fname, "a").writelines(new_lines)
    else:
        sys.stdout.writelines(new_lines)


def cleanup_radius(fname):
    keywords = [fname.split("/")[-1], "#RADIUS"]
    new_lines = []
    lines = open(fname).readlines()
    for line in lines:
        fields = line.split(",")
        test_word = fields[0].split("/")[-1].strip()
        #print(test_word, keywords)
        if test_word not in keywords:
            new_lines.append(line)
    
    open(fname, "w").writelines(new_lines)

if __name__ == "__main__":
    options = vdw_tpl2ftp_Options()
    folder = options.args.param_folder[0]
    
    vdw_db_00always = read_00always(folder)
    tpl = TPL()
    tpl.read_ftpl_folder(folder)

    files = glob.glob("%s/*.ftpl" % folder)
    files.sort()
    for fname in files:
        update_radius(fname, db=vdw_db_00always, tpl=tpl, inplace=options.args.i)
