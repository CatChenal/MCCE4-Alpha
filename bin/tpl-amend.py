#!/usr/bin/env python

"""
Amend a ftpl file with some entries from an existing tpl file.
"""

import argparse


def load_db(tpl_file):
    """
    Load the RADIUS and CHARGE record from tpl file and return a dictionary of entries.
    """
    db = {}
    with open(tpl_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.strip() == "":
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            if fields[0] == "RADIUS":
                resname = fields[1]
                atomname_stripped = fields[2]
                vdw_radius = float(fields[3])
                key = ("RADIUS", resname, atomname_stripped)
                db[key] = vdw_radius
            elif fields[0] == "CHARGE":
                confname = fields[1]
                atomname_stripped = fields[2]
                charge = float(fields[3])
                key = ("CHARGE", confname, atomname_stripped)
                db[key] = charge
    return db


def ammend_ftpl(ftpl_file, db_tpl, r=False, c=False):
    """
    Amend the ftpl file with the RADIUS and CHARGE entries from the tpl file.
    """
    new_lines = []
    with open(ftpl_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                new_lines.append(line)
                continue
            if line.strip() == "":
                new_lines.append(line)
                continue
            lfields = [x.strip(",:") for x in line.split(":")]
            if len(lfields) != 2:
                new_lines.append(line)
                continue
            fields = [x.strip() for x in lfields[0].split(",")] + [x.strip() for x in lfields[1].split(",")]
            #print(fields)
            new_line = line # default
            if r and fields[0] == "RADIUS":
                resname = fields[1][:3]
                atomname_stripped = fields[2].strip('"').strip()
                key = ("RADIUS", resname, atomname_stripped)
                if key in db_tpl:
                    vdw_radius = db_tpl[key]
                    new_line = "RADIUS, %s, %s: %s, %6.3f, %s # amended\n" % (fields[1], fields[2], fields[3], vdw_radius, fields[5])
            elif c and fields[0] == "CHARGE":
                confname = fields[1]
                #print("===", fields, "===")
                atomname_stripped = fields[2].strip('"').strip()
                key = ("CHARGE", confname, atomname_stripped)
                #print(key)
            
                if key in db_tpl:
                    charge = db_tpl[key]
                    new_line = "CHARGE, %s, %s: %6.3f # amended\n" % (fields[1], fields[2], charge)

            new_lines.append(new_line)
    return new_lines


def print_conf_charges(lines):
    """
    Print the charges of each confomrer in the ftpl file.
    """
    # get conflist
    for line in lines:
        fields =line.split(",")
        if len(fields) > 2:
            key1 = fields[0].strip()
            if key1 == "CONFLIST":
                fields = line.strip().split(":")
                confs = fields[1].split(",")
                break

    # get charges for each confomer
    conf_charges = {conf.strip(): 0.0 for conf in confs}
    for line in lines:
        line = line.strip().split("#")[0]  # cut off comments
        fields = line.split(",")
        if len(fields) > 2:
            key1 = fields[0].strip()
            if key1 == "CHARGE":
                fields = line.strip().split(":")
                confname = fields[0].strip().split(",")[1].strip()
                charge_str = fields[1].strip()
                # if charge_str is a number, add it to the charge, else issue a warning
                try:
                    charge = float(charge_str)
                except ValueError:
                    print("Warning: Expect a floating point number -> %s" % line)
                    continue
                conf_charges[confname] += charge
    
    # print the charges
    print("\nCharges for each conformer:")
    for confname, charge in conf_charges.items():
        # Test if the charge is close to an integer
        if abs(charge - round(charge)) < 0.001:
            msg = ""
        else:
            msg = "Warning: Not an integer"
        print("%s: %6.3f %s" % (confname, charge, msg))


if __name__ == "__main__":
    helpmsg = "Amend a ftpl file with some entries from an existing tpl file. Atom names in tpl file are stripped names to relax the spacing rules."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--ftpl", metavar="input.ftpl", type=str, nargs=1, help="Input ftpl file.")
    parser.add_argument("--tpl", metavar="input.tpl", type=str, nargs=1, help="Input mcce tpl file.")
    parser.add_argument("--out", metavar="output.ftpl", type=str, nargs="?", default="output.ftpl", help="Output ftpl file. Default: output.ftpl")
    parser.add_argument("-r", default=False, action="store_true", help="Get RADIUS from tpl file. Default: False")
    parser.add_argument("-c", default=False, action="store_true", help="Get CHARGE from tpl file. Default: False")
    args = parser.parse_args()

    if not args.ftpl or not args.tpl:
        print("Please specify the input ftpl and tpl files or use -h for command help.")
        print("Example: tpl-amend.py --ftpl input.ftpl --tpl input.tpl --out output.ftpl")
        print("Example: tpl-amend.py --ftpl input.ftpl --tpl input.tpl --out output.ftpl -r -c")
        exit()

    # load relevant entries from the tpl file
    db_tpl = load_db(args.tpl[0])
    #print(db_tpl)

    # open the ftpl file and ammend the lines
    new_lines = ammend_ftpl(args.ftpl[0], db_tpl, r=args.r, c=args.c)

    # Check if the charges are integers
    print_conf_charges(new_lines)

    # write the new lines to the output file
    with open(args.out, "w") as f:
        for line in new_lines:
            f.write(line)
    print("\nWrote the amended ftpl file in %s" % args.out)
    print("Amended entries have comment # amended at the end of the line.")

