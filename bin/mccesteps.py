#!/usr/bin/env python

import os
from datetime import datetime

record_file = "run.prm.record"
backup_file = "run.prm.backup~"

# write the entries to file run.prm
def export_runprm(runprm):

    lines = ["# WARRPER GENERATED run.prm at %s\n" % datetime.now().strftime("%Y/%m%d, %H:%M:%S")]
    for key in runprm:
        line = "%-20s    (%s)\n" % (runprm[key], key)
        lines.append(line)
    open("run.prm", "w").writelines(lines)
    return

# update the step section with new entries to run.prm.recorded
def record_runprm(runprm, section_id):
    toplines = []
    newid = "%s: %s\n" % (section_id, datetime.now().strftime("%Y/%m%d, %H:%M:%S"))
    bodylines = [newid]
    buttomlines = []
    for key in runprm:
        line = "%-20s    (%s)\n" % (runprm[key], key)
        bodylines.append(line)


    found = False
    buttom = False
    if os.path.exists(record_file):
        lines = open(record_file).readlines()
        # skip to the first line that starts with "STEP":
        new_lines = []
        istart = 0       
        for i in range(len(lines)):
            if lines[i][:5] == "#STEP":
                istart = i
                break
        for line in lines[istart:]:
            if not found:
                if line.find(section_id) == 0:
                    found = True
                else:
                    toplines.append(line)
            else:
                if buttom:
                    buttomlines.append(line)
                elif line.find("#STEP") == 0:
                    buttom = True
                    buttomlines.append(line)

    open(record_file, "w").writelines(toplines+bodylines+buttomlines)
    #print(buttomlines)

    return

def detect_runprm():
    msg = """!!! Detected an existing run.prm.
If you intend to use this run.prm, either use mcce directly or rename run.prm and use "-load_runprm" to load it.
For now, run.prm is overwritten and will be restored at the end of this program. 
"""
    detected = False
    if os.path.exists("run.prm"):
        lines = open("run.prm").readlines()
        if not lines[0].startswith("# WARRPER GENERATED"):
            print(msg)
            open(backup_file, "w").writelines(lines)
            detected = True
    
    return detected


def restore_runprm():
    if os.path.exists(backup_file):
        lines = open(backup_file).readlines()
        open("run.prm", "w").writelines(lines)


    
