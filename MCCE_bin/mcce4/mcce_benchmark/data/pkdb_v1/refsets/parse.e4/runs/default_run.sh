#!/bin/bash

step1.py --dry prot.pdb
step2.py -d 4
step3.py -d 4
step4.py --xts

sleep 10
