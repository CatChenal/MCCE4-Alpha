#!/usr/bin/env python

from mcce4.geom import *

p1 = (0, 0, 0)
p2 = (9, 9, 9)

line = LINE()
line.from2p(p1, p2)
line.print_me()