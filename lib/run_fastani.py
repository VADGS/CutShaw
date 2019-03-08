#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import shutil
import argparse
import re
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from CutShaw.core import fileparser
from CutShaw.core import calldocker

dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4]
print(dirpath)

for root,dirs,files in os.walk(dirpath+"/PT_genomes"):
    for file in files:
        print(file)

