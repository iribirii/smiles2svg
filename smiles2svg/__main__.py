#!/usr/bin/env python

#########################################
#       __main__ file for the code      #
#########################################

import sys
from smiles2svg import smiles2svg

if __package__ != 'smiles2svg':
    print('smiles2svg is not installed! Use: pip install smiles2svg or python setup.py install')

if __name__ == '__main__':
    smiles2svg.main()
    sys.exit()

