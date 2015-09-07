"""
pycaldb is a set of python functions which can be used to do useful things with a HEASARC-style calibration database

Version 0.1 (20131122)

pycaldb contains the following functions:
	get_cif(telescope,instrument,caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb',caldbver='',cifout="/tmp/cif")
	quizcif(telescope, instrument, cal_cnam, detnam='',cal_cbd=['','','','','','','','',''],
                          filter="none", cal_vsd="today", cal_vst="now", cal_qual=0, caldb='ftp://heasarc.gsfc.nasa.gov/caldb', caldbver='')
               get_cbds(telescope, instrument, caldb='http://heasarc.gsfc.nasa.gov/FTP/caldb')
               cmp_cbd(cbd, pname, pval, punit="", chatter=0)
               parse_cbd(bound, chatter=0)
               test_pycaldb(dummy, caldb='ftp://heasarc.gsfc.nasa.gov/caldb')
               ck_file_existence(cif,caldb="/FTP/caldb")
"""

import pycaldb
from caldb_supported_missions import *
from nustar_clockcor import *