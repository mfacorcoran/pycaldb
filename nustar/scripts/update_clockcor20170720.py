#!/usr/bin/env python
from heasarc.pycaldb import nustar_caldb as nc
caldb = '/web_chroot/FTP/caldb'
nc.update_clockcor('20170720', 'nuCclock20100101v073.fits', caldb=caldb)
