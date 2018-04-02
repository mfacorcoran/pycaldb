#!/usr/bin/env python
from heasarc.pycaldb import nustar_caldb as nc
caldb = '/web_chroot/FTP/caldb'
nc.update_clockcor('20180312', 'nuCclock20100101v079.fits', caldb=caldb)
