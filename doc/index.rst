.. pycaldb documentation master file, created by
   sphinx-quickstart on Tue May  9 13:03:35 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pycaldb: Documentation
===================================

.. toctree::
      :maxdepth: 4

      Introduction: Developing pycaldb documentation with Sphinx <introduction>
      Updating the CALDB <tutorial>
      Creating CIF summary web pages using Pycaldb <pycaldb_mkhtml>
      Pycaldb List of functions <pycaldb_functions>

Pycaldb is a python module which can be used to access, check and perform other useful things on the HEASARC
Calibration DataBase (CalDB).

it should be imported as::

 In [1]: from heasarc import *
 importing pycaldb as pc
 importing chandra_caldb_update
 importing nustar_clockcor as nuc
 importing nustar_caldb_patch as nup

Note that for some reason (12/15/17) compiling the doc with ``sphinx`` in ``PyCharm`` doesn't seem to be working.  So you can create the doc from the command line by::

    % cd /software/github/heasarc/pycaldb/doc
    % make html

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
