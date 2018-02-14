How to Create an HTML Summary of a caldb.indx file
===================================================

1) create a cif object for the specified mission/instrument:

``In[1] cif = pc.Caldbindx('nustar','fpm')``

2) get the available versions of the cif:

``In [2]: versions = cif.get_versions()``


3) set the version of the cif object to the cif you want (in this case the last cif in the list, usually the most recent):

``In [3]: cif.set_version(versions[-1])``

4) call the cif's html_summary method specifying the mission_url (and other info as needed):

``In [4]: mu = 'https://heasarc.gsfc.nasa.gov/docs/nustar'``

``In [5]: cif.html_summary(missionurl=mu)``