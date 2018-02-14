

How to Update the Caldb
=======================

The Old Way
############

(Example: an update of the Swift Caldb for a new clock correction file)

Start with a update notification e-mail which has a subject line of the format ``caldb update <mission> <version>`` where ``<mission>``
is the mission name (the value of the TELESCOP keyword in the caldb.indx file) and ``<version>`` is the update date (in ``YYYYMMDD`` format).

1. Save the update e-mail (usually from Lorella) as a text file to the directory ``/FTP/caldb/staging/data/<mission>/to_ingest``, for example::

    /FTP/caldb/staging/data/swift/to_ingest

2. In ``/FTP/caldb/staging/data/<mission>/to_ingest``, make a symbolic link from the e-mail text file to the file ``caldb.update.txt``,
for example::

    caldb_update.txt -> caldb_update_swift_20170306.txt

3. In ``/FTP/caldb/local/scripts/DATA_DELIVERIES``, run the appropriate shell script to do the ingest::

    % source swift_caldbupdate.csh

4. Create the html version of the cif via ``make_missioncif_page()``::

    make_missioncif_page(telescope='Swift',instrument='SC', version='20170505',missionurl='http://swift.gsfc.nasa.gov')

5. Now update the caldb web pages on the development side (``/www/htdocs/docs/heasarc/caldb``).  Generally the html files
that need to be updated are::

    caldb_wotsnew.html
    caldb.rss
    caldb_supported_missions.html
    the astro-update associated data file (inc/associated_data.html)

But for some missions (NuSTAR, for example) there are other files that need to be updated.  For NuSTAR these files are::

    release_<VERSION>.txt
    nustar_caldbhistory.html
    nustar_caldb.html
    nustar.rss

6. Sync the web pages from the develop side to the public side using ``sync_caldbweb()``.  ``sync_caldbweb()`` also
creates the rss feeds and the "What's New" items for the caldb and HEASARC home pages, on both the public and develop side.


The New Way
###########

NuSTAR
++++++

This section describes how to update the caldb using the pycaldb modules, using the ``NuSTAR`` clock correction file update as an example.

1. the ``NuSTAR`` team will send out a notification that a new clock correction file is available. Once that's received, from the ``caldbmgr`` account, run ``update_clockcor``::

    update_clockcor('20170614', 'nuCclock20100101v072.fits', caldb=caldb)
    
I've been updating the following lines in ``nustar_clockcor.py``::

    if __name__ == "__main__":
        caldb = '/web_chroot/FTP/caldb' # appropriate for running as caldbmgr on heasarcdev
        # nucc = nustar_get_clockcor()
        # nucckeys = nucc.keys()
        # for i in range(len(nucc[nucckeys[0]])):
        #     print nucc[nucckeys[0]][i], nucc[nucckeys[1]][i]
        update_clockcor('20170614', 'nuCclock20100101v072.fits', caldb=caldb)

then running via::

    % python nustar_clockcor

1b. **A better way to do the clock correction update**:

Create a new python script in ``/Home/lhea3/caldbmgr/update_nustar_clockcor`` containing the following::

    heasarcdev-sl [32] [caldbmgr]: more update_clockcor20170720.py
        #!/usr/bin/env python
        from heasarc.pycaldb import nustar_caldb as nc
        caldb = '/web_chroot/FTP/caldb'
        nc.update_clockcor('20170720', 'nuCclock20100101v073.fits', caldb=caldb)
        #END

make sure the file is executable, then execute it::
    % ./update_clockcor20170720.py > logs/update_clockcor20170720.log

This downloads the new clock correction file, adds it to the ``CALDB`` then creates the tar files and soft link to ``caldb.indx``

2. The onerous part is updating all the web pages on the development side (/www/htdocs).  To do this:

    * edit the CALDB what's new, rss feed and supported missions pages
    * update `nustar_caldb.html <http://heasarcdev.gsfc.nas .gov/docs/heasarc/caldb/nustar/nustar_caldb.html>`_
    * Create /www/htdocs/docs/heasarc/caldb/nustar/docs/release*.txt file for current release
    * make a new entry in the `nustar_caldbhistory.html <http://heasarcdev.gsfc.nasa.gov/docs/heasarc/caldb/nustar/docs/nustar_caldbhistory.html>`_ table
    * add rss item to /www/htdocs/docs/nustar/news/nustar.rss (use RESPONSIBLE_PARTY = HEASARC)
    * add rss item to  astro-update.rss (use RESPONSIBLE_PARTY = NuSTAR for astroupdate.rss)
    * update the astro-update `inc/associated_data.html` table
    * Change the dates on the nustar caldb tar files
    * Create the html version of the caldb.indx file using the pycaldb function `make_missioncif_page('NuSTAR','FPM','20170614')` (for example).

3. Once the web updates on the develop side have been done, then run::

        sync_caldbweb(mission,instrument, version, DoCopy=False)

Set `DoCopy=False` to see the copy commands that will be run without actually executing the commands; `DoCopy=True` will actually run the commands. `sync_caldbweb` also updates the  HEASARC rss news feeds.

