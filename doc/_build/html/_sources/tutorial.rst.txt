

Tutorial
===================================

How to Update the Caldb (Using Swift as an example)

The old way:
============

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


