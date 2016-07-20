def patch_nustar_caldb(version,
                       nupatchserver = "hassif.caltech.edu",
                       nupatchworkdir = "/pub/nustar/pipeline"
                       caldbstage = "/FTP/caldb/staging")
    """
    The following is from /web_chroot/FTP/.caldb_staging/data/nustar/README_rev1.txt

    20131028
    ========

    installing patch version 20131007 - STEPS:

    a) first created a full version of the existing nustar caldb in the
    /FTP/caldb/data/nustar in the nustar/versions are via rsync:

       [caldbmgr]: pwd
       /web_chroot/FTP/.caldb_staging/data/nustar/versions
       [caldbmgr]: mkdir 20131007
       [caldbmgr]: cd 20131007
       [caldbmgr]: rsync -avz 20130509/CALDB/data/nustar/fpm/bcf .
       [caldbmgr]: rsync -avz 20130509/CALDB/data/nustar/fpm/cpf .
       [caldbmgr]: rsync -avz 20130509/CALDB/data/nustar/fpm/index .

    b) wget the patch (using the location given in Brian Grefenstette's e-mail):

       [caldbmgr]: cd /FTP/caldb/staging/data/nustar/versions/20131007/
       [caldbmgr]: wget ftp://hassif.srl.caltech.edu/pub/nustar/pipeline/20131007/NUSTARCALDB_patch20131007.tgz


    c) Untar the patch:

      [caldbmgr]: pwd
      /web_chroot/FTP/.caldb_staging/data/nustar/versions/20131007

       [caldbmgr]: ls
       fpm  NUSTARCALDB_patch20131007.tgz
       [caldbmgr]: tar -zxvf NUSTARCALDB_patch20131007.tgz
         fpm/caldb.indx
         fpm/index/caldb.indx20131007
         fpm/bcf/arf/nuA20100101v006.arf
         fpm/bcf/arf/nuB20100101v006.arf
         ... stuff deleted ...
         fpm/bcf/instrument/nuAhkrange05_20100101v005.fits
         fpm/bcf/instrument/nuAhkrange04_20100101v005.fits
         fpm/bcf/instrument/nuAhkrange03_20100101v005.fits
         fpm/bcf/instrument/nuAhkrange02_20100101v005.fits

    c) then link the (now full) 20131007 version to staging/data/nustar/fpm:

       [caldbmgr]: pwd
            /web_chroot/FTP/.caldb_staging/data/nustar
       [caldbmgr]: ln -nfs versions/20131007/fpm .
       [caldbmgr]: ls -l
            total 12
            lrwxrwxrwx 1 caldbmgr caldb         21 Oct 28 15:31 fpm -> versions/20131007/fpm
            -rw-r--r-- 1 caldbmgr caldb       1909 Oct 28 15:27 README.txt
            drwxr-xr-x 2 caldbmgr caldb       4096 Oct 28 15:09 to_ingest
            drwxr-xr-x 7 caldbmgr caldb_stage 4096 Oct 28 15:14 versions

    d) from Brian's e-mail list of new files in the patch, created caldb_update_nustar_20131007.txt in the "to_ingest" subdirectory:

         [caldbmgr]: pwd
             /web_chroot/FTP/.caldb_staging/data/nustar/to_ingest

         [caldbmgr]: ls -l
             total 16
             -rw-r--r-- 1 caldbmgr caldb 8839 Aug 14 11:27 caldb_update_nustar_20130509.txt
             -rw-r--r-- 1 caldbmgr caldb 1636 Oct 28 15:09 caldb_update_nustar_20131007.txt
             lrwxrwxrwx 1 caldbmgr caldb   32 Aug 14 11:02 caldb_update.txt -> caldb_update_nustar_20130509.txt

    e) to ingest the data (as in other missions) run /web_chroot/FTP/caldb/local/scripts/DATA_DELIVERIES/nustar_caldbupdate.csh

    :param version: string in 'YYYYMMDD' format, from nustar patch notification (from BG)
    :param nupatchurl: URL of directory where patch tar file is located, from BG
    :return:
    """

    import os
    import ftputil
    import tarfile
    import datetime

    #
    # create patch install directory
    #
    patchinstalldir = caldbstage+"/nustar/data/versions/"+version.strip().lower()+'/fpm'
    os.mkdir(patchinstalldir)
    #
    #create a full version of the existing nustar caldb in the
    #    /FTP/caldb/data/nustar in the nustar/versions are via rsync:
    cdirs=['bcf','cpf','index']
    for c in cdirs:
        cmd="rsync -avz /FTP/caldb/data/nustar/fpm/"+c+" "+patchinstalldir+"/"
        print (cmd)
        os.system(cmd)
    #
    # get the patch file
    #
    host = ftputil.FTPHost(nupatchserver,'anonymous','caldbmgr@nasa.gov')
    patchfiledir = nupatchworkdir+'/'+version.strip()
    patchfile = 'NUSTARCALDB_patch'+version.strip()+'.tgz'
    print('Downloading {0}...'.format(patchfile))
    host.download(patchfiledir+'/'+patchfile,patchinstalldir+'/'+patchfile)
    #
    # untar the file
    #
    print('Extracting {0} to {1}'.format(patchfile, patchinstalldir))
    ptf = tarfile.open(patchinstalldir+'/'+patchfile)
    ptf.extractall(patchinstalldir)
    #
    # Link the updated version to staging/data/nustar/fpm:
    #
    #
    dst = caldbstage+'/data/nustar/fpm'
    if os.path.islink(dst):
        os.remove(dst)
    print "Creating symbolic link from {0} to {1}".format(patchinstalldir, dst)
    os.symlink(patchinstalldir, dst)
    #
    # create update file in staging/nustar/to_ingest directory
    #
    update_txt='from: caldbmgr@nasa.gov\n' \
               'Subject: caldb update nustar '+version.strip()+'\n' \
               'Date: '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M") \
               'To: caldbingest@bigbang.gsfc.nasa.gov' \
               'instrument=fpm' \
               'caldbindexfile=caldb.indx'+version.strip()+'\n\n' \
               'newfile\n\n' \
               '<ENTER LIST OF NEWFILES HERE>\n\n' \
               'endnewfile'
    ingestfile = caldbstaging+'/nustar/to_ingest/caldb_update_nustar_+'version.strip()+'.txt'
    print('Creating {0}'.format(ingestfile))
    fingest = open(ingestfile,'w')
    fingest.write(update_txt)
    fingest.close()
    print('Edit {0} and enter list of new files from the update e-mail from BG'.format(ingestfile))
    return




